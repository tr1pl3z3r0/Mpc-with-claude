function T = comparar_con_plecs(arch_plecs, arch_matlab, t1, t2)
%COMPARAR_CON_PLECS Compara los CSV de PLECS con los del modelo MATLAB.
%   T = comparar_con_plecs(arch_plecs, arch_matlab) recibe dos cell arrays de
%   5 rutas cada uno, en este orden (el mismo de los exportes de PLECS):
%     1) Gain9, Gain14            2) pαΣ, pβΣ, p0Σ, pαΔ, pβΔ, p0Δ
%     3) 6 tensiones de celda     4) iΣα, iΣβ (3ph->srf3)
%     5) estados vC: Σαβ0 y Δαβ0 (3ph->srf5)
%   Si arch_matlab es una carpeta, usa los nombres de exportar_csv.
%   t1, t2: ventana de analisis en segundos (por defecto 5 a 10 s).
%
%   Las metricas no dependen de la fase de la red ni del signo de v0, asi que
%   sirven aunque las dos simulaciones no esten sincronizadas.

if nargin < 3, t1 = 5; end
if nargin < 4, t2 = 10; end
if ischar(arch_matlab)
    c = arch_matlab;
    arch_matlab = {fullfile(c,'1_gain9_gain14.csv'), fullfile(c,'2_potencias.csv'), ...
        fullfile(c,'3_tensiones_celda.csv'), fullfile(c,'4_corriente_circ.csv'), fullfile(c,'5_estados_vc.csv')};
end
Ts = 100e-6;
P = cargar(arch_plecs, t1, t2, Ts);
M = cargar(arch_matlab, t1, t2, Ts);
mp = metricas(P); mm = metricas(M);

nombres = fieldnames(mp);
fprintf('\n%-44s %12s %12s %9s\n', 'Metrica (ventana de analisis)', 'PLECS', 'MATLAB', 'Dif. %');
for q = 1:numel(nombres)
    a = mp.(nombres{q}); b = mm.(nombres{q});
    fprintf('%-44s %12.4g %12.4g %9.1f\n', strrep(nombres{q}, '_', ' '), a, b, 100*(b - a)/max(abs(a), eps));
end
T = struct('metrica', {nombres}, 'plecs', {struct2cell(mp)}, 'matlab', {struct2cell(mm)});

% ---------------- graficos: 30 ms y 0,5 s al inicio de la ventana
try
    graficar(P, M);
catch err
    fprintf('No se pudieron dibujar los graficos: %s\n', err.message);
end
end

% =====================================================================
function graficar(P, M)
for v = 1:2
    if v == 1, dur = 0.03; else, dur = 0.5; end
    figure('Name', sprintf('PLECS vs MATLAB (%g s)', dur));
    sen = {'pαΔ [W]', @(D) D.p(:,4); 'pαΣ [W]', @(D) D.p(:,1); 'iΣα [A]', @(D) D.iS(:,1); 'VcαΔ [V]', @(D) D.xD(:,1)};
    for s = 1:size(sen,1)
        subplot(size(sen,1), 1, s); hold on; grid on;
        kp = P.t < P.t(1) + dur; km = M.t < M.t(1) + dur;
        f = sen{s,2};
        yp = f(P); ym = f(M);
        plot(P.t(kp) - P.t(1), yp(kp), 'b');
        plot(M.t(km) - M.t(1), ym(km), 'r');
        ylabel(sen{s,1});
        if s == 1, legend('PLECS', 'MATLAB'); end
    end
    xlabel('t - t_1 [s]');
end
end

% =====================================================================
function D = cargar(arch, t1, t2, Ts)
tg = (t1:Ts:t2-Ts)';
D.t = tg;
x = leer(arch{1}, tg); D.g9 = x(:,1); D.g14 = x(:,2);
D.p  = leer(arch{2}, tg);
D.vC = leer(arch{3}, tg);
D.iS = leer(arch{4}, tg);
x = leer(arch{5}, tg); D.xS = x(:,1:3); D.xD = x(:,4:6);
end

function y = leer(archivo, tg)
d = dlmread(archivo, ',', 1, 0);
[t, idx] = unique(d(:,1), 'last');      % PLECS puede repetir instantes (paso variable)
y = interp1(t, d(idx, 2:end), tg, 'linear');
end

function m = metricas(D)
Ts = D.t(2) - D.t(1); n10 = round(0.01/Ts);
amp = @(x, f) ajuste(D.t, x, f);
m.amplitud_mitad_Vdc_i_alfa_W = max(abs(D.g9));
pDa_lenta = movmean_simple(D.p(:,4), n10);
m.residuo_4Hz_p_alfa_Delta_pct = 100*amp(pDa_lenta, 4)/max(abs(D.g9));
m.pico_i_Sigma_alfa_A = max(abs(D.iS(:,1)));
m.rms_i_Sigma_alfa_A = sqrt(mean(D.iS(:,1).^2));
m.componente_8Hz_i_Sigma_alfa_A = amp(movmean_simple(D.iS(:,1), n10), 8);
m.Vc_alfa_Delta_4Hz_V = amp(D.xD(:,1), 4);
m.Vc_alfa_Delta_pp_V = max(D.xD(:,1)) - min(D.xD(:,1));
m.Vc_alfa_Sigma_8Hz_V = amp(D.xS(:,1), 8);
m.Vc_alfa_Sigma_pp_V = max(D.xS(:,1)) - min(D.xS(:,1));
m.Vc_0_Delta_medio_V = mean(D.xD(:,3));
m.Vc_0_Sigma_medio_V = mean(D.xS(:,3));
m.p_0_Delta_nivel_W = mean(abs(D.p(:,6)));
m.p_0_Sigma_medio_W = mean(D.p(:,3));
m.p_alfa_Sigma_lenta_8Hz_W = amp(movmean_simple(D.p(:,1), n10), 8);
m.rizado_celda_Pa_pp_V = max(D.vC(:,1)) - min(D.vC(:,1));
end

function a = ajuste(t, x, f)
w = 2*pi*f;
A = [cos(w*t) sin(w*t) ones(size(t))];
ok = ~isnan(x);
c = A(ok,:)\x(ok);
a = hypot(c(1), c(2));
end

function y = movmean_simple(x, n)
% Promedio movil centrado de n muestras (un periodo de v0), con NaN en los bordes.
c = cumsum([0; x]);
y = nan(size(x));
h = floor(n/2);
for k = h+1:numel(x)-n+h+1
    y(k) = (c(k-h+n) - c(k-h))/n;
end
end

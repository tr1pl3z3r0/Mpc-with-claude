function out = mmc_mpc_paper(p)
%MMC_MPC_PAPER Simula el MMC con el CCS-MPC de dos etapas conectado a la red.
%   out = mmc_mpc_paper(p) usa los parametros p (por defecto parametros_mmc()).
%
%   Todo el modelo sale del paper (Arias-Esquivel et al., IEEE TIE 2023):
%   - Planta: transformacion Sigma-Delta-alfa-beta-0 (ecs. (1)-(2)), dinamica de
%     corrientes de la ec. (3)/(4) con la red (fuente + Rl + Ll) en el puerto AC,
%     y capacitores por la potencia de cluster p = v*i (ec. (5) sin linealizar).
%     Cada cluster se modela promediado, con sus n celdas balanceadas (el paper
%     usa IPD-PWM con ordenamiento) y tension limitada a [0, n*vC] (medio puente).
%   - MPC externo: ecs. (7), (9)-(14), con restricciones de corriente (15)-(17).
%   - MPC interno: ecs. (18)-(21), con restricciones de tension (22)-(24).
%   - Puerto AC: control de corriente dq con prealimentacion (Fig. 3a, con la red
%     en lugar de la maquina). Puerto DC: control de energia total (Fig. 3b).
%   - Modo comun trapezoidal: Fig. 3a y ec. (8).
%   Pesos fijos: no se usa la adaptacion de lambda_Delta.
%
%   Convenciones: celdas y clusters en orden [Pa Pb Pc Na Nb Nc]; tensiones de
%   capacitor por submodulo; Sigma = (P+N)/2, Delta = P-N; Clarke invariante en
%   amplitud.

if nargin < 1, p = parametros_mmc(); end

Ts = p.Ts; N = round(p.Tsim/Ts); h = Ts/p.nsub;
K  = Ts/(p.n*p.C*p.vCref);                      % K de la ec. (9)
we = 2*pi*p.fe; Le = p.L + 2*p.Ll;
V0 = 0.8*p.Vdc/2*(1 - p.fe/p.fn);               % ec. (8)

% ------------------------------------------------ estado de la planta
% z = [iSigma_alfa; iSigma_beta; idc; i_alfa; i_beta; vC(6)]
z = [0; 0; 0; 0; 0; p.vC_init*ones(6,1)];

% ------------------------------------------------ estado del control
x_ac = [0; 0]; x_i = 0; x_v = 0;                % integradores de los PI
u_ext_ant = [0; 0]; u_int_ant = [0; 0];         % acciones aplicadas (retardo)
vref_apl = zeros(6,1);                          % tensiones de cluster aplicadas

% ------------------------------------------------ registro
reg = struct();
reg.t     = (0:N-1)'*Ts;
reg.p     = zeros(N,6);   % [pSa pSb pS0 pDa pDb pD0] con las ecs. (7a)-(7d)
reg.g9    = zeros(N,1);   % 1/2*Vdc*i_alfa
reg.g14   = zeros(N,1);   % -iSigma_alfa*v0Delta
reg.vC    = zeros(N,6);   % tension de submodulo [Pa Pb Pc Na Nb Nc]
reg.iS    = zeros(N,2);   % corriente circulante medida
reg.iS_ref= zeros(N,2);   % referencia del MPC externo
reg.xS    = zeros(N,3);   % vC Sigma alfa, beta, 0
reg.xD    = zeros(N,3);   % vC Delta alfa, beta, 0
reg.i     = zeros(N,2);   % corriente AC alfa, beta
reg.idq   = zeros(N,2);
reg.idc   = zeros(N,1);
reg.idc_ref = zeros(N,1);
reg.vS    = zeros(N,3);   % vSigma alfa, beta (MPC interno) y vSigma_0 (DC Port)
reg.vD    = zeros(N,3);   % vDelta alfa, beta (AC Port) y vDelta_0 (modo comun)
reg.vref  = zeros(N,6);   % tensiones de cluster pedidas, ec. (22)
reg.nact  = zeros(N,2);   % restricciones activas [externo interno] (-1: infactible)

for k = 1:N
    t = (k-1)*Ts;
    iSab = z(1:2); idc = z(3); iab = z(4:5); vC = z(6:11);

    % ---------------- mediciones
    e   = red(t, p);                                    % tension de red alfa-beta
    th  = atan2(e(2), e(1));                            % orientacion con la tension de red
    Rot = [cos(th) sin(th); -sin(th) cos(th)];          % alfa-beta -> dq
    idq = Rot*iab; edq = Rot*e;
    [xS, xD] = sigma_delta(vC(1:3), vC(4:6));
    x   = [xD(1); xD(2); xD(3); xS(1); xS(2)];          % estados de la ec. (9)

    % ---------------- puerto DC: control de energia total (Fig. 3b)
    ev  = p.vCref - xS(3);
    x_v = x_v + Ts*ev;
    idc_ref = p.kp_v*ev + p.ki_v*x_v;
    ei  = idc_ref - idc;
    x_i = x_i + Ts*ei;
    vS0 = p.Vdc/2 - (p.kp_dc*ei + p.ki_dc*x_i);         % (L/3) didc/dt = Vdc/2 - vSigma_0

    % ---------------- puerto AC: control de corriente dq (Fig. 3a con la red)
    err  = [p.id_ref; p.iq_ref] - idq;
    x_ac = x_ac + Ts*err;
    uac  = p.kp_ac*err + p.ki_ac*x_ac;
    vDdq = [-uac(1) - 2*edq(1) + we*Le*idq(2);          % (L+2Ll) di/dt = -vDelta - 2e - 2Rl i
            -uac(2) - 2*edq(2) - we*Le*idq(1)];
    vDab = Rot'*vDdq;

    % ---------------- modo comun trapezoidal (Fig. 3a, ec. (8))
    v0  = V0*trapecio(t, p);
    vD0 = -2*v0;

    % ---------------- MPC externo: capacitores (ecs. (9)-(17))
    v_ref = -vDab/2;                                    % v del puerto AC, ec. (4) sin L di/dt
    i_ref = Rot'*[p.id_ref; p.iq_ref];                  % el paper usa las referencias (Fig. 2)
    if p.retardo
        % compensacion de un paso: se predice x(k+1) con la accion que se esta aplicando
        [Bk, dk] = matrices_externo(v_ref, v0, i_ref, idc_ref, p.Vdc, K);
        xpred = x + Bk*u_ext_ant + dk;
        v0s = V0*trapecio(t + Ts, p);
        [u_ext, na_ext] = mpc_externo(xpred, v_ref, v0s, i_ref, idc_ref, p, K);
    else
        [u_ext, na_ext] = mpc_externo(x, v_ref, v0, i_ref, idc_ref, p, K);
    end

    % ---------------- MPC interno: corrientes circulantes (ecs. (18)-(24))
    vCP = p.n*vC(1:3); vCN = p.n*vC(4:6);               % tension total de capacitores por cluster
    if p.retardo
        ipred = iSab - Ts/p.L*u_int_ant;
        [u_int, na_int] = mpc_interno(ipred, u_ext, vS0, vD0, vDab, vCP, vCN, p);
    else
        [u_int, na_int] = mpc_interno(iSab, u_ext, vS0, vD0, vDab, vCP, vCN, p);
    end

    % ---------------- tensiones de cluster, ec. (22)
    vSx = clarke_inv([u_int; vS0]);
    vDx = clarke_inv([vDab; vD0]);
    vref = [vSx + vDx/2; vSx - vDx/2];

    if p.retardo
        aplicar = vref_apl; vref_apl = vref;            % lo calculado en k se aplica en k+1
        u_ext_ant = u_ext; u_int_ant = u_int;
    else
        aplicar = vref;
    end

    % ---------------- registro (valores en t_k)
    vpot = -vDab/2; v0pot = -vD0/2;                     % mismas senales que el bloque de potencias
    reg.p(k,:)  = potencias(p.Vdc, idc, iab, iSab, vpot, v0pot);
    reg.g9(k)   = 0.5*p.Vdc*iab(1);
    reg.g14(k)  = -iSab(1)*vD0;
    reg.vC(k,:) = vC';
    reg.iS(k,:) = iSab'; reg.iS_ref(k,:) = u_ext';
    reg.xS(k,:) = xS'; reg.xD(k,:) = xD';
    reg.i(k,:)  = iab'; reg.idq(k,:) = idq';
    reg.idc(k)  = idc; reg.idc_ref(k) = idc_ref;
    reg.vS(k,:) = [u_int' vS0]; reg.vD(k,:) = [vDab' vD0];
    reg.vref(k,:) = vref';
    reg.nact(k,:) = [na_ext na_int];

    % ---------------- planta: RK4 durante un periodo con tensiones de cluster constantes
    for j = 1:p.nsub
        tj = t + (j-1)*h;
        k1 = planta(tj,       z,          aplicar, p);
        k2 = planta(tj + h/2, z + h/2*k1, aplicar, p);
        k3 = planta(tj + h/2, z + h/2*k2, aplicar, p);
        k4 = planta(tj + h,   z + h*k3,   aplicar, p);
        z  = z + h/6*(k1 + 2*k2 + 2*k3 + k4);
    end
end
out = reg;
out.param = p;
end

% =====================================================================
function dz = planta(t, z, vref, p)
% Ecs. (3)-(5) con la red en el puerto AC.
vC   = z(6:11);
vapl = min(max(vref, 0), p.n*vC);                     % cluster de medio puente: 0 <= v <= n*vC
vP = vapl(1:3); vN = vapl(4:6);
s = (vP + vN)/2; d = vP - vN;
vSa = (2*s(1) - s(2) - s(3))/3; vSb = (s(2) - s(3))/sqrt(3); vS0 = (s(1) + s(2) + s(3))/3;
vDa = (2*d(1) - d(2) - d(3))/3; vDb = (d(2) - d(3))/sqrt(3);
e  = red(t, p);
Le = p.L + 2*p.Ll;
% ec. (4): vSigma + L diSigma/dt = [0 0 Vdc/2], vDelta + L di/dt + 2v = 0, v = e + Rl i + Ll di/dt
diSa = -vSa/p.L;
diSb = -vSb/p.L;
didc = 3*(p.Vdc/2 - vS0)/p.L;
dia  = (-vDa - 2*e(1) - 2*p.Rl*z(4))/Le;
dib  = (-vDb - 2*e(2) - 2*p.Rl*z(5))/Le;
% corrientes de cluster: iP = iSigma + i/2, iN = iSigma - i/2 (ec. (15))
iSx = [z(1) + z(3)/3; -z(1)/2 + sqrt(3)/2*z(2) + z(3)/3; -z(1)/2 - sqrt(3)/2*z(2) + z(3)/3];
ix  = [z(4); -z(4)/2 + sqrt(3)/2*z(5); -z(4)/2 - sqrt(3)/2*z(5)];
iarm = [iSx + ix/2; iSx - ix/2];
% ec. (5) sin linealizar: n*C*vC*dvC/dt = v*i
dvC = vapl.*iarm./(p.n*p.C*vC);
dz = [diSa; diSb; didc; dia; dib; dvC];
end

function e = red(t, p)
% Tension de red en alfa-beta, secuencia positiva.
th = 2*pi*p.fe*t;
e = p.Eg*[cos(th); sin(th)];
end

function s = trapecio(t, p)
% v0/V0: senal trapezoidal de frecuencia f0 con rampas de duracion Tr.
s = sin(2*pi*p.f0*t)/sin(pi*p.f0*p.Tr);
s = min(1, max(-1, s));
end

function [xS, xD] = sigma_delta(xP, xN)
% Ecs. (1)-(2): [Sigma; Delta] alfa-beta-0 a partir de las variables P y N (a, b, c).
xS = clarke((xP + xN)/2);
xD = clarke(xP - xN);
end

function y = clarke(x)
y = [(2*x(1) - x(2) - x(3))/3; (x(2) - x(3))/sqrt(3); (x(1) + x(2) + x(3))/3];
end

function x = clarke_inv(y)
x = [y(1) + y(3); -y(1)/2 + sqrt(3)/2*y(2) + y(3); -y(1)/2 - sqrt(3)/2*y(2) + y(3)];
end

function pw = potencias(Vdc, idc, i, iS, v, v0)
% Ecs. (7a)-(7d). (a)^c es el conjugado del producto complejo, o es el producto punto.
pSa = 0.5*Vdc*iS(1) - 0.25*(i(1)*v(1) - i(2)*v(2)) - 0.5*v0*i(1);
pSb = 0.5*Vdc*iS(2) + 0.25*(i(1)*v(2) + i(2)*v(1)) - 0.5*v0*i(2);
pS0 = Vdc*idc/6 - 0.25*(v(1)*i(1) + v(2)*i(2));
pDa = 0.5*Vdc*i(1) - 2/3*idc*v(1) - (v(1)*iS(1) - v(2)*iS(2)) - 2*v0*iS(1);
pDb = 0.5*Vdc*i(2) - 2/3*idc*v(2) + (v(1)*iS(2) + v(2)*iS(1)) - 2*v0*iS(2);
pD0 = -(v(1)*iS(1) + v(2)*iS(2)) - 2/3*idc*v0;
pw = [pSa pSb pS0 pDa pDb pD0];
end

function [B, d] = matrices_externo(v, v0, i, idc, Vdc, K)
% Ec. (9).
B = K*[-v(1) - 2*v0,  v(2);
        v(2),         v(1) - 2*v0;
       -v(1),        -v(2);
        Vdc/2,        0;
        0,            Vdc/2];
d = K*[0.5*Vdc*i(1) - 2/3*idc*v(1);
       0.5*Vdc*i(2) - 2/3*idc*v(2);
       -2/3*idc*v0;
       0.25*i(2)*v(2) - 0.25*i(1)*v(1) - 0.5*i(1)*v0;
       0.25*i(2)*v(1) + 0.25*i(1)*v(2) - 0.5*i(2)*v0];
end

function [u, nact] = mpc_externo(x, v, v0, i, idc, p, K)
% Ecs. (11)-(14): min u'Hu + 2[(Ax + d)'F]u  s.a. Gu >= W, con A = I.
[B, d] = matrices_externo(v, v0, i, idc, p.Vdc, K);
H = B'*p.Q*B + p.R;
g = B'*p.Q*(x + d);                    % (F'(Ax + d)) con F = QB
if isinf(p.imax)
    G = zeros(0,2); W = zeros(0,1);
else
    % Ecs. (16)-(17): limites de corriente de cluster
    r3 = sqrt(3);
    i1 = -idc/3 - i(1)/2;  i2 = -idc/3 + i(1)/4 - r3/4*i(2);  i3 = -idc/3 + i(1)/4 + r3/4*i(2);
    i4 = -idc/3 + i(1)/2;  i5 = -idc/3 - i(1)/4 + r3/4*i(2);  i6 = -idc/3 - i(1)/4 - r3/4*i(2);
    G = [-1 0; 0.5 -r3/2; 0.5 r3/2; 1 0; -0.5 r3/2; -0.5 -r3/2];
    W = [-p.imax - min(i1,i4); -p.imax - min(i2,i5); -p.imax - min(i3,i6);
         -p.imax + max(i1,i4); -p.imax + max(i2,i5); -p.imax + max(i3,i6)];
end
[u, nact] = qp_dos_variables(H, g, G, W);
end

function [u, nact] = mpc_interno(x, xref, vS0, vD0, vDab, vCP, vCN, p)
% Ecs. (18)-(21): x(k+1) = x(k) - Ts/L*u, min u'Hu + 2[(Ax - x*)'F]u  s.a. Gu >= W.
Bi = -p.Ts/p.L*eye(2);
H  = Bi'*p.Qi*Bi + p.Ri;
g  = Bi'*p.Qi*(x - xref);              % u* = 0 porque A = I, ec. (20)
if p.restricciones_cluster
    % Ecs. (23)-(24): 0 <= tension de cada cluster <= suma de sus tensiones de capacitor
    r3 = sqrt(3); vDa = vDab(1); vDb = vDab(2);
    v1 = -vS0 - vD0/2 - vDa/2;
    v2 = -vS0 - vD0/2 + vDa/4 - r3/4*vDb;
    v3 = -vS0 - vD0/2 + vDa/4 + r3/4*vDb;
    v4 = -vS0 + vD0/2 + vDa/2;
    v5 = -vS0 + vD0/2 - vDa/4 + r3/4*vDb;
    v6 = -vS0 + vD0/2 - vDa/4 - r3/4*vDb;
    G = [-1 0; 0.5 -r3/2; 0.5 r3/2; 1 0; -0.5 r3/2; -0.5 -r3/2];
    W = [-min(v1 + vCP(1), v4 + vCN(1)); -min(v2 + vCP(2), v5 + vCN(2)); -min(v3 + vCP(3), v6 + vCN(3));
          max(v1, v4); max(v2, v5); max(v3, v6)];
else
    G = zeros(0,2); W = zeros(0,1);
end
[u, nact] = qp_dos_variables(H, g, G, W);
end

function [u, nact] = qp_dos_variables(H, g, G, W)
% Optimo exacto de min u'Hu + 2g'u s.a. Gu >= W con u de 2 componentes.
% Con 2 variables el optimo tiene como mucho 2 restricciones activas: se resuelve
% el sistema KKT de cada conjunto activo posible y se elige el factible de menor
% costo (problema estrictamente convexo => optimo global). Si no hay ninguno
% factible, se usa el de menor violacion total y nact = -1.
m = size(G,1);
u = -H\g; nact = 0;
if m == 0 || all(G*u >= W - 1e-9), return; end
conj = {};
for a = 1:m, conj{end+1} = a; end %#ok<AGROW>
for a = 1:m, for b = a+1:m, conj{end+1} = [a b]; end, end %#ok<AGROW>
mejor = []; mejor_c = inf; mejor_v = inf; fact_ok = false;
for q = 1:numel(conj)
    S = conj{q}; ns = numel(S);
    M = [2*H, -G(S,:)'; G(S,:), zeros(ns)];
    if rcond(M) < 1e-14, continue; end
    sol = M\[-2*g; W(S)];
    uc = sol(1:2);
    viol = sum(max(W - G*uc - 1e-9, 0));
    costo = uc'*H*uc + 2*g'*uc;
    if viol == 0
        if ~fact_ok || costo < mejor_c
            mejor = uc; mejor_c = costo; fact_ok = true; nact = ns;
        end
    elseif ~fact_ok && (viol < mejor_v || (viol == mejor_v && costo < mejor_c))
        mejor = uc; mejor_v = viol; mejor_c = costo; nact = -1;
    end
end
if ~isempty(mejor), u = mejor; end
end

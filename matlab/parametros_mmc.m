function p = parametros_mmc()
%PARAMETROS_MMC Parametros del modelo MATLAB del MMC con CCS-MPC.
%   Los valores numericos son los mismos de la simulacion PLECS ("trapezoide y
%   con 2V0"); la estructura del modelo y del control sale solo del paper:
%   Arias-Esquivel et al., "Continuous Control Set Model Predictive Control of
%   a Modular Multilevel Converter for Drive Applications", IEEE TIE, 2023.

% ---------------------------------------------------------------- Convertidor
p.Vdc     = 450;        % tension del puerto DC [V]
p.n       = 3;          % submodulos por cluster
p.C       = 1;          % capacitancia de cada submodulo [F]
p.vCref   = 150;        % referencia de tension de cada submodulo, vC* [V]
p.vC_init = 150;        % tension inicial de cada submodulo [V]
p.L       = 2.5e-3;     % inductancia de cluster [H]

% ---------------------------------------------------------------- Red (puerto AC)
p.Eg = 20;              % amplitud de la tension de fase de la red [V]
p.fe = 4;               % frecuencia de la red [Hz]
p.Rl = 1;               % resistencia de linea [ohm]
p.Ll = 1e-3;            % inductancia de linea [H]
p.id_ref = 16;          % referencia de corriente d (orientada con la tension de red) [A]
p.iq_ref = 2;           % referencia de corriente q [A]

% ---------------------------------------------------------------- Simulacion
p.Ts   = 100e-6;        % periodo de muestreo del control [s]
p.Tsim = 10;            % tiempo simulado [s]
p.nsub = 5;             % subpasos RK4 de la planta por periodo de muestreo

% ---------------------------------------------------------------- Modo comun, ec. (8) y Fig. 3a
p.f0 = 100;             % frecuencia de v0 [Hz]
p.Tr = 1e-3;            % duracion de las rampas del trapecio [s]
p.fn = 100;             % frecuencia nominal en la ec. (8): V0 = 0.8*Vdc/2*(1 - fe/fn) = 172.8 V

% ---------------------------------------------------------------- MPC externo, ecs. (11)-(12) (pesos fijos)
p.Q = diag([25 25 1 1 1]);   % [lambda_alfa^Delta, lambda_beta^Delta, lambda_0^Delta, lambda_alfa^Sigma, lambda_beta^Sigma]
p.R = diag([1e-5 1e-5]);     % [lambda_alfa^U, lambda_beta^U]
p.imax = inf;                % limite de corriente de cluster, ecs. (15)-(17); inf = sin limite

% ---------------------------------------------------------------- MPC interno, ec. (19)
p.Qi = diag([1 1]);          % [lambda_alfa^Sigma', lambda_beta^Sigma']
p.Ri = diag([1e-3 1e-3]);    % [lambda_alfa^U', lambda_beta^U']
p.restricciones_cluster = true;   % ecs. (22)-(24)

% ---------------------------------------------------------------- Retardo de computo
p.retardo = false;      % true: el control calculado en k se aplica en k+1, con
                        % compensacion de un paso en ambos MPC (ref. [33] del paper)

% ---------------------------------------------------------------- Puertos AC y DC (Fig. 3)
% Ganancias disenadas desde el modelo de la ec. (4), no copiadas de PLECS.
Le = p.L + 2*p.Ll;                 % inductancia vista por la corriente AC
wc = 2*pi*100;                     % ancho de banda del lazo de corriente AC [rad/s]
p.kp_ac = Le*wc;                   % cancela el polo (L+2Ll)s + 2Rl
p.ki_ac = 2*p.Rl*wc;
wi = 2*pi*200;                     % ancho de banda del lazo de idc [rad/s]
p.kp_dc = (p.L/3)*wi;              % planta (L/3) s
p.ki_dc = p.kp_dc*wi/5;
Kv = p.Vdc/(6*p.n*p.C*p.vCref);    % d(vC0Sigma)/dt = Kv*idc + ... (ecs. (6) y (7d))
wv = 10; zv = 0.7;                 % lazo de energia: wn = 10 rad/s, zeta = 0.7
p.kp_v = 2*zv*wv/Kv;
p.ki_v = wv^2/Kv;
end

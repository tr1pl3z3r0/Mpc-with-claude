% EJECUTAR_MODELO Simula el modelo del paper y lo compara con PLECS.
%   1) Ajusta parametros si hace falta (ver parametros_mmc.m).
%   2) Pon en 'plecs' las rutas de tus 5 CSV exportados de PLECS, en el orden
%      Gain9/Gain14, potencias, tensiones de celda, iΣ (srf3), estados vC (srf5).
%   Funciona en MATLAB y en GNU Octave.

p = parametros_mmc();
% p.restricciones_cluster = false;   % para comparar con un PLECS sin restricciones
% p.imax = 18;                       % limites de corriente de cluster, ecs. (15)-(17)
% p.retardo = true;                  % retardo de computo de un paso con compensacion

out = mmc_mpc_paper(p);
exportar_csv(out, 'resultados_matlab');

plecs = {'1.csv', '2.csv', '3.csv', '4.csv', '5.csv'};   % <- tus archivos de PLECS
if all(cellfun(@(f) exist(f, 'file') == 2, plecs))
    comparar_con_plecs(plecs, 'resultados_matlab', 5, 10);
else
    disp('Pon las rutas de los CSV de PLECS en la variable plecs para comparar.');
end

function exportar_csv(out, carpeta)
%EXPORTAR_CSV Guarda los resultados con el mismo formato que los CSV de PLECS.
%   exportar_csv(out, carpeta) escribe en la carpeta:
%     1_gain9_gain14.csv     "Time / s","Gain9","Gain14"
%     2_potencias.csv        "Time / s","pαΣ","pβΣ","p0Σ","pαΔ","pβΔ","p0Δ"
%     3_tensiones_celda.csv  "Time / s", 6 x "Vm:Measured voltage" (Pa Pb Pc Na Nb Nc)
%     4_corriente_circ.csv   "Time / s","3ph->srf3:1:1","3ph->srf3:1:2"
%     5_estados_vc.csv       "Time / s","3ph->srf5:1:1",... (Sigma alfa,beta,0 y Delta alfa,beta,0)

if nargin < 2, carpeta = 'resultados_matlab'; end
if ~exist(carpeta, 'dir'), mkdir(carpeta); end
t = out.t;
escribir(fullfile(carpeta, '1_gain9_gain14.csv'), {'Gain9','Gain14'}, [t out.g9 out.g14]);
escribir(fullfile(carpeta, '2_potencias.csv'), ...
    {'pαΣ','pβΣ','p0Σ','pαΔ','pβΔ','p0Δ'}, [t out.p]);
escribir(fullfile(carpeta, '3_tensiones_celda.csv'), repmat({'Vm:Measured voltage'}, 1, 6), [t out.vC]);
escribir(fullfile(carpeta, '4_corriente_circ.csv'), {'3ph->srf3:1:1','3ph->srf3:1:2'}, [t out.iS]);
escribir(fullfile(carpeta, '5_estados_vc.csv'), ...
    {'3ph->srf5:1:1','3ph->srf5:1:2','3ph->srf5:1:3','3ph->srf5:3:1','3ph->srf5:3:2','3ph->srf5:3:3'}, ...
    [t out.xS out.xD]);
end

function escribir(archivo, nombres, datos)
f = fopen(archivo, 'w');
fwrite(f, uint8(unicode2native(['"Time / s"' sprintf(',"%s"', nombres{:}) sprintf('\n')], 'UTF-8')));
fmt = ['%.10g' repmat(',%.10g', 1, size(datos,2)-1) '\n'];
fprintf(f, fmt, datos');
fclose(f);
end

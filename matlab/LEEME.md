# Modelo MATLAB del MMC con CCS-MPC

Modelo independiente de PLECS, construido solo con las ecuaciones de:

> Arias-Esquivel, Cárdenas, Urrutia, Díaz, Tarisciotti y Clare, «Continuous Control Set Model Predictive Control of a Modular Multilevel Converter for Drive Applications», IEEE TIE, vol. 70, n.º 9, 2023.

Usa la red y los mismos valores numéricos que la simulación de PLECS, con pesos fijos (sin λΔ adaptativo). Funciona en MATLAB y en GNU Octave.

## Archivo

Todo está en `mmc_paper.m`: la función principal y, como subfunciones, los parámetros (`parametros`), la simulación (`simular`), la exportación a CSV con los encabezados de PLECS (`exportar_csv`) y la comparación (`comparar`).

## Qué sale del paper

| Parte | Ecuaciones |
|---|---|
| Transformación ΣΔαβ0 | (1)–(2) |
| Corrientes de brazo, con la red (fuente + Rl + Ll) en el puerto AC | (3)–(4) |
| Capacitores: n·C·vC·dvC/dt = v·i por cluster | (5), sin linealizar |
| Potencias | (7a)–(7d) |
| MPC externo (capacitores) y límites de corriente de cluster | (9)–(14) y (15)–(17) |
| MPC interno (corrientes circulantes) y límites de tensión de cluster | (18)–(21) y (22)–(24) |
| Modo común trapezoidal | Fig. 3a y ec. (8) |
| Puerto AC: control de corriente dq con prealimentación | Fig. 3a, con la red en lugar de la máquina |
| Puerto DC: control de energía total | Fig. 3b |

Cada cluster se modela promediado, con sus celdas balanceadas, como supone el paper al usar IPD-PWM con ordenamiento. Su tensión queda limitada a [0, n·vC] (medio puente), de forma que la energía se conserva.

## Lo que el paper no fija y se eligió aquí

- **Ganancias de los PI:**
  - AC: 100 Hz de ancho de banda, cancelando el polo de (L + 2Ll)s + 2Rl.
  - idc: 200 Hz.
  - Energía: ωn = 10 rad/s, ζ = 0,7.
- **Valores de la simulación:**
  - rampa del trapecio Tr = 1 ms y f0 = 100 Hz;
  - V0 = 0,8·Vdc/2·(1 − fe/fn) con fn = 100 Hz, que da 172,8 V;
  - referencias id = 16 A e iq = 2 A.
- **Entradas del MPC externo:** usa las referencias de corriente AC, de idc y de tensión AC, como en la Fig. 2.
- **Retardo de cómputo:** por defecto el control se aplica en el mismo período, como en PLECS. Con `p.retardo = true` se aplica en el siguiente, con la compensación de un paso de la ref. [33].
- **Solución de los QP:** exacta, enumerando los conjuntos activos posibles. Da el mismo óptimo que cualquier método de conjunto activo.

## Uso

```matlab
out = mmc_paper();                                  % simula y guarda los CSV en resultados_matlab
out = mmc_paper({'1.csv','2.csv','3.csv','4.csv','5.csv'});   % además compara con PLECS (5 a 10 s)
out = mmc_paper({}, 'restricciones_cluster', false); % cualquier campo de parametros, en pares nombre, valor
out = mmc_paper(plecs, 'imax', 18, 't1', 5, 't2', 10, 'carpeta', 'otra_carpeta');
```

Los cinco CSV de PLECS van en este orden:
1. Gain9 y Gain14.
2. Las seis potencias.
3. Las seis tensiones de celda.
4. iΣαβ (3ph->srf3).
5. Estados vC (3ph->srf5).

Para cambiar un valor de forma permanente, edítalo en la subfunción `parametros`.

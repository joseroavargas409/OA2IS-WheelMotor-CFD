# Sensores de aceite/aire

La simulación incluye **20 sensores virtuales capacitivos** que leen el campo de aceite de la
solución CFD en cada paso de tiempo:

* **S1–S8**: montajes OA²IS a 5–10 mm de paredes **fijas** (spindle, coronas), donde se podría
  instalar un sensor real (`scenarios/_sensors.yaml`).
* **V1–V12**: las 12 sondas de `komatsu-oil-sensors-v4-1` (`dashboard/constants.js`), llevadas al
  marco del modelo CFD (`scenarios/_sensors_v41.yaml`):
  x_cfd = x_v41, y_cfd = z_v41 (vertical), z_cfd = −1,3135 + 0,9356·y_v41 (axial, +y_v41 = etapa alta).
  V1, V5 y V9 caían dentro de piezas fijas de esta geometría y se reubicaron al fluido más cercano
  (15, 15 y 62 mm).

## Modelo del sensor

`oa2is_cfd/oil_sensor.py` es una traducción a Python del modelo de
`komatsu-oil-sensors-v4-1/dashboard/sensor-model.js` (licencia MIT; copia sin cambios en
`third_party/`). Se verificó que reproduce exactamente la capacitancia y el generador de ruido del
original (mismos valores en 4 casos de prueba y en la secuencia aleatoria).

| Elemento | Valor (igual a V4.1) |
|---|---|
| Electrodo coaxial | a = 4 mm, b = 6 mm, L = 25 mm, vertical; ε_aceite 2,2, ε_aire 1,0; 8 pF parásitos |
| Volumen sensible | 9 puntos: centro, ±6 mm en x y en z, ±6,25 y ±12,5 mm a lo largo del electrodo |
| Puntos dentro de un sólido | se excluyen; si los 9 están en sólido → **no disponible** (nunca «aire») |
| Película residual | objetivo 50 µm × fracción de aceite; se forma en 0,1 s y drena en 0,7 s |
| Electrónica | primer orden, τ = 0,08 s; ruido gaussiano 0,01 pF (xorshift32, semillas de V4.1) |
| Señal y estado | s = (C − C_seco)/(C_mojado − C_seco); seco < 0,2 ≤ mixto < 0,8 ≤ mojado; «mojado» si s > 0,6 |

Diferencia con V4.1: allí el aceite eran partículas SPH; aquí la fracción de aceite sale de la
solución VOF de OpenFOAM (sondas `probes` de `alpha.oil` y de la fracción de sólido `chi` en los
9 puntos de cada sensor = 180 puntos registrados en cada paso). No se calculan la separación
«aceite conectado / salpicado» ni la columna de nivel de V4.1, que dependen de las partículas.

## Cómo agregar o mover sensores

1. Agregar líneas a `scenarios/_sensors.yaml` (o a un archivo propio) con `id`, `name`, `position`.
2. Preparar el caso (los 9 puntos por sensor se escriben en `system/controlDict` y en
   `sensor_layout.json`) y correrlo.
3. `python3 tools/oil_sensors.py <caso> [--history <caso_anterior>]` → `oil_sensors.csv/.parquet`,
   `oil_sensors_summary.csv`, gráficos de capacitancia y de estados.
4. `python3 tools/render_with_sensors.py results/<caso>/oil_sensors.csv <dir_frames> <caso_anterior> <caso>`
   → video con cada sensor dibujado con el color de su estado.

Sin volver a simular, `--history` también lee sensores nuevos sobre los campos ya guardados
(cada 0,01 s); la resolución temporal queda en 0,01 s en ese tramo.

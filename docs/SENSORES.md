# Sensores de aceite/aire

La simulación incluye **20 sensores virtuales capacitivos fijos** (más 8 en la tapa de maza en los escenarios 06–09) que leen el campo de aceite de la
solución CFD en cada paso de tiempo:

* **S1–S8**: montajes OA²IS a 5–10 mm de paredes **fijas** (spindle, coronas), donde se podría
  instalar un sensor real (`scenarios/_sensors.yaml`).
* **V1–V12**: las 12 sondas de `komatsu-oil-sensors-v4-1` (`dashboard/constants.js`), llevadas al
  marco del modelo CFD (`scenarios/_sensors_v41.yaml`):
  x_cfd = x_v41, y_cfd = z_v41 (vertical), z_cfd = −1,3135 + 0,9356·y_v41 (axial, +y_v41 = etapa alta).
  V1, V5 y V9 caían dentro de piezas fijas de esta geometría y se reubicaron al fluido más cercano
  (15, 15 y 62 mm).

* **P1–P8** (escenarios 06–09): un sensor en cada **puerto de tapón magnético de la tapa de maza**,
  el lugar donde se controla el nivel en la máquina real (guía de Komatsu). Ver la sección siguiente.

## Sensores en los tapones de la tapa de maza (giran con la maza)

| id | ángulo en reposo | rol en la máquina real |
|---|---|---|
| P1, P2, P3, P4, P5, P8 | 30°, 60°, 120°, 150°, 210°, 330° | tapones magnéticos |
| **P6, P7** | **240°, 300°** | **par de control de nivel** (alineados horizontalmente con el camión detenido) |

Círculo de r = 0,270 m (medido de las fotografías de la guía), 15 mm dentro de la cara interior de la
tapa (z = −1,7195 m), electrodo a lo largo del eje de la rueda (como un sensor roscado en el puerto
del tapón). Ángulos medidos desde +x (adelante) hacia +y (arriba); 270° = abajo.

La tapa gira con la maza, así que los sensores también. La solución CFD registra un **anillo fijo de
120 posiciones** (cada 3°, 14 mm de arco, ≤ 1 celda) con los 9 puntos del electrodo en cada una
(1 080 puntos, función `ring_hub_cover_plugs`). En cada paso, `tools/plug_sensors.py` calcula el giro
de la maza θ(t) = ∫ 0,030667·ω_rotor dt (la misma cinemática que el solver aplica a las paredes de
la maza), ubica cada tapón en φ₀ + θ(t), interpola el anillo en ese ángulo y aplica el mismo modelo
capacitivo. Con el camión detenido, P6 y P7 quedan justo en la superficie libre (lectura «mixto»):
es la condición de nivel correcto de la guía.

```bash
python3 tools/plug_sensors.py <caso>      # plug_sensors.csv/.parquet, plug_sensors_summary.csv,
                                          # plug_ring.csv, plots/plug_sensors.png, plots/plug_ring.png
python3 tools/render_with_sensors.py results/<caso>/oil_sensors.csv <frames> <caso> \
        --plugs results/<caso>/plug_sensors.csv   # video con un 4.º panel en la tapa
```

`plug_ring.csv` guarda además la fracción de aceite a lo largo de todo el círculo de tapones: sirve
para evaluar cualquier otra posición de montaje en la tapa sin volver a simular.

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

# Informe de validación

Todo lo que sigue sale de ejecuciones reales del solver `interIBFoamOA2IS` (OpenFOAM v1912
parcheado) en este repositorio; los números se regeneran con `tools/postprocess.py` (por caso) y
`tools/compare.py` (comparaciones), y están en `results/*/summary.json` y
`results/validation/validation_summary.json`. Máquina: 2 núcleos, 8 GB, 1 proceso por caso.

## 1. Corridas ejecutadas

| Caso | Malla (celdas) | Escenario | t simulado | Pasos | Δt medio | Co máx / p95 | Deriva máx. del aceite |
|---|---|---|---|---|---|---|---|
| `test_01_reposo` | ib_test (145 934) | detenido, 95 L | 0,25 s | 137 | 1,8 ms | 0,23 / 0,23 | 4 × 10⁻⁶ % |
| `test_02_operacion_normal` | ib_test | 20 km/h, 887 rpm | 0,25 s | 1 725 | 0,149 ms | 1,74 / 1,02 | 7 × 10⁻⁵ % |
| `test_03_operacion_inclinada` | ib_test | + pendiente 10 %, peralte 4° | 0,25 s | 1 563 | 0,160 ms | 1,62 / 1,03 | 8 × 10⁻⁵ % |
| `test_04_transitorio` | ib_test | +1,5 / −3,0 m/s², rpm de v(t) | 0,40 s | 2 777 | 0,144 ms | 1,85 / 1,02 | 1,4 × 10⁻⁴ % |
| `test_05_bajo_inventario` | ib_test | 55 L | 0,25 s | 1 443 | 0,173 ms | 1,53 / 1,05 | 8 × 10⁻⁵ % |
| `dt_02_Co05` | ib_test | nominal con Co_max 0,5 | 0,15 s | 1 286 | 0,117 ms | 1,29 / 0,54 | 1,6 × 10⁻⁴ % |
| `mid_02_operacion_normal` | ib_mid (477 386) | nominal | 0,15 s | 1 484 | 0,103 ms | 1,49 / 0,97 | 1,2 × 10⁻⁴ % |

Residuo final de `p_rgh` ≤ 1 × 10⁻⁷ en todos los pasos de todas las corridas; error de continuidad
local ≤ 10⁻⁷. Costo: ~1,5–2 s por paso en `ib_test` y ~7–8 s en `ib_mid` (1 núcleo).
Los videos (`results/*/animation.mp4`, cortes 2D; `animation_3d.mp4`, vista 3D con corte) están
hechos con los campos escritos cada 0,01 s por el solver.

## 2. Geometría y movimiento

* **Engranajes sin interferencia** (`tests/test_gears.py`, 60 poses en un periodo de engrane, 3
  satélites por etapa): solape 0,000 mm² solar–satélite y satélite–corona; holgura mínima
  1,88 mm (alta) y 1,41 mm (baja). PASSED.
* **Cinemática en OpenFOAM** (`tests/test_kinematics_openfoam.py`, malla dinámica con órbita+giro):
  error máximo 1,5 × 10⁻⁸ m. PASSED.
* **Posición de los cuerpos dentro de la solución CFD** (`tests/test_ib_motion.py`): se muestrea el
  campo χ escrito por el solver sobre el círculo de órbita de cada etapa y se comparan los ángulos
  de los 3 satélites con la integral analítica de ω(t):

  | Caso | Giro del porta alto al final | Error máximo | Tamaño de celda | Resultado |
  |---|---|---|---|---|
  | test_02 (ω constante) | −239,5° | 3,9 mm | 18 mm | PASSED |
  | test_04 (ω variable con v(t)) | — | 6,7 mm | 18 mm | PASSED |
  | mid_02 | — | 4,3 mm | 12 mm | PASSED |

  Los engranajes están, en cada instante, donde dice la cinemática (error < ½ celda).

## 3. Conservación del aceite

El volumen integrado de aceite (∫α dV, cada paso) cambia como máximo 1,6 × 10⁻⁴ % (0,00015 L de
95 L) en todas las corridas, incluidas las de engranajes en movimiento y el transitorio
(`results/validation/inventory_all.png`). El desplazamiento de aceite por los cuerpos es
conservativo por construcción.

Matiz importante (`oil_excess_in_bodies_L`): parte del aceite queda en la fracción de celda
ocupada por un cuerpo (en celdas parcialmente cubiertas, o dentro de cuerpos sumergidos donde no
hay espacio libre al lado para devolverlo). Es aceite contado en el inventario pero no "libre":

| | t = 0,15 s | final |
|---|---|---|
| reposo (01) | 1,16 L (condición inicial en celdas parciales) | 1,29 L |
| nominal ib_test (02) | 2,60 L | 3,18 L (0,25 s) |
| nominal ib_mid | **1,73 L** | — |
| transitorio (04) | — | 4,32 L (0,40 s) |

Disminuye con la malla (2,60 → 1,73 L al pasar de 18 a 12 mm), como se espera de un error de
discretización de la frontera inmersa. Equivale a bajar el nivel libre ~5 mm. Es la principal
limitación cuantitativa del método en estas mallas.

## 4. Mecanismos detenidos frente a en movimiento (01 vs 02)

Promedios en 0,05–0,25 s (`results/validation/stopped_vs_moving.png`):

| Magnitud | Detenido | En movimiento |
|---|---|---|
| Velocidad media del aceite | 0,007 m/s | 0,74 m/s (≈ 100×) |
| Energía cinética del aceite (0,10–0,25 s) | 0,014 J | 117 J |
| Aceite por encima del eje | 0 | 3,5 % (7,2 % al final, creciendo) |
| Centroide y del aceite | −0,347 m | −0,322 m (sube 25 mm) |
| Fracción mojada, etapa alta (0,10–0,25 s) | 0,177 | 0,094 |
| Fracción mojada, etapa baja | 0,190 | 0,222 |
| Ángulo del pozo de aceite | +2,0° (fijo) | −5,3° (arrastrado en el sentido de giro) |

Con los mecanismos detenidos el aceite queda en reposo hidrostático (las velocidades espurias en el
aire junto a la superficie libre llegan a ~1 m/s, típicas de VOF con razón de densidades 850 en una
malla no alineada; en el aceite la velocidad media es 7 mm/s). Con los engranajes en movimiento el
aceite es arrastrado, lanzado por encima del eje y redistribuido: los satélites de la etapa alta lo
expulsan hacia el anillo exterior (la zona de engranajes de la etapa alta queda menos mojada),
mientras que en la etapa baja (más lenta, 0,03 × rotor) aumenta el mojado. Sensores
(fracción del tiempo con aceite, t > 0,1 s):

| Sensor | Detenido | En movimiento |
|---|---|---|
| S1 fondo etapa alta | 100 % | 53 % (salpicadura intermitente) |
| S4 superior etapa alta | 0 % | 4 % (gotas lanzadas) |
| |U| medio en S4 | 0,003 m/s | 5,1 m/s |

## 5. Sensibilidad de malla y de paso de tiempo (caso nominal, 0–0,15 s)

`results/validation/mesh_sensitivity.png`, `dt_sensitivity.png`; diferencias entre medias en
0,05–0,15 s:

| Magnitud global | ib_test vs ib_mid (18 → 12 mm, 146 k → 477 k celdas) | Co_max 1,0 vs 0,5 |
|---|---|---|
| Fracción mojada zona engranajes | 0,9 % | 1,1 % |
| Aceite sobre el eje | 3,7 % | 9,8 % |
| Velocidad media del aceite | 2,3 % | 2,0 % |
| Centroide y del aceite | 0,2 % | 0,2 % |

Las magnitudes **globales** de distribución del aceite son insensibles a la malla y al paso de
tiempo en este rango (≤ 4 %, salvo la pequeña fracción de aceite sobre el eje, que es la más
sensible). Las magnitudes **locales** del engrane no lo son: el máximo de velocidad del fluido
(chorros de aire/aceite donde un diente de satélite entra en la corona) es 29 m/s en ib_test y
45 m/s en ib_mid a t = 0,15 s, y controla el paso de tiempo. Con 1,5–2 celdas por diente el
bombeo y el atrapamiento de aceite entre dientes no están resueltos; para eso harían falta las
mallas `ib_medium`/`ib_fine` (definidas, no corridas por costo: ~10–25× el costo de ib_mid).

## 6. Respuesta a inclinación y aceleración

`results/validation/inertial_response.png`: ángulo del centroide del aceite alrededor del eje
(desde −y hacia +x) frente a la dirección de la gravedad efectiva g − a(t).

* **Pendiente 10 % + peralte 4° (03)**: el pozo de aceite se desplaza −5,2° respecto del caso
  nominal (pendiente: −5,7° cuasi-estático; nariz arriba → aceite hacia atrás) y el centroide axial
  se mueve 16 mm hacia afuera (lado bajo del peralte). Sensores: S5 (delantero) queda seco y S6
  (trasero) mojado. Cualitativa y cuantitativamente coherente.
* **Transitorio (04)**: durante la aceleración (+1,5 m/s²) el pozo sigue al nominal; desde el frenado
  (−3 m/s², g_eff girada +17°) el pozo gira hacia adelante: +4,7° respecto del nominal a los 0,25 s de
  frenado, todavía subiendo. El aceite responde con la inercia de una masa de ~95 kg que oscila con
  periodo del orden de 1 s, por lo que en 0,25 s no alcanza el ángulo cuasi-estático: es la respuesta
  dinámica esperada, no un error. S6 (trasero) pasa de 100 % a 70 % de tiempo mojado.

## 7. Bajo inventario (05, 55 L)

El nivel inicial baja 90 mm. Los sensores laterales S5 y S6 (y = −0,25 m) pasan a estar secos todo
el tiempo (100 % mojados con 95 L), S1 y S2 se secan parcialmente y aparecen estados ALERTA (9 % del
tiempo). La fracción mojada de la zona de engranajes cae a 0,04 (alta) y 0,12 (baja) frente a 0,09 y
0,22 con 95 L, y la energía cinética del aceite baja de 117 a 43 J. El conjunto S1/S2/S5/S6 distingue
claramente ambos inventarios; la regla de estado demostrativa (solo sensores de fondo) es poco
sensible y debería incluir S5/S6.

## 8. Integridad de los repositorios originales

`tools/verify_originals.sh` (solo lectura: `git --no-optional-locks`, `sha256sum -c`), contra la línea
base registrada antes de empezar (`docs/integrity/`):

| Repositorio | HEAD | `git status` | SHA-256 de archivos | HEAD en GitHub |
|---|---|---|---|---|
| wheel-motor-inspector | 64bd32e = línea base | limpio | 54/54 idénticos | = línea base |
| komatsu-oil-sensors-v4 | 675d924 = línea base | limpio | 2 592/2 592 idénticos | = línea base |
| komatsu-oil-sensors-v4-1 | 8455afb = línea base | limpio | 2 920/2 920 idénticos | = línea base |

Ningún repositorio original fue modificado (ni localmente ni en GitHub).

## 9. Corrida representativa con sensores de aceite/aire (`rep_02_test`)

Continuación del caso nominal (02, malla ib_test, 20 km/h) desde 0,25 s hasta **0,554 s**
(2 procesos MPI, 52 min de cálculo, 2 778 pasos) con los **20 sensores capacitivos** de
`docs/SENSORES.md` registrados en cada paso (180 puntos). Se eligió la malla ib_test porque las
magnitudes globales que leen los sensores no cambian con la malla (§5) y así se alcanzan 0,55 s de
tiempo físico en una hora; la malla media llegaba a ~0,2 s.

| t [s] | aceite sobre el eje | fracción mojada engranajes | velocidad media del aceite | aceite en celdas de cuerpos |
|---|---|---|---|---|
| 0,25 | 7,2 % | 0,154 | 0,75 m/s | 3,2 L |
| 0,35 | 7,9 % | 0,162 | 0,79 m/s | 2,7 L |
| 0,45 | 9,3 % | 0,156 | 0,85 m/s | 3,8 L |
| 0,554 | 11,3 % | 0,141 | 0,89 m/s | 4,9 L |

Inventario: deriva máxima 4,9 × 10⁻⁴ % en 0,554 s. El aceite sigue siendo levantado por los
engranajes (todavía no hay régimen estacionario a 0,55 s).

Lo que leen los sensores (fracción del tiempo entre 0,05 y 0,554 s; `oil_sensors_summary.csv`):
* siempre mojados: S3 (fondo bolsillo tapa), S7 (fondo anillo corona alta), V3 (fondo del cárter),
  V7 (canal central inferior), V10 (zona de salida);
* pasan de mojado a mixto a medida que el aceite es arrastrado: S1 (fondo etapa alta, desde 0,18 s),
  S2, S5 y S6 (laterales entre etapas, desde ~0,40 s) — el nivel junto a las paredes baja porque el
  aceite queda girando con los engranajes;
* salpicadura intermitente: S4 (superior etapa alta: seco → mixto desde 0,40 s), V5, V9;
* siempre secos (aire): S8, V2, V4, V6, V8, V11, V12.
Ninguno quedó «no disponible»: los volúmenes sensibles de V1–V12 no fueron cubiertos por un
engranaje en esta corrida.

## 10. Escenarios del ciclo real con sensores en los tapones (06–09)

Escenarios calibrados con datos reales de operación (`docs/ESCENARIOS.md`): rpm medida
(48,43 rpm/km/h), aceite al nivel de los tapones de la tapa de maza (94,5 L en el modelo), aire a
0,844 kg/m³ y 28 sensores capacitivos (S1–S8, V1–V12 y P1–P8 en los tapones, que giran con la maza).
Malla ib_test; números de `tools/compare_real.py` → `results/validation/real_cycle_summary.csv`
y `real_cycle.png` (promedios entre 0,15 s y el final).

| | 06 reposo | 07 acarreo | 08 subida cargado | 09 bajada retardo |
|---|---|---|---|---|
| entrada | 0 km/h, plano | 14,07 km/h, +1,6 % | 10,27 km/h, +7,8 % | 20,47 km/h, −6 %, −0,5 m/s² |
| motor | 0 rpm | 681 rpm | 497 rpm | 991 → 965 rpm |
| t simulado / pasos | 0,15 s / 87 | 0,35 s / 1 782 | 0,35 s / 1 187 | 0,35 s / 2 777 |
| cálculo | 2 min (1 núcleo) | 52 min (1) | 33 min (1) | 56 min (1, luego 2 desde 0,125 s) |
| Co máx / p95 | 0,24 / 0,24 | 2,31 / 1,03 | 1,42 / 1,05 | 1,77 / 1,02 |
| deriva máx. del aceite | 4 × 10⁻⁶ % | 9 × 10⁻⁵ % | 4 × 10⁻⁵ % | 2,6 × 10⁻⁴ % |
| velocidad media del aceite | 0,007 m/s | 0,60 m/s | 0,46 m/s | 0,82 m/s |
| aceite sobre el eje | 0 % | 5,4 % | 3,6 % | 7,1 % |
| fracción mojada engranajes (alta / baja) | 0,18 (0,18 / 0,19) | 0,16 (0,09 / 0,22) | 0,17 (0,12 / 0,22) | 0,16 (0,08 / 0,23) |
| inclinación del pozo de aceite (+ adelante) | +2,0° | −5,8° | −8,3° | −0,9° |
| nivel inferido por el anillo de tapones | −0,230 m | −0,234 m | −0,233 m | −0,238 m |

Lectura:
* **Velocidad**: el aceite levantado sobre el eje y su velocidad crecen con la velocidad del camión
  (08 < 07 < 09); la etapa alta se moja menos al girar más rápido (el aceite es expulsado hacia la
  corona), la baja casi no cambia.
* **Pendiente y retardo**: el arrastre de los engranajes inclina el pozo hacia atrás (07: −7,8°
  respecto del reposo); en 08 la subida suma −4,5° (−10,4° en total); en 09 la bajada (+3,4°) y la
  desaceleración (+2,9° cuasi-estático) lo compensan casi por completo.
* **Inventario**: conservado (deriva ≤ 2,6 × 10⁻⁴ %), también al pasar 09 de 1 a 2 procesos.
* **Sensores en los tapones**: en reposo P6 y P7 (240°/300°) quedan justo en la superficie libre
  («mixto», fracción de aceite 0,61–0,62), que es la condición de nivel correcto de la guía; los
  demás tapones están secos. Con la maza girando, cada tapón barre el círculo: P7 baja hacia el pozo
  (mojado 54–99 % del tiempo, el resto mixto), P6 sube y se seca (mojado/mixto 33 % en 07, 54 % en 08, 0 % en 09, donde
  la maza gira más rápido) y P8 entra al pozo al final de 09. Un tapón solo indica si está bajo el
  nivel en ese instante; el nivel se recupera del arco mojado que ve en una vuelta de la rueda.
* **Nivel desde los tapones**: el arco mojado del círculo de tapones (57–63°) da el nivel con
  ±4 mm (resolución del anillo de 3°) en los cuatro escenarios: el aceite del bolsillo de la tapa
  (24 L) casi no es perturbado a estas velocidades de maza (15–30 rpm). Una vuelta de rueda tarda
  2–4 s, más que lo simulado (0,35 s); el anillo fijo de 120 posiciones da lo que vería un tapón
  en cualquier ángulo sin simular la vuelta completa.
* **No disponible**: ningún sensor quedó dentro de un sólido.

Videos con todos los sensores y un panel en la tapa de maza: `results/0[6-9]_real_*/animation_sensores.mp4`;
3D del acarreo típico: `results/07_real_acarreo/animation_3d.mp4`.

## 11. Limitaciones

1. **Resolución del engrane**: 1,5–2 celdas por diente (ib_test/ib_mid). El arrastre global del aceite
   está convergido (§5), el flujo dentro del engrane no. Velocidades locales de hasta 30–150 m/s en el
   aire del engrane son artefactos de resolución; limitan Δt.
2. **Aceite en celdas de cuerpos** (§3): 1,7–4,3 L de 95 L según malla y duración.
3. **Duración**: 0,25–0,40 s de tiempo físico (≈ 3,3 vueltas del rotor; 0,67 vueltas del porta alto;
   0,1 del cubo). Suficiente para el arrastre y la salpicadura de la etapa alta; la etapa baja y el
   bolsillo de la tapa no alcanzan régimen periódico.
4. **Geometría aproximada** (modelo reconstruido, no CAD OEM) y sin pasajes de lubricación;
   rodamientos y sellos simplificados (`docs/GEOMETRIA.md`).
5. **Laminar**, sin modelo de turbulencia ni de cavitación/aireación; temperatura uniforme (60 °C):
   los datos reales no traen temperatura del aceite del motor-rueda.
6. **Corridas en serie** de 1 núcleo: el paralelo funciona (`tools/run_case.sh <caso> N`), pero el
   desplazamiento de aceite no cruza fronteras de procesador (sigue siendo conservativo).
7. El costo lo domina la ecuación de presión (acelerar el cálculo de χ redujo el tiempo por paso
   solo ~6 %). A velocidad máxima del camión (2 840 rpm) el Δt baja ~3×.
8. **Datos reales**: el camión registrado lleva un 5GDY106B y las guías describen el GDY108C; la
   geometría del modelo no es la de ninguno de los dos en detalle. El círculo de tapones (r = 0,270 m)
   se midió en fotografías, con incertidumbre de algunos milímetros a centímetros.

# Escenarios

Los YAML de `scenarios/` se generan con `python3 tools/make_scenarios.py` (no editarlos a mano) a
partir de una sola entrada física por escenario, para que velocidad del motor, aceleración del
camión e inclinación sean coherentes entre sí.

Supuestos comunes: rueda derecha; neumático clase 53/80R63 con radio cargado 1,95 m; razón total
32,6087; aceite 95 L a 60 °C (ρ = 900 kg/m³, ν = 264 cSt); aire 1,06 kg/m³; σ = 0,030 N/m;
superficie libre inicial plana y normal a la gravedad efectiva en t = 0, con el volumen exacto
calculado sobre la malla. Los primeros 0,05 s llevan el tren de engranajes desde el reposo hasta la
velocidad inicial: es un **arranque numérico** del campo de flujo, no una aceleración del camión, y
no se usa en las estadísticas.

| Escenario | Entrada | Motor | Cargas | Duración |
|---|---|---|---|---|
| `01_reposo` | mecanismos detenidos | 0 rpm | g | 0,25 s |
| `02_operacion_normal` | camión a 20 km/h constante | 887 rpm (rueda 27,2 rpm) | g | 0,25 s |
| `03_operacion_inclinada` | 20 km/h subiendo pendiente 10 % (5,71°) con peralte 4° | 887 rpm | g girada | 0,25 s |
| `04_transitorio` | 20 km/h; +1,5 m/s² entre 0,08–0,15 s; −3,0 m/s² entre 0,15–0,40 s | 887 → 904 → 784 rpm (de v(t)) | g − a(t) | 0,40 s |
| `05_bajo_inventario` | como 02 con 55 L | 887 rpm | g | 0,25 s |

## Escenarios del ciclo real (06–09)

Calibrados con datos reales de operación de la mina (carpeta `hack` del Hackathon; ver «Uso de los
datos reales» más abajo):

| Escenario | Entrada (dato real) | Motor | Cargas | Duración |
|---|---|---|---|---|
| `06_real_reposo_nivel` | detenido en plano; aceite justo en los tapones de nivel | 0 rpm | g | 0,15 s |
| `07_real_acarreo` | 14,07 km/h (mediana en movimiento), pendiente +1,6 % (mediana) | 681 rpm | g girada 0,92° | 0,35 s |
| `08_real_subida_cargado` | cargado (> 100 t) en subida > 5 %: 10,27 km/h y +7,8 % (medianas) | 497 rpm | g girada 4,46° | 0,35 s |
| `09_real_bajada_retardo` | bajada con retardo eléctrico: 20,47 km/h, −6,0 % (medianas); −0,5 m/s² (cuartil 25 %) desde 0,05 s | 991 → 965 rpm (de v(t)) | g girada −3,43° y − a(t) | 0,35 s |

Cambios respecto de 01–05 (que se conservan tal cual, con sus resultados):

1. **Velocidad → rpm medida**: 48,43 rpm de motor por km/h (mediana; p10 45,7, p90 50,4), en lugar
   del radio supuesto de 1,95 m. Con la razón 32,6 del modelo equivale a un radio de rodadura de 1,79 m.
2. **Nivel de aceite por los tapones**: la guía de Komatsu controla el nivel con el camión detenido
   en plano y dos tapones magnéticos de la tapa de maza alineados horizontalmente; el aceite debe
   llegar a ellos (máximo 90 L, GDY108C). De las fotografías de la guía: 8 tapones en un círculo de
   r ≈ 0,270 m (0,479 × círculo de pernos de la tapa, ≈ 0,565 m) a 30, 60, 120, 150, 210, 240, 300 y
   330°; el par de nivel está a 240° y 300° → nivel y = −0,234 m. El inventario se calcula sobre la
   malla como el fluido bajo ese nivel: **94,5 L** en el modelo (5 % sobre los 90 L de la guía; la
   cavidad del modelo es una simplificación de la real). Con pendiente, ese volumen se pone con la
   superficie libre normal a la gravedad efectiva.
3. **Aire a la altitud de la faena**: 81 kPa medidos (11,7 psi) → 0,844 kg/m³ a 60 °C.
4. **Sensores en los puertos reales** (`rotating_sensors`): un sensor capacitivo en cada uno de los 8
   tapones P1–P8, 15 mm dentro de la cara interior de la tapa (z = −1,7195 m), electrodo a lo largo
   del eje de la rueda; **giran con la maza** (ver `docs/SENSORES.md`). Se mantienen S1–S8 y V1–V12.

El aceite sigue a 60 °C (264 cSt): los datos no traen temperatura de aceite del motor-rueda (la del
estator tiene mediana 64 °C, p95 87 °C) y del aceite solo se conoce la viscosidad a 60 °C.

### Uso de los datos reales

Los archivos de la carpeta `hack` (registros del camión 58, 930E-5 con motor-rueda 5GDY106B, y las
guías de Komatsu del GDY108C) contienen información operacional (coordenadas, números de serie,
fallas). **No se suben al repositorio**: solo se usaron estadísticas agregadas (medianas y
percentiles de 63 720 tramas de 531 paquetes de datos), que son las que aparecen aquí. Las guías
describen el GDY108C (razón 35,03:1); el camión registrado lleva un 5GDY106B, por eso la relación
rpm/velocidad se toma medida y no de la razón nominal.

Comparación «mecanismos detenidos vs en movimiento»: `01_reposo` frente a `02_operacion_normal`
(mismo aceite, misma malla).

Por qué 20 km/h: es una velocidad típica de un camión cargado en rampa de mina y mantiene el costo
del caso nominal en ~1 h por escenario en la malla de prueba (2 núcleos compartidos). A velocidad máxima (~64 km/h,
≈ 2 840 rpm) el paso de tiempo baja ~3 veces; el flujo de trabajo es el mismo
(`motor_rpm` en el YAML).

Sensores virtuales (`scenarios/_sensors.yaml`), todos a 5–10 mm de paredes **fijas** (donde se
podría montar un sensor real), coordenadas del modelo (z = eje, y = arriba):

| id | ubicación | posición [m] |
|---|---|---|
| S1 | fondo etapa alta, cara de placa spindle | (0, −0,400, −1,025) |
| S2 | fondo entre etapas | (0, −0,400, −1,300) |
| S3 | fondo bolsillo exterior (tapa) | (0, −0,450, −1,680) |
| S4 | superior etapa alta (aire esperado) | (0, 0,400, −1,025) |
| S5 | lateral delantero entre etapas | (0,400, −0,250, −1,300) |
| S6 | lateral trasero entre etapas | (−0,400, −0,250, −1,300) |
| S7 | anillo corona alta / tubo spindle, fondo | (0, −0,500, −1,130) |
| S8 | parte superior del bolsillo exterior (aire esperado) | (0, 0,450, −1,680) |

Cada sensor registra en cada paso: fracción de aceite α, presencia (α > 0,5), velocidad, p, p_rgh y χ.
El estado del sistema se arma solo con lo que vería una instalación real (sensores de fondo
S1, S2, S3, S7): NORMAL (≥ 3 mojados), ALERTA (1–2), CRÍTICO (0). Es una regla demostrativa para
mostrar el flujo de datos hacia OA²IS, no un criterio de mantenimiento calibrado.

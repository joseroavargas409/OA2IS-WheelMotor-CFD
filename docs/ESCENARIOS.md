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
| `01_reposo` | mecanismos detenidos | 0 rpm | g | 0,3 s |
| `02_operacion_normal` | camión a 20 km/h constante | 887 rpm (rueda 27,2 rpm) | g | 0,4 s |
| `03_operacion_inclinada` | 20 km/h subiendo pendiente 10 % (5,71°) con peralte 4° | 887 rpm | g girada | 0,4 s |
| `04_transitorio` | 20 km/h; +1,5 m/s² entre 0,10–0,25 s; −3,0 m/s² entre 0,25–0,60 s | 887 → 923 → 755 rpm (de v(t)) | g − a(t) | 0,6 s |
| `05_bajo_inventario` | como 02 con 55 L | 887 rpm | g | 0,4 s |

Comparación «mecanismos detenidos vs en movimiento»: `01_reposo` frente a `02_operacion_normal`
(mismo aceite, misma malla).

Por qué 20 km/h: es una velocidad típica de un camión cargado en rampa de mina y mantiene el costo
del caso nominal en ~30 min por escenario en la malla de prueba. A velocidad máxima (~64 km/h,
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

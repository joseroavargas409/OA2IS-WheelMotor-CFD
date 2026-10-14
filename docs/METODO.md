# Método numérico

Resumen: **VOF bifásico aceite/aire (OpenFOAM v1912, MULES/PIMPLE) + frontera inmersa con
penalización de Brinkman para los engranajes y porta-satélites, que se mueven con su cinemática
exacta (giro propio + órbita) en cada paso de tiempo**, sobre una malla fija 2,5D (capas axiales)
ajustada a las paredes fijas y a las paredes de revolución que giran. Las cargas de inclinación y
aceleración del camión entran como gravedad efectiva g − a(t).

## 1. Ecuaciones

`solvers/interIBFoamOA2IS` es `interFoam` (v1912) con dos agregados:

* **Campo de presencia de sólido** χ ∈ [0, 1] y velocidad de sólido **U_s**, recalculados al inicio de
  cada paso para la pose de cada cuerpo en ese instante (`immersedBodies.H`):
  * la pose sale de integrar la tabla de velocidades angulares (órbita ω_o y giro propio ω_s) del
    escenario; para un satélite: centro c(t) = R(θ_o) c₀, ángulo propio θ_s;
  * χ de cada celda = fracción de 2×2×2 sub-puntos dentro de algún cuerpo (pruebas exactas
    punto-en-polígono sobre el perfil de evolvente regenerado, en coordenadas del cuerpo);
  * U_s = ω_o ẑ × p + ω_s ẑ × (p − c(t)) (cuerpo rígido; div U_s = 0).
* **Término de penalización** en la ecuación de momento (implícito):

  ∂(ρU)/∂t + ∇·(ρφU) − ∇·(μ∇U) = −∇p_rgh − g·h∇ρ + σκ∇α + **ρ χ (K/Δt)(U_s − U)**, con K = 1000.

  En celdas sólidas el fluido queda forzado a la velocidad del cuerpo en un paso; fuera (χ = 0) la
  ecuación es la de `interFoam`. Es un acoplamiento **mecánico real**: el aceite solo se mueve porque
  los cuerpos lo empujan y lo arrastran a través de la presión y la viscosidad (no hay fuerzas
  artificiales sobre el aceite, ni MRF, ni rotación "visual").
* **Desplazamiento conservativo del aceite**: en cada paso, el aceite del volumen que un cuerpo acaba
  de cubrir (Δχ > 0) y cualquier aceite que exceda la parte libre (1 − χ) de una celda se reparte entre
  las celdas vecinas con capacidad libre (1 − χ − α), en hasta `displacePasses` pasadas. Así α queda
  acotado, el volumen total de aceite se conserva exactamente y los cuerpos no "absorben" aceite.
  Lo que la advección VOF de un paso alcanza a meter en celdas parcialmente cubiertas se devuelve en
  el paso siguiente (≈ 1 L de los 95 L en la malla `ib_test`, ver `oil_excess_in_bodies_L`).
  El inventario se verifica cada paso (`postProcessing/oilVolume`).

Paredes ajustadas a malla: spindle y cuerpos de las coronas (fijas, `noSlip`); cubo, tapa, anillo del
porta bajo, eje del rotor y manguito intermedio (superficies de revolución que giran con su
velocidad, `rotatingWallVelocity` con tabla ω(t)). Presión: `fixedFluxPressure`. α: gradiente nulo.

## 2. Fluidos

| | ρ [kg/m³] | ν [m²/s] | fuente |
|---|---|---|---|
| Aceite (60 °C) | 900 | 2,642 × 10⁻⁴ (264 cSt) | `komatsu-oil-sensors-v4-1/dashboard/constants.js` |
| Aire | 1,06 | 1,9 × 10⁻⁵ | aire a 60 °C |
| σ | 0,030 N/m | | ídem V4.1 |

`oa2is_cfd/scenario.py` permite en cambio dar dos puntos de referencia de viscosidad (ecuación de
Walther, ASTM D341) y la temperatura.

## 3. Esquemas y control

Euler implícito; `div(rhoPhi,U)` linearUpwind; α con vanLeer + compresión (cAlpha = 1, MULES
semi-implícito, 2 correcciones); laplaciano/snGrad `limited 0.5` con 1 corrección no ortogonal;
PISO con 3 correctores; `p_rgh` GAMG (tol 10⁻⁷). Paso adaptativo con Co_max = 1 (producción) o 0,5
(sensibilidad), Co de interfase ≤ 0,5, Δt ≤ 2 ms.

## 4. Inclinación y aceleraciones (sin doble conteo)

La malla está fija en el sistema de la carcasa (que se traslada con el camión y no gira). Por lo tanto:

* pendiente (cabeceo) y peralte (alabeo) giran el vector gravedad: g₀ = R(pitch, roll)·(0, −9,81, 0);
* la aceleración longitudinal/lateral del camión a(t) entra **una sola vez**, como fuerza aparente
  −a(t) mediante `tabulatedAccelerationSource` (recalcula g·h con g_eff = g₀ − a(t));
* la rotación de los engranajes **no** se modela como sistema rotatorio (eso sí sería doble conteo
  con la frontera inmersa): los cuerpos se mueven en el sistema inercial de la carcasa.
* En el escenario transitorio la velocidad del camión v(t) es la única entrada: de ella salen la
  aceleración a = dv/dt y las rpm del motor = v/r_neumático · 32,6087 (`tools/make_scenarios.py`).

## 5. Malla

`oa2is_cfd/layered.py` + `oa2is_cfd/ib.py`: malla 2D (gmsh, triángulos/cuadriláteros) del plano de la
sección que **sigue todos los contornos fijos** (spindle, coronas, cubo, ejes) y se refina en la zona
de engranajes; luego se extruye en capas axiales cuyos cortes coinciden con todas las cotas z de las
piezas (prismas exactos). Las bolsas de fluido completamente encerradas se eliminan (si no, la
ecuación de presión es singular). Se escribe directamente en formato `polyMesh` y pasa `checkMesh`
(«Mesh OK»).

| Preset | celdas | tamaño en zona móvil | uso |
|---|---|---|---|
| `ib_test` | 145 934 | 18 mm, capas 25 mm | 5 escenarios + sensibilidad Δt |
| `ib_mid` | 477 386 | 12 mm, capas 18 mm | sensibilidad de malla (nominal) |
| `ib_coarse` | 824 852 | 9 mm, capas 13 mm | generada y probada; corrida descartada por costo (~7 h) |
| `ib_medium` / `ib_fine` | ~3 M / ~8 M (estimado) | 6 / 4 mm | definidos, no corridos (costo) |

## 6. Por qué frontera inmersa (alternativas investigadas)

Se probaron con ejecuciones reales:

1. **Malla dinámica overset** (`overInterDyMFoam`, una malla por engranaje/porta, `multiSolidBodyMotionSolver`
   con órbita+giro verificado). Se construyó y corrió (510 671 celdas, Mesh OK), pero: el corte de
   huecos en el contacto entre dientes (juego de 1,4–1,9 mm) colapsaba en cascada (se parcheó
   `inverseDistance`), la interpolación entre mallas no conserva masa en VOF, y el costo era alto.
   Además, parte de las velocidades espurias que se atribuían al overset venían del error de
   compilación descrito abajo. Queda en el repositorio (`oa2is_cfd/meshing.py`, `case.py`,
   `solvers/overInterDyMFoamOA2IS`) como alternativa documentada.
2. **Malla deslizante (AMI)**: no aplicable de forma general: los satélites giran y orbitan a la vez y
   sus volúmenes se cruzan con los de la solar y la corona.
3. **ALE / deformación de malla**: imposible con dientes que engranan (la malla se invierte).
4. **MRF estacionario**: descartado por el propio enunciado (y no representa la superficie libre).

La frontera inmersa permite mover todos los cuerpos con su cinemática exacta sobre una malla fija,
conserva el aceite y su costo es el de `interFoam`. Su limitación es la resolución: el diente se
representa con 1–2 celdas en las mallas corridas, por lo que el flujo dentro del engrane (bombeo y
atrapamiento entre dientes) está sub-resuelto; ver `docs/VALIDACION.md`.

## 7. Error de compilación de OpenFOAM v1912 en Ubuntu 24.04 (encontrado y corregido)

Síntoma: `interFoam` (sin modificar, paquete `openfoam 1912.200626-2build3`) **explota en 3–5 pasos
en cualquier malla no ortogonal** con superficie libre, incluso una caja de hexaedros con nodos
perturbados (no ortogonalidad 30°), con gravedad o solo tensión superficial, con razón de
densidades 1,1, con cualquier esquema. En una caja ortogonal funciona bien.

Causa (aislada imprimiendo los términos de `pEqn.H` celda a celda): `fvc::reconstruct()` calcula
`inv(T) & S` reutilizando el almacenamiento temporal del operando `S` para el resultado. OpenFOAM
declara los punteros de lista con `__restrict__`; con GCC 13 y -O3 el compilador escribe la
componente x del resultado antes de leer la componente x de `S` para calcular y y z. Resultado
medido en una celda: reconstrucción correcta (−1208, 625, 0,006) frente a la del binario
(−1208, **2,24 × 10⁶**, 0,006). El error solo existe si el tensor tiene términos fuera de la
diagonal (celdas no ortogonales), por eso las cajas ortogonales funcionaban.

Corrección (`patches/openfoam-v1912/oa2is-gcc13-restrict.patch`): quitar `__restrict__` de
`UList`, `ListLoopM.H` y `Matrix`, y proteger `inv()` contra alias. `tools/build_patched_openfoam.sh`
recompila `libOpenFOAM` y los solvers OA²IS en `/opt/oa2is-openfoam-patches` **sin tocar la
instalación del sistema**; los solvers se ejecutan con `LD_LIBRARY_PATH` apuntando ahí.
Verificación: la caja perturbada pasa de explotar en el paso 4 a quedar estable (Co ≈ 4 × 10⁻⁴), y el
caso real corre estable.

Riesgo residual documentado: código de plantillas instanciado **dentro** de `libfiniteVolume` del
sistema conserva el `__restrict__`. Las operaciones críticas del solver (reconstrucción, inversa,
transpuesta, producto tensor·vector) se compilan en nuestro solver o en `libOpenFOAM` reconstruida;
`FULL=1 tools/build_patched_openfoam.sh` recompila también `libfiniteVolume`.

# OA2IS-WheelMotor-CFD

Simulación CFD 3D transitoria y bifásica (aceite/aire, VOF) del sistema de lubricación de la
reducción planetaria de un motor-rueda Komatsu 930E/980E, con los engranajes, satélites y
porta-satélites **moviendo físicamente el aceite** (giro propio + órbita, cinemática exacta), hecha
solo con herramientas de código abierto (OpenFOAM, Python, gmsh, VTK/PyVista, FFmpeg).
Proyecto OA²IS — Hack the Challenge 2026 (NTT DATA × Komatsu).

![geometría](docs/img/geometria_cfd.png)

**Estado: funcional.** Se ejecutaron 5 escenarios + estudios de sensibilidad de malla y de paso de
tiempo con el solver real; resultados, sensores y videos en `results/`. Lo que el modelo resuelve
bien y lo que no está en [`docs/VALIDACION.md`](docs/VALIDACION.md).

## Qué hay aquí

| Carpeta | Contenido |
|---|---|
| `third_party/` | Copias (sin modificar) de los archivos necesarios de los repositorios originales, con procedencia y licencias (`PROVENANCE.md`) |
| `geometry/` | Exportador del modelo 3D original a STL (`export_parts.mjs`), inventario de piezas y cinemática, geometría CFD procesada (`cfd/*.stl`) |
| `oa2is_cfd/` | Biblioteca Python: engranajes (`gears.py`), cinemática (`kinematics.py`), geometría del dominio y reparaciones (`model.py`), mallador 2,5D (`layered.py`), casos de frontera inmersa (`ib.py`), escenarios (`scenario.py`), post-proceso (`ibpost.py`), render (`render.py`); y la vía overset investigada (`meshing.py`, `case.py`, `post.py`) |
| `solvers/` | `interIBFoamOA2IS` (producción: VOF + cuerpos inmersos en movimiento) y `overInterDyMFoamOA2IS` (alternativa overset) |
| `patches/` | Parches a OpenFOAM v1912 (functionObjects con GCC 13, overset, **error de compilación GCC 13 `__restrict__`**) |
| `scenarios/` | 5 escenarios YAML (generados por `tools/make_scenarios.py`) y sensores virtuales |
| `tools/` | Compilar OpenFOAM parcheado, construir casos, correr, post-procesar, comparar, exportar geometría |
| `tests/` | Engranajes sin interferencia, cinemática en OpenFOAM, posición de los cuerpos dentro de la solución CFD |
| `cases/` | Diccionarios exactos de un caso corrido (referencia) |
| `results/` | Por caso: CSV/Parquet de sensores e inventario, historia numérica, métricas de campo, gráficos, MP4; `validation/` |
| `docs/` | Método, geometría, escenarios, validación, dependencias y licencias |

## Reproducir

Probado en Ubuntu 24.04 (2 núcleos, 8 GB RAM).

```bash
# 1. dependencias (todas libres)
sudo apt-get install openfoam openmpi-bin ffmpeg xvfb
pip install --break-system-packages numpy scipy shapely gmsh pyvista trimesh pandas pyarrow matplotlib pyyaml

# 2. OpenFOAM v1912 parcheado + solvers OA2IS (no modifica la instalación del sistema;
#    instala en /opt/oa2is-openfoam-patches, ~1 h con 2 núcleos)
tools/build_patched_openfoam.sh

# 3. (opcional) geometría: exportar el modelo 3D original y la geometría CFD procesada
PATH=/opt/node22/bin:$PATH node --import ./geometry/register.mjs geometry/export_parts.mjs geometry/raw
python3 tools/export_cfd_geometry.py
python3 tests/test_gears.py                      # engranajes sin interferencia en un ciclo de engrane

# 4. escenarios (ya generados en scenarios/)
python3 tools/make_scenarios.py

# 5. malla + caso (malla ib_test: 146 k celdas, ~1 min; ib_mid: 477 k)
python3 tools/build_ib_case.py scenarios/02_operacion_normal.yaml runs/c02 --mesh runs/mesh_ib_test --res ib_test

# 6. correr (1 proceso; N procesos con descomposición jerárquica: último argumento)
SOLVER=/opt/oa2is-openfoam-patches/bin/interIBFoamOA2IS tools/run_case.sh runs/c02 1
#    o varios casos en cola con post-proceso automático:
tools/run_queue.sh 1 runs/c01 runs/c02 ...

# 7. post-proceso: sensores CSV/Parquet, inventario, gráficos, fotogramas y MP4
python3 tools/postprocess.py runs/c02 --mp4
python3 tests/test_ib_motion.py runs/c02         # cuerpos en la pose analítica dentro de la solución
python3 tools/compare.py                          # figuras y números de validación
```

Las variables de entorno de OpenFOAM que necesita el paquete Debian están en
`tools/run_case.sh` (`FOAM_ETC`, `WM_PROJECT_DIR`, ... y `LD_LIBRARY_PATH=/opt/oa2is-openfoam-patches/lib`).

## Resultados

Ver [`docs/VALIDACION.md`](docs/VALIDACION.md) (tabla de corridas, conservación, sensibilidad,
detenido vs en movimiento, respuesta a inclinación y frenado) y los videos `results/*/animation.mp4`.

## Repositorios originales

`wheel-motor-inspector`, `komatsu-oil-sensors-v4` y `komatsu-oil-sensors-v4-1` se usaron solo como
referencia de lectura y **no fueron modificados** (verificación de HEAD, `git status` y SHA-256 de
todos los archivos antes y después: `docs/VALIDACION.md` §8). Este repositorio no depende de ellos
en tiempo de ejecución.

## Licencia

GPL-3.0-or-later (derivado de solvers de OpenFOAM). `third_party/` conserva sus licencias (MIT).
Ver `docs/DEPENDENCIAS_Y_LICENCIAS.md`.

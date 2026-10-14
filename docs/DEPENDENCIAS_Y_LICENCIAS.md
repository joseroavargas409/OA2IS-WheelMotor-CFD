# Dependencias y licencias

Todo el flujo usa **solo software libre / de código abierto**. No se usó ANSYS, Fluent, STAR-CCM+,
COMSOL, SolidWorks, SimScale ni ninguna API propietaria.

## Software usado (versiones con que se ejecutó)

| Herramienta | Versión | Licencia | Para qué |
|---|---|---|---|
| Ubuntu | 24.04.5 LTS | varias (libres) | sistema base |
| OpenFOAM (ESI) | v1912, paquete Debian `openfoam 1912.200626-2build3` | GPL-3.0-or-later | solver VOF, utilidades (checkMesh, decomposePar, postProcess) |
| GCC | 13.3.0 | GPL-3.0 con excepción de runtime | compilar bibliotecas parcheadas y solvers OA²IS |
| Open MPI | 4.1.6 | BSD-3-Clause | ejecución paralela |
| Python | 3.13 | PSF | todo el pipeline |
| gmsh (PyPI) | 4.15.2 | GPL-2.0-or-later (con excepción para enlazar) | malla 2D de la sección |
| NumPy / SciPy | 2.5 / 1.18 | BSD-3-Clause | cálculo |
| Shapely | 2.2 | BSD-3-Clause | geometría 2D, perfiles, recortes |
| trimesh | 5.1 | MIT | exportar STL de la geometría procesada |
| PyVista / VTK | 0.49 / 9.7 | MIT / BSD-3-Clause | lectura de resultados OpenFOAM, cortes, render 3D |
| pandas / PyArrow | 3.0 / 26.0 | BSD-3-Clause / Apache-2.0 | tablas de sensores, CSV/Parquet |
| Matplotlib | 3.11 | PSF-compatible (Matplotlib License) | gráficos y fotogramas |
| PyYAML | 6.0 | MIT | escenarios |
| FFmpeg (libx264) | 6.1.1 | LGPL-2.1+/GPL-2+ (x264 GPL) | video MP4 |
| Xvfb | 21.1 | MIT/X11 | render 3D sin pantalla |
| Node.js | 22.22 | MIT | ejecutar el generador del modelo 3D para exportar STL |
| three.js (copiado) | la del repositorio original | MIT | generador del modelo 3D |

Instalación de las dependencias Python: `pip install --break-system-packages numpy scipy shapely gmsh
pyvista trimesh pandas pyarrow matplotlib pyyaml` (más `apt-get install openfoam openmpi-bin ffmpeg xvfb`).

## Licencia de este repositorio

**GPL-3.0-or-later** (archivo `LICENSE`). Motivo: `solvers/interIBFoamOA2IS` y
`solvers/overInterDyMFoamOA2IS` derivan de solvers de OpenFOAM (GPL-3.0-or-later) y los parches de
`patches/openfoam-v1912/` modifican OpenFOAM; usar la misma licencia en todo el repositorio evita
incompatibilidades. Los componentes de `third_party/` conservan su licencia original (MIT para
three.js y para `komatsu-oil-sensors-v4-1`, ambas compatibles con GPL-3.0); ver
`third_party/PROVENANCE.md`.

Nota sobre `wheel-motor-inspector`: el repositorio original no trae archivo LICENSE. Sus archivos
`dist/model.js`, `assembly-v2.js` y `assembly-v3.js` se reutilizan aquí porque el autor del
proyecto OA²IS (dueño de ese repositorio) lo pidió expresamente para este desarrollo. Antes de una
publicación abierta conviene que el dueño agregue una licencia explícita a ese repositorio. No se
copió ninguna imagen, figura de patente ni documento de terceros.

## Componentes que NO se redistribuyen

* La instalación de OpenFOAM del sistema (se instala con `apt`); el árbol de compilación con los
  parches queda en `/opt/oa2is-openfoam-build` y se regenera con `tools/build_patched_openfoam.sh`.
* Las mallas y los resultados de campo completos (cientos de MB): se regeneran con los scripts; el
  repositorio guarda los resultados procesados (CSV/Parquet, gráficos, MP4, resúmenes).

# Casos

Los casos de OpenFOAM se **generan** (malla + diccionarios + condición inicial) con
`tools/build_ib_case.py`; no se versionan completos porque la malla y los campos pesan cientos de MB.

`example_02_operacion_normal_ib_test/` guarda, como referencia, los diccionarios exactos del caso
nominal que se corrió (malla `ib_test`):

* `system/controlDict`, `fvSchemes`, `fvSolution` — control, esquemas y solvers lineales;
* `constant/immersedBodiesDict` — los 12 cuerpos inmersos (solares, satélites, porta-satélites,
  bandas dentadas de las coronas) con sus perfiles y tablas ω(t) de órbita y giro;
* `constant/fvOptions` + `acceleration-terms.dat` — aceleración del vehículo (cero en este caso);
* `constant/transportProperties`, `g`, `boundary` (lista de parches de la malla);
* `0/U`, `0/p_rgh` — condiciones de borde (paredes giratorias con ω(t), `fixedFluxPressure`);
* `scenario.json`, `initial_oil.json` (nivel y volumen inicial), `mesh_stats.json`.

Regenerar el caso completo:

```
python3 tools/build_ib_case.py scenarios/02_operacion_normal.yaml cases/c02 --mesh cases/mesh_ib_test --res ib_test
SOLVER=/opt/oa2is-openfoam-patches/bin/interIBFoamOA2IS tools/run_case.sh cases/c02 1
```

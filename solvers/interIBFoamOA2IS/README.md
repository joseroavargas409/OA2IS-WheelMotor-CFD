# interIBFoamOA2IS

Solver de producción de OA2IS-WheelMotor-CFD. Es una copia de
`applications/solvers/multiphase/interFoam` de OpenFOAM v1912 (GPL-3.0-or-later, © OpenCFD Ltd.)
con un agregado de frontera inmersa para cuerpos rígidos en movimiento prescrito
(engranajes, satélites, porta-satélites):

* `immersedBodies.H` — lee `constant/immersedBodiesDict`; en cada paso calcula la pose de cada
  cuerpo (integral exacta de las tablas ω(t) de órbita y giro propio), la fracción de sólido χ por
  celda (2×2×2 sub-puntos, perfiles exactos de evolvente en tabla polar), la velocidad de cuerpo
  rígido U_s y desplaza de forma conservativa el aceite de las celdas que un cuerpo acaba de cubrir.
* `createImmersedBodies.H` — campos `chi` y `Us` (se escriben con los resultados).
* `interIBFoamOA2IS.C` — llama `ib.update(t, chi, Us, alpha1)` después de avanzar el tiempo.
* `UEqn.H` — término implícito de penalización de Brinkman `ρ χ K/Δt (U_s − U)` (K = `penalty`).

Formato de `immersedBodiesDict` (generado por `oa2is_cfd/ib.py::bodies_dict`):

```
penalty 1000;        // K
displaceOil true;    // empujar el aceite fuera de las celdas que se cubren
subSamples 2;        // n -> n^3 sub-puntos por celda para chi
bodies
{
    planetHigh0
    {
        center (0.2656 0);                       // centro en t = 0 (plano xy)
        orbit ((0 0) (0.05 -18.58) (100 -18.58)); // omega de órbita [rad/s] vs t (porta-satélites)
        spin  ((0 0) (0.05 49.55) (100 49.55));   // omega de giro propio relativo [rad/s]
        primitives
        {
            gear { type polar; z (-1.181 -1.079); r ( ...radio vs ángulo... ); }
            core { type annulus; z (-1.195 -1.065); rIn 0; rOut 0.044; }
        }
    }
    ...
}
```

Tipos de primitiva: `polar` (sólido dentro de r(θ)), `polarOut` (sólido fuera de r(θ) hasta `rOut`,
para coronas), `annulus` (`rIn`–`rOut`, opcional centro `c`), `hole` (vacía lo que cubre).

Compilar: `tools/build_patched_openfoam.sh` (requiere el parche
`patches/openfoam-v1912/oa2is-gcc13-restrict.patch`; ver `docs/METODO.md` §7).

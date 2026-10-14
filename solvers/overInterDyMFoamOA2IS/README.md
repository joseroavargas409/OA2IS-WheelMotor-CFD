# overInterDyMFoamOA2IS

Copy of OpenFOAM v1912 `applications/solvers/multiphase/interFoam/overInterDyMFoam`
(GPL-3.0-or-later, (C) OpenCFD Ltd.) with one change: the maximum Courant and
interface-Courant numbers used for adaptive time stepping exclude overset
*interpolated* cells (CourantNo.H, alphaCourantNo.H). Those cells are not
solved; tiny interpolated cells at component boundaries otherwise dictate
dt ~1e-5 s. Build: tools/build_patched_openfoam.sh.

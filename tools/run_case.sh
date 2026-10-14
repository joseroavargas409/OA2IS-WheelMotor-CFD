#!/bin/bash
# Run an OA2IS case with the patched OpenFOAM libraries.
# usage: tools/run_case.sh <case_dir> [nprocs]
set -e
CASE=$(cd "$1" && pwd); NP=${2:-2}
export FOAM_ETC=/usr/share/openfoam/etc WM_PROJECT_DIR=/usr/share/openfoam WM_PROJECT_VERSION=v1912 WM_PROJECT=OpenFOAM
export LD_LIBRARY_PATH=/opt/oa2is-openfoam-patches/lib:$LD_LIBRARY_PATH
cd "$CASE"
if [ "$NP" -gt 1 ]; then
  cat > system/decomposeParDict <<EOD
FoamFile { version 2.0; format ascii; class dictionary; object decomposeParDict; }
numberOfSubdomains $NP;
method hierarchical;
coeffs { n ($NP 1 1); order xyz; }
EOD
  [ -d processor0 ] || decomposePar -force > log.decomposePar 2>&1
  mpirun --allow-run-as-root --oversubscribe -np $NP -x LD_LIBRARY_PATH -x FOAM_ETC -x WM_PROJECT_DIR \
         -x WM_PROJECT_VERSION -x WM_PROJECT ${SOLVER:-/opt/oa2is-openfoam-patches/bin/overInterDyMFoamOA2IS} -parallel >> log.solver 2>&1
  reconstructPar -newTimes > log.reconstructPar 2>&1 || true
else
  ${SOLVER:-/opt/oa2is-openfoam-patches/bin/overInterDyMFoamOA2IS} >> log.solver 2>&1
fi

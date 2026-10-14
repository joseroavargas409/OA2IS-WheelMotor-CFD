#!/bin/bash
# Build patched OpenFOAM v1912 libraries (libOpenFOAM, liboverset) for the
# Ubuntu 24.04 'openfoam' package (1912.200626-2build3) WITHOUT touching the
# system installation. Result: $PREFIX/lib/{libOpenFOAM.so,liboverset.so};
# load them with: export LD_LIBRARY_PATH=$PREFIX/lib:$LD_LIBRARY_PATH
#
# Patches (patches/openfoam-v1912/oa2is-overset-sha1.patch, GPL-3.0-or-later):
#  1. OSHA1stream: add overflow() so functionObjects work with GCC >= 13.
#  2. inverseDistance walkFront: donor-less cells next to holes stay calculated
#     instead of cascading holes through the mesh (gear-tooth contact).
#  3. inverseDistance: interpolated cells with an empty donor stencil -> calculated.
# patches/openfoam-v1912/oa2is-gcc13-restrict.patch (GPL-3.0-or-later):
#  4. Remove __restrict__ from UList/ListLoopM/Matrix. With GCC >= 13 at -O3 the
#     restrict promise is broken by OpenFOAM's tmp re-use (result field aliasing
#     an operand) and in-place tensor & vector, inv, T(), ... produce WRONG values
#     (e.g. fvc::reconstruct on non-orthogonal cells -> interFoam explodes in a
#     few steps on any non-orthogonal mesh). Libraries/solvers compiled here get
#     the fix; libOpenFOAM is rebuilt completely (header change).
#     Set FULL=1 to also rebuild libfiniteVolume (long) for template code that
#     was instantiated inside the Debian library.
set -euo pipefail
HERE=$(cd "$(dirname "$0")/.." && pwd)
PREFIX=${PREFIX:-/opt/oa2is-openfoam-patches}
WORK=${WORK:-/opt/oa2is-openfoam-build}
mkdir -p "$WORK" "$PREFIX/lib"
cd "$WORK"
if [ ! -d src-v1912 ]; then
  apt-get source -q --download-only openfoam=1912.200626-2build3 >/dev/null
  mkdir src-v1912
  tar xjf openfoam_1912.200626.orig.tar.bz2 --strip-components=1 -C src-v1912
  tar xJf openfoam_1912.200626-2build3.debian.tar.xz
  (cd src-v1912 && for p in $(cat ../debian/patches/series); do patch -p1 -N < ../debian/patches/$p; done)
  (cd src-v1912 && patch -p1 < "$HERE/patches/openfoam-v1912/oa2is-overset-sha1.patch")
  (cd src-v1912 && patch -p1 < "$HERE/patches/openfoam-v1912/oa2is-gcc13-restrict.patch")
fi
C=$WORK/src-v1912
export WM_PROJECT_DIR=$C WM_PROJECT=OpenFOAM WM_PROJECT_VERSION=v1912 FOAM_SRC=$C/src FOAM_ETC=$C/etc
export WM_ARCH=linux64 WM_LABEL_SIZE=32 WM_LABEL_OPTION=Int32 WM_COMPILER=Gcc WM_COMPILER_LIB_ARCH=64
export WM_COMPILE_OPTION=Opt WM_DIR=$C/wmake WM_OPTIONS=linux64Gcc WM_PRECISION_OPTION=DP
export WM_OSTYPE=POSIX WM_MPLIB=SYSTEMOPENMPI FOAM_MPI=openmpi-system MPI_ARCH_PATH=/usr/lib/x86_64-linux-gnu/openmpi
export FOAM_LIBBIN=$PREFIX/lib FOAM_EXT_LIBBIN=$PREFIX/lib FOAM_USER_LIBBIN=$PREFIX/lib
export FOAM_APPBIN=$WORK/bin FOAM_USER_APPBIN=$WORK/bin WM_NCOMPPROCS=${NPROC:-2} WM_COMPILER_TYPE=system
export PATH=$C/wmake:$PATH
(cd $C/wmake/src && make -s)
# lnInclude directories for every library source tree (headers only)
for mk in $(find $C/src -name Make -type d); do
  wmakeLnInclude -u "$(dirname "$mk")" >/dev/null 2>&1 || true
done
wmakeLnInclude -u $C/src/OpenFOAM
(cd $C/src/OSspecific/POSIX && ./Allwmake)
(cd $C/src/OpenFOAM && wmake -j ${NPROC:-2} libso)
[ "${FULL:-0}" = 1 ] && (cd $C/src/finiteVolume && wmake -j ${NPROC:-2} libso)
# liboverset links against the system libraries in /usr/lib
(cd $C/src/overset && wmake -j ${NPROC:-2} libso)
# OA2IS solver variant (Courant numbers exclude interpolated overset cells)
export FOAM_SOLVERS=$C/applications/solvers FOAM_USER_APPBIN=$PREFIX/bin FOAM_APPBIN=$PREFIX/bin
mkdir -p $PREFIX/bin
rm -rf "$WORK/overInterDyMFoamOA2IS" && cp -r "$HERE/solvers/overInterDyMFoamOA2IS" "$WORK/"
(cd "$WORK/overInterDyMFoamOA2IS" && wmake)
# OA2IS immersed-boundary VOF solver (production path)
rm -rf "$WORK/interIBFoamOA2IS" && cp -r "$HERE/solvers/interIBFoamOA2IS" "$WORK/"
(cd "$WORK/interIBFoamOA2IS" && wmake)
ls -la $PREFIX/lib $PREFIX/bin
echo BUILD_OK

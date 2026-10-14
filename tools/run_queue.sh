#!/bin/bash
# Run several prepared cases one after another (serial or parallel each) and
# post-process each one when it finishes.
# usage: tools/run_queue.sh <nprocs> <case_dir> [<case_dir> ...]
# env:   SOLVER (default /opt/oa2is-openfoam-patches/bin/interIBFoamOA2IS)
#        POST=0 to skip post-processing, FRAMES=1 to render frames + MP4
HERE=$(cd "$(dirname "$0")/.." && pwd)
NP=$1; shift
export SOLVER=${SOLVER:-/opt/oa2is-openfoam-patches/bin/interIBFoamOA2IS}
for c in "$@"; do
  echo "$(date -Is) START $c"
  "$HERE/tools/run_case.sh" "$c" "$NP" || echo "$(date -Is) SOLVER FAILED $c"
  echo "$(date -Is) END $c"
  if [ "${POST:-1}" = 1 ]; then
    if [ "${FRAMES:-1}" = 1 ]; then
      python3 "$HERE/tools/postprocess.py" "$c" --mp4 > "$c/log.postprocess" 2>&1 || echo "POST FAILED $c"
    else
      python3 "$HERE/tools/postprocess.py" "$c" > "$c/log.postprocess" 2>&1 || echo "POST FAILED $c"
    fi
  fi
done
echo "$(date -Is) QUEUE DONE"

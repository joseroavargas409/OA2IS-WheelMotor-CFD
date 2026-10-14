#!/bin/bash
# Read-only integrity check of the ORIGINAL repositories against the baseline
# recorded before any work started (docs/integrity/). Nothing is written inside
# the original repositories: git runs with --no-optional-locks (no index
# refresh) and files are only read by sha256sum.
# usage: tools/verify_originals.sh <dir_with_clones> [baseline_dir]
REF=${1:?dir with the three clones}
BASE=${2:-$(cd "$(dirname "$0")/.." && pwd)/docs/integrity}
ok=1
for r in wheel-motor-inspector komatsu-oil-sensors-v4 komatsu-oil-sensors-v4-1; do
  d="$REF/$r"
  head=$(git -C "$d" --no-optional-locks rev-parse HEAD)
  want=$(cat "$BASE/$r.HEAD")
  st=$(git -C "$d" --no-optional-locks status --porcelain --untracked-files=all)
  sums=$(cd "$d" && sha256sum --quiet -c "$BASE/$r.sha256" 2>&1)
  nfiles=$(git -C "$d" --no-optional-locks ls-files | wc -l)
  nbase=$(wc -l < "$BASE/$r.sha256")
  remote=$(git -C "$d" --no-optional-locks ls-remote origin HEAD 2>/dev/null | cut -f1)
  echo "== $r"
  echo "   HEAD      $head  (baseline $want)  $( [ "$head" = "$want" ] && echo OK || { echo CHANGED; ok=0; })"
  echo "   status    $( [ -z "$st" ] && echo 'clean (no modified/untracked files) OK' || { echo "DIRTY: $st"; ok=0; })"
  echo "   sha256    $nbase files checked: $( [ -z "$sums" ] && echo 'all identical OK' || { echo "MISMATCH: $sums"; ok=0; })"
  echo "   tracked   $nfiles files in git (baseline list $nbase)"
  [ -n "$remote" ] && echo "   GitHub    origin HEAD $remote $( [ "$remote" = "$want" ] && echo '(= baseline, nothing pushed) OK' || echo '(differs from baseline: changed upstream by its owner?)')"
done
[ $ok = 1 ] && echo "RESULT: the original repositories are UNMODIFIED" || echo "RESULT: DIFFERENCES FOUND"

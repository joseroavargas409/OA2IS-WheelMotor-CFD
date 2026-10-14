#!/usr/bin/env python3
"""Build the multi-zone gear-case mesh.

usage: python3 tools/build_mesh.py <case_dir> [--res coarse|fine] [--backlash-high 0.004] [--backlash-low 0.003]
"""
import argparse
import sys
import time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from oa2is_cfd.meshing import PRESETS, build_all          # noqa: E402
from oa2is_cfd.model import default_case                    # noqa: E402

ap = argparse.ArgumentParser()
ap.add_argument('case')
ap.add_argument('--res', default='coarse', choices=sorted(PRESETS))
ap.add_argument('--backlash-high', type=float, default=0.004)
ap.add_argument('--backlash-low', type=float, default=0.003)
a = ap.parse_args()
t = time.time()
st = build_all(Path(a.case), default_case((a.backlash_high, a.backlash_low)), PRESETS[a.res],
               log=lambda s: print(s, flush=True))
print(f"total cells {st['total_cells']}  ({time.time() - t:.0f} s)", flush=True)

"""Study 4 pre-registration pilot: the UNMUTATED subject on Study 1 base inputs (seeds 0-1) only,
to check the implementation and the harness. No mutant is generated or run here."""
import json, os, sys, pickle
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.insert(0, os.path.join(ROOT, "code"))
import numpy as np, shapely
from shapely.geometry import Polygon
import vor_common as vc
from vor_jobs import build_input, REFDIR, _fn
import hp_voronoi as H
out = {}
for cls in vc.CLASSES:
    for n in (12, 40):
        for seed in (0, 1):
            iid = f"{cls}|{n}|{seed}|base"
            P, _, _ = build_input(iid); bx = vc.clip_box(P)
            ref = pickle.load(open(os.path.join(REFDIR, _fn(iid) + ".pkl"), "rb"))["cells"]
            cells = [None if c is None else Polygon(c) for c in H.voronoi_cells(P, bx)]
            r = vc.compare_to_reference(P, cells, ref, bx)
            ok = not r["missing"] and not r["misassign"] and r["worst"] <= 1e-9
            out.setdefault(cls, []).append((n, seed, ok, r["worst"]))
for k, v in out.items():
    print(k, sum(x[2] for x in v), "/", len(v), "max worst %.1e" % max(x[3] for x in v))
json.dump(out, open(os.path.join(ROOT, "study4", "results", "pilot_unmutated.json"), "w"), indent=1)

"""Job functions for the Voronoi study (run inside crash-tolerant workers)."""
import os, pickle, time, warnings
import numpy as np
import shapely
import vor_common as vc

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
REFDIR = os.path.join(ROOT, "results", "vor_ref")
SUTDIR = os.path.join(ROOT, "results", "vor_sut")
os.makedirs(REFDIR, exist_ok=True)
os.makedirs(SUTDIR, exist_ok=True)
warnings.filterwarnings("ignore")


def build_input(iid):
    cls, n, seed, mr = iid.split("|")
    P = vc.generate(cls, int(n), int(seed))
    if mr == "base":
        return P, None, None
    Q, back, perm = vc.apply_mr(mr, P, int(seed))
    return Q, back, perm


def _fn(iid):
    return iid.replace("|", "__")


def job_ref(job):
    iid = job["iid"]
    P, _, _ = build_input(iid)
    bx = vc.clip_box(P)
    t = time.time()
    cells, areas, _ = vc.reference_cells(P, bx)
    dt = time.time() - t
    with open(os.path.join(REFDIR, _fn(iid) + ".pkl"), "wb") as f:
        pickle.dump({"cells": cells, "bx": bx}, f)
    return {"secs": dt, "nuniq": len(cells)}


def job_sut(job):
    iid, sut = job["iid"], job["sut"]
    P, _, _ = build_input(iid)
    bx = vc.clip_box(P)
    with open(os.path.join(REFDIR, _fn(iid) + ".pkl"), "rb") as f:
        ref = pickle.load(f)
    t = time.time()
    cells = vc.SUTS[sut](P, bx)
    dt = time.time() - t
    cells = [None if c is None else c for c in cells]
    with open(os.path.join(SUTDIR, _fn(iid) + "__" + sut + ".pkl"), "wb") as f:
        pickle.dump([None if c is None else shapely.to_wkb(c) for c in cells], f)
    r = vc.compare_to_reference(P, cells, ref["cells"], bx)
    r["secs"] = dt
    r["ncells"] = sum(c is not None and not c.is_empty for c in cells)
    r["n"] = len(P)
    return r

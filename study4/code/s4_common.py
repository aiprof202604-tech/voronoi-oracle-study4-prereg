"""Study 4 shared definitions (fixed before registration)."""
import os, sys, pickle, signal
from fractions import Fraction
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(os.path.dirname(HERE))
sys.path.insert(0, os.path.join(ROOT, "code")); sys.path.insert(0, HERE)
import vor_common as vc
from shapely.geometry import Polygon
RES = os.path.join(ROOT, "study4", "results")
CLASSES = ["G1_uniform", "G2_cluster", "G3_grid", "G4_collinear", "G5_cocircular", "G9_tinyscale"]
NS = (12, 40)
SEEDS = (300, 301, 302, 303)
MRS = ["MR-T", "MR-R", "MR-S", "MR-P", "MR-F"]
EXACT_MRS = ["MR-R", "MR-S", "MR-P", "MR-F"]
TAU = 1e-9
RUN_TIMEOUT_S = 5


def bases():
    return [f"{c}|{n}|{s}" for c in CLASSES for n in NS for s in SEEDS]


def build(bid, mr):
    cls, n, seed = bid.split("|"); n, seed = int(n), int(seed)
    P = vc.generate(cls, n, seed)
    if mr == "base":
        return P, (lambda g: g), None
    return vc.apply_mr(mr, P, seed)


def to_polys(cells):
    out = []
    for c in cells:
        if c is None or len(c) < 3:
            out.append(None); continue
        try:
            g = Polygon(c)
            out.append(g if g.is_valid and not g.is_empty else g.buffer(0))
        except Exception:
            out.append(None)
    return out


def def_checks(P, cells, bx):
    """Definition checks P1-P3 of Study 1 (code/property_oracles.py), applied to shapely polygons."""
    Abox = Fraction(bx[2] - bx[0]) * Fraction(bx[3] - bx[1])
    F = [(Fraction(float(x)), Fraction(float(y))) for x, y in P]
    groups = {}
    for i, k in enumerate(F):
        groups.setdefault(k, []).append(i)
    uniq = list(groups)
    Pf = np.array([[float(a), float(b)] for a, b in uniq])

    def area(poly):
        A = Fraction(0)
        for k in range(len(poly)):
            u, v = poly[k], poly[(k + 1) % len(poly)]
            A += u[0] * v[1] - v[0] * u[1]
        return abs(A) / 2

    def poly_of(g):
        if g is None or g.is_empty:
            return None
        if g.geom_type != "Polygon":
            g = max(g.geoms, key=lambda h: h.area)
        c = list(g.exterior.coords)[:-1]
        return [(Fraction(x), Fraction(y)) for x, y in c] if len(c) >= 3 else None
    p1 = 0; p3 = 0; tot = Fraction(0)
    for k, idxs in groups.items():
        cand = [poly_of(cells[i]) if i < len(cells) else None for i in idxs]
        cand = [c for c in cand if c]
        if not cand:
            p1 += 1; continue
        tot += area(cand[0])
        for poly in cand:
            m = len(poly)
            cx = sum(v[0] for v in poly) / m; cy = sum(v[1] for v in poly) / m
            pts = [(cx, cy)] + [((cx + v[0]) / 2, (cy + v[1]) / 2) for v in poly]
            bad = False
            for (x, y) in pts:
                d0 = (x - k[0]) ** 2 + (y - k[1]) ** 2
                fx, fy = float(x), float(y)
                df = (Pf[:, 0] - fx) ** 2 + (Pf[:, 1] - fy) ** 2
                for j in np.nonzero(df <= float(d0) * (1 + 1e-6) + 1e-300)[0]:
                    q = uniq[j]
                    if q != k and (x - q[0]) ** 2 + (y - q[1]) ** 2 < d0:
                        bad = True; break
                if bad: break
            if bad:
                p3 += 1; break
    p2 = float(abs(tot - Abox) / Abox)
    return {"P1": p1, "P2": p2, "P3": p3, "flag": bool(p1 > 0 or p2 > TAU or p3 > 0)}


class _TO(Exception):
    pass


def _alarm(signum, frame):
    raise _TO()


def run_subject(fn, P, bx):
    """Run subject function fn(points, box) with a timeout; returns (status, cells)."""
    signal.signal(signal.SIGALRM, _alarm); signal.alarm(RUN_TIMEOUT_S)
    try:
        cells = fn(P, bx)
        signal.alarm(0)
        if not isinstance(cells, list) or len(cells) != len(P):
            return "BADSHAPE", None
        return "OK", cells
    except _TO:
        return "TIMEOUT", None
    except Exception as e:
        signal.alarm(0)
        return "ERR:" + type(e).__name__, None

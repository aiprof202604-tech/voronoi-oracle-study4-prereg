"""Voronoi study: input generators, exact reference, SUT adapters, comparison.

Every SUT adapter returns a list `cells` aligned with the input order: cells[i] is a
shapely Polygon (possibly empty) for site i clipped to the box, or None if the SUT
produced no cell for that site. Exceptions propagate and are recorded as type E.
"""
from __future__ import annotations
from fractions import Fraction
import math
import zlib
import numpy as np
import shapely
from shapely.geometry import Polygon, MultiPoint, box as sbox

# ---------------------------------------------------------------- generators
CLASSES = ["G1_uniform", "G2_cluster", "G3_grid", "G4_collinear", "G5_cocircular",
           "G6_duplicate", "G7_nearduplicate", "G8_largeoffset", "G9_tinyscale"]


def _pyth_points(R):
    """Integer points on the circle x^2+y^2=R^2 (exactly cocircular)."""
    pts = []
    for x in range(-R, R + 1):
        y2 = R * R - x * x
        y = math.isqrt(y2)
        if y * y == y2:
            pts.append((x, y))
            if y:
                pts.append((x, -y))
    return sorted(set(pts))


def generate(cls: str, n: int, seed: int) -> np.ndarray:
    rng = np.random.default_rng([zlib.crc32(cls.encode()), n, seed])
    if cls == "G1_uniform":
        P = rng.random((n, 2))
    elif cls == "G2_cluster":
        k = max(2, n // 10)
        C = rng.random((k, 2))
        P = C[rng.integers(0, k, n)] + rng.normal(0, 1e-3, (n, 2))
    elif cls == "G3_grid":
        s = math.ceil(math.sqrt(n))
        g = np.array([(i, j) for i in range(s) for j in range(s)], float)
        P = g[rng.permutation(len(g))[:n]]
    elif cls == "G4_collinear":
        t = rng.random(n)
        a, b = rng.random(2), rng.random(2)
        P = np.outer(1 - t, a) + np.outer(t, b)
        P[: max(1, n // 4)] = rng.random((max(1, n // 4), 2))  # mostly collinear + some free
    elif cls == "G5_cocircular":
        circ = np.array(_pyth_points(5525), float)  # 5525 has many representations
        P = circ[rng.permutation(len(circ))[:n]]
        if len(P) < n:
            P = np.vstack([P, rng.random((n - len(P), 2)) * 2000 - 1000])
    elif cls == "G6_duplicate":
        m = n - max(1, n // 10)
        U = rng.random((m, 2))
        P = np.vstack([U, U[rng.integers(0, m, n - m)]])
        P = P[rng.permutation(n)]
    elif cls == "G7_nearduplicate":
        m = n - max(1, n // 10)
        U = rng.random((m, 2))
        D = U[rng.integers(0, m, n - m)] + rng.choice([-1, 1], (n - m, 2)) * 1e-9
        P = np.vstack([U, D])[rng.permutation(n)]
    elif cls == "G8_largeoffset":
        P = rng.random((n, 2)) * 1000.0 + np.array([512345.0, 4012345.0])
    elif cls == "G9_tinyscale":
        P = rng.random((n, 2)) * 1e-7
    else:
        raise ValueError(cls)
    return np.ascontiguousarray(P, dtype=float)


def clip_box(P: np.ndarray):
    lo, hi = P.min(0), P.max(0)
    ext = float(max(hi - lo)) or 1.0
    m = 0.25 * ext
    return (float(lo[0] - m), float(lo[1] - m), float(hi[0] + m), float(hi[1] + m))


# ---------------------------------------------------------------- metamorphic transforms
# each returns (P', forward map for polygons, inverse map for polygons, permutation idx or None)
MRS = ["MR-T", "MR-R", "MR-S", "MR-P", "MR-F", "MR-D"]


def apply_mr(mr: str, P: np.ndarray, seed: int):
    rng = np.random.default_rng([seed, 99])
    ext = float(max(P.max(0) - P.min(0))) or 1.0
    if mr == "MR-T":
        d = np.array([1000.0 * ext, -700.0 * ext])
        return P + d, (lambda g: shapely.affinity.translate(g, -d[0], -d[1])), None
    if mr == "MR-R":   # exact 90-degree rotation
        Q = np.column_stack([-P[:, 1], P[:, 0]])
        return Q, (lambda g: shapely.affinity.affine_transform(g, [0, 1, -1, 0, 0, 0])), None
    if mr == "MR-S":
        return P * 1024.0, (lambda g: shapely.affinity.scale(g, 1 / 1024, 1 / 1024, origin=(0, 0))), None
    if mr == "MR-P":
        perm = rng.permutation(len(P))
        return P[perm], (lambda g: g), perm
    if mr == "MR-F":
        Q = P.copy(); Q[:, 0] = -Q[:, 0]
        return Q, (lambda g: shapely.affinity.scale(g, -1, 1, origin=(0, 0))), None
    if mr == "MR-D":
        j = int(rng.integers(0, len(P)))
        return np.vstack([P, P[j:j + 1]]), (lambda g: g), np.r_[np.arange(len(P)), j]
    raise ValueError(mr)


# ---------------------------------------------------------------- exact reference
def _clip(poly, a, b, c):
    """Keep part of polygon (list of Fraction pairs) with a*x+b*y <= c."""
    out = []
    n = len(poly)
    for k in range(n):
        p, q = poly[k], poly[(k + 1) % n]
        fp = a * p[0] + b * p[1] - c
        fq = a * q[0] + b * q[1] - c
        if fp <= 0:
            out.append(p)
        if (fp < 0 < fq) or (fq < 0 < fp):
            t = fp / (fp - fq)
            out.append((p[0] + t * (q[0] - p[0]), p[1] + t * (q[1] - p[1])))
    return out


def reference_cells(P: np.ndarray, bx):
    """Exact Voronoi cells of the *float* input values, clipped to box bx.

    Returns dict: unique coordinate tuple -> list of float vertex tuples, and exact areas.
    """
    F = [(Fraction(float(x)), Fraction(float(y))) for x, y in P]
    uniq = sorted(set(F))
    x0, y0, x1, y1 = (Fraction(v) for v in bx)
    base = [(x0, y0), (x1, y0), (x1, y1), (x0, y1)]
    cells, areas = {}, {}
    # sort others by distance for early shrink (performance only)
    for p in uniq:
        poly = base
        others = sorted((q for q in uniq if q != p), key=lambda q: (q[0] - p[0]) ** 2 + (q[1] - p[1]) ** 2)
        for q in others:
            # |x-p|^2 <= |x-q|^2  <=>  2(q-p).x <= |q|^2-|p|^2
            a, b = 2 * (q[0] - p[0]), 2 * (q[1] - p[1])
            c = q[0] ** 2 + q[1] ** 2 - p[0] ** 2 - p[1] ** 2
            # quick reject: if all vertices satisfy, skip
            if all(a * v[0] + b * v[1] <= c for v in poly):
                continue
            poly = _clip(poly, a, b, c)
            if not poly:
                break
        A = Fraction(0)
        for k in range(len(poly)):
            u, v = poly[k], poly[(k + 1) % len(poly)]
            A += u[0] * v[1] - v[0] * u[1]
        areas[p] = abs(A) / 2
        cells[p] = [(float(u), float(v)) for u, v in poly]
    return cells, areas, F


# ---------------------------------------------------------------- SUT adapters
def _ghosts(P, bx):
    cx, cy = (bx[0] + bx[2]) / 2, (bx[1] + bx[3]) / 2
    D = 8.0 * math.hypot(bx[2] - bx[0], bx[3] - bx[1])
    return np.array([[cx - D, cy - D], [cx + D, cy - D], [cx + D, cy + D], [cx - D, cy + D]])


def sut_scipy(P, bx):
    from scipy.spatial import Voronoi
    G = np.vstack([P, _ghosts(P, bx)])
    vor = Voronoi(G)
    B = sbox(*bx)
    cells = []
    for i in range(len(P)):
        reg = vor.regions[vor.point_region[i]]
        if not reg or -1 in reg:
            cells.append(None); continue
        poly = Polygon(vor.vertices[reg])
        cells.append(poly.buffer(0).intersection(B) if not poly.is_valid else poly.intersection(B))
    return cells


def sut_shapely(P, bx):
    B = sbox(*bx)
    out = shapely.voronoi_polygons(MultiPoint(P), extend_to=B, ordered=True)
    geoms = list(out.geoms)
    # ordered=True: polygons follow input point order (duplicates collapse -> fewer polygons)
    cells = [None] * len(P)
    if len(geoms) == len(P):
        for i, g in enumerate(geoms):
            cells[i] = g.intersection(B)
    else:
        # documented behaviour for repeated points: one polygon per distinct point, order of first
        # occurrence; map by containment of the site (adapter cannot do better without the spec)
        seen = {}
        k = 0
        for i, p in enumerate(map(tuple, P)):
            if p in seen:
                continue
            seen[p] = k; k += 1
        for i, p in enumerate(map(tuple, P)):
            j = seen[p]
            if j < len(geoms) and cells[i] is None:
                cells[i] = geoms[j].intersection(B)
    return cells


def sut_pyvoronoi(P, bx):
    import pyvoronoi
    G = np.vstack([P, _ghosts(P, bx)])
    c = G.mean(0)
    R = float(np.abs(G - c).max())
    factor = (2 ** 30) / R           # map into the int32 range required by Boost.Polygon
    pv = pyvoronoi.Pyvoronoi(factor)
    Gc = G - c
    for x, y in Gc:
        pv.AddPoint([float(x), float(y)])
    pv.Construct()
    verts = pv.GetVertices()
    B = sbox(*bx)
    cells = [None] * len(P)
    for cell in pv.GetCells():
        if not cell.contains_point or cell.is_open:
            continue
        i = cell.site
        if i >= len(P):
            continue
        vs = [(verts[v].X + c[0], verts[v].Y + c[1]) for v in cell.vertices if v != -1]
        if len(vs) >= 3:
            cells[i] = Polygon(vs).intersection(B)
    return cells


def sut_freud(P, bx):
    import freud
    c = np.array([(bx[0] + bx[2]) / 2, (bx[1] + bx[3]) / 2])
    L = 8.0 * math.hypot(bx[2] - bx[0], bx[3] - bx[1])
    fbox = freud.box.Box.square(L)
    X = np.zeros((len(P), 3))
    X[:, :2] = P - c
    v = freud.locality.Voronoi()
    v.compute((fbox, X))
    B = sbox(*bx)
    cells = []
    for poly in v.polytopes:
        cells.append(Polygon(poly[:, :2] + c).intersection(B))
    return cells


def sut_geovoronoi(P, bx):
    from geovoronoi import voronoi_regions_from_coords
    B = sbox(*bx)
    region_polys, region_pts = voronoi_regions_from_coords(P, B)
    cells = [None] * len(P)
    for r, idxs in region_pts.items():
        for i in idxs:
            cells[i] = region_polys[r]
    return cells


SUTS = {"scipy": sut_scipy, "shapely": sut_shapely, "pyvoronoi": sut_pyvoronoi,
        "freud": sut_freud, "geovoronoi": sut_geovoronoi}


# ---------------------------------------------------------------- comparison
def compare_to_reference(P, cells, ref_cells, bx, taus=(1e-12, 1e-9, 1e-6)):
    """Return per-tau failure info: worst normalised symmetric difference and failure type."""
    Abox = (bx[2] - bx[0]) * (bx[3] - bx[1])
    keys = [(Fraction(float(x)), Fraction(float(y))) for x, y in P]
    groups = {}
    for i, k in enumerate(keys):
        groups.setdefault(k, []).append(i)
    worst = 0.0
    missing = 0
    misassign = 0
    for k, idxs in groups.items():
        R = Polygon(ref_cells[k]) if len(ref_cells[k]) >= 3 else Polygon()
        ds = []
        for i in idxs:
            g = cells[i]
            if g is None or g.is_empty:
                ds.append(None); continue
            try:
                d = g.symmetric_difference(R).area / Abox
            except Exception:
                d = abs(g.area - R.area) / Abox
            ds.append(d)
        valid = [d for d in ds if d is not None]
        if not valid:
            missing += 1
            worst = max(worst, R.area / Abox)
            continue
        best = min(valid)
        if len(idxs) > 1:
            # duplicates: accept if one copy matches and others are same/empty
            worst = max(worst, best)
            bad_dups = [d for d in valid if d > 1e-6 and d != best]
            if bad_dups:
                misassign += 1
                worst = max(worst, max(bad_dups))
        else:
            worst = max(worst, best)
            if best > 0.5 * max(R.area / Abox, 1e-300) and best > 1e-6:
                misassign += 1
    return {"worst": worst, "missing": missing, "misassign": misassign}


def cells_distance(cA, cB, bx, pairs):
    """Max normalised symmetric difference between cells cA[i] and cB[j] for (i,j) in pairs."""
    Abox = (bx[2] - bx[0]) * (bx[3] - bx[1])
    worst = 0.0
    for i, j in pairs:
        a, b = cA[i], cB[j]
        ea = a is None or a.is_empty
        eb = b is None or b.is_empty
        if ea and eb:
            continue
        if ea != eb:
            worst = max(worst, ((b if ea else a).area) / Abox); continue
        try:
            d = a.symmetric_difference(b).area / Abox
        except Exception:
            d = abs(a.area - b.area) / Abox
        worst = max(worst, d)
    return worst

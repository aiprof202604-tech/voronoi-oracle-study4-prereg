"""Study 4 subject: a small, readable Voronoi implementation by half-plane clipping in double
precision (pure Python). Cell i is the clip box intersected with the half-planes
{x : |x - p_i| <= |x - p_j|} for all other distinct sites p_j. It is the subject into which
first-order faults are seeded by mutate.py. Returns one polygon (list of (x, y)) or None per site."""


def clip(poly, a, b, c):
    """Keep the part of convex polygon `poly` with a*x + b*y <= c (Sutherland-Hodgman, one edge)."""
    out = []
    n = len(poly)
    for k in range(n):
        ux, uy = poly[k]
        vx, vy = poly[(k + 1) % n]
        fu = a * ux + b * uy - c
        fv = a * vx + b * vy - c
        if fu <= 0:
            out.append((ux, uy))
        if (fu < 0 and fv > 0) or (fu > 0 and fv < 0):
            t = fu / (fu - fv)
            out.append((ux + t * (vx - ux), uy + t * (vy - uy)))
    return out


def voronoi_cells(points, box):
    x0, y0, x1, y1 = box
    pts = [(float(p[0]), float(p[1])) for p in points]
    cells = []
    for i in range(len(pts)):
        px, py = pts[i]
        poly = [(x0, y0), (x1, y0), (x1, y1), (x0, y1)]
        for j in range(len(pts)):
            if j == i:
                continue
            qx, qy = pts[j]
            if qx == px and qy == py:
                continue
            a = 2.0 * (qx - px)
            b = 2.0 * (qy - py)
            c = qx * qx + qy * qy - px * px - py * py
            poly = clip(poly, a, b, c)
            if len(poly) < 3:
                break
        cells.append(poly if len(poly) >= 3 else None)
    return cells

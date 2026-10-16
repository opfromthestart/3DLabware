#!/usr/bin/env python3
"""Self-tests for threelib: every primitive is checked for closedness,
outward winding, and volume against a closed-form value."""

import math
import sys

import threelib as t


def _tri_areas(mesh):
    out = []
    for tri in mesh.tris:
        a, b, c = (mesh.verts[j] for j in tri)
        u = (b[0] - a[0], b[1] - a[1], b[2] - a[2])
        w = (c[0] - a[0], c[1] - a[1], c[2] - a[2])
        n = (u[1] * w[2] - u[2] * w[1], u[2] * w[0] - u[0] * w[2],
             u[0] * w[1] - u[1] * w[0])
        out.append(0.5 * math.sqrt(sum(x * x for x in n)))
    return out


fails = []


def check(label, ok, detail=""):
    print(f"  [{'PASS' if ok else 'FAIL'}] {label}" + (f" -- {detail}" if detail else ""))
    if not ok:
        fails.append(label)


def closed_shells(m):
    """Every vertex must be used by an even number of edges (manifold-ish)."""
    from collections import defaultdict
    edges = defaultdict(int)
    for a, b, c in m.tris:
        for u, v in ((a, b), (b, c), (c, a)):
            edges[(min(u, v), max(u, v))] += 1
    return all(n == 2 for n in edges.values())


def near(label, got, want, rel=2e-2):
    check(label, abs(got - want) <= abs(want) * rel,
          f"{got:.4f} vs {want:.4f} (rel tol {rel:g})")


print("== threelib primitives ==")

m = t.Mesh()
t.box(m, 0, 0, 0, 3, 4, 5)
check("box closed", closed_shells(m))
near("box volume", m.signed_volume(), 60.0, 1e-9)

m = t.Mesh()
t.cylinder(m, 10, 0, 7, 128)
check("cylinder closed", closed_shells(m))
check("cylinder outward", m.signed_volume() > 0)
near("cylinder volume", m.signed_volume(), math.pi * 100 * 7, 5e-3)

m = t.Mesh()
t.tube(m, 20, 25, 0, 10, 128)
check("tube closed", closed_shells(m))
check("tube outward", m.signed_volume() > 0)
near("tube volume", m.signed_volume(), math.pi * (625 - 400) * 10, 5e-3)

m = t.Mesh()
t.cone(m, 10, 0, 0, 12, 128)
check("cone closed", closed_shells(m))
near("cone volume", m.signed_volume(), math.pi * 100 * 12 / 3, 5e-3)

m = t.Mesh()
t.cone(m, 10, 4, 0, 12, 128)
near("frustum volume", m.signed_volume(),
     math.pi * 12 / 3 * (100 + 40 + 16), 5e-3)

def rho(k):
    """Area ratio of an inscribed k-gon to the circle it approximates."""
    return (k / 2.0) * math.sin(2 * math.pi / k) / math.pi


m = t.Mesh()
t.sphere(m, 10, 0, seg=96, rings=48)
check("sphere closed", closed_shells(m))
near("sphere volume", m.signed_volume(),
     4 / 3 * math.pi * 1000 * rho(96), 1e-2)

# A torus is faceted in BOTH directions, so compare against the product of
# the two inscribed-polygon ratios rather than a blanket tolerance.
m = t.Mesh()
t.torus(m, 20, 5, seg=96, sides=32)
check("torus closed", closed_shells(m))
near("torus volume", m.signed_volume(),
     2 * math.pi ** 2 * 20 * 25 * rho(96) * rho(32), 1e-3)

for ax in ("x", "y", "z"):
    m = t.Mesh()
    t.wedge(m, 0, 0, 0, 10, 10, 10, axis=ax)
    check(f"wedge({ax}) closed", closed_shells(m))
    check(f"wedge({ax}) outward", m.signed_volume() > 0)
    check(f"wedge({ax}) volume is half the box",
          abs(m.signed_volume() - 500.0) < 1e-9, f"{m.signed_volume():.2f} mm^3")

# --- Tessellation --------------------------------------------------------
for r, tol in ((10, 0.01), (55, 0.01), (2, 0.001)):
    n = t.auto_segments(r, tol)
    sag = r * (1 - math.cos(math.pi / n))
    check(f"auto_segments(r={r}, tol={tol})", sag <= tol,
          f"n={n}, sagitta={sag:.6f} <= {tol}")

# --- Clipping ------------------------------------------------------------
disc = t.disc_polygon(20, 180)
sq = t.rect_polygon(-40, -40, 40, 40)
clipped = t.clean_polygon(t.clip_convex(sq, disc))
check("clip keeps a polygon", len(clipped) >= 3, f"{len(clipped)} points")
check("clip is inside the disc", all(
    abs(math.hypot(x, y)) <= 20 + 1e-9 for x, y in clipped))
n_gon = 180
expect = 0.5 * n_gon * 400 * math.sin(2 * math.pi / n_gon)
near("clipped area == inscribed n-gon", abs(t.polygon_area(clipped)), expect, 1e-3)

check("clipping a far-away box returns nothing",
      t.clip_convex(t.rect_polygon(100, 100, 110, 110), disc) == [])

# --- Rotations preserve volume -------------------------------------------
m = t.Mesh()
t.box(m, -5, -5, -5, 5, 5, 5)
v0 = m.signed_volume()
t.rotate_x(m, 37)
t.rotate_y(m, 11)
t.rotate_z(m, 90)
near("rotation preserves volume", m.signed_volume(), v0, 1e-9)
check("rotation preserves outward winding", m.signed_volume() > 0)

# --- revolve: a solid of revolution against closed forms ------------------
m = t.Mesh()
t.revolve(m, [(0.0, 0.0), (5.0, 10.0)], 512)          # cone on a point
near("revolve makes a cone", m.signed_volume(), math.pi * 25 * 10 / 3, 0.01)
check("revolve winds outward", m.signed_volume() > 0)

m = t.Mesh()
t.revolve(m, [(0.0, 0.0), (2.0, 0.0), (2.0, 5.0), (0.0, 5.0)], 512)
near("revolve caps a cylinder", m.signed_volume(), math.pi * 4 * 5, 0.01)

m = t.Mesh()
t.revolve(m, [(0.0, 0.0), (4.0, 3.0), (0.0, 6.0)], 512)   # bicone, two poles
near("revolve handles two poles", m.signed_volume(), 2 * math.pi * 16 * 3 / 3,
     0.01)

m = t.Mesh()
t.revolve(m, [(0.0, 0.0), (6.0, 0.0), (5.0, 20.0), (2.0, 20.0), (2.0, 40.0),
              (0.0, 40.0)], 512)
want = (math.pi * 20 / 3 * (36 + 30 + 25) + math.pi * 4 * 20)
near("revolve makes a stepped solid", m.signed_volume(), want, 0.05)
check("revolve leaves no degenerate triangles",
      all(area > 1e-9 for area in _tri_areas(m)), "every triangle has area")

# --- Mesh plumbing -------------------------------------------------------
m = t.Mesh()
t.box(m, 0, 0, 0, 1, 1, 1)
m2 = t.Mesh()
t.box(m2, 5, 5, 5, 6, 6, 6)
m.merge(m2)
check("merge concatenates", len(m.verts) == 16 and len(m.tris) == 24)
check("merge keeps shells separate", len(m.shells()) == 2)

t.write_stl("/tmp/_threelib_test.stl", m)
import os
check("STL written", os.path.getsize("/tmp/_threelib_test.stl") == 84 + 50 * 24)
os.remove("/tmp/_threelib_test.stl")

print("\n" + ("ALL CHECKS PASSED" if not fails else f"FAILURES: {fails}"))
sys.exit(1 if fails else 0)

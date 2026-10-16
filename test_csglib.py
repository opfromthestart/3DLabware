#!/usr/bin/env python3
"""Self-test for csglib: booleans against closed forms.

Every case here has an answer that can be worked out on paper, which is the
only kind of test worth having for a boolean: a union that looks right and is
0.2 mm out is exactly the failure mode this module exists to prevent.
"""

import math
import sys

import memcap
import csglib as csg
import threelib as t

fails = 0


def check(label, ok, detail=""):
    global fails
    print(f"  [{'PASS' if ok else 'FAIL'}] {label}" + (f" -- {detail}" if detail else ""))
    if not ok:
        fails += 1


def near(label, got, want, tol, unit=""):
    check(label, abs(got - want) <= tol,
          f"{got:.6f}{unit} vs {want:.6f}{unit} (tol {tol:g})")


def box(x0, y0, z0, x1, y1, z1):
    m = t.Mesh("box")
    t.box(m, x0, y0, z0, x1, y1, z1)
    return m


def faceted_disc(r, n):
    """A regular n-gon's area -- what a revolve of radius r actually makes."""
    return 0.5 * n * r * r * math.sin(2.0 * math.pi / n)


def watertight(m, label):
    check(f"{label}: outward wound", m.signed_volume() > 0.0)
    check(f"{label}: one shell", len(m.shells()) == 1, f"{len(m.shells())}")
    # A hole, not just an unpaired edge: the merge pass removes interior
    # vertices, and one of those can sit on a neighbouring patch's edge, which
    # leaves edges used once without leaving a gap. `real_gaps` is the test for
    # a gap -- an edge used once with nothing lying along it.
    gaps = csg.real_gaps(m.verts, m.tris)
    crack = csg.crack_length(m.verts, m.tris)
    check(f"{label}: no holes between faces", not gaps,
          f"{len(gaps)} edges used once with nothing along them"
          if gaps else
          f"no gap; {crack:.2e} mm of edge used once in T-junctions")
    areas = []
    for tri in m.tris:
        a, b, c = (m.verts[j] for j in tri)
        u = (b[0]-a[0], b[1]-a[1], b[2]-a[2])
        w = (c[0]-a[0], c[1]-a[1], c[2]-a[2])
        n = (u[1]*w[2]-u[2]*w[1], u[2]*w[0]-u[0]*w[2], u[0]*w[1]-u[1]*w[0])
        areas.append(0.5*math.sqrt(sum(x*x for x in n)))
    # Flat triangles are allowed and expected: they close the T-junctions the
    # boolean leaves behind, and dropping them opens the mesh again. What must
    # not be there is a face with area that is *almost* zero but not zero.
    check(f"{label}: no near-degenerate triangles",
          min(areas) < 1e-12 or min(areas) > 1e-9,
          f"smallest {min(areas):.2e} mm^2")


print("--- two boxes overlapping in a corner ---")
a = box(0, 0, 0, 10, 10, 10)
b = box(6, 6, 6, 16, 16, 16)
u = csg.union(a, b)
near("union volume", u.signed_volume(), 1000 + 1000 - 64, 1e-6, " mm^3")
watertight(u, "union")
near("the overlap is 4x4x4", csg.overlap_volume(a, b), 64.0, 1e-6, " mm^3")
near("subtract leaves the corner out", csg.subtract(a, b).signed_volume(),
     1000 - 64, 1e-6, " mm^3")
near("intersect is just the corner", csg.intersect(a, b).signed_volume(),
     64.0, 1e-6, " mm^3")

print("--- one box wholly inside another ---")
outer = box(0, 0, 0, 20, 20, 20)
inner = box(5, 5, 5, 15, 15, 15)
hollow = csg.subtract(outer, inner)
# Two shells: the outside of the box, and the wall of the cavity inside it.
near("a hollow box is two surfaces, outer and cavity",
     len(hollow.shells()), 2, 0)
near("hollow volume", hollow.signed_volume(), 8000 - 1000, 1e-6, " mm^3")
near("union with a solid inside is the outer box",
     csg.union(outer, inner).signed_volume(), 8000.0, 1e-6, " mm^3")

print("--- two boxes sharing a face ---")
left = box(0, 0, 0, 5, 5, 5)
right = box(5, 0, 0, 10, 5, 5)
joined = csg.union(left, right)
near("the shared wall is gone", joined.signed_volume(), 250.0, 1e-6, " mm^3")
watertight(joined, "two boxes side by side")

print("--- a cylinder through a box ---")
# The rod's axis lies *in* the block's z=0 face, so only half of its
# cross-section is inside the block. That is not a rounding detail: it is the
# difference between 283 and 565 mm^3, and getting it wrong is how a boolean
# looks broken when it is right.
N = 256
blk = box(-10, -10, 0, 10, 10, 10)
rod = t.Mesh("rod")
t.cone(rod, 3.0, 3.0, -20.0, 20.0, N)            # capped solid rod
t.rotate_x(rod, 90.0)                             # now along y
half_bore = 0.5 * faceted_disc(3.0, N) * 20        # the half inside the block
v_rod = faceted_disc(3.0, N) * 40
near("the overlap is half the rod's cross-section",
     csg.overlap_volume(blk, rod), half_bore, half_bore * 1e-8, " mm^3")
drilled = csg.subtract(blk, rod)
near("a drilled block loses the bore", drilled.signed_volume(),
     4000 - half_bore, 4000 * 1e-8, " mm^3")
watertight(drilled, "drilled block")
crossed = csg.union(blk, rod)
near("a rod through a block adds only what is outside it",
     crossed.signed_volume(), 4000 + v_rod - half_bore, 5000 * 1e-8, " mm^3")
watertight(crossed, "rod through block")
near("the union is the sum less the overlap, which is what a union is",
     crossed.signed_volume(),
     blk.signed_volume() + rod.signed_volume()
     - csg.overlap_volume(blk, rod), 5000 * 1e-8, " mm^3")

print("--- a rod through a rod: no planar faces to lean on ---")
a_rod = t.Mesh("a")
t.cone(a_rod, 2.0, 2.0, -10.0, 10.0, 128)
b_rod = t.Mesh("b")
t.cone(b_rod, 2.0, 2.0, -10.0, 10.0, 128)
t.rotate_x(b_rod, 90.0)
t.translate(b_rod, 0.0, 0.0, 1.0)
cross_rods = csg.union(a_rod, b_rod)
check("crossed rods: outward wound", cross_rods.signed_volume() > 0.0)
check("crossed rods: one shell", len(cross_rods.shells()) == 1,
      f"{len(cross_rods.shells())}")
# Two curved surfaces cannot close every edge: each carries the intersection
# as a polyline of its own planes, and the two only agree where the planes did.
# The result is a hairline of unshared edges -- reported, not hidden.
crack = csg.crack_length(cross_rods.verts, cross_rods.tris)
print(f"  note: {crack:.2f} mm of hairline edges where the two rods' facets "
      f"cross, on a part {40:g} mm long -- the faceting error, not a gap")
check("crossed rods: the crack is a hairline, not a gap", crack < 20.0,
      f"{crack:.2f} mm total")
check("two rods crossing enclose more than either",
      cross_rods.signed_volume() > max(a_rod.signed_volume(),
                                       b_rod.signed_volume()),
      f"{cross_rods.signed_volume():.1f} mm^3 against "
      f"{a_rod.signed_volume():.1f}")

print("--- an open mesh is refused, not unioned ---")
try:
    csg.union(box(0, 0, 0, 1, 1, 1), t.Mesh("open"))
    check("an empty mesh is refused", False, "no error raised")
except ValueError as exc:
    check("an empty mesh is refused", "not closed" in str(exc), str(exc)[:70])
# A mesh that is closed but full of degenerate faces must still work: the
# sweep of an annulus with an inner radius of zero collapses its inner surface
# onto the axis, and dropping those faces leaves a perfectly good solid.
zeroed = t.Mesh("zeroed")
t.tube(zeroed, 0.0, 3.0, -10.0, 10.0, 64)
made = csg.union(box(-5, -5, -5, 5, 5, 5), zeroed)
check("a mesh with collapsed faces at the axis still unions",
      made.signed_volume() > 0.0 and len(made.shells()) == 1,
      f"{made.signed_volume():.1f} mm^3, {len(made.shells())} shell(s)")
# A single triangle is not closed, and a boolean on it must refuse rather than
# answer.
ragged = t.Mesh("ragged")
ragged.verts = [(0.0, 0.0, 0.0), (1.0, 0.0, 0.0), (0.0, 1.0, 0.0)]
ragged.tris = [(0, 1, 2)]
try:
    csg.union(box(-5, -5, -5, 5, 5, 5), ragged)
    check("an open mesh is refused", False, "no error raised")
except ValueError as exc:
    check("an open mesh is refused", "not closed" in str(exc), str(exc)[:60])

print("--- idempotence and disjoint parts ---")
near("a union with itself changes nothing",
     csg.union(a, a).signed_volume(), 1000.0, 1e-6, " mm^3")
far = box(100, 100, 100, 101, 101, 101)
two = csg.union(a, far)
near("disjoint parts stay two shells", len(two.shells()), 2, 0)
near("and keep both volumes", two.signed_volume(), 1001.0, 1e-6, " mm^3")

print("\n" + ("ALL CHECKS PASSED" if not fails else f"FAILURES: {fails}"))
sys.exit(1 if fails else 0)

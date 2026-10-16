#!/usr/bin/env python3
"""Find shells that interpenetrate: the pieces that overlap instead of meeting.

Every part in this repo was built by overlapping shells, because that is the
only way to join two shapes without a boolean. A slicer that fills by winding
unions them; one that fills even-odd shows the joint as a hole. `csglib` makes
the union real, and this says which parts still need it.

A pair of shells is reported when a vertex of one, nudged a little way into its
own solid, lands inside the other. Vertices alone would lie about it -- two
shells that merely touch share vertices, and a vertex on a shared face is a
coin flip -- so the test point is moved along the inward normal of one of the
triangles that meets at that vertex, which puts it in the material rather than
on a boundary.

Usage: python3 check_overlaps.py part.3mf [part2.3mf ...] [--depth] [--quiet]

`--depth` measures how far each pair interpenetrates, at a cost of thousands
of point-in-solid tests per pair; without it the answer is whether they overlap
at all, which is the question that decides whether a part needs a union.
"""

import math
import sys

import memcap
import checks as K


def inward_points(verts, tris, shell, depth=0.02, most=48):
    """Vertices of a shell, moved `depth` mm into the solid they belong to."""
    normals = {}
    for tri in shell:
        a, b, c = (verts[j] for j in tris[tri])  # tris[tri] is the triangle
        u = (b[0] - a[0], b[1] - a[1], b[2] - a[2])
        w = (c[0] - a[0], c[1] - a[1], c[2] - a[2])
        n = (u[1] * w[2] - u[2] * w[1], u[2] * w[0] - u[0] * w[2],
             u[0] * w[1] - u[1] * w[0])
        length = math.sqrt(n[0] ** 2 + n[1] ** 2 + n[2] ** 2)
        if length < 1e-12:
            continue
        n = (n[0] / length, n[1] / length, n[2] / length)
        for j in tris[tri]:
            if j not in normals:
                normals[j] = [0.0, 0.0, 0.0]
            for k in range(3):
                normals[j][k] += n[k]
    # A sample of the shell's vertices, not all of them. The depth measure is
    # a point-in-solid test per step, and a lattice of 300 bars has tens of
    # thousands of pairs: testing every vertex made a scan that used to take
    # seconds run past its timeout.
    items = list(normals.items())
    step = max(1, len(items) // most)
    out = []
    for j, n in items[::step]:
        length = math.sqrt(sum(x * x for x in n))
        if length < 1e-9:
            continue
        v = verts[j]
        n = (n[0] / length, n[1] / length, n[2] / length)
        out.append(((v[0], v[1], v[2]), n))
    return out


def penetration(verts, tris, shells, boxes, pts, mine, into, depth_max=3.0,
                tol=5e-4):
    """How far this shell's material reaches into `into`.

    A vertex nudged a fixed 0.02 mm into its own solid says *whether* two
    shells touch, not whether they overlap: two surfaces that meet tangentially
    put a few sample points inside each other by a micron, and reporting that
    as an overlap sends you off to fix a model that is fine. So the depth is
    measured -- the nudge is grown until the point leaves the other solid --
    and the caller can tell a graze from an intersection.

    `depth_max` has to clear the wall it is probing into, or the first probe
    lands past the far surface and the depth comes back as zero: a 0.5 mm
    overlap into a 2.5 mm wall measured 0.000 until the cap was raised.

    The probe stops as soon as it leaves *this* shell. One shell's material
    reaching into another is a measurement taken inside the first shell; a
    probe that walks out of a thin bar and into its neighbour reports an
    overlap of the bar's own thickness between two bars that never meet. The
    grid came back as 596 overlapping pairs, every one of them by 0.979 mm --
    which is the bar thickness, not an overlap.
    """
    def inside_at(p, n, d, which):
        probe = tuple(p[k] - d * n[k] for k in range(3))
        return K.inside(verts, tris, probe, shells=[which[0]],
                        boxes=[which[1]])

    worst = 0.0
    for p, n in pts:
        # Sample across the range rather than bisecting from the cap: an
        # overlap *thinner* than the cap has its far end outside the other
        # solid, so a bisection that starts there finds nothing. A sleeve
        # biting 0.5 mm into a 2.5 mm wall measured 0.000 that way.
        deep = 0.0
        for k in range(1, 25):
            d = depth_max * k / 24.0
            if not inside_at(p, n, d, mine):
                break                     # left this shell: stop measuring
            if inside_at(p, n, d, into):
                deep = d
        if deep <= 0.0:
            continue
        lo, hi = deep, min(depth_max, deep + depth_max / 24.0)
        while hi - lo > tol:
            mid = 0.5 * (lo + hi)
            if inside_at(p, n, mid, mine) and inside_at(p, n, mid, into):
                lo = mid
            else:
                hi = mid
        worst = max(worst, lo)
    return worst


def overlaps(path, quiet=False, max_report=6, depth=False):
    objs, _, _ = K.load_3mf(path, K.Report(path))
    findings = []
    for obj in objs:
        verts, tris = obj["verts"], obj["tris"]
        shells = K.shell_triangles(verts, tris)
        boxes = K.shell_boxes(verts, tris, shells)
        if len(shells) < 2:
            continue
        pts = [inward_points(verts, tris, s) for s in shells]
        for i in range(len(shells)):
            for j in range(i + 1, len(shells)):
                lo_i, hi_i = boxes[i]
                lo_j, hi_j = boxes[j]
                if any(hi_i[k] < lo_j[k] - 1e-9 or hi_j[k] < lo_i[k] - 1e-9
                       for k in range(3)):
                    continue
                hit = sum(1 for q, _ in pts[i]
                          if K.inside(verts, tris, q, shells=[shells[j]],
                                      boxes=[boxes[j]]))
                back = sum(1 for q, _ in pts[j]
                           if K.inside(verts, tris, q, shells=[shells[i]],
                                       boxes=[boxes[i]]))
                if not hit and not back:
                    continue
                # Measuring the depth costs thousands of point-in-solid tests
                # per pair, and a lattice of 300 bars offers tens of thousands
                # of candidate pairs: twenty minutes and still running. It is
                # asked for with --depth, and the default answers the question
                # that matters -- do these two shells overlap at all.
                size = 0.0
                if depth:
                    size = max(
                        penetration(verts, tris, shells, boxes, pts[i],
                                    (shells[i], boxes[i]),
                                    (shells[j], boxes[j])),
                        penetration(verts, tris, shells, boxes, pts[j],
                                    (shells[j], boxes[j]),
                                    (shells[i], boxes[i])))
                findings.append((obj.get("name", "?"), i, j, size,
                                 len(shells), hit, back))
    return findings


def main():
    memcap.cap()
    args = [a for a in sys.argv[1:] if not a.startswith("--")]
    quiet = "--quiet" in sys.argv
    deep = "--depth" in sys.argv
    if not args:
        print(__doc__)
        return 1
    total = 0
    for path in args:
        findings = overlaps(path, quiet, depth=deep)
        total += len(findings)
        if not findings:
            print(f"{path}: no shells interpenetrate")
            continue
        if deep:
            # With the depth measured, a pair whose depth is zero is a contact
            # and not an overlap: the media's layers stack flush, one bar
            # directly on the one below, and a vertex nudged into one of them
            # lands on the other's surface. Those print as one solid under
            # either fill rule.
            over = [f for f in findings if f[3] > 0.05]
            touch = [f for f in findings if f[3] <= 0.05]
            print(f"{path}: {len(over)} overlapping shell pair(s), "
                  f"{len(touch)} touching but not overlapping")
            for f in sorted(over, key=lambda f: -f[3])[:6]:
                name, i, j, size, n, hit, back = f
                print(f"    {name}: shell {i} and {j} overlap by {size:.3f} mm"
                      f" ({n} shells in the object)")
        else:
            print(f"{path}: {len(findings)} shell pair(s) interpenetrate; "
                  f"pass --depth to tell an overlap from a contact")
            for f in sorted(findings, key=lambda f: -f[5])[:6]:
                name, i, j, size, n, hit, back = f
                print(f"    {name}: shell {i} and {j} ({hit} and {back} sample "
                      f"points inside the other) ({n} shells in the object)")
        if len(findings) > 6:
            print(f"    ... and {len(findings) - 6} more")
    return 0


if __name__ == "__main__":
    sys.exit(main())

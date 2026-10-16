#!/usr/bin/env python3
"""Verify a weigh boat 3MF by re-reading it, knowing only what was asked for.

The part is a shallow square dish with a flat floor and walls that lean outward.
What has to hold: the floor is the size and thickness asked for; the rim is
where the lean puts it; the wall is the thickness asked for *across* it, which
is not the same as the thickness across the floor once the wall leans; the
cavity is the depth asked for and is open at the top; and it prints flat without
support.

Everything is measured off the mesh. The wall is measured across itself, by
taking the material's span along x at a height and multiplying by the cosine of
the lean, because a leaning wall is thicker than its horizontal span.

Usage: python3 verify_weighboat.py weighboat_xxxxxxxx.3mf
"""

import math
import sys

import checks as K
import csglib as csg
import memcap
import printability as PR
import threelib as t

# --- What was asked for ----------------------------------------------------
BOAT = 74.0
DEPTH = 15.0
FLARE = 6.0
WALL = 1.2
TOL = 0.05

LEAN = math.atan2(FLARE, DEPTH)


def main():
    memcap.cap()
    path = sys.argv[1] if len(sys.argv) > 1 else None
    if path is None:
        print(__doc__)
        return 2

    objects, _, _ = K.load_3mf(path, K.Report("weighboat"))
    obj = objects[0]
    verts, tris = obj["verts"], obj["tris"]
    rep = K.Report("weighboat")
    shells = K.shell_triangles(verts, tris)

    xs = [v[0] for v in verts]
    zs = [v[2] for v in verts]
    size = (max(xs) - min(xs), max(v[1] for v in verts) - min(v[1] for v in verts),
            max(zs) - min(zs))

    def solid(x, y, z):
        return K.inside(verts, tris, (x, y, z), shells=shells)

    # --- one solid, and sound ---------------------------------------------
    rep.check("the boat is one solid", len(shells) == 1,
              f"{len(shells)} shells in the mesh")
    gaps = csg.real_gaps(verts, tris)
    rep.check("no hole in the surface", not gaps,
              f"{len(gaps)} edges are used once with nothing along them"
              if gaps else
              f"{csg.crack_length(verts, tris):.2e} mm of edge used once, all "
              f"of it in T-junctions")

    # --- the floor ---------------------------------------------------------
    floor_x = [abs(v[0]) for v in verts if v[2] < 0.05]
    rep.near("the floor is the size asked for", 2.0 * max(floor_x), BOAT, TOL)

    # Its thickness, on the axis: where the material starts and stops.
    span = [z for z in [0.02 + 0.02 * i for i in range(400)] if solid(0.0, 0.0, z)]
    if span:
        rep.near("the floor is the thickness asked for",
                 max(span) - min(span), WALL, TOL)
        rep.near("the floor sits on the bed", min(span), 0.0, TOL)
        # The cavity above it, up to the rim. On the axis everything above the
        # floor is empty, so the cavity's depth is the part's height less the
        # floor rather than anything measured on the axis.
        rep.near("the cavity is the depth asked for",
                 size[2] - max(span), DEPTH - WALL, TOL)
    else:
        rep.check("the boat has a floor", False, "no material on the axis")

    # --- the walls ---------------------------------------------------------
    rep.near("the boat is as deep as asked for", size[2], DEPTH, TOL)
    rim_x = [abs(v[0]) for v in verts if v[2] > DEPTH - 0.3]
    rep.near("the rim is where the lean puts it",
             2.0 * max(rim_x), BOAT + 2.0 * FLARE, TOL)

    # Across the wall, at a height where the cavity is wide open. The material
    # spans from the cavity's surface out to the outside; the span along x is
    # the wall divided by the cosine of the lean.
    z_mid = 0.5 * DEPTH
    reach = int(100 * (0.5 * BOAT + FLARE + 2.0))
    scan = [x for x in [0.01 * i for i in range(reach)] if solid(x, 0.0, z_mid)]
    if scan:
        rep.near("the wall is the thickness asked for",
                 (max(scan) - min(scan)) * math.cos(LEAN), WALL, TOL)
    else:
        rep.check("the boat has a wall", False, "no material at mid height")

    # The wall leans out, not in: the outside at the rim is wider than at the
    # floor. A wall leaning in would be an overhang.
    rep.check("the walls lean outward", max(rim_x) > max(floor_x),
              f"rim half-width {max(rim_x):.1f} mm against the floor's "
              f"{max(floor_x):.1f} mm")

    # --- open at the top ---------------------------------------------------
    rep.check("the cavity is open at the top",
              not solid(0.0, 0.0, DEPTH + 0.5),
              "nothing found on the axis above the rim")

    # --- printing ----------------------------------------------------------
    bad = PR.overhangs(verts, tris, shells=shells, exposed_only=True)
    rep.check("no overhang needs support, at 45 deg from horizontal", not bad,
              f"{len(bad)} faces slope less than 45 deg from horizontal"
              if bad else
              f"the walls lean {math.degrees(LEAN):.1f} deg out, so the cavity's "
              f"surface is {90.0 - math.degrees(LEAN):.1f} deg from horizontal")
    bed = PR.bed_contact(verts, tris)
    rep.check("it stands on the bed", bed > 1000.0,
              f"{bed:.0f} mm^2 of contact with the bed")

    print(f"  part size         : {size[0]:.1f} x {size[1]:.1f} x {size[2]:.1f} mm")
    print(f"  vertices / tris   : {len(verts)} / {len(tris)}")
    return 0 if rep.finish() else 1


if __name__ == "__main__":
    sys.exit(main())

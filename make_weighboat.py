#!/usr/bin/env python3
"""Generate a weigh boat: a shallow flared dish for weighing powder into.

A square dish with a flat floor and walls that lean outward as they rise. The
lean is what makes it printable: the outside of a leaning wall faces up and out,
and the inside -- the surface of the cavity -- is steeper than 45 degrees from
horizontal for any lean under 45 degrees, so neither needs support. Printed flat
on its bed.

    BOAT = 74        floor size, mm, square
    DEPTH = 15       inside depth, mm
    FLARE = 6        how far each wall leans out over that depth, mm
    WALL = 1.2       wall and floor thickness, mm

The body is one solid: the outside is lofted from the floor's rectangle to a
wider one at the rim, the cavity is lofted the same way, inset by the wall
thickness, and subtracted. The cavity's loft runs past the rim so the cut is a
crossing at the top rather than two faces meeting in the same plane.

Usage: python3 make_weighboat.py [output.3mf]
"""

import math
import sys

import csglib as csg
import memcap
import printability as PR
import threelib as t

# --- Parameters ------------------------------------------------------------
BOAT = 44.0          # floor size (square), mm
DEPTH = 15.0         # inside depth, mm
FLARE = 6.0          # outward lean of each wall over that depth, mm
WALL = 0.8           # wall and floor thickness, mm
TOL = 0.01           # max chord deviation on round surfaces

PARAMS = dict(boat=BOAT, depth=DEPTH, flare=FLARE, wall=WALL, tol=TOL)


def default_out():
    return t.stamped_filename("weighboat", PARAMS)


def half_at(z, **p):
    """Half-width of the *outside* at height z, and of the cavity inside it.

    The inset is the wall measured *across* it, not across the floor: a wall
    leaning over by `lean` has to be inset by `wall / cos(lean)` horizontally to
    come out `wall` thick. Insetting by `wall` alone leaves a 1.11 mm wall where
    1.2 mm was asked for.
    """
    out = 0.5 * p["boat"] + p["flare"] * z / p["depth"]
    return out, out - p["wall"] / math.cos(lean(p))


def lean(p):
    """How far the walls lean from vertical, in radians."""
    return math.atan2(p["flare"], p["depth"])


def section(half, z):
    """A square ring at height z, wound counter-clockwise seen from +Z."""
    return [(half, -half, z), (half, half, z), (-half, half, z), (-half, -half, z)]


def build(**p):
    if p["wall"] <= 0.0:
        raise ValueError("a boat with no wall is not a boat")
    if p["flare"] < 0.0:
        raise ValueError(f"a {p['flare']:g} mm flare leans inward, which is an "
                         f"overhang; a weigh boat's walls lean out")

    outer = t.Mesh("outer")
    t.loft(outer, [section(half_at(0.0, **p)[0], 0.0),
                   section(half_at(p["depth"], **p)[0], p["depth"])])

    # The cavity, inset by the wall at every height so the wall is uniform, and
    # run past the rim so the subtraction crosses the top edge.
    inner = t.Mesh("inner")
    t.loft(inner, [section(half_at(p["wall"], **p)[1], p["wall"]),
                   section(half_at(p["depth"] + 2.0, **p)[1],
                           p["depth"] + 2.0)])

    mesh = csg.subtract(outer, inner, name="weighboat")
    return mesh, {"shells": len(mesh.shells()),
                  "cracks": csg.crack_length(mesh.verts, mesh.tris),
                  "gaps": len(csg.real_gaps(mesh.verts, mesh.tris)),
                  "bbox": mesh.bbox(),
                  "bed": PR.bed_contact(mesh.verts, mesh.tris)}


def main():
    memcap.cap()
    out = sys.argv[1] if len(sys.argv) > 1 else default_out()
    mesh, s = build(**PARAMS)
    desc = (f"Weigh boat: {PARAMS['boat']:g} mm square floor, "
            f"{PARAMS['depth']:g} mm deep, walls leaning {PARAMS['flare']:g} mm "
            f"out, {PARAMS['wall']:g} mm wall.")
    t.write_3mf(out, (mesh.name, mesh), title="Weigh boat",
                description=desc, params=PARAMS)

    bb = s["bbox"]
    top = PARAMS["boat"] + 2.0 * PARAMS["flare"]
    print(f"wrote {out}")
    print(f"  params            : {t.describe(PARAMS)}")
    print(f"  param hash        : {t.params_hash(PARAMS)}")
    print(f"  floor             : {PARAMS['boat']:g} mm square, "
          f"{PARAMS['wall']:g} mm thick")
    print(f"  rim               : {top:g} mm square at {PARAMS['depth']:g} mm up, "
          f"so each wall leans {PARAMS['flare']:g} mm out "
          f"({_lean(PARAMS):.1f} deg from vertical)")
    print(f"  size              : {bb[3] - bb[0]:.1f} x {bb[4] - bb[1]:.1f} x "
          f"{bb[5] - bb[2]:.1f} mm")
    print(f"  volume            : {mesh.signed_volume():.1f} mm^3 of material")
    print(f"  capacity          : {_capacity(PARAMS):.1f} ml to the rim")
    print(f"  bed contact       : {s['bed']:.1f} mm^2")
    print(f"  print orientation : flat on its floor. The walls lean out, so no "
          f"surface needs support")
    print(f"  vertices / tris   : {len(mesh.verts)} / {len(mesh.tris)}")
    print(f"  shells            : {s['shells']}")
    print(f"  gaps              : {s['gaps']} real, {s['cracks']:.2e} mm of "
          f"edge used once")


def _lean(p):
    return math.degrees(lean(p))


def _capacity(p):
    """Volume of the cavity, as a frustum, in millilitres."""
    lo = p["boat"]
    hi = p["boat"] + 2.0 * p["flare"]
    h = p["depth"] - p["wall"]
    # Square frustum: V = h/3 * (A1 + A2 + sqrt(A1*A2))
    a1, a2 = lo * lo, hi * hi
    return h / 3.0 * (a1 + a2 + (a1 * a2) ** 0.5) / 1000.0


if __name__ == "__main__":
    main()

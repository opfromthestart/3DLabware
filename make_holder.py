#!/usr/bin/env python3
"""Generate the filter holder / funnel cradle as a 3MF.

A standalone part: a sheet with a central hole and a conical cradle rising
from it. It does not import the filter and derives nothing from it -- the two
parts are designed separately, and `verify_parts.py` is what checks they
agree. Nothing here changes when the filter changes; a mismatch shows up in
the verifier, which is where a mismatch belongs.

    z = 0.0 .. 0.4     sheet, 90 mm disc, 18 mm central hole
    z = 0.4 .. 10.4    conical cradle, bore r 11.5 -> 38.0

Two jobs, done by two different features:

* The sheet's *top face* carries the filter. The hole is deliberately
  smaller than the filter's outlet outer radius, so the flat annulus at the
  bottom of the filter lands on the sheet instead of dropping through. It is
  still wider than the filter's bore, so it does not restrict the flow. A
  positive seat has to work this way round: the filter's outside is one
  continuous taper with no horizontal face anywhere, so there is nothing
  further out for a rim to press against.
* The cradle's *bore* keeps the filter concentric. It is a straight cone
  `clearance` larger than the filter's funnel, which is also a straight cone
  once the spout is gone, so the gap stays constant over the whole height.

For that gap to stay constant the tapers must match, which is a relationship
between the two parts rather than a property of either:
`(cradle_r1 - cradle_r0) / cradle_h` must equal the funnel's taper.

Usage: python3 make_holder.py [output.3mf]
"""

import math
import sys

import csglib as csg
import threelib as t

# --- Parameters (mm) -------------------------------------------------------
PLATE_DIA = 90.0     # outer diameter of the sheet
PLATE_H = 0.4        # sheet thickness
HOLE_R = 9.0         # central hole; see the seating note above
CRADLE_R0 = 11.5     # cradle bore at the sheet's top face
CRADLE_R1 = 38.0     # cradle bore at the rim
CRADLE_H = 12.5      # cradle height above the sheet; must span the
                     # filter's funnel, so (cradle_r1 - cradle_r0) /
                     # cradle_h equals the funnel's own taper
WALL = 0.8           # radial thickness of the cradle wall
TOL = 0.01           # max chord deviation on round surfaces

# The filter this is shaped for. Recorded for the record only -- nothing below
# is computed from it. verify_parts.py checks these against the real filter.
FITS_FILTER = dict(outlet_r=10.0, bore_r=8.8, funnel_slope=2.65)

PARAMS = dict(plate_dia=PLATE_DIA, plate_h=PLATE_H, hole_r=HOLE_R,
              cradle_r0=CRADLE_R0, cradle_r1=CRADLE_R1, cradle_h=CRADLE_H,
              wall=WALL, tol=TOL)


def default_out():
    return t.stamped_filename("holder", PARAMS)


def build(plate_dia=PLATE_DIA, plate_h=PLATE_H, hole_r=HOLE_R,
          cradle_r0=CRADLE_R0, cradle_r1=CRADLE_R1, cradle_h=CRADLE_H,
          wall=WALL, tol=TOL):
    plate_r = plate_dia / 2.0
    n = t.auto_segments(plate_r, tol)

    plate = t.Mesh("plate")
    t.tube(plate, hole_r, plate_r, 0.0, plate_h, n)
    # Start the wall at z=0, below the sheet's top face, so the two solids
    # interpenetrate instead of merely touching on a shared plane -- and then
    # union them, because interpenetrating is the same thing as overlapping,
    # and a slicer that fills even-odd prints the overlap as a hole.
    z_top = plate_h + cradle_h
    # cradle_r0 is the bore at the sheet's *top face*, because that is where
    # the filter begins. The wall is built from z=0 so it interpenetrates the
    # sheet, so extrapolate the cone down rather than starting it at
    # cradle_r0: measuring from the underside instead adds a whole taper's
    # worth of gap before the filter starts, and the clearance grows from
    # 1.5 mm to 2.5 mm across the cradle.
    slope = (cradle_r1 - cradle_r0) / cradle_h
    r_base = cradle_r0 - slope * plate_h
    cradle = t.Mesh("cradle")
    t.tapered_tube(cradle, r_base, r_base + wall,
                   cradle_r1, cradle_r1 + wall, 0.0, z_top, n)
    mesh = csg.union(plate, cradle, name="holder")

    sagitta = plate_r * (1 - math.cos(math.pi / n))
    return mesh, {"n": n, "sagitta": sagitta, "hole_r": hole_r,
                  "shells": len(mesh.shells()),
                  "cracks": csg.crack_length(mesh.verts, mesh.tris),
                  "cradle_r0": cradle_r0, "cradle_r1": cradle_r1,
                  "z_top": z_top,
                  "seat_w": FITS_FILTER["outlet_r"] - hole_r,
                  "slope": (cradle_r1 - cradle_r0) / cradle_h}


def main():
    out = sys.argv[1] if len(sys.argv) > 1 else default_out()
    mesh, s = build(**PARAMS)

    desc = (f"Filter holder. {PARAMS['plate_dia']:g} mm sheet "
            f"{PARAMS['plate_h']:g} mm thick with a {HOLE_R * 2:g} mm hole, "
            f"and a {PARAMS['wall']:g} mm cradle wall rising "
            f"{s['z_top'] - PARAMS['plate_h']:g} mm, bore r "
            f"{s['cradle_r0']:g} -> {s['cradle_r1']:g}.")
    t.write_3mf(out, (mesh.name, mesh), title="Filter holder",
                description=desc, params=PARAMS)

    bb = mesh.bbox()
    print(f"wrote {out}")
    print(f"  params            : {t.describe(PARAMS)}")
    print(f"  param hash        : {t.params_hash(PARAMS)}")
    print(f"  sheet             : {PARAMS['plate_dia']:g} mm dia x "
          f"{PARAMS['plate_h']:g} mm, hole {s['hole_r'] * 2:g} mm")
    print(f"  cradle            : bore r {s['cradle_r0']:g} -> {s['cradle_r1']:g} "
          f"over z {PARAMS['plate_h']:g}..{s['z_top']:g} "
          f"(taper {s['slope']:g} mm/mm)")
    print(f"  carries filter on : the sheet's top face, r {s['hole_r']:g}.."
          f"{FITS_FILTER['outlet_r']:g} -> {s['seat_w']:g} mm of seat "
          f"(for a filter with a {FITS_FILTER['outlet_r']:g} mm outlet)")
    print(f"  round facets      : {s['n']} segments, chord {s['sagitta']:.5f} mm")
    print(f"  bounding box      : {bb[3] - bb[0]:g} x {bb[4] - bb[1]:g}"
          f" x {bb[5] - bb[2]:g} mm")
    print(f"  vertices / tris   : {len(mesh.verts)} / {len(mesh.tris)}")
    print(f"  volume            : {mesh.signed_volume():.2f} mm^3")
    print(f"  shells            : {len(mesh.shells())}")


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
"""Generate a stopper: a solid ground male joint with a flat disc handle.

A ground joint seals by friction on a taper, so a stopper is the same cone the
rest of this project's joints are, with nothing bored through it. The handle is
a flat disc standing on edge, its plane through the axis, like the head of a
glass stopper: the rim is a long lever for very little material, and the two
flat faces are vertical walls that print as themselves.

    JOINT = "14/25"   the ground joint it plugs
    DISC_R = 6.6      disc radius, mm (13.2 mm across)
    DISC_T = 5        disc thickness, mm
    DISC_DROP = 2.5   how far the disc's rim dips into the top of the cone, mm

The disc stays inside the cone the fitting surface implies. That cone is
`rho = r_mouth + TAPER (z - glen)` continued past the mouth, and a disc of
radius R centred at `glen + R - d` is inside it only while

    R <= (r_mouth - TAPER d) / (sqrt(1 + TAPER^2) - TAPER)

which for a 14 mm mouth is 6.7 mm at a 2.5 mm drop -- 13.4 mm across,
measured to the corners of the rim rather than to the rim, because those
corners are the furthest the material gets. A wider disc is a wider
lever, and it protrudes past the taper's cone -- at which point its underside
leaves the cone where the rim is nearly horizontal, and the part needs support
it was designed not to need. The build reports the limit and refuses to go
past it.

Staying inside the cone also settles the printing, for the same reason: the
rim crosses the cone's top face at `sqrt(2 R d - d^2)` from the axis and leaves
it rising at `atan(reach / (R - d))`, which is 50 degrees here against a
printer's 45. The cone itself is 87 degrees from the horizontal and the disc's
faces are vertical, so the whole part prints standing on its tip with nothing
under it.

It is solid on purpose, which has a consequence worth
stating: a solid stopper on a flask that is heated or that holds something
volatile will pressurise, and the pressure goes into the glass or into the
stopper until it pops. This is a plug for a cold vessel, or one vented
elsewhere.

Usage: python3 make_stopper.py [output.3mf]
"""

import math
import sys

import csglib as csg
import junctionlib as J
import printability as PR
import threelib as t

# --- Parameters ------------------------------------------------------------
JOINT = "15/25"      # the ground joint it plugs: "14/20", "19/26", "24/40"...
DISC_R = 6.6         # handle disc radius, mm
DISC_T = 5.0         # handle disc thickness, mm
DISC_DROP = 2.5      # how far the disc dips into the top of the ground cone
TOL = 0.01           # max chord deviation on round surfaces

PARAMS = dict(joint=JOINT, disc_r=DISC_R, disc_t=DISC_T, disc_drop=DISC_DROP,
              tol=TOL)


def default_out():
    return t.stamped_filename("stopper", PARAMS)


def dims(p):
    """Everything derived from the parameter set, in one place."""
    r_mouth, r_tip, glen = J.joint_dims(p["joint"])
    r, d = p["disc_r"], p["disc_drop"]
    # Half the width of the disc where it crosses the cone's top face. The
    # disc's centre is above that face, so this is its widest point at or below
    # it, and therefore the whole test.
    reach = math.sqrt(max(0.0, 2.0 * r * d - d * d))
    # The fitting surface is a cone of slope TAPER on the radius, and continued
    # past the mouth it is the envelope the handle has to live inside. The
    # disc's furthest material is not its rim but the corners where the rim
    # meets the flat faces, at `hypot(R, T/2)` from the axis, and they sit at
    # the disc's centre height -- so that is where the envelope is tested.
    root = math.sqrt(1.0 + J.TAPER ** 2)
    return {"r_mouth": r_mouth, "r_tip": r_tip, "glen": glen, "reach": reach,
            "rim": math.hypot(r, p["disc_t"] / 2.0),
            "cone_at_rim": (r_mouth + J.TAPER * (r - d)) / root,
            "r_max": max_radius(p),
            "z_centre": glen + r - d, "z_top": glen + r - d + r}


def max_radius(p, lo=0.5, hi=20.0):
    """The largest disc that stays inside the cone the fitting surface implies.

    `hypot(R, T/2) <= (r_mouth + TAPER (R - d)) / sqrt(1 + TAPER^2)`, solved by
    bisection: the inequality is monotone in R on this side of the tangent, and
    it is not worth an algebra slip to write it out.
    """
    d = p["disc_drop"]
    root = math.sqrt(1.0 + J.TAPER ** 2)
    r_mouth = J.joint_dims(p["joint"])[0]

    def inside(r):
        return (math.hypot(r, p["disc_t"] / 2.0) * root
                <= r_mouth + J.TAPER * (r - d))

    if not inside(lo):
        return lo
    while hi - lo > 1e-6:
        mid = 0.5 * (lo + hi)
        lo, hi = (mid, hi) if inside(mid) else (lo, mid)
    return 0.5 * (lo + hi)


def fused_area(p):
    """The area of the disc that ends up inside the cone -- a circular segment.

    This is the cross-section that carries the twist, so it is worth a number:
    a disc 7 mm in radius dropping 2.5 mm into a 7 mm mouth fuses over about
    35 mm^2, which is far more than a hand can shear off.
    """
    r, drop = p["disc_r"], p["disc_drop"]
    if drop <= 0.0:
        return 0.0
    drop = min(drop, r)
    return (r * r * math.acos(1.0 - drop / r)
            - (r - drop) * math.sqrt(max(0.0, 2.0 * r * drop - drop * drop)))


def body_profile(p):
    """The stopper's cross-section: the ground cone, flat at both ends."""
    d = dims(p)
    return [(0.0, 0.0), (d["r_tip"], 0.0), (d["r_mouth"], d["glen"]),
            (0.0, d["glen"])]


def _facet(k):
    """What fraction of a circle a k-gon is: the volume scaling of a revolve."""
    return k * math.sin(2.0 * math.pi / k) / (2.0 * math.pi)


def build(**p):
    d = dims(p)
    for key in ("joint", "disc_r", "disc_t", "disc_drop"):
        if p[key] is None:
            raise ValueError(f"{key} is not set")
    if p["disc_r"] <= 0.0 or p["disc_t"] <= 0.0:
        raise ValueError("a disc needs a radius and a thickness")
    if p["disc_drop"] < 0.0:
        raise ValueError("the disc cannot float above the cone")
    if p["disc_drop"] >= p["disc_r"]:
        raise ValueError("the disc would be buried past its own centre")
    if d["rim"] > d["cone_at_rim"] + 1e-9:
        raise ValueError(
            f"a {p['disc_r'] * 2:g} mm disc {p['disc_t']:g} mm thick reaches "
            f"{d['rim']:.2f} mm from the axis, and the cone the fitting "
            f"surface implies is only {d['cone_at_rim']:.2f} mm there: the "
            f"handle pokes out through the taper's cone by "
            f"{d['rim'] - d['cone_at_rim']:.2f} mm, which leaves the cone's "
            f"top face where the rim is nearly horizontal and needs support. "
            f"The largest disc that stays inside is {d['r_max'] * 2:.2f} mm "
            f"across at a {p['disc_drop']:g} mm drop")

    # The body: one solid of revolution, the ground cone with flat ends.
    n_body = t.auto_segments(d["r_mouth"], p["tol"])
    mesh = t.Mesh("body")
    t.revolve(mesh, body_profile(p), n_body)

    # The handle: a flat disc standing on edge, its plane through the axis. A
    # prism is built flat in XY and extruded up Z, so it is built that way and
    # then turned: rotate_x(90) sends the extrusion axis to -Y and the
    # polygon's own Y to +Z, which leaves the disc's plane the XZ plane -- the
    # one containing the axis -- with its thickness across Y.
    n_disc = t.auto_segments(p["disc_r"], p["tol"])
    disc = t.Mesh("disc")
    t.prism(disc, t.disc_polygon(p["disc_r"], n_disc),
            -p["disc_t"] / 2.0, p["disc_t"] / 2.0)
    t.rotate_x(disc, 90.0)
    t.translate(disc, 0.0, 0.0, d["z_centre"])

    # The disc dips into the cone, so the two surfaces cross. Union them: a
    # slicer that fills even-odd counts the crossing region twice and prints a
    # hole where the handle joins the body.
    mesh = csg.union(mesh, disc, name="stopper")

    return mesh, {"n_body": n_body, "n_disc": n_disc, "dims": d,
                  "shells": len(mesh.shells()),
                  "cracks": csg.crack_length(mesh.verts, mesh.tris),
                  "facet_body": _facet(n_body), "facet_disc": _facet(n_disc),
                  "fused": fused_area(p), "bbox": mesh.bbox(),
                  "bed": math.pi * d["r_tip"] ** 2,
                  "v_body": (math.pi * d["glen"] / 3.0
                             * (d["r_tip"] ** 2 + d["r_tip"] * d["r_mouth"]
                                + d["r_mouth"] ** 2)),
                  "v_disc": math.pi * p["disc_r"] ** 2 * p["disc_t"]}


def main():
    out = sys.argv[1] if len(sys.argv) > 1 else default_out()
    mesh, s = build(**PARAMS)
    d = s["dims"]

    desc = (f"Ground stopper, {PARAMS['joint']} male, solid, with a "
            f"{PARAMS['disc_r'] * 2:g} mm disc handle {PARAMS['disc_t']:g} mm "
            f"thick standing on edge.")
    t.write_3mf(out, (mesh.name, mesh), title="Ground joint stopper",
                description=desc, params=PARAMS)

    bb = s["bbox"]
    # The rim crosses the top face `reach` from the axis and is `R - drop`
    # below the disc's centre, so the slope it leaves at is the ratio of those.
    leave = math.degrees(math.atan2(d["reach"],
                                    PARAMS["disc_r"] - PARAMS["disc_drop"]))
    print(f"wrote {out}")
    print(f"  params            : {t.describe(PARAMS)}")
    print(f"  param hash        : {t.params_hash(PARAMS)}")
    print(f"  joint             : {PARAMS['joint']} male, tip "
          f"{d['r_tip'] * 2:g} mm, mouth {d['r_mouth'] * 2:g} mm, ground "
          f"{d['glen']:g} mm at 1:{1.0 / J.TAPER:g} on the diameter")
    print(f"  handle            : disc {PARAMS['disc_r'] * 2:g} mm across, "
          f"{PARAMS['disc_t']:g} mm thick, standing on edge, plane through the "
          f"axis")
    print(f"  fusion            : disc dips {PARAMS['disc_drop']:g} mm into "
          f"the cone, crossing the top face {d['reach'] * 2:.2f} mm wide, "
          f"fusing {s['fused']:.1f} mm^2")
    print(f"  inside the cone   : the disc's furthest corner is "
          f"{d['rim']:.2f} mm from the axis where the fitting surface's cone "
          f"is {d['cone_at_rim']:.2f} mm, so the handle stays inside it "
          f"({d['cone_at_rim'] - d['rim']:+.2f} mm clear); the largest disc "
          f"that would is {d['r_max'] * 2:.2f} mm across at this drop")
    print(f"  height            : {bb[5] - bb[2]:.1f} mm, widest "
          f"{max(bb[3] - bb[0], bb[4] - bb[1]):.1f} mm")
    print(f"  round facets      : {s['n_body']} on the cone, {s['n_disc']} on "
          f"the disc")
    print(f"  volume            : cone {s['v_body'] * s['facet_body']:.2f} + "
          f"disc {s['v_disc'] * s['facet_disc']:.2f} mm^3, less the "
          f"{s['fused'] * PARAMS['disc_t']:.1f} mm^3 they share")
    print(f"  bed contact       : {s['bed']:.1f} mm^2, the tip's disc")
    print(f"  print orientation : standing on the tip; the cone is 87 deg from "
          f"the horizontal, the disc's faces are vertical, and its edge leaves "
          f"the cone at {leave:.0f} deg, so it beats the "
          f"{PR.OVERHANG_DEG:g} deg the printer can hold")
    print(f"  vertices / tris   : {len(mesh.verts)} / {len(mesh.tris)}")
    print(f"  shells            : {s['shells']} -- the cone and the disc "
          f"unioned, so no surface lies inside another")
    print(f"  cracks            : {s['cracks']:.2e} mm of edge used once")
    print(f"  note              : solid, so a heated or volatile flask "
          f"pressurises behind it")


if __name__ == "__main__":
    main()

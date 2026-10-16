#!/usr/bin/env python3
"""Generate a fraction collection adapter: one inlet, three flask outlets.

Rotating fraction collection. The inlet is a female ground joint standing on
the axis of rotation, so the whole adapter turns on the male taper above the
condenser while the vacuum holds it seated. The outlets hang below it, one per
receiving flask, so turning the adapter routes the distillate to the next
flask -- which is how you take a fraction when the thermometer says the next
cut has started.

    IN_OUT_ANGLE = 155     inlet axis to an outlet axis
    OUT_OUT_ANGLE = 43     outlet axis to its neighbour's
    COUNT = 3              receiving flasks

Both angles are yours to set, and they decide the third thing between them:
how far around the inlet's axis neighbouring outlets sit. Two arms, each
IN_OUT_ANGLE from the same axis and `step` apart around it, meet at

    cos(OUT_OUT_ANGLE) = cos^2(IN_OUT_ANGLE) + sin^2(IN_OUT_ANGLE) * cos(step)

so the step is not a free parameter, it is whatever those two angles imply:

    step = acos( (cos(OUT_OUT_ANGLE) - cos^2(IN_OUT_ANGLE)) / sin^2(IN_OUT_ANGLE) )

Two angles do not always go together. With the inlet at 155 deg each outlet is
26 deg off straight down, and its neighbours can be anywhere from coincident to
50 deg apart; ask for 80 and there is no arrangement at all, and it is refused
rather than approximated.

What the step buys you is whether the outlets close a ring. When `COUNT * step`
comes to exactly 360 the outlets are the classic carousel, evenly spaced around
the axis, and one index step puts every outlet where its neighbour was -- three
flasks in a circle, a third of a turn each time. Anything less than that spreads
them into a fan: still a rotation of one another, still indexed by turning
through `step`, but the flasks then sit along an arc rather than a circle, with
a gap left over of `360 - COUNT * step` degrees. Both are usable; the printout
says which one you have built, and the sum, so you can see how near the ring it
is.

The step is also what sets the flask bench: neighbouring mouths stand
`2 * OUT_LEN * sin(IN_OUT_ANGLE) * sin(step/2)` apart, and that is the space a
flask has to fit in.

Everything the arms meet in is built the way `junctionlib` describes: the arms'
outer solids and a blend ball, minus the bores taken out last, so the passage
is one open cavity instead of tubes plugging each other. Each arm's wall stops
short of the junction by its trim, so that the inlet's socket and the outlets'
ground cones are the only things a joint ever sees.

Printing. The inlet points up and the outlets hang below, so every wall of the
passage prints clean: the worst of them is `IN_OUT_ANGLE - 90` degrees off
vertical, which at the defaults is 25. The down-facing surfaces are all on the
*outside* of the junction -- the crotch under and between the outlets, where a
branch always has one and where support is ordinary and comes off. This part
does not print support-free; it does print with nothing inside a passage.
`verify_collector.py` measures both.

Usage: python3 make_collector.py [output.3mf]
"""

import math
import sys

import csglib as csg
import junctionlib as J
import threelib as t

# --- Parameters ------------------------------------------------------------
INLET_JOINT, INLET_GENDER, INLET_LEN = "14/20", "female", 36.0
OUT_JOINT, OUT_GENDER, OUT_LEN = "14/20", "male", 70.0
IN_OUT_ANGLE = 120.0    # inlet axis to an outlet axis, degrees
OUT_OUT_ANGLE = 43.0    # outlet axis to its neighbour's, degrees
COUNT = 3               # receiving flasks
HANDLE = 0.0            # >0 adds a knob above the inlet, that far wide (mm)
WALL = 2.5              # tube wall thickness
TIP_WALL = 1.0          # wall at a male tip
FLARE = 0.0             # length over which a joint's boss tapers to the tube
TOL = 0.01              # max chord deviation on round surfaces

PARAMS = dict(inlet_joint=INLET_JOINT, inlet_gender=INLET_GENDER,
              inlet_arm=INLET_LEN, out_joint=OUT_JOINT, out_gender=OUT_GENDER,
              out_arm=OUT_LEN, in_out_angle=IN_OUT_ANGLE,
              out_out_angle=OUT_OUT_ANGLE, count=COUNT, handle=HANDLE,
              wall=WALL, tip_wall=TIP_WALL, flare=FLARE, tol=TOL)


def azimuth_step(in_out_angle, out_out_angle, count):
    """How far around the inlet's axis neighbouring outlets sit, in degrees.

    The two angles the part is asked for fix this one. Spherical law of
    cosines, on the triangle between the inlet's axis and two outlet axes.
    """
    count = int(round(count))
    if count < 2:
        raise ValueError(f"{count} outlets is not a fraction collector")
    if not 90.0 < in_out_angle < 180.0:
        raise ValueError(
            f"an outlet {in_out_angle:g} deg from the inlet is not hanging "
            f"below it; the angle has to be over 90 for the outlets to drain")
    si = math.sin(math.radians(in_out_angle))
    ci = math.cos(math.radians(in_out_angle))
    c = (math.cos(math.radians(out_out_angle)) - ci * ci) / (si * si)
    if not -1.0 <= c <= 1.0:
        widest = math.degrees(math.acos(max(-1.0, min(1.0,
                                                     2.0 * ci * ci - 1.0))))
        raise ValueError(
            f"no arrangement has outlets {out_out_angle:g} deg apart and "
            f"{in_out_angle:g} deg off the inlet: at that inlet angle they can "
            f"be anywhere from 0 up to {widest:g} deg apart")
    step = math.degrees(math.acos(c))
    if (count - 1) * step >= 360.0:
        raise ValueError(
            f"outlets {out_out_angle:g} deg apart need {step:g} deg each "
            f"around the axis, and {count} of them would run past a full turn "
            f"and overlap")
    return step


def ring_gap(in_out_angle, out_out_angle, count):
    """Degrees left over when the outlets do not close a ring.

    Zero means the outlets are evenly spaced around the axis: one index step
    puts every outlet where its neighbour was. Any more than that and they
    stand in a fan, with this much of the circle unused.
    """
    return 360.0 - int(round(count)) * azimuth_step(in_out_angle,
                                                    out_out_angle, count)


def outlet_dirs(in_out_angle, out_out_angle, count):
    """Unit vectors for the outlets, `step` apart about the inlet's axis.

    Straight out from the junction like every other arm here, so an outlet
    hanging below has an angle to +z of `in_out_angle`, i.e. more than 90.
    """
    count = int(round(count))
    th = math.radians(in_out_angle)
    step = math.radians(azimuth_step(in_out_angle, out_out_angle, count))
    out = []
    for k in range(count):
        ph = k * step
        out.append((math.sin(th) * math.cos(ph),
                    math.sin(th) * math.sin(ph),
                    math.cos(th)))
    return out


def spacing(out_len, in_out_angle, out_out_angle, count):
    """How far apart neighbouring mouths stand, and so the flask that fits."""
    r = out_len * math.sin(math.radians(in_out_angle))
    return 2.0 * r * math.sin(
        math.radians(azimuth_step(in_out_angle, out_out_angle, count)) / 2.0)


def default_out():
    return t.stamped_filename("collector", PARAMS)


def build(**p):
    wall, tip_wall, flare = p["wall"], p["tip_wall"], p["flare"]
    count = p["count"]
    dirs = outlet_dirs(p["in_out_angle"], p["out_out_angle"], count)
    # The inlet stands on +z: that axis is what the adapter turns about, and
    # what the outlets are spaced around, so it is the part's own axis too.
    arms = [J.arm("in", p["inlet_joint"], p["inlet_gender"], p["inlet_arm"],
                  (0.0, 0.0, 1.0), wall, tip_wall)]
    for k, vec in enumerate(dirs):
        arms.append(J.arm(f"out{k + 1}", p["out_joint"], p["out_gender"],
                          p["out_arm"], vec, wall, tip_wall))

    bore_j = max(a["bore"] for a in arms)
    for a in arms:
        if a["glen"] + flare > a["len"]:
            raise ValueError(f"arm {a['tag']} is {a['len']:g} mm long but its "
                             f"{a['spec']} joint needs {a['glen'] + flare:g} mm")
        if a["bore"] <= 0.0:
            raise ValueError(f"tip wall {tip_wall:g} mm leaves no bore at a "
                             f"{a['spec']} joint")

    cuts, r_cap, tube_r, n = J.plan(arms, wall, flare, bore_j, p["tol"])
    J.check_arms(arms, cuts, r_cap, wall, flare, tube_r, bore_j)

    # The arms and the junction are built into one mesh, as separate shells
    # that interpenetrate: measured, up to 2.8 mm deep at the outlets.
    #
    # They were unioned and it did not finish: at this part's facet count --
    # about eleven thousand triangles across four arms and the junction -- a
    # BSP boolean was still running after fifty minutes, and with better
    # splitting-plane selection it ran the machine out of memory instead. Every
    # facet gets split by every plane of the other solid and the count
    # multiplies; the same union on the holder, at a fifth of the size, takes
    # seconds. The junction's blend is *meant* to overlap the arms, so this
    # part keeps its overlaps until the boolean is replaced with something that
    # does not explode, and a slicer that fills by winding unions it.
    # One mesh per arm, unioned one at a time with the junction, which is what
    # `csglib` requires: a mesh of all the arms at once intersects itself and a
    # boolean over it returns a wrong solid. Measured: the result is one shell
    # with no overlaps and no gaps -- 0 of each, against 4 and 6 overlapping
    # pairs up to 3 mm deep before -- in about fifteen seconds, and the mouths
    # come through untouched (material at 5.00..6.05 mm at the collector's mouth,
    # identical before and after).
    #
    # The verifier still reports a missing mouth annulus: a boolean splits that
    # face into a few hundred coplanar triangles, and `checks.annuli`'s ring
    # test -- the largest gap between consecutive vertex radii -- depends on how
    # a face was triangulated. The ring test needs fixing; the geometry does not.
    arm_meshes = []
    for a, cut in zip(arms, cuts):
        one = t.Mesh("arm")
        J.sweep_arm(one, a, cut, n, wall, flare, tube_r, bore_j)
        arm_meshes.append(one)

    # The junction. All four bores open into one cavity -- that is the whole
    # point of the part -- so it is meshed from the bores outwards, never
    # filled in.
    seg = n
    rings = max(6, n // 4)
    junction_mesh = t.Mesh("junction")
    jv, jf = J.junction(junction_mesh, arms, seg, rings, r_cap, tube_r, wall,
                        flare, bore_j)
    import time as _time
    _t0 = _time.time()
    mesh = csg.union(*arm_meshes, junction_mesh, name="collector")
    print(f"  union took {_time.time() - _t0:.1f} s")
    if p["handle"] > 0.0:
        t.torus(mesh, p["handle"], 1.5, seg=n, sides=8, cz=p["inlet_arm"] + 3.0)

    sagitta = tube_r * (1.0 - math.cos(math.pi / n))
    return mesh, {"n": n, "sagitta": sagitta, "arms": arms, "tube_r": tube_r,
                  "bore_junction": bore_j, "bbox": mesh.bbox(),
                  "cuts": cuts, "clip": r_cap, "junction_verts": jv,
                  "junction_faces": jf, "segments": seg, "rings": rings,
                  "spacing": spacing(p["out_arm"], p["in_out_angle"],
                                     p["out_out_angle"], count),
                  "step": azimuth_step(p["in_out_angle"], p["out_out_angle"],
                                       count),
                  "ring_gap": ring_gap(p["in_out_angle"], p["out_out_angle"],
                                       count)}


def main():
    out = sys.argv[1] if len(sys.argv) > 1 else default_out()
    mesh, s = build(**PARAMS)

    desc = (f"Fraction collector. Inlet {PARAMS['inlet_joint']} "
            f"{PARAMS['inlet_gender']} on the rotation axis; "
            + f"{PARAMS['count']} outlets {PARAMS['out_joint']} "
              f"{PARAMS['out_gender']} {PARAMS['out_out_angle']:g} deg apart, "
              f"{PARAMS['in_out_angle']:g} deg off the inlet, "
              f"{PARAMS['out_arm']:g} mm, mouths {s['spacing']:.1f} mm apart. "
            + f"{PARAMS['wall']:g} mm wall, {s['bore_junction'] * 2:g} mm bore.")
    t.write_3mf(out, (mesh.name, mesh), title="Fraction collector",
                description=desc, params=PARAMS)

    bb = s["bbox"]
    print(f"wrote {out}")
    print(f"  params            : {t.describe(PARAMS)}")
    print(f"  param hash        : {t.params_hash(PARAMS)}")
    for a, cut in zip(s["arms"], s["cuts"]):
        free, base = ((a["r_tip"], a["r_mouth"]) if a["gender"] == "male"
                      else (a["r_mouth"], a["r_tip"]))
        d = a["dir"]
        print(f"  arm {a['tag']:5}       : {a['gender']:6} {a['spec']:6} "
              f"free end {free * 2:g} -> {a['glen']:g} mm in {base * 2:g} mm, "
              f"{a['len']:g} mm to the junction")
        print(f"            dir    : ({d[0]:7.4f}, {d[1]:7.4f}, {d[2]:7.4f})  "
              f"bore {a['bore'] * 2:g} mm, wall stops {cut:g} mm short")
    axis = s["arms"][0]["dir"]
    angles = [math.degrees(math.acos(max(-1, min(1, sum(
        axis[k] * a["dir"][k] for k in range(3)))))) for a in s["arms"][1:]]
    print(f"  inlet to outlets  : " + "/".join(f"{v:.3f}" for v in angles)
          + f" deg (asked for {PARAMS['in_out_angle']:g})")
    out_dirs = [a["dir"] for a in s["arms"][1:]]
    pairs = [math.degrees(math.acos(max(-1, min(1, sum(
        out_dirs[i][k] * out_dirs[j][k] for k in range(3))))))
        for i, j in ((0, 1), (1, 2))]
    print(f"  outlet to outlet  : " + "/".join(f"{v:.3f}" for v in pairs)
          + f" deg between neighbours (asked for "
            f"{PARAMS['out_out_angle']:g})")
    even = 360.0 / int(round(PARAMS["count"]))
    off = s["step"] - even
    closes = abs(off) < 1e-6
    print(f"  around the axis   : {s['step']:.3f} deg per step"
          + (f", the ring closes exactly (even spacing is {even:.3f})"
             if closes else
             f", {abs(off) * int(round(PARAMS['count'])):.3f} deg "
             f"{'over' if off > 0 else 'short of'} a full turn at "
             f"{even:.3f} deg even, so a fan"))
    if not closes:
        ci = math.cos(math.radians(PARAMS["in_out_angle"]))
        si = math.sin(math.radians(PARAMS["in_out_angle"]))
        ring = math.degrees(math.acos(max(-1.0, min(1.0,
            ci * ci + si * si * math.cos(math.radians(even))))))
        # Invert the same relation the other way: what inlet angle closes the
        # ring at the outlet angle asked for? (The acos gives the tilt off
        # straight down, so it comes back off 180.)
        tilt = math.degrees(math.acos(math.sqrt(max(0.0, (
            math.cos(math.radians(PARAMS["out_out_angle"])) + 0.5) / 1.5))))
        print(f"  for a closed ring : OUT_OUT_ANGLE = {ring:.3f} at this "
              f"IN_OUT_ANGLE, or IN_OUT_ANGLE = {180.0 - tilt:.3f} at this "
              f"OUT_OUT_ANGLE")
    print(f"  outlet spacing    : {s['spacing']:.1f} mm between neighbouring "
          f"mouths (what sets the largest flask)")
    print(f"  tube              : {s['tube_r'] * 2:g} mm OD, "
          f"junction ball r {s['tube_r']:g} mm, bore {s['bore_junction'] * 2:g} mm")
    print(f"  junction          : open -- the four bores meet, junction piece "
          f"{s['junction_verts']} verts / {s['junction_faces']} faces "
          f"out to r {s['clip']:g} mm")
    print(f"  round facets      : {s['n']} segments, chord {s['sagitta']:.5f} mm"
          f"; junction {s['segments']}x{s['rings']} radial grid")
    print(f"  bounding box      : {bb[3] - bb[0]:.1f} x {bb[4] - bb[1]:.1f}"
          f" x {bb[5] - bb[2]:.1f} mm")
    print(f"  vertices / tris   : {len(mesh.verts)} / {len(mesh.tris)}")
    print(f"  volume            : {mesh.signed_volume():.2f} mm^3 "
          f"(shells summed; they overlap at the junction)")
    print(f"  shells            : {len(mesh.shells())} "
          f"(the arms + the junction; they overlap, so it is one solid)")


if __name__ == "__main__":
    main()

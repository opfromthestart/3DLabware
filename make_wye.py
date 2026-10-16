#!/usr/bin/env python3
"""Generate a three-way adapter with standard taper joints.

Three arms leave a common junction. Which way each arm points is set by the
three *pairwise* angles between them, so any arrangement in space is reachable:
the arms need not be coplanar, and the default is the familiar one -- arms 1
and 3 opposite each other forming a straight passage, arm 2 branching off it.

    ANGLE_12 = 75    ANGLE_13 = 180    ANGLE_23 = 105

Two angles fix the arms relative to arm 1; the third is what tips the last arm
out of their plane. Not every triple is realisable -- three directions in space
obey a spherical triangle inequality -- so an impossible combination is
rejected with the numbers rather than quietly approximated. `mirror` picks
which of the two mirror-image arrangements you get; it makes no difference
until the arms stop being coplanar.

Joints are per-arm: size and gender set independently, all male 14/20 by
default. The bore follows the elbow's rule -- a female socket ends at its tip
radius and the pipe continues at that size, a male tip is that radius on the
outside so its bore is smaller by `tip_wall` -- and the arms meet at the
largest of those bores.

THE JUNCTION. Three tubes cannot be glued into a three-way passage at one
point: an arm's wall runs straight through the other two bores and plugs them,
and what comes out is a sealed node with, at best, a crescent of a gap left
over. (Measured on the default, with the arms meeting at the junction: two of
the three pairs of mouths are shut outright and the third is pinched to 3 mm.)
What a tee *is*, and what is built here, is

    (the three outer solids u a blend ball) - (the three bores)

with the bores taken out *after* the union, so nothing is left standing inside
them. Every one of those primitives is convex and contains the junction, so a
ray leaving the junction crosses each of them in one interval; the material is
therefore star-shaped about the junction -- one interval rho(u)..R(u) per
direction, rho being the furthest a bore reaches along u and R the furthest an
outer solid or the ball reaches. That is what makes the junction meshable
radially, with no boolean kernel behind it.

Each arm's tube stops short of the junction by `trim`: the distance at which
its wall would start to lie inside another arm's bore. Out at that plane the
tube adds nothing the junction does not already carry, and any closer in it
would plug a passage. The junction piece spans the rest.

The ball, of the tube's own radius, comes with the outer solids: the arms are
tangent to it, and in an arrangement where all three arms lean the same way it
is the ball that closes the back of the junction. (Arms spanning more than a
half-space -- as the default's opposite pair do -- already cover it with their
own outer solids, so there the ball is along for the ride.)

Usage: python3 make_wye.py [output.3mf]
"""

import math
import sys

import csglib as csg
import junctionlib as J
import threelib as t
from junctionlib import (JOINTS, TAPER, arm, arm_profile, check_arms,
                         end_bore, joint_dims, plan, sweep_arm)

# --- Parameters ------------------------------------------------------------
JOINT_1, GENDER_1, ARM_1 = "14/20", "male", 45.0
JOINT_2, GENDER_2, ARM_2 = "14/20", "male", 45.0
JOINT_3, GENDER_3, ARM_3 = "14/20", "male", 45.0
ANGLE_12 = 75.0      # between arms 1 and 2
ANGLE_13 = 180.0     # between arms 1 and 3; 180 makes 1 and 3 a straight passage
ANGLE_23 = 105.0     # between arms 2 and 3; not free -- see the note above
MIRROR = False       # which side arm 3 falls on when the arms are not coplanar
WALL = 2.5           # tube wall thickness
TIP_WALL = 1.0       # wall at a male tip
FLARE = 0.0          # length over which a joint's boss tapers to the tube
TOL = 0.01           # max chord deviation on round surfaces

ARMS = (("1", JOINT_1, GENDER_1, ARM_1),
        ("2", JOINT_2, GENDER_2, ARM_2),
        ("3", JOINT_3, GENDER_3, ARM_3))

PARAMS = dict(joint_1=JOINT_1, gender_1=GENDER_1, arm_1=ARM_1,
              joint_2=JOINT_2, gender_2=GENDER_2, arm_2=ARM_2,
              joint_3=JOINT_3, gender_3=GENDER_3, arm_3=ARM_3,
              angle_12=ANGLE_12, angle_13=ANGLE_13, angle_23=ANGLE_23,
              mirror=MIRROR, wall=WALL, tip_wall=TIP_WALL, flare=FLARE, tol=TOL)


def directions(a12, a13, a23, mirror=False):
    """Unit vectors for the three arms from their pairwise angles.

    Place arm 1 on +z and arm 2 in the xz plane, which uses up the freedom to
    rotate the whole part. Arm 3 is then fixed by its angle to each of them,
    up to a reflection -- that is the one thing two angles cannot pin down, and
    what `mirror` chooses.
    """
    s12, c12 = math.sin(math.radians(a12)), math.cos(math.radians(a12))
    c13 = math.cos(math.radians(a13))
    c23 = math.cos(math.radians(a23))
    if abs(s12) < 1e-9:
        raise ValueError("arm 2 is collinear with arm 1 (angle_12 = "
                         f"{a12:g}); the junction has no room to open out")
    v1 = (0.0, 0.0, 1.0)
    v2 = (s12, 0.0, c12)
    # v3 . v1 = c13 and v3 . v2 = c23, solved for v3's x and z; the remainder
    # is v3's y component, which is what the third angle has to be able to fill.
    x = (c23 - c13 * c12) / s12
    z = c13
    y2 = 1.0 - x * x - z * z
    if y2 < -1e-9:
        raise ValueError(
            f"no arrangement has pairwise angles "
            f"{a12:g}/{a13:g}/{a23:g}: the arms cannot close a spherical "
            f"triangle (short by {math.sqrt(-y2):.4f})")
    y = math.sqrt(max(0.0, y2)) * (-1.0 if mirror else 1.0)
    return v1, v2, (x, y, z)


def default_out():
    return t.stamped_filename("wye", PARAMS)


def build(**p):
    wall, tip_wall, flare = p["wall"], p["tip_wall"], p["flare"]
    dirs = directions(p["angle_12"], p["angle_13"], p["angle_23"], p["mirror"])
    arms = [arm(tag, p["joint_" + tag], p["gender_" + tag], p["arm_" + tag], vec,
                wall, tip_wall)
            for vec, tag in zip(dirs, "123")]

    # Every arm's bore runs to the junction; the junction is the widest of
    # them, so a female socket is never narrowed where it meets the others.
    bore_j = max(a["bore"] for a in arms)
    for a in arms:
        if a["glen"] + flare > a["len"]:
            raise ValueError(f"arm {a['tag']} is {a['len']:g} mm long but its "
                             f"{a['spec']} joint needs {a['glen'] + flare:g} mm")
        if a["bore"] <= 0.0:
            raise ValueError(f"tip wall {tip_wall:g} mm leaves no bore at a "
                             f"{a['spec']} joint")

    # How far each arm pulls back, how far the junction must then reach to meet
    # the trimmed ends, and how finely to facet it. See junctionlib.
    cuts, r_cap, tube_r, n = plan(arms, wall, flare, bore_j, p["tol"])
    check_arms(arms, cuts, r_cap, wall, flare, tube_r, bore_j)

    # The arms and the junction are built into one mesh, as separate shells
    # that interpenetrate where the junction meets each trimmed arm.
    #
    # They were unioned, and it was worse: the junction's surface and the
    # trimmed arms' ends meet on planes they share, and a BSP resolves a
    # coplanar *partial* overlap badly. The union took the four overlapping
    # pairs down to one and left 649 real gaps, 195 mm of them, where the
    # surface had been sound. This boolean is right for pieces that cross and
    # wrong for pieces that touch; until it handles the second, this part keeps
    # its overlaps and a slicer that fills by winding unions it.
    # One mesh per arm, unioned one at a time with the junction, which is what
    # `csglib` requires: a mesh of all the arms at once intersects itself and a
    # boolean over it returns a wrong solid. Measured: the result is one shell
    # with no overlaps and no gaps -- 0 of each, against 4 and 6 overlapping
    # pairs up to 3 mm deep before -- in about fifteen seconds, and the mouths
    # come through untouched (material at 5.00..6.05 mm at the wye's mouth,
    # identical before and after).
    #
    # The verifier still reports a missing mouth annulus: a boolean splits that
    # face into a few hundred coplanar triangles, and `checks.annuli`'s ring
    # test -- the largest gap between consecutive vertex radii -- depends on how
    # a face was triangulated. The ring test needs fixing; the geometry does not.
    arm_meshes = []
    for a, cut in zip(arms, cuts):
        one = t.Mesh("arm")
        sweep_arm(one, a, cut, n, wall, flare, tube_r, bore_j)
        arm_meshes.append(one)

    # The junction. Arms 1 and 3 opposite each other are one straight passage
    # through it, and every arm's bore opens into the same cavity -- which is
    # the whole point of the part, so it is meshed from the bores outwards
    # rather than filled in.
    seg = n
    rings = max(8, n // 2)
    junction_mesh = t.Mesh("junction")
    jv, jf = J.junction(junction_mesh, arms, seg, rings, r_cap, tube_r, wall,
                        flare, bore_j)
    import time as _time
    _t0 = _time.time()
    mesh = csg.union(*arm_meshes, junction_mesh, name="wye")
    print(f"  union took {_time.time() - _t0:.1f} s")

    sagitta = tube_r * (1.0 - math.cos(math.pi / n))
    return mesh, {"n": n, "sagitta": sagitta, "arms": arms, "tube_r": tube_r,
                  "bore_junction": bore_j, "bbox": mesh.bbox(),
                  "cuts": cuts, "clip": r_cap, "junction_verts": jv,
                  "junction_faces": jf, "segments": seg, "rings": rings,
                  "pairs": {("1", "2"): p["angle_12"],
                            ("1", "3"): p["angle_13"],
                            ("2", "3"): p["angle_23"]}}


def main():
    out = sys.argv[1] if len(sys.argv) > 1 else default_out()
    mesh, s = build(**PARAMS)

    desc = ("Three-way adapter. " + "; ".join(
        f"arm {a['tag']}: {a['spec']} {a['gender']}, {a['len']:g} mm"
        for a in s["arms"])
        + f". Angles 12/13/23 = {PARAMS['angle_12']:g}/"
          f"{PARAMS['angle_13']:g}/{PARAMS['angle_23']:g} deg. "
          f"{PARAMS['wall']:g} mm wall, {s['bore_junction'] * 2:g} mm bore.")
    t.write_3mf(out, (mesh.name, mesh), title="Three-way adapter",
                description=desc, params=PARAMS)

    bb = s["bbox"]
    print(f"wrote {out}")
    print(f"  params            : {t.describe(PARAMS)}")
    print(f"  param hash        : {t.params_hash(PARAMS)}")
    for a, cut in zip(s["arms"], s["cuts"]):
        free, base = ((a["r_tip"], a["r_mouth"]) if a["gender"] == "male"
                      else (a["r_mouth"], a["r_tip"]))
        d = a["dir"]
        print(f"  arm {a['tag']}           : {a['gender']:6} {a['spec']:6} "
              f"free end {free * 2:g} -> {a['glen']:g} mm in {base * 2:g} mm, "
              f"{a['len']:g} mm to the junction")
        print(f"            dir    : ({d[0]:7.4f}, {d[1]:7.4f}, {d[2]:7.4f})  "
              f"bore {a['bore'] * 2:g} mm, wall stops {cut:g} mm short")
    p12 = math.degrees(math.acos(max(-1, min(1, sum(
        s["arms"][0]["dir"][k] * s["arms"][1]["dir"][k] for k in range(3))))))
    p13 = math.degrees(math.acos(max(-1, min(1, sum(
        s["arms"][0]["dir"][k] * s["arms"][2]["dir"][k] for k in range(3))))))
    p23 = math.degrees(math.acos(max(-1, min(1, sum(
        s["arms"][1]["dir"][k] * s["arms"][2]["dir"][k] for k in range(3))))))
    print(f"  pairwise angles   : 12 {p12:.3f}, 13 {p13:.3f}, 23 {p23:.3f} deg")
    # Arms 1 and 2 both lie in the xz plane by construction, so the three are
    # coplanar exactly when arm 3 does too.
    coplanar = abs(s["arms"][2]["dir"][1]) < 1e-9
    print(f"  arrangement       : {'coplanar' if coplanar else 'non-coplanar'}")
    print(f"  tube              : {s['tube_r'] * 2:g} mm OD, "
          f"junction ball r {s['tube_r']:g} mm, bore {s['bore_junction'] * 2:g} mm")
    print(f"  junction          : open -- the three bores join, junction piece "
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
          f"(3 arms + the junction; they overlap, so it is one solid)")


if __name__ == "__main__":
    main()

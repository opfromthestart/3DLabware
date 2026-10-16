#!/usr/bin/env python3
"""Generate a distillation elbow / bent adapter with standard taper joints.

A bent tube with a standard taper ground joint at each end, each independently
sized and independently male or female.

STANDARD TAPER JOINTS. A designation like `14/20` gives the mouth diameter in
millimetres and the length of the ground zone in millimetres. Every size shares
the same 1:10 taper, measured on the *diameter*, so the diameter falls by B/10
over the ground length whatever the size:

    14/20   14 mm at the mouth -> 12 mm at 20 mm depth
    24/40   24 mm at the mouth -> 20 mm at 40 mm depth
    45/50   45 mm at the mouth -> 40 mm at 50 mm depth

The wide end is at the mouth for both genders: a male joint is that cone as an
outer surface, a female joint is the same cone as a socket bore. They mate
because they are the same cone.

ANGLE. `angle` is the included angle between the two arm axes, measured
through the bend. A listing that says "75/105 degrees" is describing one elbow
two ways: an included 105 degrees leaves the second arm 75 degrees off the
line of the first. Both numbers are printed below so you can match whichever
the listing used.

SIZING. The tube is as fat as the smallest joint's tip allows, so a male
joint's ground cone runs straight into the tube with no step at its base. The
bore follows from the wall thickness. A female joint's socket is wider than
the tube by construction -- its mouth is the joint's diameter -- so its outer
wall flares back down to the tube over `flare` millimetres. Set `flare` to 0
for a hard shoulder, which is what glassware actually has and which prints as
an overhang.

Usage: python3 make_elbow.py [output.3mf]
"""

import math
import sys

import threelib as t

# Standard taper joints: mouth diameter / ground length, in mm.
JOINTS = {
    "5/13": (5.0, 13.0), "7/16": (7.0, 16.0), "10/18": (10.0, 18.0),
    "12/21": (12.0, 21.0), "14/20": (14.0, 20.0), "19/22": (19.0, 22.0),
    "24/40": (24.0, 40.0), "29/42": (29.0, 42.0), "34/45": (34.0, 45.0),
    "40/50": (40.0, 50.0), "45/50": (45.0, 50.0), "50/50": (50.0, 50.0),
    "55/50": (55.0, 50.0), "60/50": (60.0, 50.0),
}
TAPER = 1.0 / 20.0        # radius change per unit length: 1:10 on the diameter

# --- Parameters ------------------------------------------------------------
JOINT_A = "14/20"    # mouth diameter / ground length, or a (dia, len) pair
GENDER_A = "female"  # "male" or "female"; female end takes a male joint
JOINT_B = "14/20"
GENDER_B = "male"
ARM_A = 45.0        # straight length from joint A's mouth to the turn
ARM_B = 45.0        # straight length from the turn to joint B's mouth
ANGLE = 75.0        # included angle between the two arm axes
WALL = 2.5           # tube wall thickness
TIP_WALL = 1.0       # wall at a male joint's tip; the bore shrinks by
                     # 2*TIP_WALL between a female and a male end
BEND_R = 25.0        # centreline radius of the turn
FLARE = 0.0          # length over which a female socket's boss tapers to the tube
BEND_SEG = 24        # segments across the turn
TOL = 0.01           # max chord deviation on round surfaces


def parse_joint(spec):
    """'14/20' or (14, 20) -> (mouth diameter, ground length) in mm."""
    if isinstance(spec, str):
        dia, length = spec.split("/")
        return float(dia), float(length)
    return float(spec[0]), float(spec[1])


def joint_dims(spec):
    """(mouth radius, tip radius, ground length) for a joint designation."""
    dia, length = parse_joint(spec)
    r_mouth = dia / 2.0
    return r_mouth, r_mouth - length * TAPER, length


JA, LA = parse_joint(JOINT_A)
JB, LB = parse_joint(JOINT_B)
R_TIP_A = JA / 2.0 - LA * TAPER
R_TIP_B = JB / 2.0 - LB * TAPER
TUBE_R = min(R_TIP_A, R_TIP_B)      # the tube is as fat as the tightest tip allows
BORE_R = TUBE_R - WALL
DELTA = 180.0 - ANGLE               # how far the centreline actually turns

PARAMS = dict(joint_a=JOINT_A, gender_a=GENDER_A, joint_b=JOINT_B,
              gender_b=GENDER_B, arm_a=ARM_A, arm_b=ARM_B, angle=ANGLE,
              wall=WALL, tip_wall=TIP_WALL, bend_r=BEND_R, flare=FLARE,
              bend_seg=BEND_SEG, tol=TOL)


def default_out():
    return t.stamped_filename("elbow", PARAMS)


def _centreline(s, arm_a, bend_r, delta):
    """Point and unit tangent at arc length s along the bent centreline."""
    arc = bend_r * delta
    if s <= arm_a:
        return (0.0, 0.0, s), (0.0, 0.0, 1.0)
    if s <= arm_a + arc:
        phi = (s - arm_a) / bend_r
        return ((bend_r * (1.0 - math.cos(phi)), 0.0,
                 arm_a + bend_r * math.sin(phi)),
                (math.sin(phi), 0.0, math.cos(phi)))
    rem = s - arm_a - arc
    px = bend_r * (1.0 - math.cos(delta))
    pz = arm_a + bend_r * math.sin(delta)
    tan = (math.sin(delta), 0.0, math.cos(delta))
    return (px + rem * tan[0], 0.0, pz + rem * tan[2]), tan


def build(**p):
    r_mouth_a, r_tip_a, len_a = joint_dims(p["joint_a"])
    r_mouth_b, r_tip_b, len_b = joint_dims(p["joint_b"])
    wall = p["wall"]
    tip_wall = p["tip_wall"]
    fem_a = p["gender_a"] == "female"
    fem_b = p["gender_b"] == "female"

    # The bore is not one number. A female joint's socket ends at its tip
    # radius and the pipe has to carry straight on at that size, or there is a
    # shoulder for solids to sit on and for an incoming male joint to bottom
    # out against. A male joint's tip is that same radius on its *outside*, so
    # its bore has to be smaller or the tip comes to a knife edge. Letting the
    # bore shrink along the elbow lets both ends be right: it starts at the
    # socket's own diameter and loses 2*TIP_WALL across the length.
    bore_a = r_tip_a if fem_a else r_tip_a - tip_wall
    bore_b = r_tip_b if fem_b else r_tip_b - tip_wall
    if min(bore_a, bore_b) <= 0.0:
        raise ValueError(f"tip wall {tip_wall:g} mm leaves no bore at a "
                         f"{min(r_tip_a, r_tip_b) * 2:g} mm joint tip")

    # The plain tube has to be at least as wide as the fattest of the two
    # joints where they meet it, or the transition would step *inwards*.
    def inner_end_outer(spec, gender):
        r_mouth, r_tip, _ = joint_dims(spec)
        return r_tip + wall if gender == "female" else r_mouth

    tube_r = max(inner_end_outer(p["joint_a"], p["gender_a"]),
                 inner_end_outer(p["joint_b"], p["gender_b"]),
                 max(bore_a, bore_b) + wall)
    flare = p["flare"]
    if len_a + flare > p["arm_a"] or len_b + flare > p["arm_b"]:
        raise ValueError("arm is shorter than its joint plus flare")

    delta = math.radians(180.0 - p["angle"])
    arc = p["bend_r"] * delta
    total = p["arm_a"] + arc + p["arm_b"]
    n = t.auto_segments(tube_r + wall, p["tol"])

    # Stations, in order from joint A's mouth to joint B's mouth. Each entry is
    # (arc length, bore radius, outer radius). A female joint needs two entries
    # at its ground depth -- the socket ends there and the bore steps in -- and
    # the sweep lofts the flat shoulder between them.
    # The two ends taper in *opposite* directions, and that is the whole point.
    # A male joint's ground cone has its narrow tip at the free end and widens
    # to the joint's nominal diameter at its base; a female joint's socket has
    # the nominal diameter at its mouth and narrows going in. Put the male's
    # wide end at the free end and the two will not mate -- the cone is simply
    # the wrong way round.
    def a_out(s):
        if p["gender_a"] == "male":
            return r_tip_a + TAPER * s
        return r_mouth_a - TAPER * s + wall

    def b_out(s):
        depth = total - s
        if p["gender_b"] == "male":
            return r_tip_b + TAPER * depth
        return r_mouth_b - TAPER * depth + wall

    def bore_at(s):
        """Bore radius at arc length s: the socket cone, then the shrink."""
        if s <= len_a:
            return r_mouth_a - TAPER * s if fem_a else bore_a
        if s >= total - len_b:
            depth = s - (total - len_b)
            return r_tip_b + TAPER * depth if fem_b else bore_b
        # Linear between the two ends. At the female end this is tangent to the
        # socket cone, which is what removes the shoulder: the socket's last
        # millimetre and the pipe's first are the same diameter.
        u = (s - len_a) / ((total - len_b) - len_a)
        return bore_a + (bore_b - bore_a) * u

    st = [(0.0, bore_at(0.0), a_out(0.0)),
          (len_a, bore_at(len_a), a_out(len_a))]
    if flare > 0.0 and abs(a_out(len_a) - tube_r) > 1e-9:
        st.append((len_a + flare, bore_at(len_a + flare), tube_r))
    st.append((p["arm_a"], bore_at(p["arm_a"]), tube_r))
    for k in range(1, p["bend_seg"] + 1):
        sk = p["arm_a"] + arc * k / p["bend_seg"]
        st.append((sk, bore_at(sk), tube_r))
    s_taper = total - len_b
    if flare > 0.0 and abs(b_out(s_taper) - tube_r) > 1e-9:
        st.append((s_taper - flare, bore_at(s_taper - flare), tube_r))
    st.append((s_taper, bore_at(s_taper), b_out(s_taper)))
    st.append((total, bore_at(total), b_out(total)))
    st.sort(key=lambda e: e[0])

    # Dedupe coincident stations that carry the same radii (the bend's first
    # station lands on the arm's last one).
    clean = []
    for s, ri, ro in st:
        if clean and abs(s - clean[-1][0]) < 1e-12 \
                and abs(ri - clean[-1][1]) < 1e-12 and abs(ro - clean[-1][2]) < 1e-12:
            continue
        clean.append((s, ri, ro))

    stations = []
    for s, r_in, r_out in clean:
        pos, tan = _centreline(s, p["arm_a"], p["bend_r"], delta)
        stations.append((pos, tan, r_in, r_out))

    mesh = t.Mesh("elbow")
    t.sweep(mesh, stations, n)

    sagitta = tube_r * (1.0 - math.cos(math.pi / n))
    return mesh, {"n": n, "sagitta": sagitta, "stations": stations,
                  "arc": arc, "total": total, "delta": delta,
                  "tube_r": tube_r, "bore_r": bore_a,
                  "bore_a": bore_a, "bore_b": bore_b,
                  "r_mouth_a": r_mouth_a, "r_tip_a": r_tip_a, "len_a": len_a,
                  "r_mouth_b": r_mouth_b, "r_tip_b": r_tip_b, "len_b": len_b,
                  "bbox": mesh.bbox()}


def main():
    out = sys.argv[1] if len(sys.argv) > 1 else default_out()
    mesh, s = build(**PARAMS)

    desc = (f"Distillation elbow. {PARAMS['joint_a']} {PARAMS['gender_a']} "
            f"joint to {PARAMS['joint_b']} {PARAMS['gender_b']} joint, "
            f"{PARAMS['arm_a']:g} + {PARAMS['arm_b']:g} mm arms, "
            f"{PARAMS['angle']:g} deg included "
            f"({180.0 - PARAMS['angle']:g} deg off straight), "
            f"{PARAMS['bend_r']:g} mm bend, {PARAMS['wall']:g} mm wall, "
            f"{s['bore_r'] * 2:g} mm bore.")
    t.write_3mf(out, (mesh.name, mesh), title="Distillation elbow",
                description=desc, params=PARAMS)

    bb = s["bbox"]
    print(f"wrote {out}")
    print(f"  params            : {t.describe(PARAMS)}")
    print(f"  param hash        : {t.params_hash(PARAMS)}")
    for tag, spec, gender, rm, rt, ln in (
            ("A", PARAMS["joint_a"], PARAMS["gender_a"],
             s["r_mouth_a"], s["r_tip_a"], s["len_a"]),
            ("B", PARAMS["joint_b"], PARAMS["gender_b"],
             s["r_mouth_b"], s["r_tip_b"], s["len_b"])):
        ratio = ln / ((rm - rt) * 2.0) if rm > rt else float("inf")
        # A male's free end is its narrow tip; a female's is its wide mouth.
        free, base = (rt, rm) if gender == "male" else (rm, rt)
        print(f"  joint {tag}           : {gender:6} {spec:6} "
              f"free end {free * 2:g} mm -> {ln:g} mm in {base * 2:g} mm "
              f"(1:{ratio:g} taper on the diameter)")
    print(f"  angle             : included {PARAMS['angle']:g} deg = "
          f"{180.0 - PARAMS['angle']:g} deg off straight "
          f"(centreline turns {math.degrees(s['delta']):g} deg)")
    print(f"  arms              : {PARAMS['arm_a']:g} + {PARAMS['arm_b']:g} mm "
          f"straight, {s['arc']:.2f} mm of arc, {s['total']:.2f} mm centreline")
    print(f"  tube              : {s['tube_r'] * 2:g} mm OD, "
          f"{PARAMS['wall']:g} mm wall")
    print(f"  bore              : {s['bore_a'] * 2:g} mm at the "
          f"{PARAMS['gender_a']} end -> {s['bore_b'] * 2:g} mm at the "
          f"{PARAMS['gender_b']} end ({abs(s['bore_a'] - s['bore_b']) * 2:g} mm)")
    print(f"  bend radius       : {PARAMS['bend_r']:g} mm centreline")
    print(f"  round facets      : {s['n']} segments, chord {s['sagitta']:.5f} mm")
    print(f"  bounding box      : {bb[3] - bb[0]:.1f} x {bb[4] - bb[1]:.1f}"
          f" x {bb[5] - bb[2]:.1f} mm")
    print(f"  vertices / tris   : {len(mesh.verts)} / {len(mesh.tris)}")
    print(f"  volume            : {mesh.signed_volume():.2f} mm^3")
    print(f"  shells            : {len(mesh.shells())}")


if __name__ == "__main__":
    main()

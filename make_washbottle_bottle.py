#!/usr/bin/env python3
"""Gas-washing bottle: the vessel half of a two-part wash bottle.

Two parts mate through a standard taper ground joint: the bottle has the
female (socket) neck, the head (make_head.py) is the male plug that carries
the two tubes. The joint is 24/40 -- a 24 mm mouth on a 40 mm ground length,
at 1:20 on the diameter -- so the empty bottle seals against the head well
enough to force gas through the dip tube.

    BODY_R   30 mm  (60 mm OD cylinder)
    BODY_H   95 mm  cylindrical body height
    WALL      2 mm  wall thickness (constant, inside follows outside)
    BOTTOM    2 mm  base thickness
    JOINT  24/40    female neck (socket) with 0.15 mm print clearance per radius
    TOL     0.01    chord deviation on round surfaces

Usage: python3 make_bottle.py [output.3mf]
"""
import math
import sys

import csglib as csg
import junctionlib as J
import threelib as t

# --- Parameters ------------------------------------------------------------
BODY_R = 30.0
BODY_H = 95.0
WALL = 2.0
BOTTOM = 2.0
SHOULDER_H = 20.0
JOINT = "24/40"
CLEARANCE = 0.15
TOL = 0.01

PARAMS = dict(body_r=BODY_R, body_h=BODY_H, wall=WALL, bottom=BOTTOM,
              shoulder_h=SHOULDER_H, joint=JOINT, clearance=CLEARANCE, tol=TOL)

def default_out():
    return t.stamped_filename("washbottle_bottle", PARAMS)

def _joint():
    r_mouth, r_tip, glen = J.joint_dims(JOINT)
    return r_mouth, r_tip, glen

def outer_profile():
    r_mouth, r_tip, glen = _joint()
    # exterior: body cylinder, sloped shoulder, then the socket's outer wall
    # (which follows the female cone + wall thickness).
    s_bottom = r_tip + CLEARANCE          # socket interior radius at its deep end
    s_top = r_mouth + CLEARANCE           # at the opening
    neck_out_bottom = s_bottom + WALL
    neck_out_top = s_top + WALL
    shoulder_top_z = BODY_H + SHOULDER_H
    top_z = shoulder_top_z + glen
    return [(0.0, 0.0),
            (BODY_R, 0.0),
            (BODY_R, BODY_H),
            (neck_out_bottom, shoulder_top_z),
            (neck_out_top, top_z),
            (0.0, top_z)]

def inner_profile():
    r_mouth, r_tip, glen = _joint()
    s_bottom = r_tip + CLEARANCE
    s_top = r_mouth + CLEARANCE
    body_in = BODY_R - WALL
    shoulder_top_z = BODY_H + SHOULDER_H
    top_z = shoulder_top_z + glen
    # Interior solid: closed at the base (2 mm thick), open at the top by
    # ending the cavity beyond the outer solid so its cap sits outside the part.
    return [(0.0, BOTTOM),
            (body_in, BOTTOM),
            (body_in, BODY_H),
            (s_bottom, shoulder_top_z),
            (s_top, top_z),
            (0.0, top_z + 5.0)]

def build(**p):
    r_mouth, r_tip, glen = _joint()
    n = t.auto_segments(BODY_R, TOL)
    outer = t.Mesh("outer")
    t.revolve(outer, outer_profile(), n)
    inner = t.Mesh("inner")
    t.revolve(inner, inner_profile(), n)
    mesh = csg.subtract(outer, inner, name="bottle")
    s_top = r_mouth + CLEARANCE
    return mesh, {"n": n, "top_z": outer_profile()[-2][1],
                  "socket_bottom_r": r_tip + CLEARANCE,
                  "socket_top_r": s_top, "r_mouth": r_mouth, "r_tip": r_tip,
                  "glen": glen, "bbox": mesh.bbox()}

def main():
    out = sys.argv[1] if len(sys.argv) > 1 else default_out()
    mesh, s = build(**PARAMS)
    d = s
    bb = s["bbox"]
    desc = ("Gas washing bottle (250 mL class): 60 mm OD vessel with a 24/40 "
            "female ground neck and a 2 mm wall.")
    t.write_3mf(out, (mesh.name, mesh), title="Gas washing bottle",
                description=desc, params=PARAMS)
    print(f"wrote {out}")
    print(f"  params            : {t.describe(PARAMS)}")
    print(f"  param hash        : {t.params_hash(PARAMS)}")
    print(f"  body              : {BODY_R * 2:g} mm OD x {BODY_H:g} mm, "
          f"wall {WALL:g} mm, base {BOTTOM:g} mm")
    print(f"  neck (female)     : {JOINT} socket, inner r {s['socket_bottom_r']:.2f}"
          f"..{s['socket_top_r']:.2f} mm over {s['glen']:g} mm")
    print(f"  height            : {bb[5] - bb[2]:.1f} mm")
    print(f"  bounding box      : {bb[3]-bb[0]:g} x {bb[4]-bb[1]:g} "
          f"x {bb[5]-bb[2]:g} mm")
    print(f"  volume (hollow)   : {mesh.signed_volume():.1f} mm^3")
    print(f"  vertices / tris   : {len(mesh.verts)} / {len(mesh.tris)}")
    print(f"  shells            : {len(mesh.shells())}")

if __name__ == "__main__":
    main()

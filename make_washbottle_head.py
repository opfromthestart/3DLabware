#!/usr/bin/env python3
"""Gas-washing bottle: the head (stopper + tubes) half.

Two parts mate through a standard taper ground joint: the bottle
(make_bottle.py) is the female (socket) neck, this head is the male plug.
A 24/40 joint -- a 24 mm mouth on a 40 mm ground, 1:20 on the diameter --
gives an airtight-ish seal by contact, forcing the gas down the dip tube
and up through the liquid.

The head carries two vertical passages:
  * a dip tube (axis at x = +DIP_X) that runs from the top connector down
    through the head and out below the male tip, to near the bottle bottom;
  * a vent hole (axis at x = -VENT_X) whose lower end is flush with the
    underside of the head: it has no tube sticking down into the liquid.

Both are drilled through a solid head, and the dip tube alone gets an extra
material wall so it survives below the joint.

    JOINT    24/40
    TUBE_R   3.0 mm outer radius (6 mm OD) / 2.0 mm inner (4 mm ID)
    DIP_LEN  110 mm  below the male tip (to ~10 mm above the bottle base)
    TOP_R    18 mm   grip / flange radius
    TOL      0.01
Usage: python3 make_head.py [output.3mf]
"""
import math
import sys

import csglib as csg
import junctionlib as J
import threelib as t

# --- Parameters ------------------------------------------------------------
JOINT = "25/40"
TUBE_OUT = 3.6
TUBE_IN = 2.6
DIP_X = 4.0
VENT_X = -4.0
DIP_BELOW = 110.0
TOP_R = 18.0
TOP_H = 12.0
TOL = 0.01

PARAMS = dict(joint=JOINT, tube_out=TUBE_OUT, tube_in=TUBE_IN, dip_x=DIP_X,
              vent_x=VENT_X, dip_below=DIP_BELOW, top_r=TOP_R, top_h=TOP_H,
              tol=TOL)

def default_out():
    return t.stamped_filename("washbottle_head", PARAMS)

def _joint():
    return J.joint_dims(JOINT)

def build(**p):
    r_mouth, r_tip, glen = _joint()
    n = t.auto_segments(max(r_mouth, TOP_R, TUBE_OUT), TOL)

    # Male (ground) cone, tip at z = 0, mouth at z = glen. Solid of revolution.
    male = t.Mesh("male")
    t.revolve(male, [(0.0, 0.0), (r_tip, 0.0), (r_mouth, glen), (0.0, glen)], n)

    # Grip flange / top body above the joint, plus the material for the dip
    # tube. One boolean chain, repaired once at the end, so no micro-gaps
    # remain where the facet planes meet.
    top = t.Mesh("top")
    t.cylinder(top, TOP_R, glen, glen + TOP_H, n)
    dip_outer = t.Mesh("dip_outer")
    t.cylinder(dip_outer, TUBE_OUT, -DIP_BELOW, glen + TOP_H + 2.0, n)
    t.translate(dip_outer, DIP_X, 0.0, 0.0)
    base = csg.combine([("start", male), ("union", top),
                        ("union", dip_outer)], name="base")

    # Drill voids.
    dip_void = t.Mesh("dip_void")
    t.cylinder(dip_void, TUBE_IN, -DIP_BELOW - 1.0, glen + TOP_H + 3.0, n)
    t.translate(dip_void, DIP_X, 0.0, 0.0)

    vent_void = t.Mesh("vent_void")
    t.cylinder(vent_void, TUBE_IN, -0.5, glen + TOP_H + 3.0, n)
    t.translate(vent_void, VENT_X, 0.0, 0.0)

    head = csg.combine([("start", base), ("subtract", dip_void),
                        ("subtract", vent_void)], name="head")
    return head, {"n": n, "glen": glen, "r_tip": r_tip, "r_mouth": r_mouth,
                  "bbox": head.bbox()}

def main():
    out = sys.argv[1] if len(sys.argv) > 1 else default_out()
    mesh, s = build(**PARAMS)
    bb = s["bbox"]
    desc = ("Gas washing bottle head: 24/40 male ground plug with a dip tube "
            "(4 mm bore) to near the bottle bottom and a vent hole at the top.")
    t.write_3mf(out, (mesh.name, mesh), title="Gas washing bottle head",
                description=desc, params=PARAMS)
    print(f"wrote {out}")
    print(f"  params            : {t.describe(PARAMS)}")
    print(f"  param hash        : {t.params_hash(PARAMS)}")
    print(f"  joint             : {JOINT} male, tip {s['r_tip']*2:g} mm, "
          f"mouth {s['r_mouth']*2:g} mm, ground {s['glen']:g} mm")
    print(f"  tubes             : OD {TUBE_OUT*2:g} / ID {TUBE_IN*2:g} mm, "
          f"dip {DIP_BELOW:g} mm below the tip at x={DIP_X:g}, vent at x={VENT_X:g}")
    print(f"  height            : {bb[5]-bb[2]:.1f} mm")
    print(f"  bounding box      : {bb[3]-bb[0]:g} x {bb[4]-bb[1]:g} "
          f"x {bb[5]-bb[2]:g} mm")
    print(f"  volume            : {mesh.signed_volume():.1f} mm^3")
    print(f"  vertices / tris   : {len(mesh.verts)} / {len(mesh.tris)}")
    print(f"  shells            : {len(mesh.shells())}")

if __name__ == "__main__":
    main()

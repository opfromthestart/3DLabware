#!/usr/bin/env python3
"""Powder-addition funnel: the cap.

A 51 mm-over male ground plug for the body's 50/35 mouth, bored straight
through so the auger's shaft rides in it.

It follows the repo's seal convention -- make_washbottle_head.py plugs a
24/40 neck with a "25/40" cone. The cap's cone is 1 mm over the socket's
mouth diameter, so the two tapers coincide nowhere: pushed home, they jam
and the cap stops at the depth where the interference closes -- ~21 mm of
the 35 mm ground engaged, measured from the meshes in verify_pfunnel_fit.
A lid that seals by its own taper, not by an accident of calibration.

The 7.5 mm hole through the whole cap is the auger's bushing: the shaft
is 7.0, so 0.25 mm of radial play. That is a rotating through-hole and it
breathes a little -- which is exactly why the body carries its own
equalizer bore: the gas path from flask to headspace runs through the
wall, so the shaft's leak is only the margin, never the seal.

The cap is ONE stepped revolution -- tip face, ground cone, mouth
ledge, grip, top face -- with the bushing as a single transverse
subtract; a revolve profile cannot run back down the bore, and must not. It started life as a cone unioned
with a 12-sided grip prism; the boolean left 33-12 real gaps no matter
what embed depth or facet count was tried, because the grip's flat walls
meet the cone's flat cap in coplanar seams a BSP cannot resolve. A round
grip revolved together with the cone has no seam at all, and it prints
the same: 2.9 deg taper, vertical grip, every ledge facing up.

Usage: python3 make_pfunnel_cap.py [output.3mf]
"""

import math
import sys

import csglib as csg
import junctionlib as J
import printability as PR
import threelib as t

import make_pfunnel_body as B

# --- Parameters ------------------------------------------------------------
MOUTH_JOINT = B.PARAMS["mouth_joint"]   # 50/35: the body's socket
CAP_OVER = 1.0                          # mm over the mouth diameter: the seal
SHAFT_R = 3.5                           # the auger's shaft radius, mm
BUSH_CLEAR = 0.25                       # radial play in the bushing
GRIP_R = 12.0
GRIP_H = 20.0
TOL = 0.01

PARAMS = dict(joint=MOUTH_JOINT, cap_over=CAP_OVER, shaft_r=SHAFT_R,
              bush_clear=BUSH_CLEAR, grip_r=GRIP_R, grip_h=GRIP_H,
              tol=TOL, body_hash=t.params_hash(B.PARAMS))

def default_out():
    return t.stamped_filename("pfunnel_cap", PARAMS)

def dims(p):
    m_mouth, m_tip, m_glen = J.joint_dims(p["joint"])
    c_mouth = m_mouth + p["cap_over"] / 2.0     # the plug's mouth radius
    c_tip = c_mouth - m_glen * J.TAPER          # and its tip
    bush_r = p["shaft_r"] + p["bush_clear"]
    # Seat depth is not restated here: the two cones have different slopes
    # (the cavity's rungs off a 36.4 deg cone foot), so where they jam is
    # measured off the meshes by verify_pfunnel_fit.py.
    return dict(m_mouth=m_mouth, m_tip=m_tip, m_glen=m_glen,
                c_mouth=c_mouth, c_tip=c_tip, bush_r=bush_r,
                height=m_glen + p["grip_h"])

def build(**p):
    d = dims(p)
    n = t.auto_segments(d["c_mouth"], p["tol"])
    # Solid of revolution: tip face, ground cone, mouth ledge, grip, top face.
    # Starts and ends on the axis so revolve makes no disc caps there.
    prof = [
        (0.0, 0.0),
        (d["c_tip"], 0.0),
        (d["c_mouth"], d["m_glen"]),
        (p["grip_r"], d["m_glen"]),   # ledge facing up (z flat, step inward)
        (p["grip_r"], d["height"]),
        (0.0, d["height"]),
    ]
    # a revolve profile must run in increasing z; the ledge is fine (same z,
    # radius steps in), the top face likewise.
    cap = t.Mesh("cap_solid")
    t.revolve(cap, prof, n)
    bore = t.Mesh("bore")
    t.cylinder(bore, d["bush_r"], -1.0, d["height"] + 1.0, n)
    # One subtract, crossing the tip and top faces transversally -- the only
    # boolean here, and the only pattern the seams tolerate.
    cap = csg.combine([("start", cap), ("subtract", bore)], name="pfunnel_cap")
    return cap, d

def main():
    out = sys.argv[1] if len(sys.argv) > 1 else default_out()
    mesh, d = build(**PARAMS)
    bb = mesh.bbox()
    desc = ("Powder funnel cap: 51 mm male ground plug with a 7.5 mm "
            "through-bushing for the auger shaft and a round grip.")
    t.write_3mf(out, (mesh.name, mesh), title="Powder funnel cap",
                description=desc, params=PARAMS)
    print(f"wrote {out}")
    print(f"  params            : {t.describe(PARAMS)}")
    print(f"  param hash        : {t.params_hash(PARAMS)}")
    print(f"  joint             : {MOUTH_JOINT} socket -> cap mouth "
          f"{d['c_mouth']*2:g} mm on a {d['m_glen']:g} mm ground")
    print(f"  bushing           : {d['bush_r']*2:g} mm hole, shaft "
          f"{SHAFT_R*2:g}, play {BUSH_CLEAR:g} mm radial")
    print(f"  grip              : {GRIP_R*2:g} mm round, {GRIP_H:g} mm high")
    print(f"  height            : {bb[5]-bb[2]:.2f} mm")
    print(f"  volume            : {mesh.signed_volume():.1f} mm^3")
    print(f"  vertices / tris   : {len(mesh.verts)} / {len(mesh.tris)}")
    shells = mesh.shells()
    print(f"  shells            : {len(shells)}")
    findings, summary = PR.assess(mesh.verts, mesh.tris,
                                  nozzle=0.4, layer_height=0.2)
    print(f"  printability      : {summary}")
    for f in findings:
        if f["level"] != "info":
            print(f"    [{f['level']}] {f['msg']}")

if __name__ == "__main__":
    main()

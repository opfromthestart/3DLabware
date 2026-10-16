#!/usr/bin/env python3
"""Powder-addition funnel: the body.

One printed piece that holds the powder, meters it, and equalizes the
vessel it feeds. Read it bottom-up, in the pose it prints and works in:

    z 0..42      29/42 male stem; it plugs the neck of the receiving flask
    z 0..42      powder bore, 10 mm, straight through the stem
    z 42..43     the bore converges to the throat by cone -- no shelf
    z 43..44     the throat, 7 mm: the choke. The auger above pushes powder
                 into it; when the auger stops, the plug arching over this
                 pinch holds the rest up. This is what makes "closed"
                 actually closed -- flights alone would leak.
    z 44..50     the bore opens back out, a cone (no ceiling to bridge)
    z 50..74     the screw section, 14 mm straight bore: the auger's home
    z 74..96     the hopper cone, 14 mm to 46.5 mm, 36.4 deg off vertical:
                 steeper than a powder's angle of repose, shallower than
                 the printer's overhang limit -- it empties and it prints
    z 96..131    the mouth: a 50/35 female ground joint. The fill port,
                 and sealed by the cap (make_pfunnel_cap.py) when not filling.

The equalizer is a 3.5 mm bore drilled inside the body's own wall, not a
sidearm: it starts at the stem's tip face next to the powder bore, leans
out from x=8.2 to x=11.0 crossing the ground zone inside the cone wall,
then runs straight up through the tube wall and opens into the hopper
cone's wall at z~77, well above the screw section's bore. Gas therefore
has a path tip-to-headspace that never shares a wall with the powder
path below z=74 -- the flask and the hopper see the same pressure while
the powder is only allowed out past the auger. A sidearm cannot do this:
any port through a seated ground cone's contact zone destroys the seal
it lives beside, and there is no room outside the cone anyway.

The wall's ligaments are tight on purpose and are measured in
`verify_pfunnel_body.py`: 2.25 mm between the equalizer and the screw
bore, 1.45 mm between the equalizer and the powder bore at the tip face.

Printability: the stem cone and socket cone lean 1:20 (2.9 deg), the
hopper cone 36.45 deg, the collar chamfer and both throat tapers sit
past 45 deg from horizontal -- the body carries NO horizontal
down-facing face. The first build had two: the bore's flat step into the
throat (40 mm^2) and the tube's flat collar over the stem mouth
(144 mm^2); both became cones, and the verifier now asserts the flat
ceiling total is nil.

Usage: python3 make_pfunnel_body.py [output.3mf]
"""

import math
import sys

import csglib as csg
import junctionlib as J
import memcap
import printability as PR
import threelib as t

# --- Parameters ------------------------------------------------------------
STEM_JOINT = "29/42"       # male ground plug into the receiving flask
MOUTH_JOINT = "50/35"      # female ground socket for the cap; fill port
STEM_BORE_R = 5.0          # powder bore through the stem, mm
THROAT_R = 3.5             # the choke below the screw section, mm
SCREW_BORE_R = 7.0         # the auger's bore, mm
TUBE_OD_R = 16.0           # outside of the screw-section tube, mm
CONE_D = 22.0              # hopper cone height, 14 mm -> 46.5 mm
CONE_WALL = 7.5            # hopper wall at the cone (carries the equalizer)
EQ_BORE_R = 1.75           # equalizer bore radius, mm
EQ_X0 = 8.2                # ... opens on the tip face at this x,
EQ_X1 = 11.0               # ... then runs vertical at this x,
EQ_TOP = 86.0              # ... and stops here, inside the hopper void
TOL = 0.01

PARAMS = dict(stem_joint=STEM_JOINT, mouth_joint=MOUTH_JOINT,
              stem_bore_r=STEM_BORE_R, throat_r=THROAT_R,
              screw_bore_r=SCREW_BORE_R, tube_od_r=TUBE_OD_R,
              cone_d=CONE_D, cone_wall=CONE_WALL,
              eq_bore_r=EQ_BORE_R, eq_x0=EQ_X0, eq_x1=EQ_X1,
              eq_top=EQ_TOP, tol=TOL)

def default_out():
    return t.stamped_filename("pfunnel_body", PARAMS)

def dims(p):
    """Every station and radius, derived from the stamped parameters."""
    s_mouth, s_tip, s_glen = J.joint_dims(p["stem_joint"])
    m_mouth, m_tip, m_glen = J.joint_dims(p["mouth_joint"])
    z_throat0 = s_glen + 1.0              # 43: the pinch starts
    z_throat1 = z_throat0 + 1.0           # 44: 1 mm of straight 7 mm throat
    z_screw0 = z_throat1 + 6.0            # 50: bore fully open again
    z_screw1 = z_screw0 + 24.0               # 74: top of the straight bore
    z_cone1 = z_screw1 + p["cone_d"]         # 96: hopper cone done
    z_mouth = z_cone1 + m_glen               # 131: rim plane
    # The hopper cone ends exactly at the socket's tip radius, so the whole
    # cavity -- cone plus socket -- is one 1:20 generatrix. A cap cone of the
    # same taper then rides that line: zero gap along the full engaged
    # length, preload built by over-travelting the wedge. A cone foot a few
    # thousandths off (atan 36.4 deg gives 23.198) makes the two tapers
    # disagree in slope, and the joint would touch at its tip and nowhere
    # else -- sealed on paper, ventilating in the hand.
    r_cone1 = m_tip                          # 23.25 for a 50/35 mouth
    cone_slope = (r_cone1 - p["screw_bore_r"]) / p["cone_d"]
    d = dict(s_mouth=s_mouth, s_tip=s_tip, s_glen=s_glen,
             m_mouth=m_mouth, m_tip=m_tip, m_glen=m_glen,
             z_throat0=z_throat0, z_throat1=z_throat1,
             z_screw0=z_screw0, z_screw1=z_screw1, z_cone1=z_cone1,
             z_mouth=z_mouth, r_cone1=r_cone1, cone_slope=cone_slope,
             r_od_cone1=r_cone1 + p["cone_wall"])
    return d

def build(**p):
    memcap.cap(quiet=True)
    d = dims(p)
    n = t.auto_segments(d["r_od_cone1"], p["tol"])

    # Outer solid: one closed revolve profile, tip face at z=0, rim at the top.
    # The stem's cone ends at z=s_glen; the tube stands on its mouth plane as
    # a vertical wall; the hopper cone rises at 36.4 deg; the foot of the
    # socket steps *down* (an upward-facing ledge, free to print); the socket
    # wall then runs parallel to the cavity, 2.75 mm thick.
    outer_prof = [
        (0.0, 0.0),
        (d["s_tip"], 0.0),
        (d["s_mouth"], d["s_glen"]),
        # the tube stands on the stem mouth through a 2 mm chamfer, not a
        # flat collar: a 1.5 mm downward-facing ring would otherwise be the
        # largest ceiling the body carries (and it was, at 144 mm^2)
        (p["tube_od_r"], d["s_glen"] + 2.0),
        (p["tube_od_r"], d["z_screw1"]),
        (p["tube_od_r"], d["z_screw1"] + 2.0),
        (d["r_od_cone1"], d["z_cone1"]),
        (d["r_od_cone1"], d["z_cone1"] + 1.0),
        (d["m_mouth"] + 3.0, d["z_cone1"] + 1.0),
        (d["m_mouth"] + 3.0, d["z_mouth"]),
        (0.0, d["z_mouth"]),
    ]
    outer = t.Mesh("outer")
    t.revolve(outer, outer_prof, n)

    # Powder void: bore -> throat -> cone open -> straight screw bore ->
    # hopper cone -> socket cavity. Overshoot both ends so the subtraction
    # cuts clean faces at z=0 and above the rim.
    void_prof = [
        (0.0, -2.0),
        (p["stem_bore_r"], -2.0),
        (p["stem_bore_r"], d["z_throat0"] - 2.0),   # converge to the throat
        (p["throat_r"], d["z_throat0"]),            # by cone, never by a
        (p["throat_r"], d["z_throat1"]),            # horizontal ceiling
        (p["screw_bore_r"], d["z_screw0"]),
        (p["screw_bore_r"], d["z_screw1"]),
        (d["r_cone1"], d["z_cone1"]),
        (d["m_mouth"], d["z_mouth"]),
        (d["m_mouth"], d["z_mouth"] + 2.0),
        (0.0, d["z_mouth"] + 2.0),
    ]
    void = t.Mesh("void")
    t.revolve(void, void_prof, n)

    # Equalizer: a slanted leg from the tip face out to x1, then vertical
    # up into the cone's wall. Its top cap must finish inside the hopper
    # void, which it does: the void radius passes x1+eq_bore_r at z~77.
    # Equalizer: one continuous loft, lean then straight -- two overlapping
    # drills would intersect each other and seed sliver faces along the seam.
    # The shared circle at the bend keeps the surface continuous, and the
    # bottom cap at z=-3 opens the hole through the tip face.
    eq = t.Mesh("eq")
    z_bend = d["z_throat0"] + 4.0
    secs = []
    for zc, xc in ((-3.0, p["eq_x0"]), (z_bend, p["eq_x1"]),
                   (p["eq_top"], p["eq_x1"])):
        secs.append([(x, y, zc) for (x, y) in
                     t.disc_polygon(p["eq_bore_r"], n, cx=xc, cy=0.0)])
    t.loft(eq, secs)

    body = csg.combine([("start", outer), ("subtract", void),
                        ("subtract", eq)],
                       name="pfunnel_body")
    return body, d

def main():
    out = sys.argv[1] if len(sys.argv) > 1 else default_out()
    mesh, d = build(**PARAMS)
    bb = mesh.bbox()
    desc = ("Powder-addition funnel body: hopper with a 50/35 mouth, 14 mm "
            "auger tube over a 7 mm throat, 29/42 stem, and an internal "
            "3.5 mm equalizer bore from the tip face to the headspace.")
    t.write_3mf(out, (mesh.name, mesh), title="Powder funnel body",
                description=desc, params=PARAMS)
    slope = math.degrees(math.atan(
        (d["r_cone1"] - PARAMS["screw_bore_r"]) / PARAMS["cone_d"]))
    eq_ang = math.degrees(math.atan2(
        PARAMS["eq_x1"] - PARAMS["eq_x0"],
        d["z_throat0"] + 3.0))
    print(f"wrote {out}")
    print(f"  params            : {t.describe(PARAMS)}")
    print(f"  param hash        : {t.params_hash(PARAMS)}")
    print(f"  stem joint        : {PARAMS['stem_joint']} male, mouth at "
          f"z={d['s_glen']:g}, powder bore {PARAMS['stem_bore_r']*2:g} mm")
    print(f"  throat            : {PARAMS['throat_r']*2:g} mm from z="
          f"{d['z_throat0']:g} to {d['z_throat1']:g}, under the "
          f"{PARAMS['screw_bore_r']*2:g} mm screw bore")
    print(f"  screw section     : {PARAMS['screw_bore_r']*2:g} mm bore, z="
          f"{d['z_screw0']:g}..{d['z_screw1']:g}, wall out to "
          f"{PARAMS['tube_od_r']*2:g} mm")
    print(f"  hopper cone       : {PARAMS['screw_bore_r']*2:g} to "
          f"{d['r_cone1']*2:.2g} mm over {PARAMS['cone_d']:g} mm "
          f"({slope:.1f} deg off vertical), z={d['z_screw1']:g}..{d['z_cone1']:g}")
    print(f"  mouth joint       : {PARAMS['mouth_joint']} female, rim at "
          f"z={d['z_mouth']:g}")
    print(f"  equalizer         : {PARAMS['eq_bore_r']*2:g} mm bore from "
          f"(x={PARAMS['eq_x0']:g}, z=0) leaning to x={PARAMS['eq_x1']:g} "
          f"({eq_ang:.1f} deg off vertical), up to z={PARAMS['eq_top']:g}")
    print(f"  size              : {bb[3]-bb[0]:.2f} x {bb[4]-bb[1]:.2f} x "
          f"{bb[5]-bb[2]:.2f} mm")
    print(f"  volume            : {mesh.signed_volume():.1f} mm^3")
    print(f"  vertices / tris   : {len(mesh.verts)} / {len(mesh.tris)}")
    shells = mesh.shells()
    print(f"  shells            : {len(shells)}")
    print(f"  gaps              : {len(csg.real_gaps(mesh.verts, mesh.tris))} "
          f"real, {csg.crack_length(mesh.verts, mesh.tris):.2e} mm of edge "
          f"used once")
    findings, summary = PR.assess(mesh.verts, mesh.tris,
                                  nozzle=0.4, layer_height=0.2)
    print(f"  printability      : {summary}")
    for f in findings:
        if f["level"] != "info":
            print(f"    [{f['level']}] {f['msg']}")

if __name__ == "__main__":
    main()

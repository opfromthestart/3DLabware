#!/usr/bin/env python3
"""Powder-addition funnel: the auger.

The spiral that pushes powder down the body's 14 mm tube. Read it in the
pose it *works* in (the file is written upside down: it prints standing
on its knob, so every layer is the knob disc, a slim shaft, the band
turning around it, and the shaft tip -- no overhang anywhere):

    use z 48        shaft tip, a millimetre or two above the throat ledge,
                    out of the 7 mm pinch so nothing scrapes
    use z 51..99    the flight, two arms of one right-hand helix at 42 mm
                    pitch, from r 3.3 (buried in the shaft) to r 6.5 -- the ribbon
                    corners, hypot(6.5, 0.8) = 6.55, keep 0.45 mm of
                    bore-wall gap at every azimuth
    use z 161.7+1.5 knob underside, resting clear above the seated cap's
                    top so print shrinkage of the ground seat cannot make
                    the knob carry the cap instead of the bushing

**Why the pitch is 42 mm and not "tight like a real screw".** An Archimedean
flight is a helicoid; at radius r its surface leans atan(r * 2*pi/pitch)
off horizontal. A 12 mm pitch would put the flight faces 16 deg off
level -- shelves the printer must bridge, unprintable standing up, which
is the rule `make_column.py` already refuses to break for its helix. At
42 mm the outermost fibre leans 44.2 deg, inside the 45 limit, and the
two arms still meter: each turn of the pair displaces powder down by the
flight spacing that crosses any given height -- a coarse but hand-turnable
dose. Turning clockwise (looking down) walks the right-hand helix's low
point around the tube downward, dragging powder with it.

**Why clearance, not a sealing flight.** 0.5 mm radial at the tip is too
big to seal gas and too big to meter paste -- but it is not meant to: the
powder plugs its own flights once the bed above feeds, the 7 mm throat
below arches against anything that has settled, and the body's internal
equalizer removes every pressure gradient that would push powder past a
stopped screw. The stopcock here is the screw *plus the choke*, and both
are printed at their honest sizes rather than 0.08 mm apart like a key
that would jam the first time a grain of anything got between them.

The shells overlap by intent (flight roots buried 0.2 mm in the shaft):
the printed part is a union of four closed shells, exactly the standing
convention of `make_column.py`'s staging. Its reported volume is the
sum, not the union; the union arithmetic is in the report.

Usage: python3 make_pfunnel_screw.py [output.3mf]
"""

import math
import sys

import printability as PR
import threelib as t

from make_column import twisted_band
import make_pfunnel_body as B
import make_pfunnel_cap as C

# --- Parameters ------------------------------------------------------------
PITCH = 42.0            # mm per turn; see the lead-angle argument above
ARMS = 2                # starts: two arms of one helix, 180 deg apart
FLIGHT_R = 6.5          # tip radius; the bore clearance is body-derived
FLIGHT_TH = 1.6         # ribbon thickness (four beads)
FLIGHT_Z0 = 51.0        # use-frame flight bottom, above the ramp shoulder
FLIGHT_Z1 = 99.0        # use-frame flight top, below a seated cap tip
KNOB_R = 14.0
KNOB_H = 8.0
KNOB_CLEAR = 1.5        # gap knob underside .. seated cap top, mm
ROOT_BITE = 0.2         # flight root buried in the shaft, by that much
OVERHANG = PR.OVERHANG_DEG

PARAMS = dict(pitch=PITCH, arms=ARMS, flight_r=FLIGHT_R,
              flight_th=FLIGHT_TH, flight_z0=FLIGHT_Z0, flight_z1=FLIGHT_Z1,
              knob_r=KNOB_R, knob_h=KNOB_H, knob_clear=KNOB_CLEAR,
              root_bite=ROOT_BITE, overhang=OVERHANG,
              body_hash=t.params_hash(B.PARAMS),
              cap_hash=t.params_hash(C.PARAMS))

def default_out():
    return t.stamped_filename("pfunnel_screw", PARAMS)

def dims(p):
    bd = B.dims(dict(B.PARAMS))
    cd = C.dims(dict(C.PARAMS))
    cd["shaft_r"] = float(C.PARAMS["shaft_r"])
    # Where a seated cap's tip stops: cavity cone from r_cone1 at z_cone1 to
    # m_mouth at the rim; the cap is fatter than that cone, so its tip -- the
    # first thing to jam -- stops where its radius equals the cavity's.
    slope_c = (bd["m_mouth"] - bd["r_cone1"]) / bd["m_glen"]
    z_seat_tip = bd["z_cone1"] + (cd["c_tip"] - bd["r_cone1"]) / slope_c
    cap_top = z_seat_tip + cd["height"]
    knob_bottom = cap_top + p["knob_clear"]
    # shaft tip: above the ramp's floor (z_throat1 + a hair clears the 7 mm
    # pinch at 42..44 and its 14 mm opening ramp).
    z_shaft0 = bd["z_throat0"] + 6.0
    d = dict(bd=bd, cd=cd, slope_c=slope_c, z_seat_tip=z_seat_tip,
             cap_top=cap_top, knob_bottom=knob_bottom,
             z_shaft0=z_shaft0, z_shaft1=knob_bottom + p["knob_h"],
             total=knob_bottom + p["knob_h"],
             build_flip=knob_bottom + p["knob_h"] - 0.0)
    return d

def build(**p):
    d = dims(p)
    bd, cd = d["bd"], d["cd"]
    if p["flight_r"] >= float(B.PARAMS["screw_bore_r"]):
        raise ValueError("the flight tips are inside the bore, not the bore's wall")
    n = t.auto_segments(p["knob_r"], 0.01)
    screw = t.Mesh("pfunnel_screw")

    # Everything is built in the print pose -- knob at z=0 growing up --
    # by mirroring the use frame: b = z_shaft1 - u.
    def b_of(u):
        return d["z_shaft1"] - u

    t.cylinder(screw, p["knob_r"], 0.0, p["knob_h"], n)
    t.cylinder(screw, cd["shaft_r"], p["knob_h"] - 1.0, b_of(d["z_shaft0"]), n)

    turns = (p["flight_z1"] - p["flight_z0"]) / p["pitch"]
    lean = None
    for i in range(int(p["arms"])):
        lean = twisted_band(screw, b_of(p["flight_z1"]), b_of(p["flight_z0"]),
                            cd["shaft_r"] - p["root_bite"], p["flight_r"],
                            turns, p["flight_th"], p["overhang"],
                            offset_a=i * 2.0 * math.pi / p["arms"])

    info = dict(d=d, turns=turns, lean=lean,
                clearance=float(B.PARAMS["screw_bore_r"]) - p["flight_r"])
    return screw, info

def main():
    out = sys.argv[1] if len(sys.argv) > 1 else default_out()
    mesh, info = build(**PARAMS)
    d = info["d"]
    bb = mesh.bbox()
    desc = ("Powder funnel auger: 7 mm shaft, two-start 42 mm pitch helical "
            "flight 0.45 mm clear of the 14 mm bore, round knob. Prints "
            "standing on its knob.")
    t.write_3mf(out, (mesh.name, mesh), title="Powder funnel auger",
                description=desc, params=PARAMS)
    print(f"wrote {out}")
    print(f"  params            : {t.describe(PARAMS)}")
    print(f"  param hash        : {t.params_hash(PARAMS)}")
    print(f"  flight            : {PARAMS['arms']} arms, pitch "
          f"{PARAMS['pitch']:g}, tip r {PARAMS['flight_r']:g} -> clearance "
          f"{info['clearance']:g} mm/side in the {float(B.PARAMS['screw_bore_r'])*2:g} "
          f"mm bore, lean {info['lean']:.1f} deg off vertical")
    print(f"  use frame         : shaft tip z={d['z_shaft0']:g}, flight "
          f"z={PARAMS['flight_z0']:g}..{PARAMS['flight_z1']:g}, knob "
          f"z={d['knob_bottom']:.1f}..{d['z_shaft1']:.1f}, seated cap top "
          f"{d['cap_top']:.1f}")
    print(f"  build pose        : stands on the knob; use frame mirrored "
          f"(b = {d['z_shaft1']:.1f} - u)")
    print(f"  size              : {bb[3]-bb[0]:.2f} x {bb[4]-bb[1]:.2f} x "
          f"{bb[5]-bb[2]:.2f} mm")
    print(f"  volume (sum)      : {mesh.signed_volume():.1f} mm^3 -- four "
          f"overlapping shells, union is smaller by the buried roots")
    print(f"  vertices / tris   : {len(mesh.verts)} / {len(mesh.tris)}")
    print(f"  shells            : {len(mesh.shells())}")
    findings, summary = PR.assess(mesh.verts, mesh.tris,
                                  nozzle=0.4, layer_height=0.2)
    print(f"  printability      : {summary}")
    for f in findings:
        if f["level"] != "info":
            print(f"    [{f['level']}] {f['msg']}")

if __name__ == "__main__":
    main()

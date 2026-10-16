#!/usr/bin/env python3
"""Verify the powder-funnel auger from its 3MF.

The part is four intentionally overlapping closed shells -- knob, shaft,
two flight arms -- standing on the knob in the print pose. The verifier
identifies each by its bounding signature, then measures from the mesh:

  * the flight tip radius and its clearance against the body's bore;
  * the helix itself: unwrap each arm's azimuth against height and fit
    the pitch -- right-hand (+) in the print frame, which is a proper
    rotation of the use frame, so still right-hand when it works: a
    clockwise turn (looking down) walks powder to the tip;
  * the two arms are one helix, 180 deg apart, not one arm built twice;
  * the only down-facing faces are the short ribbon-end shelves;
  * a closed-form volume sum.

Usage: python3 verify_pfunnel_screw.py [pfunnel_screw.3mf] [--angle N]
"""

import math
import sys

import checks as K
import junctionlib as J
import make_pfunnel_body as B
import make_pfunnel_cap as C
import make_pfunnel_screw as S
import printability as PR

PATH, ANGLE = K.cli(sys.argv[1:], S.default_out())
P = K.read_params(PATH)
rep = K.Report(f"Powder funnel auger: {PATH}")
objs, _, unit = K.load_3mf(PATH, rep)
rep.check("one object", len(objs) == 1, f"{len(objs)}")
V, T = objs[0]["verts"], objs[0]["tris"]

pd = {k: (float(P[k]) if isinstance(S.PARAMS[k], (int, float))
          else P.get(k, S.PARAMS[k])) for k in S.PARAMS}

import threelib as t
rep.check("body coupling hash matches make_pfunnel_body",
          str(P.get("body_hash")) == t.params_hash(B.PARAMS),
          f"file says {P.get('body_hash')}, body now {t.params_hash(B.PARAMS)}")
rep.check("cap coupling hash matches make_pfunnel_cap",
          str(P.get("cap_hash")) == t.params_hash(C.PARAMS),
          f"file says {P.get('cap_hash')}, cap now {t.params_hash(C.PARAMS)}")

rep.check("signed volume positive", K.signed_volume(V, T) > 0,
          f"{K.signed_volume(V, T):.1f} mm^3 sum of shells")
tri_shells = K.shell_triangles(V, T)
rep.check("four shells: knob, shaft, two flight arms",
          len(tri_shells) == 4, f"{len(tri_shells)}")


def sig(sh):
    vs = [V[i] for k in sh for i in T[k]]
    rs = [math.hypot(x, y) for x, y, z in vs]
    zs = [z for x, y, z in vs]
    return {"rmax": max(rs), "rmin": min(rs), "z0": min(zs), "z1": max(zs),
            "vs": vs}

sigs = [sig(s) for s in tri_shells]
knobs = [s for s in sigs if abs(s["rmax"] - pd["knob_r"]) < 0.05]
shafts = [s for s in sigs if abs(s["rmax"] - float(C.PARAMS["shaft_r"])) < 0.05]
bands = [s for s in sigs if abs(s["rmax"] - pd["flight_r"]) < 0.05
         and s["rmin"] < 0.05 + float(C.PARAMS["shaft_r"])]
rep.check("exactly one knob, one shaft, two arms",
          len(knobs) == 1 and len(shafts) == 1 and len(bands) == 2,
          f"{len(knobs)}/{len(shafts)}/{len(bands)}")

k, sh = knobs[0], shafts[0]
rep.near("knob height", k["z1"], pd["knob_h"], 1e-6, "{:.3f}")
rep.near("knob on the bed", k["z0"], 0.0, 1e-6, "{:.3f}")
d = S.dims(pd)
rep.near("shaft reaches the shaft tip (build pose)",
         sh["z1"], d["z_shaft1"] - d["z_shaft0"], 1e-6, "{:.3f}")

half_th = pd["flight_th"] / 2.0
tip_corner = math.hypot(pd["flight_r"], half_th)
root_corner = math.hypot(float(C.PARAMS["shaft_r"]) - pd["root_bite"], half_th)
for i, s in enumerate(bands):
    # the ribbon's section is a radial rectangle: its corners sit further
    # out than the nominal radius, hypot(r, th/2), and that is the real tip
    rep.near(f"arm {i} tip corner radius", s["rmax"], tip_corner, 0.01,
             "{:.4f}")
    rep.near(f"arm {i} root corner radius", s["rmin"], root_corner, 0.01,
             "{:.4f}")
    rep.near(f"arm {i} flight span", s["z1"] - s["z0"],
             pd["flight_z1"] - pd["flight_z0"], 1e-6, "{:.3f}")
    rep.near(f"arm {i} build bottom", s["z0"], d["z_shaft1"] - pd["flight_z1"],
             0.05, "{:.3f}")

# -- the helix: unwrap azimuth against height and fit rate --------------------
clear = float(B.PARAMS["screw_bore_r"]) - tip_corner
rep.check("radial clearance in the bore", clear >= 0.4,
          f"{clear:.3f} mm between flight corner and wall")

def arm_rate(s):
    vs = sorted(s["vs"], key=lambda p: p[2])
    phis, ups = [], []
    prev = None
    for x, y, z in vs:
        if z < s["z0"] + 2 or z > s["z1"] - 2:
            continue
        a = math.atan2(y, x)
        if prev is not None:
            while a - prev > math.pi:
                a -= 2 * math.pi
            while a - prev < -math.pi:
                a += 2 * math.pi
        prev = a
        phis.append(a)
        ups.append(z)
    n = len(phis)
    zm, am = sum(ups) / n, sum(phis) / n
    return (sum((ups[i] - zm) * (phis[i] - am) for i in range(n))
            / sum((z - zm) ** 2 for z in ups))

want_rate = 2 * math.pi * (pd["flight_z1"] - pd["flight_z0"]) / pd["pitch"] \
    / (pd["flight_z1"] - pd["flight_z0"])
rates = [arm_rate(s) for s in bands]
for i, r in enumerate(rates):
    rep.near(f"arm {i} helix rate dphi/dz", abs(r), want_rate, 0.01,
             "{:.5f} rad/mm")
rep.check("right-hand helix (azimuth grows with build height)",
          all(r > 0 for r in rates), f"rates {rates}")

# the arms are 180 deg apart, not the same arm twice: at a shared height
# their azimuths differ by ~pi (mod the unwrap).
zz = (bands[0]["z0"] + bands[0]["z1"]) / 2.0
angs = []
for s in bands:
    cand = [(math.atan2(y, x), z) for x, y, z in s["vs"] if abs(z - zz) < 1.0]
    angs.append(sum(a for a, z in cand) / len(cand))
diff = abs((angs[0] - angs[1] + math.pi) % (2 * math.pi) - math.pi)
rep.near("arm phase difference", diff, math.pi, 0.05, "{:.3f} rad")

# -- printability ---------------------------------------------------------------
angle = PR.OVERHANG_DEG if ANGLE is None else ANGLE
ov = PR.overhangs(V, T, max_angle=angle)
sl = [f for f in ov if f["kind"] == "sloped"]
rep.check(f"no sloped overhang at {angle:g} deg", not sl,
          f"{len(sl)} faces, {sum(f['area'] for f in sl):.2f} mm^2")
flat = [f for f in ov if f["kind"] == "flat" and f["area"] > 0.01]
rep.check("down-facing faces are only the ribbon-end shelves",
          sum(f["area"] for f in flat) < 60.0,
          f"{len(flat)} faces, {sum(f['area'] for f in flat):.1f} mm^2")
findings, summary = PR.assess(V, T, nozzle=0.4, layer_height=0.2)
rep.check("nothing narrower than the nozzle",
          summary["n_thin_impossible"] == 0, str(summary["n_thin_impossible"]))
rep.check("bed contact", summary["bed_contact"] > 500.0,
          f"{summary['bed_contact']:.1f} mm^2")

# -- closed-form sum --------------------------------------------------------------
r_mean = (float(C.PARAMS["shaft_r"]) - pd["root_bite"] + pd["flight_r"]) / 2
turns = (pd["flight_z1"] - pd["flight_z0"]) / pd["pitch"]
band_len = turns * math.hypot(2 * math.pi * r_mean, pd["pitch"])
band = pd["arms"] * (pd["flight_r"] - float(C.PARAMS["shaft_r"])
                       + pd["root_bite"]) * pd["flight_th"] * band_len
knobv = math.pi * pd["knob_r"] ** 2 * pd["knob_h"]
shaftv = math.pi * float(C.PARAMS["shaft_r"]) ** 2 * (
    d["z_shaft1"] - d["z_shaft0"] - pd["knob_h"] + 1.0)
want = knobv + shaftv + band
got = K.signed_volume(V, T)
rep.near("volume sum against closed form", got / want, 1.0, 0.05,
         f"got {got:.0f}, closed form {want:.0f} mm^3")

K.render(V, T, "pfunnel_screw_side.png", view="xz", su=4.0)
sys.exit(0 if rep.finish() else 1)

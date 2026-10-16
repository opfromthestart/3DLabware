#!/usr/bin/env python3
"""Verify the powder funnel ASSEMBLY: body + cap + auger, in one pose.

Each part passes its own verifier and still may not fit -- the class of
bug that survives to the printer -- so this reads all three 3MFs, puts
them together in the body's frame, and asks only mesh questions:

  * coupling: every file stamps the hashes of the parts it was derived
    from; a change to any other part shows up here, not in silence;
  * the cap seats where the two measured cones actually interfere. The
    cavity cone and the cap cone are fitted from the meshes (plane
    cuts), the seat is their crossing, and what gets asserted is the
    contact length and the interference depth -- the seal -- not a
    nameplate;
  * the knob hovers clear of a seated cap (the bushing guides the
    shaft, the knob must never carry the cap);
  * the flight tip never touches the bore wall, the shaft sits in its
    bushing with the stamped play, and the throat is clear of the
    shaft's tip;
  * no shaft or flight vertex is inside another part's material; the
    cap's vertices interfere only as deep as the seal allowance says;
  * the metering annulus between shaft and bore stays open around most
    of the circle where the powder travels, and the equalizer's
    channel survives its neighbours.

The auger file stands on its knob (print pose). The body-frame pose is
u = z_shaft1 - b applied to (x, y) as-is with y mirrored -- a 180 deg
rotation, a proper rigid motion, so handedness and the
clockwise-feeds-downward behaviour are rotation-invariant, as they must
be.

Usage: python3 verify_pfunnel_fit.py
"""

import math
import sys

import checks as K
import junctionlib as J
import make_pfunnel_body as B
import make_pfunnel_cap as C
import make_pfunnel_screw as S
import threelib as t


def stamp_params(mod):
    p = K.read_params(mod.default_out())
    return {k: (float(p[k]) if isinstance(v, (int, float))
                else p.get(k, v)) for k, v in mod.PARAMS.items()}


rep = K.Report("Powder funnel assembly")
bp, cp, sp = stamp_params(B), stamp_params(C), stamp_params(S)

rep.check("cap was built against the current body",
          str(cp.get("body_hash")) == t.params_hash(B.PARAMS),
          str(cp.get("body_hash")))
rep.check("auger was built against the current body",
          str(sp.get("body_hash")) == t.params_hash(B.PARAMS),
          str(sp.get("body_hash")))
rep.check("auger was built against the current cap",
          str(sp.get("cap_hash")) == t.params_hash(C.PARAMS),
          str(sp.get("cap_hash")))


def load(mod):
    objs, _, _ = K.load_3mf(mod.default_out(), rep)
    V, T = objs[0]["verts"], objs[0]["tris"]
    sh = K.shell_triangles(V, T)
    return V, T, sh, K.shell_boxes(V, T, sh)


VB, TB, SB, XB = load(B)
VC, TC, SC, XC = load(C)
VS, TS, SS, XS = load(S)

bd = B.dims(bp)
cd = C.dims(cp)
sd = S.dims(sp)
zs1 = sd["z_shaft1"]


def in_body(p):
    return K.inside(VB, TB, p, shells=SB, boxes=XB)


def in_seated_cap(p):
    """A body-frame point inside the cap, whose tip sits at z_seat."""
    return K.inside(VC, TC, (p[0], p[1], p[2] - Z_SEAT), shells=SC, boxes=XC)


def in_pose_screw(p):
    """A body-frame point inside the auger, posed for use."""
    x, y, z = p
    return K.inside(VS, TS, (x, -y, zs1 - z), shells=SS, boxes=XS)

# -- fit the two mating cones from the meshes --------------------------------
def cone_fit(V, T, z0, z1, rlo, rhi, n=7):
    xs, ys = [], []
    for k in range(n):
        z = z0 + (z1 - z0) * k / (n - 1)
        rs = [r for r in K.plane_radii(V, T, (0, 0, 0), (0, 0, 1), z)
              if rlo <= r <= rhi]
        if rs:
            xs.append(z)
            ys.append(sum(rs) / len(rs))
    m = len(xs)
    zm, rm = sum(xs) / m, sum(ys) / m
    a = sum((xs[i] - zm) * (ys[i] - rm) for i in range(m)) \
        / sum((x - zm) ** 2 for x in xs)
    return a, rm - a * zm

A_CAV, B_CAV = cone_fit(VB, TB, bd["z_cone1"] + 2, bd["z_mouth"] - 1,
                        22.5, 26.0)
A_CAP, B_CAP = cone_fit(VC, TC, 2.0, cd["m_glen"] - 2.0, 23.0, 26.0)
rep.near("cavity taper dr/dz", A_CAV, J.TAPER, 0.004, "{:.5f}")
rep.near("cap taper dr/dz", A_CAP, J.TAPER, 0.004, "{:.5f}")

# Seat: the tip ring stops exactly where the fitted cavity reaches the cap's
# tip radius -- the first interference, closed form on the two fitted lines.
Z_SEAT = (cd["c_tip"] - B_CAV) / A_CAV
cap_top = Z_SEAT + cd["height"]
rep.check("cap seats within the ground's reach",
          bd["z_cone1"] + 3 < Z_SEAT < bd["z_mouth"] - 15,
          f"tip at z={Z_SEAT:.2f}: {bd['z_mouth'] - Z_SEAT:.1f} mm of the "
          f"{cd['m_glen']:g} mm ground below the rim")

over = float(cp["cap_over"]) / 2.0
band = [math.hypot(x, y) - (A_CAV * (zc + Z_SEAT) + B_CAV)
        for (x, y, zc) in VC if zc <= cd["m_glen"] + 1e-9
        and math.hypot(x, y) > 20.0]
deep = [e for e in band if e > over + 0.06]
touch = [e for e in band if abs(e) < 0.06]
rep.check("no cap vertex interferes deeper than the seal allowance",
          not deep, f"worst {max(band):.3f} mm vs allowance {over:.2f} mm")
rep.check("contact exists -- it is a seat, not a hover", len(touch) > 50,
          f"{len(touch)} vertices within 0.06 mm of the cavity")
# A revolved cone has vertices only at its two station rings, so the
# contact length cannot be counted off vertices -- measure the surfaces
# where a plane actually cuts them. Both parts, at shared heights inside
# the overlap: cap rings at plane (0,0,1) depth z-Z_SEAT in its own verts,
# body cavity at depth z. They must agree to 0.02 mm over a real length.
zc_ = []
for k in range(25):
    z = Z_SEAT + 1.0 + k * (bd["z_mouth"] - Z_SEAT - 1.0) / 24.0
    rc = [r for r in K.plane_radii(VC, TC, (0, 0, -Z_SEAT), (0, 0, 1), z)
          if 20.0 <= r <= 26.0]
    rb = [r for r in K.plane_radii(VB, TB, (0, 0, 0), (0, 0, 1), z)
          if 20.0 <= r <= 26.0]
    if rc and rb and abs(sum(rc) / len(rc) - sum(rb) / len(rb)) < 0.02:
        zc_.append(z)
rep.check("the two surfaces ride together over a real length of ground",
          len(zc_) >= 15 and (max(zc_) - min(zc_)) >= 15.0,
          f"surfaces agree within 0.02 mm from z={min(zc_):.1f} to "
          f"{max(zc_):.1f} -- engagement {max(zc_) - Z_SEAT:.1f} mm of "
          f"{cd['m_glen']:g}")

# -- the auger, posed ---------------------------------------------------------
mv = list(VS)   # every vertex: knob, shaft tip, and 108 mm of turning flight
hits = [tuple(round(c, 2) for c in (x, y, zs1 - zb)) for (x, y, zb) in mv
        if in_body((x, y, zs1 - zb))]
rep.check("no working vertex of the auger is inside the body", not hits,
          f"{len(mv)} verts checked, hits {hits[:4]}")
hits = [tuple(round(c, 2) for c in (x, y, zs1 - zb)) for (x, y, zb) in mv
        if math.hypot(x, y) < 3.6 and in_seated_cap((x, y, zs1 - zb))]
rep.check("no shaft vertex is inside the seated cap", not hits,
          f"hits {hits[:4]}" if hits else "0.25 mm play throughout")

rs = [r for r in K.plane_radii(VB, TB, (0, 0, 0), (0, 0, 1), 62.0) if r < 9]
bore62 = max(rs)
tip_max = max(math.hypot(x, y) for (x, y, zb) in VS
              if 3.6 < math.hypot(x, y) < 10 and zs1 - zb < 130)
rep.check("flight clears the bore wall at every azimuth",
          bore62 - tip_max >= 0.4,
          f"bore {bore62:.3f}, flight corner {tip_max:.3f}, gap "
          f"{bore62 - tip_max:.3f} mm")

u_tip = min(zs1 - zb for (x, y, zb) in VS if math.hypot(x, y) <= 3.6)
rep.check("shaft tip stands clear of the throat opening",
          u_tip >= 46.5, f"shaft tip z={u_tip:.2f}")
gap_knob = (zs1 - float(sp["knob_h"])) - cap_top
# the requirement is that the knob NEVER carries the seated cap, and the
# cone-fit noise on Z_SEAT is a few hundredths: a band, not a point value.
rep.check("knob floats clear of the seated cap",
          0.8 <= gap_knob <= float(sp["knob_clear"]) + 0.15,
          f"gap {gap_knob:.3f} mm (nominal {float(sp['knob_clear']):g})")

shaft_r = max(math.hypot(x, y) for (x, y, zb) in VS if math.hypot(x, y) < 3.55)
bore_r = min(r for r in K.plane_radii(VC, TC, (0, 0, 0), (0, 0, 1), 15.0))
rep.near("bushing radial play", bore_r - shaft_r, float(cp["bush_clear"]),
         0.03, "{:.3f} mm")

# -- the powder's road when everything is together -----------------------------
zmid = (float(sp["flight_z0"]) + float(sp["flight_z1"])) / 2.0
r_mid = 0.5 * (shaft_r + 0.05 + bore62 - 0.05)
open_ang = sum(1 for k in range(72)
               if not (in_pose_screw((r_mid * math.cos(2 * math.pi * k / 72),
                                     r_mid * math.sin(2 * math.pi * k / 72),
                                     zmid))))
rep.check("metering annulus mostly open at flight mid-height",
          open_ang > 40, f"{open_ang}/72 azimuths clear at z={zmid:g}")

for z in (43.5, 46.0):     # the throat under the shaft tip must stay air
    rep.check(f"throat open at z={z:g}", not in_body((0.0, 0.0, z)),
              "axis point is air")

x1 = float(bp["eq_x1"])
for z in (55.0, 75.0, 85.0):
    ok = all(not in_body((x1 - 1.0, 0.0, z)) and not in_body((x1 + 1.0, 0.0, z))
             for _ in [0])
    rep.check(f"equalizer channel open at z={z:g}", ok,
              f"test points at x={x1 - 1:g} and {x1 + 1:g}")

print("  --- interface ledger ---")
print(f"    fill port     : ground-seated cap, interference {over:.2f} mm")
print(f"    shaft/bushing : rotating joint, play {bore_r - shaft_r:.3f} mm "
      f"-- breathes, and the equalizer bore makes that irrelevant to pressure")
print(f"    outlet        : the stem's two tip-face mouths -- the only link "
      f"to the receiving vessel")
print(f"    assembly      : {cap_top + float(sp['knob_h']) + float(sp['knob_clear']):.0f} mm tall "
      f"over the stem tip when capped and dosing")

sys.exit(0 if rep.finish() else 1)

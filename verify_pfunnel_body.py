#!/usr/bin/env python3
"""Verify the powder-funnel body from its 3MF.

Re-parsed from disk; the numbers trusted from the builder are only the
stamped parameters. What the mesh has to prove:

  * the powder passage is open on the axis from the tip face to the rim;
  * the equalizer is a separate open path from its own mouth on the tip
    face into the hopper void, and the wall between it and the powder
    path is material where they run side by side;
  * the tip face carries both mouths;
  * the stem's outside and the socket's cavity are fitted 1:20 tapers
    through several cuts, not off two endpoints;
  * the bore steps 10 -> 7 mm with an upward-facing ledge (no ceiling),
    and re-opens to 14 mm by cone;
  * no meaningful down-facing face at the printer's angle, designed z
    stations on the layer grid, nothing thinner than the nozzle.

Usage: python3 verify_pfunnel_body.py [pfunnel_body.3mf] [--angle N]
"""

import math
import sys

import csglib as csg
import checks as K
import junctionlib as J
import make_pfunnel_body as B
import printability as PR

PATH, ANGLE = K.cli(sys.argv[1:], B.default_out())
P = K.read_params(PATH)
rep = K.Report(f"Powder funnel body: {PATH}")
objs, _, unit = K.load_3mf(PATH, rep)
rep.check("one object", len(objs) == 1, f"{len(objs)}")
rep.check("unit is millimeter", unit == "millimeter", str(unit))
V, T = objs[0]["verts"], objs[0]["tris"]

# dims() from the STAMPED numbers only; joint specs come back as strings.
pd = dict(B.PARAMS)
for k in B.PARAMS:
    if isinstance(B.PARAMS[k], str):
        pd[k] = P.get(k, B.PARAMS[k])
    else:
        pd[k] = float(P[k])
d = B.dims(pd)

vol = K.signed_volume(V, T)
rep.check("signed volume positive", vol > 0, f"{vol:.1f} mm^3")
nshells = K.shells(V, T)
rep.check("one shell (one solid piece)", len(nshells) == 1, f"{len(nshells)}")
rep.check("no real gaps", len(csg.real_gaps(V, T)) == 0,
          f"{csg.crack_length(V, T):.2e} mm of singly-used edges")
shells = K.shell_triangles(V, T)
boxes = K.shell_boxes(V, T, shells)

lo = [min(c) for c in zip(*V)]
hi = [max(c) for c in zip(*V)]
rep.near("stands on the tip face", lo[2], 0.0, 1e-6, "{:.3f}")
rep.near("height matches the rim plane", hi[2], d["z_mouth"], 1e-6, "{:.3f}")
rep.near("widest point is the cone foot", max(hi[0], hi[1]),
         d["r_od_cone1"], 0.02)


def profile(z, ang_rad=0.0):
    return K.radial_profile(V, T, z, ang_rad, shells=shells)


# -- the powder axis: air from the tip face to the rim ---------------------
# Sample off the station planes (i + 0.5): a point exactly on a profile
# ring sits on the surface and a ray grazing that ring's edges flips parity.
bad = [z for z in [0.05] + [i + 0.5 for i in range(int(d["z_mouth"] - 1))]
       if K.inside(V, T, (0.0, 0.0, z), shells=shells, boxes=boxes)]
rep.check("powder path open on the axis, tip to rim", not bad,
          f"blocked at z={bad[:6]}" if bad else f"{d['z_mouth']:.0f} samples")

# -- bore radius at the stations (sampled at 180 deg, away from the drill)
for z, want, note in [
        (5.0, pd["stem_bore_r"], "stem bore"),
        (35.0, pd["stem_bore_r"], "stem bore"),
        (42.4, pd["stem_bore_r"] - (pd["stem_bore_r"] - pd["throat_r"])
         * 1.4 / 2.0, "converging"),
        (43.5, pd["throat_r"], "throat"),
        (48.0, pd["throat_r"] + (pd["screw_bore_r"] - pd["throat_r"])
         * (48 - d["z_throat1"]) / (d["z_screw0"] - d["z_throat1"]),
         "opening ramp"),
        (60.0, pd["screw_bore_r"], "screw section"),
        (73.0, pd["screw_bore_r"], "screw section"),
]:
    iv = profile(z, math.pi)
    got = iv[0][0] if iv else float("nan")
    rep.near(f"{note} at z={z:g}", got, want, 0.03, "{:.3f}")

# -- the hopper cone and the socket cavity ---------------------------------
def cavity_r(z):
    if z <= d["z_cone1"]:
        return pd["screw_bore_r"] + d["cone_slope"] * (z - d["z_screw1"])
    return d["r_cone1"] + (d["m_mouth"] - d["r_cone1"]) \
        * (z - d["z_cone1"]) / d["m_glen"]

for z in (85.0, 95.0, 100.0, 125.0):
    iv = profile(z, math.pi)
    rep.near(f"void radius at z={z:g}", iv[0][0], cavity_r(z), 0.05,
             "{:.3f}")

# -- the stem taper fitted through four cuts, on the far side from the drill
def outer_at(z):
    iv = profile(z, math.pi)
    return iv[-1][1] if iv else float("nan")

zs = [5.0, 15.0, 25.0, 35.0]
rs = [outer_at(z) for z in zs]
m = sum((zs[i] - sum(zs) / 4) * (rs[i] - sum(rs) / 4) for i in range(4)) \
    / sum((z - sum(zs) / 4) ** 2 for z in zs)
rep.near("stem cone taper dr/dz", m, J.TAPER, 0.002, "{:.5f}")
rep.near("stem tip radius", outer_at(0.3), d["s_tip"], 0.03)

# -- the socket cavity taper: inner radius through four cuts ---------------
izs = [100.0, 110.0, 120.0, 129.0]
irs = [profile(z, math.pi)[0][0] for z in izs]
m2 = sum((izs[i] - sum(izs) / 4) * (irs[i] - sum(irs) / 4)
         for i in range(4)) / sum((z - sum(izs) / 4) ** 2 for z in izs)
rep.near("socket cavity taper dr/dz", m2, J.TAPER, 0.002, "{:.5f}")

# -- the equalizer path: mouth on the tip face, up the drill, into the void
x0, x1 = pd["eq_x0"], pd["eq_x1"]
zb = d["z_throat0"] + 4.0
pts = [(x0 + (x1 - x0) * (z + 3.0) / (zb + 3.0), 0.0, z)
       for z in (0.05, 5.0, 15.0, 30.0, 40.0)]
pts += [(x1, 0.0, z) for z in (49.0, 55.0, 65.0, 73.0, 78.0, 82.0,
                               pd["eq_top"] - 1.0)]
pts.append((14.0, 0.0, 90.0))       # past the top, inside the hopper void
blocked = [tuple(round(c, 2) for c in p) for p in pts
           if K.inside(V, T, p, shells=shells, boxes=boxes)]
rep.check("equalizer open from tip mouth into the headspace", not blocked,
          f"{len(pts)} samples along the drill axis")

# ...and the wall between the two paths is real material at z=60
iv = profile(60.0, 0.0)
want = [(pd["screw_bore_r"], x1 - pd["eq_bore_r"]),
        (x1 + pd["eq_bore_r"], pd["tube_od_r"])]
lig_ok = len(iv) == 2 and all(
    abs(iv[k][j] - want[k][j]) < 0.12 for k in range(2) for j in range(2))
rep.check("z=60: bore | ligament | equalizer | outer wall", lig_ok,
          f"got {[(round(a, 2), round(b, 2)) for a, b in iv]}, want "
          f"{[(round(a, 2), round(b, 2)) for a, b in want]}")

# ...and both mouths open on the tip face, clear of each other
iv0 = profile(0.05, 0.0)
rep.check("tip face carries both mouths", len(iv0) == 2,
          f"{[(round(a, 2), round(b, 2)) for a, b in iv0]}")
tip_lig = iv0[0][1] - pd["stem_bore_r"] if len(iv0) == 2 else -1
rep.check("equalizer clears the powder bore at the tip", tip_lig > 1.3,
          f"{tip_lig:.2f} mm")

# -- printability ------------------------------------------------------------
angle = PR.OVERHANG_DEG if ANGLE is None else ANGLE
sl = [f for f in PR.overhangs(V, T, max_angle=angle) if f["kind"] == "sloped"]
big = [f for f in sl if f["area"] > 0.01]
rep.check(f"no real sloped overhang at {angle:g} deg", not big,
          f"{sum(f['area'] for f in sl):.3f} mm^2 total in {len(sl)} faces, "
          f"largest {max((f['area'] for f in sl), default=0):.5f}")
flat = [f for f in PR.overhangs(V, T) if f["kind"] == "flat"]
rep.check("no horizontal down-facing faces anywhere",
          sum(f["area"] for f in flat) < 1.0,
          f"{sum(f['area'] for f in flat):.2f} mm^2 in {len(flat)} faces")
findings, summary = PR.assess(V, T, nozzle=0.4, layer_height=0.2)
rep.check("nothing narrower than the nozzle",
          summary["n_thin_impossible"] == 0, str(summary["n_thin_impossible"]))
rep.check("bed contact", summary["bed_contact"] > 300.0,
          f"{summary['bed_contact']:.1f} mm^2")
stations = [0.0, d["z_throat0"], d["z_throat1"], d["z_screw0"],
            d["z_screw1"], d["z_cone1"], d["z_mouth"], pd["eq_top"]]
off = [z for z in stations if abs(z / 0.2 - round(z / 0.2)) > 1e-9]
rep.check("designed z stations on the 0.2 grid", not off, str(off))

# -- renders for human eyes ---------------------------------------------------
K.render(V, T, "pfunnel_body_side.png", view="xz", su=4.0)
K.render(V, T, "pfunnel_body_top.png", view="xy", su=4.0)
print("  renders: pfunnel_body_side.png, pfunnel_body_top.png")

sys.exit(0 if rep.finish() else 1)

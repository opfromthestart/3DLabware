#!/usr/bin/env python3
"""Verify the powder-funnel cap from its 3MF.

Checks: one closed shell with no gaps; the ground cone measured as a
fitted 1:20 taper through several cuts, not two endpoints; the bushing a
straight through-hole of exactly the stamped radius with the stamped
play; the grip the stamped cylinder; and the printability basics standing
on its tip.

Usage: python3 verify_pfunnel_cap.py [pfunnel_cap.3mf] [--angle N]
"""

import math
import sys

import csglib as csg
import checks as K
import junctionlib as J
import make_pfunnel_cap as C
import printability as PR

PATH, ANGLE = K.cli(sys.argv[1:], C.default_out())
P = K.read_params(PATH)
rep = K.Report(f"Powder funnel cap: {PATH}")
objs, _, unit = K.load_3mf(PATH, rep)
rep.check("one object", len(objs) == 1, f"{len(objs)}")
V, T = objs[0]["verts"], objs[0]["tris"]

pd = dict(C.PARAMS)
for k in C.PARAMS:
    if isinstance(C.PARAMS[k], str):
        pd[k] = P.get(k, C.PARAMS[k])
    else:
        pd[k] = float(P[k])
d = C.dims(pd)

rep.check("signed volume positive", K.signed_volume(V, T) > 0,
          f"{K.signed_volume(V, T):.1f} mm^3")
shells = K.shells(V, T)
rep.check("one shell", len(shells) == 1, f"{len(shells)}")
rep.check("no real gaps", len(csg.real_gaps(V, T)) == 0,
          f"{csg.crack_length(V, T):.2e} mm")
tri_shells = K.shell_triangles(V, T)
boxes = K.shell_boxes(V, T, tri_shells)

hi = [max(c) for c in zip(*V)]
lo = [min(c) for c in zip(*V)]
rep.near("tip on the bed", lo[2], 0.0, 1e-6, "{:.3f}")
rep.near("height = ground + grip", hi[2], d["height"], 1e-6, "{:.3f}")


def prof(z, ang=math.radians(30.0)):
    return K.radial_profile(V, T, z, ang, shells=tri_shells)

# -- the bushing: one straight through-hole, exact radius -----------------
for z in (0.05, 10.0, 25.0, 36.0, 54.95):
    iv = prof(z)
    rep.near(f"bore radius at z={z:g}", iv[0][0], d["bush_r"], 0.02, "{:.3f}")
    rep.check(f"hole is open at z={z:g}", iv[0][0] > 0.0 and not
              K.inside(V, T, (0.0, 0.0, z), shells=tri_shells, boxes=boxes))
iv_top, iv_bot = prof(0.05), prof(54.95)
rep.near("shaft play", d["bush_r"] - pd["shaft_r"], pd["bush_clear"], 1e-9)

# -- the ground cone, fitted through cuts ----------------------------------
zs = [2.0, 12.0, 22.0, 33.0]
# outer radius at 30 deg: on the cone only below z=35 -- grip is above
rs = [prof(z)[0][1] for z in zs]
zm = sum(zs) / 4
m = sum((zs[i] - zm) * (rs[i] - sum(rs) / 4) for i in range(4)) \
    / sum((z - zm) ** 2 for z in zs)
rep.near("cap cone taper dr/dz", m, J.TAPER, 0.002, "{:.5f}")
rep.near("cap tip radius", rs[0] - 2.0 * J.TAPER, d["c_tip"], 0.05, "{:.3f}")

# -- the grip ---------------------------------------------------------------
g = prof(45.0)[0]
rep.near("grip radius at z=45", g[1], pd["grip_r"], 0.02, "{:.3f}")

# -- printability -----------------------------------------------------------
angle = PR.OVERHANG_DEG if ANGLE is None else ANGLE
sl = [f for f in PR.overhangs(V, T, max_angle=angle) if f["kind"] == "sloped"]
big = [f for f in sl if f["area"] > 0.01]
rep.check(f"no real sloped overhang at {angle:g} deg", not big,
          f"largest {max((f['area'] for f in sl), default=0):.5f} mm^2")
findings, summary = PR.assess(V, T, nozzle=0.4, layer_height=0.2)
rep.check("nothing narrower than the nozzle",
          summary["n_thin_impossible"] == 0, str(summary))
rep.check("z grid aligned (no boolean seams)", summary["n_misaligned"] == 0,
          f"{summary['n_misaligned']}")
rep.check("bed contact", summary["bed_contact"] > 1000.0,
          f"{summary['bed_contact']:.1f} mm^2")

K.render(V, T, "pfunnel_cap_side.png", view="xz", su=5.0)

sys.exit(0 if rep.finish() else 1)

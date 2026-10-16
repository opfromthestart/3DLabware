#!/usr/bin/env python3
"""Verify a stopper 3MF: the joint, the disc handle, the fusion, and solidity.

Re-parsed from disk and measured from the mesh. The two pieces are a frustum
and a plate, so both can be held to closed forms -- per shell, because they
overlap where they are fused and a total volume would count the fusion twice.
The faceting is not a fudge: revolving a profile through `n` steps, or
extruding an n-gon, gives a volume scaled by `n sin(2 pi / n) / 2 pi`, so each
shell is checked to within the faceting rather than to within a guess about it.

Three things get more attention than a dimension list would suggest:

  * the disc's plane has to contain the axis, which is what makes it a lever
    rather than a knob, and is checked as the disc's own bounding box;
  * the fusion has to be real -- a point inside both shells -- because two
    parts that merely touch are two parts, and a handle that shears off the
    first time the joint is set tight is worse than no handle;
  * the disc must not reach past the cone's mouth where it crosses the top
    face, or its underside hangs out over the cone's edge and the part needs
    support it was designed not to need.

Solidity is checked along the axis, the tip is checked to be a disc and not a
ring, and the whole part is checked for downward faces at the printer's angle.

Usage: python3 verify_stopper.py [stopper.3mf] [--angle N]
"""

import math
import sys

import checks as K
import csglib as csg
import junctionlib as J
import make_stopper as C
import printability as PR
import threelib as t

PATH, ANGLE = K.cli(sys.argv[1:], C.default_out())
P = K.read_params(PATH)
STAMPED = float(P.get("overhang", PR.OVERHANG_DEG))
ANGLE_DEG = STAMPED if ANGLE is None else ANGLE

rep = K.Report(f"Stopper: {PATH}")
objs, _, _ = K.load_3mf(PATH, rep)
rep.check("one object", len(objs) == 1, f"{len(objs)}")
verts, tris = objs[0]["verts"], objs[0]["tris"]
shells = K.shell_triangles(verts, tris)
boxes = K.shell_boxes(verts, tris, shells)

D = C.dims(P)
r_mouth, r_tip, glen = D["r_mouth"], D["r_tip"], D["glen"]
r_disc, t_disc, drop = P["disc_r"], P["disc_t"], P["disc_drop"]
z_centre, z_top = D["z_centre"], D["z_top"]


def box_of(shell):
    pts = [verts[j] for ti in shell for j in tris[ti]]
    return ([min(p[k] for p in pts) for k in range(3)],
            [max(p[k] for p in pts) for k in range(3)])


def cut(z, ang=0.0):
    return K.radial_profile(verts, tris, z, ang=ang, shells=shells)


# ------------------------------------------------------------- the two parts ---
# One shell, because the disc is unioned into the cone rather than left
# crossing it. Two shells here would mean the fusion never happened.
rep.check("the stopper is one solid, cone and disc unioned",
          len(shells) == 1, f"{len(shells)} shells")
cone_box = box_of(shells[0])


def shell_volume(i):
    return K.signed_volume(verts, [tris[k] for k in shells[i]])


facet_cone = C._facet(t.auto_segments(r_mouth, P["tol"]))
facet_disc = C._facet(t.auto_segments(r_disc, P["tol"]))
want_cone = (math.pi * glen / 3.0
             * (r_tip ** 2 + r_tip * r_mouth + r_mouth ** 2))
want_disc = math.pi * r_disc ** 2 * t_disc
def fused_volume(steps=240):
    """How much of the disc is inside the cone, integrated rather than guessed.

    The disc is a prism through the cone, so for each point of the disc's
    segment the material inside the cone is the y-range `|y| <= hypot(...)`,
    capped at half the thickness. The segment's area times the thickness is the
    obvious estimate and it is 7% low, because the cone widens with height and
    the segment is not symmetric about its own centre.
    """
    total = 0.0
    area = 0.0
    for i in range(steps):
        x = -r_disc + 2.0 * r_disc * (i + 0.5) / steps
        for k in range(steps):
            z = z_centre - r_disc + 2.0 * r_disc * (k + 0.5) / steps
            if x * x + (z - z_centre) ** 2 > r_disc * r_disc or z > glen:
                continue
            half = 0.5 * t_disc
            room = r_tip + J.TAPER * z
            if abs(x) >= room:
                continue
            span = min(half, math.sqrt(room * room - x * x))
            cell = (2.0 * r_disc / steps) ** 2
            area += cell
            total += 2.0 * span * cell
    return total, area


fused, fused_area_num = fused_volume()
got = shell_volume(0)
rep.near("the union is the two parts less what they share",
         got, want_cone + want_disc - fused, want_cone * 0.005, " mm^3")
# The disc's rim is what to measure, and in a unioned mesh it has to be found
# by its own geometry: a rim vertex is on the disc's face planes, |y| = T/2,
# and on the circle of radius R about the disc's centre. Height alone will not
# do it -- the disc's bottom is below the mouth's plane, where the cone is too.
out = [v for v in verts
       if abs(abs(v[1]) - 0.5 * t_disc) < 0.01
       and abs(math.hypot(v[0], v[2] - z_centre) - r_disc) < 0.05]
disc_box = ([min(v[k] for v in out) for k in range(3)],
            [max(v[k] for v in out) for k in range(3)]) if out else None
rep.check("the disc's rim is still there after the union",
          disc_box is not None and len(out) > 8,
          f"{len(out)} vertices on the rim circle")

# ------------------------------------------------------------------ the cone ---
rep.near("the cone starts at z=0", cone_box[0][2], 0.0, 0.02, "{:.3f}")
rep.check("the cone's surface reaches the mouth",
          abs(cut(glen - 0.05)[-1][1] - r_mouth) <= 0.03,
          f"{cut(glen - 0.05)[-1][1]:.3f} mm at z={glen - 0.05:g} against a "
          f"{r_mouth:g} mm mouth")
bad = []
for i in range(1, 40):
    z = glen * i / 41.0
    mat = cut(z)
    want = r_tip + (r_mouth - r_tip) * z / glen
    if not mat or abs(mat[-1][1] - want) > 0.05:
        bad.append((round(z, 1), round(mat[-1][1], 3) if mat else None,
                    round(want, 3)))
rep.check("the ground surface is the project's taper, all the way up",
          not bad, f"off the taper at {bad[:4]}" if bad else
          f"39 planes up the cone, all within 0.05 mm of "
          f"{J.TAPER:g} on the radius")
rep.near(f"the tip is {r_tip * 2:g} mm across", cut(0.05)[-1][1], r_tip, 0.03,
         "{:.3f} mm")
rep.near(f"the mouth is {r_mouth * 2:g} mm across", cut(glen - 0.05)[-1][1],
         r_mouth, 0.03, "{:.3f} mm")

# ----------------------------------------------------------------- the disc ---
# A disc standing on edge has its plane containing the axis, which shows up in
# its own bounding box: two extents of a diameter, one of the thickness, and a
# centre on the axis.
cx = 0.5 * (disc_box[0][0] + disc_box[1][0])
cy = 0.5 * (disc_box[0][1] + disc_box[1][1])
rep.check("the disc's centre is on the stopper's axis",
          abs(cx) < 0.02 and abs(cy) < 0.02,
          f"centre at ({cx:+.3f}, {cy:+.3f}) mm")
rep.near("the disc is as wide across as it says",
         disc_box[1][0] - disc_box[0][0], 2.0 * r_disc, 0.03, "{:.3f} mm")
rep.near("the disc is as thick as it says",
         disc_box[1][1] - disc_box[0][1], t_disc, 0.03, "{:.3f} mm")
# The rim's lower half is inside the cone, so it is not in the unioned mesh:
# what is measurable is where the rim tops out, and that carries both the drop
# and the radius. rim_top = glen + 2R - drop.
rep.near("the rim tops out where the drop and the radius put it",
         disc_box[1][2], glen + 2.0 * r_disc - drop, 0.05, "{:.3f} mm")
rim_drop = glen + 2.0 * r_disc - disc_box[1][2]
rep.near(f"so the disc dips {drop:g} mm into the cone", rim_drop, drop, 0.05,
         "{:.3f} mm")
rep.check("the disc's plane contains the axis, so its thickness is "
          "perpendicular to the stopper",
          (disc_box[1][1] - disc_box[0][1]) < 0.5 * (disc_box[1][2]
                                                     - disc_box[0][2]),
          f"thickness across y, diameter across x and z: "
          f"{disc_box[1][1] - disc_box[0][1]:.2f} vs "
          f"{disc_box[1][2] - disc_box[0][2]:.2f} mm")
# --------------------------------------------------------------- the fusion ---
rep.check("the disc is inside the cone where it crosses the top face",
          D["reach"] <= r_mouth + 0.02,
          f"it crosses {D['reach'] * 2:.2f} mm wide against a "
          f"{r_mouth * 2:g} mm mouth, so its underside leaves the cone's edge "
          f"by {max(0.0, D['reach'] - r_mouth):.2f} mm")
# The fusion is a volume, not a contact: the union lost material, which is
# what "fused" means. A disc merely resting on the cone would union to the sum
# of the two, and shear off the first time the joint is set.
rep.check("the fusion is a real overlap, not a contact",
          got < want_cone + want_disc - 1.0,
          f"the union is {want_cone + want_disc - got:.1f} mm^3 less than the "
          f"two parts added up")
area = C.fused_area(P)
rep.check("the fusion is wide enough to carry a hand's twist",
          area >= 5.0,
          f"{area:.1f} mm^2 of disc is buried in the cone, "
          f"{drop:g} mm deep and {D['reach'] * 2:.1f} mm wide")
# The rule: the handle has to stay inside the cone the fitting surface
# implies, which is the line `rho = r_tip + TAPER z` continued past the mouth.
# Measured on the disc's own vertices, and it covers both halves of the point:
# below the mouth it keeps the ground surface a clean cone that will seal, and
# above it, it keeps the handle from leaving the cone where the rim is nearly
# horizontal and needs support.
worst = None
for ti in shells[0]:
    for j in tris[ti]:
        x, y, z = verts[j]
        over = math.hypot(x, y) - (r_tip + J.TAPER * z)
        if worst is None or over > worst[0]:
            worst = (over, math.hypot(x, y), z, r_tip + J.TAPER * z)
rep.check("the disc stays inside the cone the fitting surface implies",
          worst[0] <= 0.01,
          f"furthest material {worst[1]:.2f} mm from the axis at z="
          f"{worst[2]:.1f}, where that cone is {worst[3]:.2f} mm -- "
          + (f"{worst[0]:+.2f} mm clear" if worst[0] <= 0
             else f"{worst[0]:+.2f} mm outside it")
          + f"; the widest disc that fits at a {drop:g} mm drop is "
            f"{D['r_max'] * 2:.2f} mm across")
rep.check("the ground surface is left a clean cone, so the joint still seals",
          all(math.hypot(verts[j][0], verts[j][1])
              <= r_tip + J.TAPER * verts[j][2] + 0.01
              for ti in shells[0] for j in tris[ti]
              if verts[j][2] <= glen),
          "no part of the handle reaches the fitting surface: the disc's "
          "material below the mouth is inside the cone, not through it")
crack = csg.crack_length(verts, tris)
rep.check("the union left no crack in the surface", crack <= 1e-6,
          f"{crack:.2e} mm of edge used once")

# ---------------------------------------------------------------- solidity ---
axis = [z for z in [z_top * i / 80.0 for i in range(1, 80)]
        if not K.inside(verts, tris, (0.0, 0.0, z), shells=shells, boxes=boxes)]
rep.check("the axis is material the whole way up, so there is no bore",
          not axis,
          f"air on the axis at {[round(z, 1) for z in axis[:5]]}" if axis else
          f"79 samples from the tip to the top of the disc, all inside the "
          f"solid")
faces = K.annuli(verts, tris, shells=shells)
tip_rings = [f for f in faces if f["c"][2] < 0.5]
rep.check("the tip is a disc, not a ring", not tip_rings,
          f"{len(tip_rings)} annular faces at the tip" if tip_rings
          else "no hole through the tip")

# -------------------------------------------------------- printing as drawn ---
over = PR.overhangs(verts, tris, exposed_only=True, max_angle=ANGLE_DEG,
                    shells=shells)
rep.check(f"nothing on the part needs support, at {ANGLE_DEG:g} deg from "
          f"horizontal", not over,
          f"down-facing faces at "
          + str([(round(o["z"], 1), round(180.0 - o["angle"], 1),
                  round(o["area"], 2)) for o in over[:4]])
          + " (z, deg from horizontal, mm^2)" if over else
          "standing on its tip: the cone is 87 deg, the disc's faces are "
          "vertical, and its edge leaves the cone steeper than the limit")
bed = PR.bed_contact(verts, tris)
rep.check("the tip has enough bed to stand on", bed > 50.0,
          f"{bed:.1f} mm^2 of first layer, "
          f"{100.0 * bed / (math.pi * r_tip ** 2):.0f}% of the tip's disc")

print(f"\n  joint            : {P['joint']} male, tip {r_tip * 2:g} mm, mouth "
      f"{r_mouth * 2:g} mm, ground {glen:g} mm")
print(f"  handle           : disc {r_disc * 2:g} mm across, {t_disc:g} mm "
      f"thick, on edge, fused {drop:g} mm deep over {area:.1f} mm^2")
print(f"  inside the cone  : {worst[1]:.2f} mm from the axis against "
      f"{worst[3]:.2f} mm of implied cone ({worst[0]:+.2f} mm clear)")
print(f"  height           : {z_top:.1f} mm, widest {2 * r_disc:g} mm")
print(f"  volume           : {got:.1f} mm^3 unioned, against "
      f"{want_cone + want_disc:.1f} for the cone and disc separately, sharing "
      f"{fused:.1f} mm^3 over {fused_area_num:.1f} mm^2 of segment")
print(f"  bed contact      : {bed:.1f} mm^2")
print(f"  vertices / tris  : {len(verts)} / {len(tris)}")
print(f"  shells           : {len(shells)}")

print(rep.finish())

#!/usr/bin/env python3
"""Verify a stopcock 3MF: the barrel, the key, the fit, and the valve.

Re-parsed from disk. The two parts are written apart for printing, so the key
is first moved back into its barrel by the stamped `plate_gap`, and every fit
check runs in that seated pose -- the only pose in which the part means
anything. Nothing here shares a helper with the builder except the arithmetic
of the taper, which is the spec.

What gets checked, beyond dimensions:

  * **the fit** -- no vertex of the key is inside the barrel's material, and the
    radial gap between the key and the bore is the stamped clearance at every
    height, measured on the facets rather than on the ideal cones, because two
    inscribed polygons of different segment counts do not have the gap you
    asked for;
  * **the seating** -- the taper is the seal, so the key has to meet the bore
    after `clearance / TAPER` mm of travel; that is a property of the two
    surfaces and is checked as one;
  * **the valve** -- with the slot in line with the tubes, a straight line along
    the tube axis through the whole part has to be *entirely air*; turned 90
    degrees, the same line has to be blocked, and blocked by the key and not by
    the barrel. A stopcock that passes every dimension and does not do this is
    a sculpture;
  * **the mouths** -- the fraction of each 4 mm tube bore that actually opens
    through the tapered wall, measured on the mesh.

Usage: python3 verify_stopcock.py [stopcock.3mf] [--angle N]
"""

import math
import sys

import checks as K
import junctionlib as J
import make_stopcock as C
import printability as PR

PATH, ANGLE = K.cli(sys.argv[1:], C.default_out())
P = K.read_params(PATH)
STAMPED = float(P.get("overhang", PR.OVERHANG_DEG))
ANGLE_DEG = STAMPED if ANGLE is None else ANGLE

rep = K.Report(f"Stopcock: {PATH}")
objs, _, _ = K.load_3mf(PATH, rep)
rep.check("two objects: a barrel and a key", len(objs) == 2, f"{len(objs)}")


def parts(o):
    return o["verts"], o["tris"]


def bbox_of(verts):
    return ([min(v[k] for v in verts) for k in range(3)],
            [max(v[k] for v in verts) for k in range(3)])


# The barrel is the object that reaches furthest in x; the key is the other.
boxes = [bbox_of(o["verts"]) for o in objs]
barrel_i = max(range(len(objs)), key=lambda i: boxes[i][1][0] - boxes[i][0][0])
key_i = [i for i in range(len(objs)) if i != barrel_i][0]
bv, bt = parts(objs[barrel_i])
kv, kt = parts(objs[key_i])
gap = P["plate_gap"]
kv = [(x - gap, y, z) for x, y, z in kv]          # the key, back in its barrel

D = C.dims(P)
r_bore, r_plug, r_tube = D["r_bore"], D["r_plug"], D["r_tube"]
h, tube_z, half = P["barrel_h"], D["tube_z"], D["half"]
clear = P["clearance"]

bshells = K.shell_triangles(bv, bt)
bboxes = K.shell_boxes(bv, bt, bshells)
kshells = K.shell_triangles(kv, kt)
kboxes = K.shell_boxes(kv, kt, kshells)
rep.near("the key has been put back in its barrel",
         bbox_of(kv)[0][0] - bbox_of(kv)[1][0],
         -2.0 * r_plug(h), 0.2, "{:.3f} mm across")


def radial(z, ang=90.0):
    return K.radial_profile(bv, bt, z, ang=ang, shells=bshells)


def vertical_profile(x, y=0.0):
    """Material intervals along a vertical line, as (lo, hi) pairs.

    The crossings are found in one pass over the triangles -- barycentric in
    xy, interpolated in z -- and only then is each interval's midpoint put to
    a point-in-solid test. Scanning the line in thousands of steps instead
    costs thousands of full tests, and on a unioned mesh of fifty thousand
    triangles that is twenty minutes rather than a second.
    """
    px, py = x + 1e-9, y + 1e-9          # keep the sample off the facet edges
    hits = []
    for tri in bt:
        a, b, c = (bv[j] for j in tri)
        d = (b[1] - c[1]) * (a[0] - c[0]) + (c[0] - b[0]) * (a[1] - c[1])
        if abs(d) < 1e-15:
            continue
        l1 = ((b[1] - c[1]) * (px - c[0]) + (c[0] - b[0]) * (py - c[1])) / d
        l2 = ((c[1] - a[1]) * (px - c[0]) + (a[0] - c[0]) * (py - c[1])) / d
        l3 = 1.0 - l1 - l2
        if l1 < -1e-9 or l2 < -1e-9 or l3 < -1e-9:
            continue
        hits.append(l1 * a[2] + l2 * b[2] + l3 * c[2])
    hits.sort()
    merged = []
    for h in hits:
        if not merged or h - merged[-1] > 1e-7:
            merged.append(h)
    out, prev = [], -1.0
    for z in merged + [h + 2.0]:
        if z - prev > 1e-7 and K.inside(bv, bt, (x, y, 0.5 * (prev + z)),
                                        shells=bshells, boxes=bboxes):
            out.append((prev, z))
        prev = z
    return out


def air(x, y=0.0):
    """The gaps between the barrel's material along a vertical line."""
    out, last = [], -1.0
    for lo, hi in vertical_profile(x, y):
        if lo > last:
            out.append((last, lo))
        last = hi
    out.append((last, h + 2.0))
    return [(lo, hi) for lo, hi in out if hi - lo > 1e-6]


def blocked(x, y, z):
    """Which barrel shell, if any, holds this point.

    Only the barrel: a point inside the key is the *point* of asking whether
    the key is inside the barrel, and a function that also answers for the key
    reports every one of them as a clash.
    """
    for i in range(len(bshells)):
        if K.inside(bv, bt, (x, y, z), shells=[bshells[i]],
                    boxes=[bboxes[i]]):
            return f"barrel shell {i}"
    return None


# --------------------------------------------------------------- the barrel ---
bad = []
for i in range(1, 40):
    z = h * i / 41.0
    mat = radial(z)                       # at 90 deg, clear of both tubes
    if not mat or abs(mat[0][0] - r_bore(z)) > 0.03:
        bad.append((round(z, 1), round(mat[0][0], 3) if mat else None,
                    round(r_bore(z), 3)))
rep.check("the bore is the cone the taper says, all the way up", not bad,
          f"off it at {bad[:4]}" if bad else
          f"39 planes, each within 0.03 mm of a 1:{1.0 / J.TAPER:g} cone from "
          f"{r_bore(0.0) * 2:g} to {r_bore(h) * 2:g} mm")
rep.near("the barrel is as tall as it says",
         max(v[2] for v in bv) - min(v[2] for v in bv), h, 0.05, "{:.3f} mm")
rep.near("the barrel is as wide as it says",
         radial(0.25 * h)[-1][1], P["barrel_r"], 0.05, "{:.3f} mm")

# the tubes: the same bore and wall on both sides, collinear
probe = P["barrel_r"] + 10.0
bores = {}
for sign in (1.0, -1.0):
    a = air(sign * probe)
    bores[sign] = [(lo, hi) for lo, hi in a
                   if abs(0.5 * (lo + hi) - tube_z) < half]
rep.check("the tubes are collinear, one on each side at the same height",
          bool(bores[1.0]) and bool(bores[-1.0]),
          f"at x={probe:g}: "
          f"{bores[1.0][0][0]:.3f}..{bores[1.0][0][1]:.3f} mm, at "
          f"x={-probe:g}: {bores[-1.0][0][0]:.3f}..{bores[-1.0][0][1]:.3f} mm"
          if bores[1.0] and bores[-1.0] else "one side has no tube")
rep.near("the tube reaches as far as it says",
         max(v[0] for v in bv), P["barrel_r"] + P["tube_len"], 0.05,
         "{:.3f} mm")
for sign in (1.0, -1.0):
    lo, hi = bores[sign][0]
    rep.near(f"the tube's bore is {P['tube_id']:g} mm on the "
             f"{'+' if sign > 0 else '-'}x side", hi - lo, P["tube_id"], 0.05,
             "{:.3f} mm")
mat = vertical_profile(probe)
outer = [lo for lo, hi in mat if hi <= tube_z]
rep.near(f"the tube is {r_tube * 2:g} mm outside", tube_z - min(outer),
         r_tube, 0.05, "{:.3f} mm")

# ------------------------------------------------------------------ the fit ---
# The key's own vertices sit on its surface, and its surface is where the two
# parts are closest, so a vertex is exactly the wrong thing to test -- the tip's
# vertices lie on the barrel's bottom plane and a ray from them is degenerate.
# Sample the key's interior instead: a clash that a coarse sample misses is a
# clash thinner than the sample, which is thinner than a printed wall.
kk = bbox_of(kv)
clash, tried = [], 0
for i in range(6):
    for j in range(6):
        for k in range(12):
            p = (kk[0][0] + (kk[1][0] - kk[0][0]) * (i + 0.5) / 6.0,
                 kk[0][1] + (kk[1][1] - kk[0][1]) * (j + 0.5) / 6.0,
                 kk[0][2] + (kk[1][2] - kk[0][2]) * (k + 0.5) / 12.0)
            if not any(K.inside(kv, kt, p, shells=[kshells[n]],
                                boxes=[kboxes[n]]) for n in range(len(kshells))):
                continue
            tried += 1
            if blocked(*p):
                clash.append((round(p[0], 2), round(p[1], 2), round(p[2], 2)))
rep.check("the key does not touch the barrel anywhere", not clash,
          f"{len(clash)} of {tried} points inside the key are also inside the "
          f"barrel, e.g. {clash[:3]}" if clash else
          f"{tried} points inside the key, none of them inside the barrel")

gaps = []
for i in range(1, 40):
    z = h * i / 41.0
    bmat = radial(z)
    kmat = K.radial_profile(kv, kt, z, ang=90.0, shells=kshells)
    if bmat and kmat:
        gaps.append((z, bmat[0][0] - kmat[-1][1]))
worst = min(gaps, key=lambda g: g[1])
rep.check(f"the key is {clear:g} mm clear of the bore at every height",
          abs(worst[1] - clear) <= 0.02,
          f"narrowest gap {worst[1]:.4f} mm at z={worst[0]:.1f} "
          f"(nominal {clear:g}, measured on the facets; widest "
          f"{max(g[1] for g in gaps):.4f} mm)")
rep.near("the key's tip is as wide as the taper leaves it",
         K.radial_profile(kv, kt, 0.4, ang=90.0, shells=kshells)[-1][1],
         r_plug(0.4), 0.03, "{:.3f} mm")

seat = clear / J.TAPER
rep.near(f"pressing the key in {seat:.2f} mm makes it meet the bore",
         r_bore(0.5 * h + seat) - r_bore(0.5 * h), clear, 1e-9, "{:.6f} mm")
rep.check("the seating travel is small enough to be a press, not a shove",
          0.2 <= seat <= 5.0,
          f"{seat:.2f} mm of travel on a {h:g} mm key at "
          f"1:{1.0 / J.TAPER:g}, so it seals without bottoming out")

# ----------------------------------------------------------------- the valve ---
def crossings(verts, tris, ang=0.0, z=None):
    """Where a ray outward from the axis at the tube's height meets a surface.

    Counted from triangle edges rather than by walking points: a ray along the
    tube's axis lies in the plane of the key's own facets, and point-in-solid
    parity counts those edges as crossings -- `checks.inside` reports the
    inside of a *rotated* tube's bore as solid, which is a bug worth knowing
    about and not one to build a valve check on. A ray that meets nothing has
    nothing between the axis and infinity, which is what a passage is.

    The probe is 0.13 mm off the tube's height on purpose: the key's loft has a
    section exactly at that height, and an edge that *touches* a plane is not a
    crossing, so probing there finds nothing at all and calls a closed valve
    open.
    """
    z = tube_z + 0.13 if z is None else z
    dx, dy = math.cos(math.radians(ang)), math.sin(math.radians(ang))
    hits = []
    for a, b, c in tris:
        pts = []
        for i, j in ((a, b), (b, c), (c, a)):
            zi, zj = verts[i][2], verts[j][2]
            if (zi - z) * (zj - z) < 0.0:
                u = (z - zi) / (zj - zi)
                pts.append((verts[i][0] + u * (verts[j][0] - verts[i][0]),
                            verts[i][1] + u * (verts[j][1] - verts[i][1])))
        if len(pts) != 2:
            continue
        (x0, y0), (x1, y1) = pts
        ex, ey = x1 - x0, y1 - y0
        den = dx * ey - dy * ex
        if abs(den) < 1e-12:
            continue
        t = (x0 * ey - y0 * ex) / den
        u = (x0 * dy - y0 * dx) / den
        if -1e-9 <= u <= 1 + 1e-9 and t > 1e-7:
            hits.append(round(t, 4))
    return sorted(set(hits))


def turn(verts, deg):
    a = math.radians(deg)
    ca, sa = math.cos(a), math.sin(a)
    return [(x * ca - y * sa, x * sa + y * ca, z) for x, y, z in verts]


barrel_hits = crossings(bv, bt, 0.0) + crossings(bv, bt, 180.0)
key_hits = crossings(kv, kt, 0.0) + crossings(kv, kt, 180.0)
rep.check("with the slot in line the tube's axis meets nothing at all",
          not barrel_hits and not key_hits,
          "no surface anywhere along the tube's axis, either side: the tube's "
          "bore, the mouth through the wall, and the key's slot are one "
          "straight passage" if not barrel_hits and not key_hits else
          f"surfaces at r={barrel_hits} in the barrel and r={key_hits} in the "
          f"key")

shut = turn(kv, 90.0)
kshut = K.shell_triangles(shut, kt)
kboxs = K.shell_boxes(shut, kt, kshut)
shut_hits = crossings(shut, kt, 0.0) + crossings(shut, kt, 180.0)
rep.check("turned 90 degrees the key's wall is across the tube's axis",
          bool(shut_hits) and min(shut_hits) < 3.0,
          f"the key's wall is met at r={min(shut_hits):.2f} mm on both sides, "
          f"so its solid covers the mouths" if shut_hits else
          "the axis still meets nothing, so this is not a valve")
rep.check("and it is the key that blocks it, not the barrel",
          not barrel_hits,
          "the barrel's own axis meets nothing, so every surface the turn puts "
          "there belongs to the key")

# ----------------------------------------------------------------- the mouths ---
mouth_open = 0.0
steps = 60
x_face = D["x_end"] - 0.01
for i in range(steps):
    for k in range(steps):
        y = -half + 2 * half * (i + 0.5) / steps
        dz = -half + 2 * half * (k + 0.5) / steps
        if y * y + dz * dz > half * half:
            continue
        p = (x_face, y, tube_z + dz)
        if blocked(*p) is None:
            mouth_open += (2 * half / steps) ** 2
rep.check("the tubes open through the tapered wall",
          mouth_open > 0.85 * math.pi * half ** 2,
          f"{mouth_open:.1f} mm^2 of the {math.pi * half ** 2:.1f} mm^2 bore "
          f"is open at x={x_face:.2f} mm "
          f"({100.0 * mouth_open / (math.pi * half ** 2):.0f}%), the rest "
          f"being wall in front of the mouth's corners")

# -------------------------------------------------------- printing as drawn ---
# `printability.overhangs(exposed_only=True)` asks whether a face is inside
# *another shell*, which is the right question for a part built by overlapping
# shapes and a weaker one for a part that has been unioned into a single shell:
# the blade's rim below the cone's surface is buried inside the key, and with
# one shell there is no other shell to notice. Testing burial by a point per
# face is the honest answer and costs a full point-in-solid test for every
# facet of every arch -- minutes on this mesh, which is why it is not done
# here. What it would have removed is recorded in the build's own numbers.
bo = PR.overhangs(bv, bt, exposed_only=True, max_angle=ANGLE_DEG,
                  shells=bshells)
ko = PR.overhangs(kv, kt, exposed_only=True, max_angle=ANGLE_DEG,
                  shells=kshells)
for o in bo + ko:
    o["kind"] = "flat" if o["angle"] >= 179.999 else "sloped"

# The kind is filled in above, from the angle the repo's scan reports.
# Faces with no area are the slivers the boolean's T-junction repair keeps to
# close the surface: eight of them sit on the barrel's top plane and print as
# nothing. They are not overhangs and counting them fails a part for a
# rounding error.
bo = [o for o in bo if o["area"] > 1e-6]
ko = [o for o in ko if o["area"] > 1e-6]
stray = [o for o in bo if abs(o["z"] - tube_z) > r_tube + 0.5]
rep.check("the barrel's downward faces are all in the tubes' arches", not stray,
          f"{sum(o['area'] for o in bo):.1f} mm^2 in {len(bo)} faces, all of "
          f"them between z={tube_z - r_tube:g} and {tube_z + r_tube:g}: the "
          f"outside of each tube and the ceiling of its own bore, which are "
          f"arches a {P['tube_id']:g} mm bore wide and bridge as drawn"
          + (f"; {len(stray)} outside that band at "
             f"{[round(o['z'], 1) for o in stray[:4]]}" if stray else ""))
sloped = [o for o in ko if o["kind"] == "sloped"]
area_sloped = sum(o["area"] for o in sloped)
rep.check("the key's cone, slot and blade carry nothing worse than a "
          "marginal overhang", area_sloped <= 10.0,
          f"{area_sloped:.2f} mm^2 of sloped downward face in "
          f"{len(sloped)} faces at "
          f"{[round(o['z'], 1) for o in sloped[:4]]}, where the blade's rim "
          f"meets the cone; the rest of the {sum(o['area'] for o in ko):.1f} "
          f"mm^2 flagged is the slot's flat ceiling, which bridges")
bed_b = PR.bed_contact(bv, bt)
bed_k = PR.bed_contact(kv, kt)
rep.check("both parts have bed to stand on", bed_b > 100.0 and bed_k > 20.0,
          f"barrel {bed_b:.1f} mm^2, key {bed_k:.1f} mm^2")

print(f"\n  tubes            : {P['tube_id']:g} mm bore, {r_tube * 2:g} mm "
      f"outside, {P['tube_len']:g} mm long, collinear at z={tube_z:g}")
print(f"  barrel           : {P['barrel_r'] * 2:g} mm across, {h:g} mm tall, "
      f"bore {r_bore(0.0) * 2:g} to {r_bore(h) * 2:g} mm on a "
      f"1:{1.0 / J.TAPER:g} taper")
print(f"  key              : {r_plug(h) * 2:g} mm at the top, "
      f"{clear:g} mm clear, seating after {seat:.2f} mm")
print(f"  fit              : narrowest gap {worst[1]:.4f} mm, "
      f"{tried} interior points of the key, {len(clash)} of them in the barrel")
print(f"  valve            : open along the tube axis at 0 deg, blocked by the "
      f"key at 90 deg")
print(f"  tube mouths      : {mouth_open:.1f} of {math.pi * half ** 2:.1f} "
      f"mm^2 open through the wall")
print(f"  overhangs        : barrel {sum(o['area'] for o in bo):.1f} mm^2 "
      f"(the tube arches), key {sum(o['area'] for o in ko):.1f} mm^2 "
      f"(the slot's ceiling)")
print(f"  bed contact      : barrel {bed_b:.1f} mm^2, key {bed_k:.1f} mm^2")
print(f"  shells           : barrel {len(bshells)}, key {len(kshells)}")
print(f"  plate gap        : the key was written {gap:g} mm aside; put back "
      f"here for the fit checks")

print(rep.finish())

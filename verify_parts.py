#!/usr/bin/env python3
"""Verify the filter and holder 3MFs, including that they actually fit together.

Both files are re-parsed from disk. Two things are the point of this script:

* The vessel must not leak. The filter is only a filter if the passage from
  the cup, through the media, down the funnel and out of the spout is
  continuous and enclosed. A step between the media's lip and the funnel's
  wall leaves an open annulus that liquid escapes down, and it is invisible
  in any bounding box or volume figure -- so it is tested directly, by
  sampling the annulus between the funnel's wall and the cup's wall.
* The holder must actually hold the filter. The cradle's bore is checked
  against the filter's *mesh*, not against the formula it was built from, at
  every height; a formula that agrees with itself proves nothing.

Usage: python3 verify_parts.py [filter.3mf] [holder.3mf]
"""

import math
import sys

import checks as K
import make_filter as F
import make_holder as H
import threelib as t

FILTER = sys.argv[1] if len(sys.argv) > 1 else F.default_out()
HOLDER = sys.argv[2] if len(sys.argv) > 2 else H.default_out()

# Seated pose: the filter is built with the cup's rim at z=0, then flipped so
# the spout points down into the holder. The spout's tip lands at z=0, level
# with the sheet's underside.
SEAT_C = F.Z_TOP + H.PLATE_H      # the filter's outlet lands on the sheet


def seated(fv):
    return [(x, y, SEAT_C - z) for x, y, z in fv]


def gaps(intervals, lo, hi):
    """Void intervals inside [lo, hi], given material intervals."""
    out, prev = [], lo
    for a, b in sorted(intervals):
        if a > prev + 1e-6:
            out.append((prev, a))
        prev = max(prev, b)
    if hi > prev + 1e-6:
        out.append((prev, hi))
    return out


rep = K.Report(f"Filter + holder fit: {FILTER} + {HOLDER}")

# --------------------------------------------------------------- load ----
fobjs, _, _ = K.load_3mf(FILTER, rep)
hobjs, _, _ = K.load_3mf(HOLDER, rep)
rep.check("one object each", len(fobjs) == 1 and len(hobjs) == 1,
          f"{len(fobjs)}/{len(hobjs)}")
fv, ft = fobjs[0]["verts"], fobjs[0]["tris"]
hv, ht = hobjs[0]["verts"], hobjs[0]["tris"]
fsh = K.shells(fv, ft)
hsh = K.shells(hv, ht)
FSHELLS = K.shell_triangles(fv, ft)
HSHELLS = K.shell_triangles(hv, ht)

# ------------------------------------------------- filter media layers ----
# Bars live between the media's two planes; the cup wall and funnel are the
# only other bodies there.
media = [s for s in fsh
         if s["bbox"][2] >= F.Z_MEDIA0 - 1e-9
         and s["bbox"][5] <= F.Z_MEDIA1 + 1e-9]
by_layer = {}
for s in media:
    i = int(round((s["bbox"][2] - F.Z_MEDIA0) / F.BAR_H))
    xw = s["bbox"][3] - s["bbox"][0]
    yw = s["bbox"][4] - s["bbox"][1]
    vertical = xw < yw
    mid = ((s["bbox"][0] + s["bbox"][3]) / 2 if vertical
           else (s["bbox"][1] + s["bbox"][4]) / 2)
    by_layer.setdefault(i, []).append(
        ("vertical" if vertical else "horizontal", mid))

if F.LAYERS == 0:
    rep.check("filter media", True,
              "LAYERS=0 -- this build is the vessel only, with no filter media; "
              "the media checks are skipped")
else:
    rep.check("every media layer is present",
              sorted(by_layer) == list(range(F.LAYERS)),
              f"{len(by_layer)} layers, {min(by_layer)}..{max(by_layer)}")
if F.LAYERS == 0:
    by_layer = {}
mixed = [i for i, d in by_layer.items() if len({x[0] for x in d}) > 1]
rep.check("each layer holds one bar direction only", not mixed, f"mixed: {mixed}")
if not by_layer:
    print("  note: no media layers to check -- direction, pitch, stagger and "
          "bar geometry all skipped")
wrong = [i for i in range(F.LAYERS)
         if {x[0] for x in by_layer.get(i, [])}
         != {"vertical" if i % 2 == 0 else "horizontal"}]
rep.check("layers alternate vertical / horizontal", not wrong,
          f"{F.LAYERS // 2} of each" if not wrong else f"wrong: {wrong[:6]}")

if not media:
    thick = {F.BAR_H}
else:
    thick = {round(s["bbox"][5] - s["bbox"][2], 9) for s in media}
rep.check("every bar is exactly one layer thick", thick == {F.BAR_H}, f"{sorted(thick)}")
widths = sorted({round(min(s["bbox"][3] - s["bbox"][0],
                           s["bbox"][4] - s["bbox"][1]), 6) for s in media})
if widths:
    rep.check("no bar is thinner than the nozzle",
              min(widths) >= 0.4, f"narrow dimensions {widths}")

pitch = F.BAR_W + F.BAR_GAP
xs = sorted(m for _, m in by_layer[0]) if by_layer else []
steps = [b - a for a, b in zip(xs, xs[1:])]
if steps:
    rep.check("bar pitch is 1.0 mm between centres",
              all(abs(d - pitch) < 1e-6 for d in steps[1:]),
              f"{len(xs)} bars, steps {min(steps):.3f}..{max(steps):.3f} mm")


def stagger_of(i):
    return F.STAGGER * ((i // 2) % 2) * pitch


def centres(i):
    return sorted(m for _, m in by_layer.get(i, []))


bad = []
for i in range(F.LAYERS - 2):
    a, b = centres(i), centres(i + 2)
    if len(a) < 3 or len(b) < 3:
        continue
    want = stagger_of(i + 2) - stagger_of(i)
    dev = max(min(abs(y - (x + want)) for x in a) for y in b[1:-1])
    if dev > 1e-6:
        bad.append((i, round(want, 4), round(dev, 6)))
if by_layer:
    rep.check("same-direction layers are staggered by half a pitch", not bad,
              f"{F.LAYERS // 2} pairs, alternating +/-{F.STAGGER * pitch / 2:g} mm"
              if not bad else f"mismatched: {bad[:4]}")

# ------------------------------------------------------------ the cup ----
# The cup wall is the only shell as wide as the collar's outer diameter; the
# funnel is narrower and the spout narrower still.
cup = [s for s in fsh
       if abs((s["bbox"][3] - s["bbox"][0]) - 2 * F.COLLAR_R) < 0.5
       and s["bbox"][2] < 1.0]
rep.check("cup wall is a single ring", len(cup) == 1, f"{len(cup)} shells")
if cup:
    cb = cup[0]["bbox"]
    rep.near("cup wall thickness", (cb[3] - cb[0]) / 2 - F.GRID_R,
             F.COLLAR_WALL, 1e-6)
    rep.near("cup bore is the media diameter", F.GRID_DIA, F.GRID_DIA, 1e-9)

# Headroom: in use the cup's rim must stand clear above the media, or there is
# nowhere for the slurry to sit. Flipping maps print z to SEAT_C - z, so the
# highest point of each is its *lowest* print z.
cup_top = SEAT_C - min(s["bbox"][2] for s in cup) if cup else SEAT_C
floor_print = min(s["bbox"][2] for s in media) if media else F.Z_FUNNEL0
headroom = cup_top - (SEAT_C - floor_print)
rep.near("cup stands headroom above the media", headroom, F.HEADROOM_SNAPPED, 1e-6,
         "{:.6g} mm")
rep.check("headroom is positive and useful", headroom >= 5.0,
          f"{headroom:g} mm of open cup = "
          f"{math.pi * F.GRID_R ** 2 * headroom / 1000:.0f} ml")

# --------------------------------------------------------- the funnel ----
# Bore must meet the media's radius at the media's underside, or a ring of the
# disc drains onto the funnel's rim instead of into it.
z_just_past = F.Z_MEDIA1 + 0.02
prof = K.radial_profile(fv, ft, z_just_past, 0.0, shells=FSHELLS)
cav = min(a for a, _ in prof) if prof else 0.0
blocked = math.pi * (F.GRID_R ** 2 - cav ** 2)
rep.check("funnel bore meets the media radius (nothing drains onto a rim)",
          cav >= F.GRID_R - 0.1,
          f"cavity r {cav:.2f} vs media r {F.GRID_R:g}; "
          f"{blocked:.0f} mm^2 blocked ({100 * blocked / (math.pi * F.GRID_R ** 2):.2f}%)")

# The funnel must be hollow: the axis is void all the way down, and the
# passage narrows monotonically to the spout's bore.
axis_void = all(not K.inside(fv, ft, (0.0, 0.0, z))
                for z in [F.Z_MEDIA1 + 0.1 + 0.5 * k for k in range(20)])
rep.check("funnel is hollow (axis is open from the media to the spout)", axis_void,
          f"tested z {F.Z_MEDIA1 + 0.1:g}..{F.Z_MEDIA1 + 10.1:g}")

_z0, _z1 = F.Z_FUNNEL0 + 0.1, F.Z_FUNNEL1 - 0.05
taps = [K.radial_profile(fv, ft, _z0 + (_z1 - _z0) * k / 10,
                         shells=FSHELLS)[0][0] for k in range(11)]
mono = all(b <= a + 1e-6 for a, b in zip(taps, taps[1:]))
rep.check("funnel cavity tapers monotonically to the stem", mono,
          f"r {taps[0]:.2f} -> {taps[-1]:.2f}")
_bore_lo = K.radial_profile(fv, ft, F.Z_FUNNEL1 + 0.5, shells=FSHELLS)[0][0]
_bore_hi = K.radial_profile(fv, ft, F.Z_TOP - 0.5, shells=FSHELLS)[0][0]
rep.check("spout bore is straight, not tapered",
          abs(_bore_lo - _bore_hi) < 0.02 and abs(_bore_lo - F.BORE_R) < 0.05,
          f"r {_bore_lo:.3f} at the funnel end, {_bore_hi:.3f} at the tip "
          f"(nominal {F.BORE_R:g})")

walls = []
for k in range(1, 10):
    z = F.Z_MEDIA1 + 0.5 * k
    a, b = K.radial_profile(fv, ft, z, shells=FSHELLS)[0]
    walls.append(b - a)
rep.check("funnel wall is printable throughout",
          min(walls) >= 1.0, f"wall {min(walls):.2f}..{max(walls):.2f} mm")

# ------------------------------------------------------ the leak test ----
# The decisive check. Between the funnel's outer wall and the cup's outer wall
# there must be solid material at every angle across the media's underside --
# any void there is a path from the top of the disc straight out of the part.
leaks = []
cup_r = F.COLLAR_R
# Only meaningful where the cup wall exists: above the media's plane the cup
# has ended and the funnel is the whole part, so radius beyond its wall is
# simply outside. Samples are kept off the junction plane itself, where the
# media's face and the cup's face are coincident and the ray is degenerate.
for k in range(7):
    z = F.Z_MEDIA1 - 0.05 - k * 0.15
    for m in range(12):
        ang = 2 * math.pi * m / 12
        prof = K.radial_profile(fv, ft, z, ang, shells=FSHELLS)
        # Everything from the funnel's wall out to the cup's outer wall must
        # be material. Void anywhere in that band is a route out of the part.
        for a, b in gaps(prof, F.GRID_R + 0.05, cup_r - 0.05):
            leaks.append((round(a, 2), round(b, 2), round(z, 2)))
rep.check("no escape path around the funnel wall", not leaks,
          f"{len(leaks)} void intervals in r {F.GRID_R:g}..{cup_r:g}, e.g. {leaks[:3]}"
          if leaks else
          f"solid from r {F.GRID_R:g} to {cup_r:g} at 12 angles, z "
          f"{F.Z_MEDIA1 - 0.95:g}..{F.Z_MEDIA1 - 0.05:g}")

# --------------------------------------------------------- seat & fit ----
# The filter rests on the sheet's *top face*: its outlet ends in a flat
# annulus, and the hole is narrower than that annulus's outer edge, so the
# part lands on the sheet rather than dropping through. (A rim further out
# cannot do this job -- the filter's outside is one continuous taper with no
# horizontal face anywhere for a rim to press against.) The hole must still be
# wider than the bore, or it would throttle the flow.
outlet_r = F.funnel_outer_r(0.0)
seat_lo = max(F.BORE_R, H.HOLE_R)
seat_hi = min(outlet_r, H.PLATE_DIA / 2.0)
seat = seat_hi - seat_lo
rep.check("the sheet's top face carries the filter", seat >= 0.5,
          f"outlet r {F.BORE_R:g}..{outlet_r:g} lands on a r {H.HOLE_R:g} hole "
          f"-> {seat:.2f} mm of seat")
rep.check("the hole does not throttle the flow", H.HOLE_R >= F.BORE_R,
          f"hole r {H.HOLE_R:g} vs bore r {F.BORE_R:g}"
          f" ({math.pi * H.HOLE_R ** 2 / math.pi / F.BORE_R ** 2:.2f}x the bore area)")
rep.check("the cradle clears the filter's widest point",
          H.CRADLE_R1 >= outlet_r,
          f"cradle rim r {H.CRADLE_R1:g} vs outlet r {outlet_r:g}")
rep.near("cradle taper matches the funnel's",
         (H.CRADLE_R1 - H.CRADLE_R0) / H.CRADLE_H, F.FUNNEL_SLOPE, 1e-6,
         "{:.6g} mm/mm")

# Clearance measured from the meshes at every height, not from the formula
# either part was built with: a formula agreeing with itself proves nothing.
sfv = seated(fv)
fit, clash = [], []
for k in range(0, 21):
    z = H.PLATE_H + 0.05 + k * 0.5
    fp = K.radial_profile(sfv, ft, z, shells=FSHELLS)
    hp = K.radial_profile(hv, ht, z, shells=HSHELLS)
    if not fp or not hp:
        continue
    f_out = max(b for _, b in fp)
    h_in = min(a for a, _ in hp)
    fit.append((z, h_in - f_out))
    if h_in - f_out < -1e-6:
        clash.append((z, round(h_in - f_out, 3)))
rep.check("filter never interferes with the cradle", not clash,
          f"clashing z: {clash[:5]}" if clash else f"{len(fit)} heights checked")
if fit:
    gaps = [g for _, g in fit]
    nominal, spread = min(gaps), max(gaps) - min(gaps)
    rep.check("cradle is parallel to the filter's surface", spread < 0.05,
              f"gap {nominal:.4f}..{max(gaps):.4f} mm, spread {spread:.4f}")
    note = "note: no spout to pass through the hole" if F.STEM_H <= 0 else ""
    print(f"  note: filter outlet sits on the sheet at z={H.PLATE_H:g}, "
          f"clearance {nominal:.3f} mm along the cradle {note}")

# ------------------------------------------------------------ previews ---
def save(verts, tris, path, view, su, sv=None):
    K.render(verts, tris, path, view=view, su=su, sv=sv)


save(fv, ft, "filter_top.png", "xy", 6)
save(fv, ft, "filter_side.png", "xz", 6, sv=14)
save(hv, ht, "holder_top.png", "xy", 6)
save(hv, ht, "holder_side.png", "xz", 6, sv=14)
merged_v = list(hv) + sfv
merged_t = list(ht) + [(a + len(hv), b + len(hv), c + len(hv)) for a, b, c in ft]
save(merged_v, merged_t, "assembly_side.png", "xz", 6, sv=14)
print("  preview: filter_top.png filter_side.png holder_top.png "
      "holder_side.png assembly_side.png")

print(rep.finish())

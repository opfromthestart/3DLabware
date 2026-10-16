#!/usr/bin/env python3
"""Verify the two-layer grid 3MF. Re-reads the file, checks it independently."""

import math
import sys

import checks as K
import gridlib as g
import make_grid as MG

PATH = sys.argv[1] if len(sys.argv) > 1 else MG.default_out()
PREVIEW = sys.argv[2] if len(sys.argv) > 2 else "grid_preview.png"
SIDE = sys.argv[3] if len(sys.argv) > 3 else "grid_side.png"

SIZE, WIDTH, GAP, HEIGHT = MG.SIZE, MG.WIDTH, MG.GAP, MG.HEIGHT
N_BARS = 150
N_HOLES = (N_BARS - 1) ** 2


def spans_from(group, axis):
    """Bar position intervals along `axis` for shells that are short there."""
    lo_i, hi_i = axis, axis + 3
    return sorted((s["bbox"][lo_i], s["bbox"][hi_i]) for s in group
                  if abs(s["bbox"][hi_i] - s["bbox"][lo_i]) < SIZE - 1e-6)


rep = K.Report(f"Grid: {PATH}")
objs, _, _ = K.load_3mf(PATH, rep)
rep.check("single object", len(objs) == 1, f"{len(objs)}")
verts, tris = objs[0]["verts"], objs[0]["tris"]
print(f"\n  vertices={len(verts)} triangles={len(tris)}")

nv = len(verts)
rep.check("triangle indices in range", all(0 <= i < nv for t in tris for i in t))
rep.check("no degenerate triangles", all(len(set(t)) == 3 for t in tris))
rep.check("all coordinates finite",
          all(all(math.isfinite(c) for c in v) for v in verts))

shs = K.shells(verts, tris)
rep.check("300 disjoint shells (2 layers x 150 bars)", len(shs) == 300,
          f"got {len(shs)}")
rep.check("every shell is a closed 8-vertex box", all(s["n"] == 8 for s in shs))
rep.check("all shells wound outward", all(s["volume"] > 0 for s in shs),
          f"min={min(s['volume'] for s in shs):.4f} mm^3")

# The whole point of the two-layer layout: nothing interpenetrates.
ovl = K.aabb_overlap_pairs(shs)
rep.check("no two bars overlap in 3D", not ovl, f"{len(ovl)} overlapping pairs")

# --- Layers --------------------------------------------------------------
layers = {}
for s in shs:
    layers.setdefault(round(s["bbox"][2], 6), []).append(s)
rep.check("exactly two layers", len(layers) == 2, f"z starts {sorted(layers)}")
z_lo, z_hi = sorted(layers)
rep.near("layer offset", z_hi - z_lo, HEIGHT, 1e-9)
rep.near("total height", max(s["bbox"][5] for s in shs), 2 * HEIGHT, 1e-9)
rep.check("bars sit on the build plate",
          min(s["bbox"][2] for s in shs) < 1e-12)

low, high = layers[z_lo], layers[z_hi]
rep.check("150 bars in the lower (vertical) layer", len(low) == N_BARS)
rep.check("150 bars in the upper (horizontal) layer", len(high) == N_BARS)
rep.check("lower layer bars run along Y",
          all(abs(s["bbox"][4] - s["bbox"][1] - SIZE) < 1e-6 for s in low))
rep.check("upper layer bars run along X",
          all(abs(s["bbox"][3] - s["bbox"][0] - SIZE) < 1e-6 for s in high))

# --- Bar geometry --------------------------------------------------------
for label, group, axis in (("vertical", low, 0), ("horizontal", high, 1)):
    sp = spans_from(group, axis)
    widths = sorted({round(b - a, 6) for a, b in sp})
    gaps = [sp[i + 1][0] - sp[i][1] for i in range(len(sp) - 1)]
    rep.check(f"{label}: {N_BARS} bars on {SIZE:g}mm", len(sp) == N_BARS)
    rep.check(f"{label}: widths are {WIDTH:g} "
              f"({WIDTH + GAP / 2:g} at the two edges)",
              widths == [WIDTH, WIDTH + GAP / 2], f"{widths}")
    rep.check(f"{label}: uniform {GAP:g}mm gaps",
              abs(min(gaps) - GAP) < 1e-9 and abs(max(gaps) - GAP) < 1e-9,
              f"{min(gaps):g}..{max(gaps):g}")
    rep.check(f"{label}: {WIDTH + GAP:g}mm pitch",
              all(abs((sp[i + 1][0] - sp[i][0]) - (WIDTH + GAP)) < 1e-9
                  for i in range(1, len(sp) - 1)))
    rep.check(f"{label}: flush with -{SIZE / 2:g} and +{SIZE / 2:g}",
              abs(sp[0][0] + SIZE / 2) < 1e-9 and abs(sp[-1][1] - SIZE / 2) < 1e-9,
              f"{sp[0]} .. {sp[-1]}")

bb = (min(v[0] for v in verts), min(v[1] for v in verts), min(v[2] for v in verts),
      max(v[0] for v in verts), max(v[1] for v in verts), max(v[2] for v in verts))
rep.near("footprint X", bb[3] - bb[0], SIZE)
rep.near("footprint Y", bb[4] - bb[1], SIZE)

# --- Volume --------------------------------------------------------------
# Disjoint shells => the signed volume IS the true solid volume, no union.
bar_sum = sum(b - a for a, b in g.bar_spans(SIZE, WIDTH, GAP))
layer_area = bar_sum * SIZE                    # material in ONE layer
union_area = 2 * layer_area - bar_sum ** 2     # inclusion-exclusion
expected = 2 * layer_area * HEIGHT
rep.near("volume == 2 layers x one-direction area x thickness",
         K.signed_volume(verts, tris), expected, 1e-6, "{:.4f}")
rep.near("top-view union area == footprint minus holes",
         union_area, SIZE * SIZE - N_HOLES * WIDTH * WIDTH, 1e-6, "{:.2f}")
# Guards against regressing to a single layer or to both layers being grids.
rep.check("volume is NOT union area x height "
          "(each layer holds one bar direction only)",
          abs(expected - union_area * 2 * HEIGHT) > 1.0,
          f"{expected:.1f} vs {union_area * 2 * HEIGHT:.1f} mm^3")

# --- Renders -------------------------------------------------------------
SC = 10
win = (-SIZE / 2, SIZE / 2, -SIZE / 2, SIZE / 2)
img, dims = K.render(verts, tris, PREVIEW, view="xy", su=SC, window=win)
area = K.raster_area_mm2(img, SC, SC)
solid = SIZE * SIZE - N_HOLES * WIDTH * WIDTH
budget = N_HOLES * 4 * (WIDTH * SC) / 2 / (SC * SC)   # half-pixel fill bias
rep.check("top view area within rasterisation budget",
          abs(area - union_area) <= budget,
          f"{area:.1f} vs {union_area:.1f} mm^2 (budget +/-{budget:.1f})")
rep.check("top view size", dims[0] == int(SIZE * SC) + 1, f"{dims}")

# Area is badly biased by pixel quantisation here (bars are 5 px wide), so
# count features along a scanline instead. On a gap row only the vertical
# layer's bars are present, giving exactly N_BARS runs.
sp = spans_from(low, 0)
gap_mid = (sp[len(sp) // 2][1] + sp[len(sp) // 2 + 1][0]) / 2.0
runs = K.scanline_runs(img, gap_mid, -SIZE / 2, SC)
rep.check("scanline across a gap row counts 150 bars", runs == N_BARS,
          f"{runs} runs at y={gap_mid:g} mm")
bar_mid = (sp[len(sp) // 2][0] + sp[len(sp) // 2][1]) / 2.0
runs2 = K.scanline_runs(img, bar_mid, -SIZE / 2, SC)
rep.check("scanline across a bar row is continuous "
          "(upper-layer bars span the full width)", runs2 == 1,
          f"{runs2} run at y={bar_mid:g} mm")

# Side elevation, Z exaggerated 25x, so the two layers are visible.
_, sdim = K.render(verts, tris, SIDE, view="yz", su=4, sv=100,
                   window=(-SIZE / 2, SIZE / 2, 0, 2 * HEIGHT), flip_v=True)
rep.check("side view rendered", sdim[1] > 10, f"{sdim} (Z exaggerated 25x)")

print(f"  preview: {PREVIEW}  side: {SIDE}  ({N_HOLES} holes)")
sys.exit(0 if rep.finish() else 1)

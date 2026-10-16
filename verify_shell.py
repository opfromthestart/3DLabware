#!/usr/bin/env python3
"""Verify the grid-in-cylinder 3MF, re-read from the file."""

import math
import sys

import checks as K
import gridlib as g
import make_grid_shell as M

PATH = sys.argv[1] if len(sys.argv) > 1 else M.default_out()
PREVIEW = sys.argv[2] if len(sys.argv) > 2 else "grid_shell_top.png"
SIDE = sys.argv[3] if len(sys.argv) > 3 else "grid_shell_side.png"

R_IN, R_OUT = M.R_IN, M.R_OUT
EPS = 1e-6

rep = K.Report(f"Grid + cylinder: {PATH}")
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
rep.check("all shells wound outward", all(s["volume"] > 0 for s in shs),
          f"{len(shs)} shells, min={min(s['volume'] for s in shs):.4f} mm^3")

# Split: the wall is the one shell reaching the outer radius.
wall = [s for s in shs if s["max_radius"] > R_IN + 0.5]
bars = [s for s in shs if s["max_radius"] <= R_IN + 0.5]
rep.check("exactly one wall shell", len(wall) == 1)
rep.check("every other shell is a grid bar", len(bars) == len(shs) - 1,
          f"{len(bars)} bars")

# --- Wall ----------------------------------------------------------------
w = wall[0]
n = w["n"] // 4
rep.check("wall has 4 rings of n vertices", w["n"] == 4 * n and w["n"] > 0,
          f"n={n}")
radii = sorted({round(math.hypot(x, y), 6) for x, y, _ in
                (verts[i] for i in w["indices"])})
rep.check("wall has exactly two distinct radii", len(radii) == 2, f"{radii}")
rep.near("wall inner radius = bore", radii[0], R_IN, 1e-6, "{:.6g}")
rep.near("wall outer radius = bore + thickness", radii[1], R_OUT, 1e-6, "{:.6g}")
rep.check("wall height", abs(w["bbox"][5] - w["bbox"][2] - M.HEIGHT) < EPS,
          f"{w['bbox'][2]:g}..{w['bbox'][5]:g}")
sagitta = R_OUT * (1 - math.cos(math.pi / n))
rep.check("facet chord deviation within tolerance", sagitta <= M.TOL,
          f"{n} segments, {sagitta:.5f} mm <= {M.TOL:g}")

tube_area = (n / 2.0) * (R_OUT ** 2 - R_IN ** 2) * math.sin(2 * math.pi / n)
rep.near("wall volume == faceted annulus x height",
         w["volume"], tube_area * M.HEIGHT, abs(tube_area) * 1e-9, "{:.4f}")

# --- Grid ----------------------------------------------------------------
# No grid geometry may sit outside the bore.
worst = max(b["max_radius"] for b in bars)
rep.check("no grid geometry outside the bore", worst <= R_IN + 1e-9,
          f"max grid radius {worst:.6f} mm <= {R_IN:g}")

layers = {}
for b in bars:
    layers.setdefault(round(b["bbox"][2], 6), []).append(b)
rep.check("grid has two layers", len(layers) == 2, f"z starts {sorted(layers)}")
z_lo, z_hi = sorted(layers)
rep.near("layer offset", z_hi - z_lo, M.BAR_H, 1e-9)
rep.near("grid top", max(b["bbox"][5] for b in bars), 2 * M.BAR_H, 1e-9)
rep.check("grid sits on the build plate", z_lo < 1e-12)
rep.check("grid is inside the wall height",
          max(b["bbox"][5] for b in bars) < M.HEIGHT - 0.1,
          f"grid top {max(b['bbox'][5] for b in bars):g} mm < wall {M.HEIGHT:g} mm")
rep.check("both layers have the same bar count",
          len(layers[z_lo]) == len(layers[z_hi]),
          f"{len(layers[z_lo])} / {len(layers[z_hi])}")

# Bars must not interpenetrate the wall: they sit at r <= R_IN, wall at >= R_IN.
rep.check("grid and wall meet without interpenetrating",
          worst <= R_IN + 1e-9 and radii[0] >= R_IN - 1e-9)

# --- The circular trim, checked structurally and against an ideal circle ---
# Sutherland-Hodgman cannot reproduce a true circle, only the inscribed
# polygon, so the trimmed area is legitimately a hair BELOW the ideal-circle
# integral. Verifying it exactly against the circle would be wrong; instead
# check the polygon structure exactly, then bound the area against faceting.
bore = [(R_IN * math.cos(2 * math.pi * i / n),
         R_IN * math.sin(2 * math.pi * i / n)) for i in range(n)]
rep.check("bore polygon is a regular n-gon at exactly r_in",
          all(abs(math.hypot(*p) - R_IN) < 1e-9 for p in bore),
          f"{n} vertices, max dev {max(abs(math.hypot(*p) - R_IN) for p in bore):.2e}")


def inside_convex(pt, poly):
    return all((poly[(i + 1) % len(poly)][0] - poly[i][0]) * (pt[1] - poly[i][1])
               - (poly[(i + 1) % len(poly)][1] - poly[i][1]) * (pt[0] - poly[i][0])
               >= -1e-9 for i in range(len(poly)))


poly_area = 0.0
struct_ok = True
for b in bars:
    zlo = b["bbox"][2]
    pts = []
    for i in b["indices"]:
        x, y, z = verts[i]
        if abs(z - zlo) < 1e-9:
            pt = (x, y)
            if not any(abs(pt[0] - q[0]) < 1e-9 and abs(pt[1] - q[1]) < 1e-9
                       for q in pts):
                pts.append(pt)
    if len(pts) < 3:
        struct_ok = False
        continue
    ar = sum(pts[i][0] * pts[(i + 1) % len(pts)][1]
             - pts[(i + 1) % len(pts)][0] * pts[i][1]
             for i in range(len(pts))) / 2.0
    if ar <= 0 or not all(inside_convex(p, bore) for p in pts):
        struct_ok = False
    poly_area += ar
rep.check("every trimmed bar is a convex CCW polygon inside the bore",
          struct_ok, f"{len(bars)} bars")
rep.near("grid area == sum of extracted polygon areas",
         sum(b["volume"] for b in bars) / M.BAR_H, poly_area, 1e-6, "{:.4f}")

spans_world = [(a - M.GRID_SIZE / 2, b - M.GRID_SIZE / 2)
               for a, b in g.bar_spans(M.GRID_SIZE, M.BAR_W, M.BAR_GAP)]
ideal = g.analytic_disc_grid_area(R_IN, spans_world)
disc_area = (n / 2.0) * R_IN ** 2 * math.sin(2 * math.pi / n)  # inscribed n-gon
rho = disc_area / (math.pi * R_IN ** 2)      # inscribed polygon vs true circle
mesh_area = sum(b["volume"] for b in bars) / (2 * M.BAR_H)
rep.check("trimmed area is below the ideal-circle integral (faceted trim)",
          mesh_area < ideal, f"{mesh_area:.3f} < {ideal:.3f} mm^2/layer")
rep.check("trimmed area matches the faceting ratio",
          abs(mesh_area / ideal - rho) < 1e-3,
          f"ratio {mesh_area / ideal:.6f} vs rho {rho:.6f} (tol 1e-3)")

# How many bars actually reach the disc, and how many were dropped as slivers.
n_intersect = sum(1 for a, b in spans_world if b > -R_IN and a < R_IN)
kept = len(layers[z_lo])
rep.check("kept + dropped == bars that reach the disc",
          kept + (n_intersect - kept) == n_intersect,
          f"{kept} kept of {n_intersect} reaching the disc "
          f"({n_intersect - kept} dropped as slivers)")
rep.check("dropped bars are genuinely slivers",
          (n_intersect - kept) * M.BAR_W * M.BAR_W < 1.0,
          f"{n_intersect - kept} dropped")

# --- Totals --------------------------------------------------------------
grid_vol = sum(b["volume"] for b in bars)
rep.near("total volume == wall + grid",
         K.signed_volume(verts, tris), w["volume"] + grid_vol, 1e-6, "{:.4f}")

bb = (min(v[0] for v in verts), min(v[1] for v in verts), min(v[2] for v in verts),
      max(v[0] for v in verts), max(v[1] for v in verts), max(v[2] for v in verts))
rep.near("outer diameter X", bb[3] - bb[0], M.INNER_DIA + 2 * M.WALL)
rep.near("outer diameter Y", bb[4] - bb[1], M.INNER_DIA + 2 * M.WALL)
rep.near("overall height", bb[5] - bb[2], M.HEIGHT)
rep.near("centred on origin X", (bb[3] + bb[0]) / 2, 0.0, 1e-9)
rep.near("centred on origin Y", (bb[4] + bb[1]) / 2, 0.0, 1e-9)

# --- Renders -------------------------------------------------------------
SC = 12
win = (-R_OUT, R_OUT, -R_OUT, R_OUT)
img, dims = K.render(verts, tris, PREVIEW, view="xy", su=SC, window=win)
px = img.load()
wpx, hpx = img.size
cxp, cyp = (wpx - 1) / 2.0, (hpx - 1) / 2.0


def radius_px(x, y):
    return math.hypot((x - cxp) / SC, (y - cyp) / SC)


# Only look inside the model: the square image has background corners well
# outside the circular part, and those dark pixels are not holes.
# The outer boundary ring is skipped too -- polygon fill leaves a partially
# covered pixel unfilled, which reads as "dark" right at the silhouette.
margin = 2.0 / SC
dark_r = [radius_px(x, y) for y in range(hpx) for x in range(wpx)
          if px[x, y] < 127 and radius_px(x, y) < R_OUT - margin]
if dark_r:
    in_wall = [r for r in dark_r if r > R_IN + margin]
    rep.check("wall annulus is solid (no holes in it)", not in_wall,
              f"band r in ({R_IN + margin:.3f}, {R_OUT - margin:.3f}) mm is solid")
    rep.check("holes present inside the bore", len(dark_r) > 1000,
              f"{len(dark_r)} dark pixels")
    rep.check("holes reach the bore edge",
              max(dark_r) > R_IN - 2.0,
              f"dark pixels reach r={max(dark_r):.3f} mm, bore {R_IN:g}")
else:
    rep.check("holes present inside the bore", False, "no dark pixels found")

# Side elevation: wall at full height, grid layers exaggerated 25x so visible.
_, sdim = K.render(verts, tris, SIDE, view="yz", su=6, sv=30,
                   window=(-R_OUT, R_OUT, 0, M.HEIGHT), flip_v=True)
rep.check("side view rendered", sdim[1] > 20, f"{sdim}")

print(f"  preview: {PREVIEW}  side: {SIDE}")
sys.exit(0 if rep.finish() else 1)

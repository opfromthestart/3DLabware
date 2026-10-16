#!/usr/bin/env python3
"""Generate the grid floor + cylindrical shell 3MF.

A two-layer grid floor trimmed to the bore of a cylindrical wall:

    bore (inside diameter)  110 mm   -> r_in  = 55 mm
    wall thickness           3 mm    -> r_out = 58 mm
    overall height          10 mm    -> z = 0 .. 10
    grid floor               0.2 mm per layer, two layers (0 .. 0.4)

The grid is trimmed to r_in with the *same* facet polygon that builds the
tube's bore, so every trimmed bar edge lands exactly on a bore facet
instead of poking through it.

Usage: python3 make_grid_shell.py [output.3mf]
"""

import math
import sys

import gridlib as g
import threelib as t

# --- Parameters (mm) -------------------------------------------------------
INNER_DIA = 110.0   # bore diameter
WALL = 1.5          # radial wall thickness
HEIGHT = 10.0       # overall part height
GRID_SIZE = INNER_DIA + 2*WALL  # the grid square before trimming (must exceed 2*r_in)
BAR_W = 0.5         # grid bar width
BAR_GAP = 1.0       # grid bar gap
BAR_H = 0.2         # grid layer thickness
TOL = 0.01          # max chord deviation on the round surfaces

R_IN = INNER_DIA / 2.0
R_OUT = R_IN + WALL

# Single source of truth for params and output name.
PARAMS = dict(inner_dia=INNER_DIA, wall=WALL, height=HEIGHT,
              grid_size=GRID_SIZE, bar_w=BAR_W, bar_gap=BAR_GAP,
              bar_h=BAR_H, tol=TOL)


def default_out():
    return t.stamped_filename("grid_shell", PARAMS)


def build(inner_dia=INNER_DIA, wall=WALL, height=HEIGHT,
          grid_size=GRID_SIZE, bar_w=BAR_W, bar_gap=BAR_GAP, bar_h=BAR_H,
          tol=TOL):
    r_in = inner_dia / 2.0
    r_out = r_in + wall
    n = t.auto_segments(r_out, tol)

    # One polygon, two jobs: the tube's bore and the grid's clip region.
    bore = t.disc_polygon(r_in, n)

    mesh = t.Mesh("grid_shell")
    t.tube(mesh, r_in, r_out, 0.0, height, n)
    stats = g.add_grid(mesh, grid_size, bar_w, bar_gap, bar_h,
                       z0=0.0, center=(0.0, 0.0), clip_poly=bore)
    sagitta = r_out * (1 - math.cos(math.pi / n))
    return mesh, stats, n, sagitta


def main():
    params = PARAMS
    out = sys.argv[1] if len(sys.argv) > 1 else default_out()
    mesh, stats, n, sagitta = build(**params)

    desc = (f"Grid floor in a cylindrical shell. Bore {INNER_DIA:g} mm "
            f"(r={R_IN:g}), wall {WALL:g} mm (r_out={R_OUT:g}), height "
            f"{HEIGHT:g} mm. Two-layer grid {BAR_W:g} mm bars / {BAR_GAP:g} mm "
            f"gaps trimmed to the bore; horizontal bars sit "
            f"{BAR_H:g} mm above vertical bars.")
    t.write_3mf(out, (mesh.name, mesh), title="Grid floor in cylinder",
                description=desc, params=params)

    bb = mesh.bbox()
    print(f"wrote {out}")
    print(f"  params            : {t.describe(params)}")
    print(f"  param hash        : {t.params_hash(params)}")
    print(f"  bore / wall         : r {R_IN:g} -> {R_OUT:g} mm "
          f"({INNER_DIA:g} mm ID, {WALL:g} mm wall)")
    print(f"  height              : {HEIGHT:g} mm")
    print(f"  round facets        : {n} segments, chord deviation {sagitta:.5f} mm")
    for lay in stats["layers"]:
        print(f"  grid {lay['label']:<10}: {lay['bars']:>3} bars kept, "
              f"{lay['area']:>9.2f} mm^2, z {lay['z'][0]:g}..{lay['z'][1]:g}")
    print(f"  bars dropped (sliver): {stats['dropped']}")
    print(f"  bounding box        : {bb[3] - bb[0]:g} x {bb[4] - bb[1]:g}"
          f" x {bb[5] - bb[2]:g} mm")
    print(f"  vertices / tris    : {len(mesh.verts)} / {len(mesh.tris)}")
    print(f"  total volume        : {mesh.signed_volume():.3f} mm^3")
    print(f"  shells              : {len(mesh.shells())} (1 wall + {sum(l['bars'] for l in stats['layers'])} bars)")


if __name__ == "__main__":
    main()

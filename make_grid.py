#!/usr/bin/env python3
"""Generate the 150 mm grid floor as a 3MF file.

Two interleaved layers, so no two bars ever interpenetrate:

    z = 0.0 .. 0.2   vertical bars   (long axis along Y)
    z = 0.2 .. 0.4   horizontal bars (long axis along X)

Usage: python3 make_grid.py [output.3mf]
"""

import sys

import gridlib as g
import threelib as t

SIZE = 150.0    # square footprint
WIDTH = 0.5     # bar width
GAP = 0.5       # clear gap between bars (=> 1.0 mm pitch)
HEIGHT = 0.2    # thickness of one bar layer

# Single source of truth: the generator writes here and the verifier reads
# from here, so a filename can never drift away from the parameters.
PARAMS = dict(size=SIZE, width=WIDTH, gap=GAP, height=HEIGHT)


def default_out():
    return t.stamped_filename("grid", PARAMS)


def build(size=SIZE, width=WIDTH, gap=GAP, height=HEIGHT):
    mesh = t.Mesh("grid")
    stats = g.add_grid(mesh, size, width, gap, height, z0=0.0)
    # Both layers must land on the same square; a mis-centred layer shows up
    # here as a footprint larger than `size`.
    bb = mesh.bbox()
    for got, axis in ((bb[3] - bb[0], "X"), (bb[4] - bb[1], "Y")):
        if abs(got - size) > 1e-9:
            raise AssertionError(f"grid {axis} footprint {got:g} != {size:g} mm")
    return mesh, stats


def main():
    params = PARAMS
    out = sys.argv[1] if len(sys.argv) > 1 else default_out()
    mesh, stats = build(**params)

    description = (f"Two-layer grid {SIZE:g}x{SIZE:g} mm. Vertical bars at "
                   f"z=0..{HEIGHT:g} mm, horizontal bars at "
                   f"z={HEIGHT:g}..{2 * HEIGHT:g} mm. Bar width {WIDTH:g} mm, "
                   f"gap {GAP:g} mm ({WIDTH + GAP:g} mm pitch).")
    t.write_3mf(out, (mesh.name, mesh),
                title=f"Grid {SIZE:g}mm (2 layers)", description=description,
                params=params)

    bb = mesh.bbox()
    bar_sum = sum(b - a for a, b in stats["spans"])
    expected_vol = 2 * bar_sum * SIZE * HEIGHT

    print(f"wrote {out}")
    print(f"  params             : {t.describe(params)}")
    print(f"  param hash         : {t.params_hash(params)}")
    print(f"  bars per direction : {stats['n_per_direction']}"
          f"  (x2 layers = {2 * stats['n_per_direction']} boxes)")
    print(f"  pitch              : {WIDTH + GAP:g} mm")
    print("  layer z-ranges     : "
          + ", ".join(f"{l['label']} {l['z'][0]:g}..{l['z'][1]:g}"
                      for l in stats["layers"]))
    print(f"  bounding box       : {bb[3] - bb[0]:g} x {bb[4] - bb[1]:g}"
          f" x {bb[5] - bb[2]:g} mm")
    print(f"  vertices / tris    : {len(mesh.verts)} / {len(mesh.tris)}")
    print(f"  volume             : {mesh.signed_volume():.3f} mm^3"
          f" (expected {expected_vol:.3f})")
    print(f"  shells             : {len(mesh.shells())} (disjoint, no overlap)")


if __name__ == "__main__":
    main()

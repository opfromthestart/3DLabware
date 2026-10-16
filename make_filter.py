#!/usr/bin/env python3
"""Generate a Buchner-style filter cartridge as a 3MF.

The part is one continuous vessel: a cup whose floor is the filter media, a
hollow funnel below it, and a spout. Read in USE orientation, top to bottom:

    collar, no filter    open cup wall, `headroom` deep, holds the slurry
    collar, filter       the 70 mm media disc, sealed all round
    funnel, no filter    open conical cavity guiding filtrate to the stem
    spout, no filter     the outlet bore

Two things this shape has to get right, both of which a solid cone fails:

* The funnel wall must *continue* the collar -- its bore is the media's
  radius and its outer wall is the collar's outer wall. Any step between the
  two leaves an annular gap under the media lip, and liquid running off the
  outer edge of the disc escapes down the outside of the funnel instead of
  through it.
* The funnel must be *hollow*. A cone with only a stem-sized bore through it
  presents a solid ceiling to the whole disc: filtrate lands on top of it and
  pools, and only the small fraction directly over the bore drains.

The cup wall is thick (`collar_wall` 4 mm) and the funnel wall thin
(`funnel_wall` 1.5 mm) on purpose. The step between those two -- 2.5 mm of
flat annulus facing away from the funnel -- is the ledge the holder's rim
seats against. With a 1.5 mm collar the funnel's outer surface would meet the
collar's outer surface exactly, leaving nothing to rest on.

PRINT ORIENTATION. As written, z=0 is the collar's rim and the media sits
`headroom` mm above the bed, so the media bridges the open cup and needs
support; the alternative (media on the bed) puts the collar's skirt below
z=0. Either way something floats -- this shape cannot be printed flat on the
bed without support. Reposition and support it as convenient.

Usage: python3 make_filter.py [output.3mf]
"""

import math
import sys

import csglib as csg
import gridlib as g
import threelib as t

# --- Parameters (mm) -------------------------------------------------------
GRID_DIA = 90.0      # filter media diameter
LAYERS = 3          # alternating grid layers
BAR_W = 0.5          # bar width
BAR_GAP = 1.0        # clear gap between bars
BAR_H = 0.12          # thickness of one layer
GRID_SQ = 90.0       # square the grid is laid out in before trimming
HEADROOM = 60.0      # open cup depth above the media (holds the slurry)
COLLAR_WALL = 1.5    # cup wall thickness; see the seating note above
FUNNEL_WALL = 1.5    # funnel wall thickness
STEM_DIA = 20.0      # spout outer diameter
FUNNEL_H = (GRID_DIA-STEM_DIA)/3      # funnel height
STEM_H = 10.0        # spout length
BORE_R = 8.8         # spout bore radius
STAGGER = 0.0        # half-pitch offset between same-direction layers
TOL = 0.01           # max chord deviation on round surfaces

GRID_R = GRID_DIA / 2.0
STEM_R = STEM_DIA / 2.0
GRID_H = LAYERS * BAR_H
COLLAR_R = GRID_R + COLLAR_WALL
FUNNEL_TOP_OUT = GRID_R + FUNNEL_WALL     # funnel outer radius where it leaves the cup
FUNNEL_SLOPE = (FUNNEL_TOP_OUT - STEM_R) / FUNNEL_H
HEADROOM_SNAPPED = round(math.ceil(HEADROOM / BAR_H) * BAR_H, 9)

# Stack heights, z=0 at the collar's rim. Each face is a whole number of
# layers: a rim ending mid-layer is sliced across and its real height becomes
# whatever the layer height happens to be.
Z_COLLAR0 = 0.0
Z_MEDIA0 = HEADROOM_SNAPPED
Z_MEDIA1 = Z_MEDIA0 + GRID_H
Z_FUNNEL0 = Z_MEDIA1
Z_FUNNEL1 = Z_FUNNEL0 + FUNNEL_H
Z_TOP = Z_FUNNEL1 + STEM_H

PARAMS = dict(grid_dia=GRID_DIA, layers=LAYERS, bar_w=BAR_W, bar_gap=BAR_GAP,
              bar_h=BAR_H, grid_sq=GRID_SQ, headroom=HEADROOM_SNAPPED,
              collar_wall=COLLAR_WALL, funnel_wall=FUNNEL_WALL,
              funnel_h=FUNNEL_H, stem_dia=STEM_DIA, stem_h=STEM_H,
              bore_r=BORE_R, stagger=STAGGER, tol=TOL)


def funnel_outer_r(z):
    """Filter's outer radius at height z in USE orientation (spout tip at z=0).

    Single source of truth for the surface the holder has to clear: the holder
    derives its bore from this rather than restating the cone.
    """
    if z <= STEM_H:
        return STEM_R
    return min(FUNNEL_TOP_OUT, STEM_R + FUNNEL_SLOPE * (z - STEM_H))


# Height at which the funnel meets the cup: the top of the usable taper, and
# the plane the holder's rim seats the filter on.
Z_SEAT_USE = STEM_H + FUNNEL_H


def default_out():
    return t.stamped_filename("filter", PARAMS)


def build(**p):
    grid_r = p["grid_dia"] / 2.0
    stem_r = p["stem_dia"] / 2.0
    grid_h = p["layers"] * p["bar_h"]
    collar_r = grid_r + p["collar_wall"]
    funnel_top = grid_r + p["funnel_wall"]

    z_media0 = p["headroom"]
    z_media1 = z_media0 + grid_h
    z_funnel0 = z_media1
    z_funnel1 = z_funnel0 + p["funnel_h"]
    z_top = z_funnel1 + p["stem_h"]
    n = t.auto_segments(max(collar_r, funnel_top, stem_r), p["tol"])

    # The media, the cup, the funnel and the spout are built separately so the
    # three solid pieces can be unioned. Their overruns exist *for* that union
    # -- the funnel runs half a millimetre up into the cup, the spout a
    # millimetre into the funnel, both so the union has no coincident faces --
    # and the union was never done. The pieces were merged instead, which
    # leaves the funnel and the cup interpenetrating by 1.4 mm: measured, and
    # the one real overlap left in this part.
    media = t.Mesh("media")
    disc = t.disc_polygon(grid_r, n)

    # Filter media. Within one bar direction consecutive layers are staggered
    # by half a pitch, so each bar lands over the gaps of the layer below and
    # halves the span it would otherwise bridge.
    spans = g.bar_spans(p["grid_sq"], p["bar_w"], p["bar_gap"])
    pitch = p["bar_w"] + p["bar_gap"]
    layers = []
    for i in range(p["layers"]):
        along = "y" if i % 2 == 0 else "x"     # even = vertical, odd = horizontal
        stagger = p["stagger"] * ((i // 2) % 2) * pitch
        lay = g.add_bar_layer(media, spans, p["grid_sq"],
                              z_media0 + i * p["bar_h"], p["bar_h"], along,
                              clip_poly=disc, shift=stagger,
                              min_width=p["bar_w"])
        lay["index"] = i
        lay["direction"] = "vertical" if along == "y" else "horizontal"
        lay["stagger"] = stagger
        layers.append(lay)

    # Cup wall. Its bore is the same polygon the media is trimmed to, so they
    # are flush and there is no ledge for solids to settle on. It spans the
    # headroom *and* the media, so the two share a cylindrical face.
    cup = t.Mesh("cup")
    t.tube(cup, grid_r, collar_r, Z_COLLAR0, z_media1, n)

    # Funnel. Bore continues the media's radius, outer wall continues the
    # cup's, so there is no step for liquid to escape through.
    #
    # It is defined by its profile at the media's underside (bore = media
    # radius, outer = cup wall) and at the stem. It runs half a millimetre up
    # into the cup so the union has no coincident faces, so the cone is
    # *extrapolated* rather than started at those radii: a cone pinched by the
    # overrun would reach the media's underside already narrowed, leaving a
    # ring of the disc draining onto the funnel's rim instead of into it.
    in_slope = (grid_r - p["bore_r"]) / p["funnel_h"]
    out_slope = (funnel_top - stem_r) / p["funnel_h"]
    # ...but only as far up as the cup wall can hide it. The funnel flares
    # going down, so the overrun sticks out sideways: with a cup wall no
    # thicker than the funnel's own, half a millimetre of overrun puts a
    # 1.3 mm lip on the outside of the part. When the two walls are equal the
    # cross-sections at the junction already match and they can simply butt.
    ov = min(0.5, max(0.0, (collar_r - funnel_top) / out_slope))
    funnel = t.Mesh("funnel")
    t.tapered_tube(funnel, grid_r + in_slope * ov,
                   funnel_top + out_slope * ov,
                   p["bore_r"], stem_r, z_funnel0 - ov, z_funnel1, n)

    # Spout, overlapping the funnel so the union has no coincident faces.
    spout = t.Mesh("spout")
    t.tube(spout, p["bore_r"], stem_r, z_funnel1 - 1.0, z_top, n)

    # The three solid pieces cross one another, so they union; the media only
    # touches the cup's bore -- it is trimmed flush to it on purpose, so there
    # is no ledge for solids to settle on -- and is merged in as it is.
    mesh = csg.union(cup, funnel, spout, name="filter")
    mesh.merge(media)

    sagitta = grid_r * (1 - math.cos(math.pi / n))
    return mesh, {"layers": layers, "n": n, "sagitta": sagitta,
                  "collar_r": collar_r, "funnel_top": funnel_top,
                  "z": {"collar": (Z_COLLAR0, z_media1),
                        "media": (z_media0, z_media1),
                        "funnel": (z_funnel0 - 0.5, z_funnel1),
                        "spout": (z_funnel1 - 1.0, z_top),
                        "top": z_top}}


def main():
    out = sys.argv[1] if len(sys.argv) > 1 else default_out()
    mesh, s = build(**PARAMS)
    z = s["z"]

    desc = (f"Buchner-style filter. {PARAMS['grid_dia']:g} mm disc of "
            f"{PARAMS['layers']} alternating {PARAMS['bar_w']:g} mm bars / "
            f"{PARAMS['bar_gap']:g} mm gaps, {GRID_H:g} mm thick. "
            f"{PARAMS['headroom']:g} mm open cup above the media, "
            f"{PARAMS['collar_wall']:g} mm cup wall. Hollow funnel "
            f"{PARAMS['funnel_wall']:g} mm wall from a "
            f"{PARAMS['grid_dia']:g} mm bore to a {PARAMS['stem_dia']:g} mm "
            f"stem over {PARAMS['funnel_h']:g} mm, then a "
            f"{PARAMS['stem_h']:g} mm spout with a "
            f"{PARAMS['bore_r'] * 2:g} mm bore.")
    t.write_3mf(out, (mesh.name, mesh), title="Filter cartridge",
                description=desc, params=PARAMS)

    grid_area = sum(l["area"] for l in s["layers"])
    solid_vol = mesh.signed_volume() - grid_area * PARAMS["bar_h"]
    bb = mesh.bbox()
    disc_area = math.pi * GRID_R ** 2

    print(f"wrote {out}")
    print(f"  params            : {t.describe(PARAMS)}")
    print(f"  param hash        : {t.params_hash(PARAMS)}")
    print(f"  filter disc       : {PARAMS['grid_dia']:g} mm dia, "
          f"{PARAMS['layers']} layers x {PARAMS['bar_h']:g} mm = {GRID_H:g} mm")
    print(f"  layer alternation : "
          f"{sum(1 for l in s['layers'] if l['direction'] == 'vertical')} vertical / "
          f"{sum(1 for l in s['layers'] if l['direction'] == 'horizontal')} horizontal")
    if len(s['layers'])>0:
        print(f"  bars per layer    : {s['layers'][0]['bars']} "
            f"(dropped as slivers: {sum(l['dropped'] for l in s['layers'])})")
    print(f"  cup               : {s['collar_r'] * 2:g} mm OD x "
          f"{PARAMS['collar_wall']:g} mm wall, bore flush with the media")
    print(f"  headroom          : {PARAMS['headroom']:g} mm of open cup above "
          f"the media ({disc_area * PARAMS['headroom'] / 1000.0:.0f} ml)")
    print(f"  funnel            : bore r {PARAMS['grid_dia'] / 2:g} -> "
          f"{PARAMS['bore_r']:g}, outer r {s['funnel_top']:g} -> "
          f"{PARAMS['stem_dia'] / 2:g}, {PARAMS['funnel_h']:g} mm tall")
    print(f"  seating ledge     : r {s['funnel_top']:g}..{s['collar_r']:g} "
          f"({s['collar_r'] - s['funnel_top']:g} mm) facing the holder, "
          f"{Z_SEAT_USE:g} mm below the spout tip")
    print(f"  round facets      : {s['n']} segments, chord {s['sagitta']:.5f} mm")
    print(f"  z ranges          : collar {z['collar'][0]:g}..{z['collar'][1]:g}, "
          f"media {z['media'][0]:g}..{z['media'][1]:g}, "
          f"funnel {z['funnel'][0]:g}..{z['funnel'][1]:g}, "
          f"spout {z['spout'][0]:g}..{z['spout'][1]:g}")
    print(f"  bounding box      : {bb[3] - bb[0]:g} x {bb[4] - bb[1]:g}"
          f" x {bb[5] - bb[2]:g} mm")
    print(f"  vertices / tris   : {len(mesh.verts)} / {len(mesh.tris)}")
    print(f"  volume            : {mesh.signed_volume():.2f} mm^3 "
          f"(media {grid_area * PARAMS['bar_h']:.1f}, body {solid_vol:.1f})")
    print(f"  shells            : {len(mesh.shells())}")


if __name__ == "__main__":
    main()

"""gridlib -- a reusable "grid floor" pattern built on threelib.

The grid is two interleaved layers of parallel bars:

    vertical bars   (long axis along Y)   occupy z0            .. z0 + height
    horizontal bars (long axis along X)   occupy z0 + offset   .. + height

With the default `offset == height` the layers tile without overlapping, so
every bar is a separate closed box and the mesh has no interpenetrating
geometry at all. Raise `offset` above `height` if you want the two layers
to physically interleave (e.g. to weave them through each other).

Bars are laid out symmetrically and flush with the requested square, and a
circular/other convex clip region can be supplied so the grid becomes the
floor of a round vessel.
"""

from __future__ import annotations

import math

import threelib as t

__all__ = ["bar_spans", "add_bar_layer", "add_grid", "grid_volume", "analytic_disc_grid_area"]


def bar_spans(total, width, gap, edge=None):
    """Bar intervals along one axis, summing exactly to `total`.

    Returns a list of (start, end) intervals, symmetric and flush at both
    ends. An exact fit needs `total` to be a whole multiple of the pitch
    (width + gap); when it is not, the remainder is absorbed by widening the
    two end bars to `width + gap/2`. Without this the grid silently ends up
    the wrong size or leaves a ragged last cell.

    >>> spans = bar_spans(150.0, 0.5, 0.5)
    >>> len(spans), spans[0], spans[-1]
    (150, (0.0, 0.75), (149.25, 150.0))
    """
    pitch = width + gap
    edge = width + gap / 2.0 if edge is None else edge

    # Too small for the symmetric edge scheme: centre one bar instead.
    if total < 2 * edge + width + gap:
        start = (total - width) / 2.0
        return [(start, start + width)]

    spans = [(0.0, edge)]
    x = edge + gap
    while x + width + gap <= total - edge + 1e-9:
        spans.append((x, x + width))
        x += pitch
    spans.append((total - edge, total))
    return spans


def add_bar_layer(mesh, spans, size, z0, height, along="y",
                  center=(0.0, 0.0), clip_poly=None, min_area=1e-2,
                  shift=0.0, min_width=None):
    """Add ONE layer of parallel bars; returns per-layer stats.

    along="y" -> bars run along Y, positioned by spans along X ("vertical")
    along="x" -> bars run along X, positioned by spans along Y ("horizontal")

    `spans` must come from `bar_spans(size, ...)`, i.e. intervals covering
    0..size. Use this directly when you need more than two layers -- a filter
    with 30 alternating layers is not expressible as pairs.

    `shift` slides the pattern along its own span axis. Shifting by half a
    pitch staggers the layer so its bars land over the *gaps* of the layer
    below, which halves the unsupported span each bar bridges. Lay `size`
    out with margin for it: a shifted pattern runs `shift` past one edge.

    `min_width` drops bars that survived clipping but are narrower than a
    full bar along the bar's *narrow* axis. Staggering leaves bars hanging
    half off a round boundary, and the offcut is a sliver: structurally
    useless, and narrower than a nozzle, so it prints as a blob or not at
    all. The test is on the narrow axis only, so a bar clipped short along
    its length -- the normal case at the edge of a disc -- is kept.
    """
    cx, cy = center
    half = size / 2.0
    ox, oy = cx - half, cy - half
    built = area = dropped = clipped = 0
    for a, b in spans:
        a += shift
        b += shift
        if along == "y":
            rect = (ox + a, cy - half, ox + b, cy + half)
        else:
            rect = (cx - half, oy + a, cx + half, oy + b)
        poly = t.rect_polygon(*rect)
        if clip_poly is not None:
            poly = t.clean_polygon(t.clip_convex(poly, clip_poly))
            if len(poly) < 3 or abs(t.polygon_area(poly)) < min_area:
                dropped += 1
                continue
            if min_width is not None:
                narrow = ((max(q[0] for q in poly) - min(q[0] for q in poly))
                          if along == "y"
                          else (max(q[1] for q in poly) - min(q[1] for q in poly)))
                if narrow < min_width - 1e-9:
                    dropped += 1
                    continue
            clipped += 1
        t.prism(mesh, poly, z0, z0 + height)
        built += 1
        area += abs(t.polygon_area(poly))
    return {"bars": built, "area": area, "dropped": dropped,
            "clipped": clipped, "z": (z0, z0 + height)}


def add_grid(mesh, size, width=0.5, gap=0.5, height=0.2,
             z0=0.0, offset=None, center=(0.0, 0.0), clip_poly=None,
             min_area=1e-2, name="grid"):
    """Add a two-layer grid floor to `mesh`; returns a stats dict.

    size      side length of the square the grid is laid out in
    width     bar width
    gap       clear gap between neighbouring bars
    height    thickness of one bar layer
    z0        bottom of the lower (vertical-bar) layer
    offset    gap between the two layers; defaults to `height` (stacked)
    center    (cx, cy) the grid is centred on
    clip_poly convex CCW polygon the grid is trimmed to (e.g. a disc);
              bars are clipped with Sutherland-Hodgman, so pass the *same*
              polygon that builds the surrounding wall and the trimmed edge
              will coincide exactly with it
    min_area  trimmed pieces smaller than this are dropped, which removes
              slivers too thin to slice

    Bars run the full length of the square in their long direction; only
    the clipping region trims them short.
    """
    offset = height if offset is None else offset
    spans = bar_spans(size, width, gap)   # spans run 0..size

    stats = {
        "size": size, "width": width, "gap": gap, "height": height,
        "offset": offset, "spans": spans, "n_per_direction": len(spans),
        "clip": clip_poly is not None, "dropped": 0, "clipped": 0,
        "layers": [],
    }

    for label, z_lo, along in (("vertical", z0, "y"),
                               ("horizontal", z0 + offset, "x")):
        lay = add_bar_layer(mesh, spans, size, z_lo, height, along,
                            center=center, clip_poly=clip_poly,
                            min_area=min_area)
        lay["label"] = label
        stats["layers"].append(lay)
        stats["dropped"] += lay["dropped"]
        stats["clipped"] += lay["clipped"]
    return stats


def grid_volume(solid_area_per_layer, height, n_layers):
    """Analytic volume of a layered grid: sum of (bar area x height)."""
    return solid_area_per_layer * height * n_layers


def analytic_disc_grid_area(radius, spans, samples=20001):
    """Area of a set of bars trimmed to a disc, via numeric integration.

    A bar long in X and narrow in Y, trimmed to a disc of radius R, has
    area = integral over its y-range of 2*sqrt(R^2 - y^2) dy, because that is
    the chord length of the disc at height y. `spans` are the bar intervals
    along Y; bars are assumed to run the full disc width.

    This is deliberately a different code path from threelib's polygon
    clipping (Simpson's rule vs. Sutherland-Hodgman), so agreement is real
    evidence that the clip is correct.
    """
    total = 0.0
    for a, b in spans:
        lo, hi = max(-radius, a), min(radius, b)
        if hi <= lo:
            continue
        n = samples if (samples % 2 == 0) else samples + 1
        h = (hi - lo) / n
        acc = 0.0
        for i in range(n + 1):
            y = lo + i * h
            f = 2.0 * math.sqrt(max(0.0, radius * radius - y * y))
            w = 1 if i in (0, n) else (4 if i % 2 else 2)
            acc += w * f
        total += acc * h / 3.0
    return total

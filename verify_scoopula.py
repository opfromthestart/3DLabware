#!/usr/bin/env python3
"""Verify a scoopula 3MF by re-reading it, knowing only what was asked for.

The part is a long, narrow trough whose long axis runs *in line* with a flat bar
handle: closed at the heel, running out to a thin edge at the toe, open along the
top, printed lying on its back. Its channel is a **U** -- a round bottom tangent
to two flanks that lean outward -- and its wall is thin, because that is what a
drawn steel scoopula is. Nothing here is read from the builder: every number is
measured off the mesh by intersecting test lines with its triangles, not by a
point-in-solid parity test, which this part's boolean slivers can flip.

Four checks are here because earlier versions had neither, or had a weaker one.

* **The cavity is open.** The old verifier said "the cavity is open above the
  floor" because one probe point on the axis read as void -- which is exactly
  what a *sealed* void reads as. That build's "trough" was a hollow shell with a
  0.04 mm roof over it: one shell, zero gaps, right bounding box, right volume,
  and it would never have held powder. Here every column inside the measured
  mouth is fired upward through the whole height of the part and its crossings
  are counted: one band, from the bed to the floor, and nothing over it.
* **The scoop points the right way.** The old handle was joined *across* the
  trough's axis, so the part was a trowel with a sideways blade, and the old
  verifier agreed because it restated the builder's geometry instead of
  measuring the artifact. So the mouth's length is compared against its width,
  the handle is required on the axis beyond the heel, and the whole part's
  length, width and height are compared against what was asked for.
* **The heel is closed at the height the powder is carried at, not just at the
  floor.** An open-backed scoop is the same mistake as a roofed one, in the
  opposite direction, and a cut across the section at 2.4 mm cannot see it: the
  channel has to be solid in a row taken three quarters of the way up the part.
* **The bottom is round, not a V.** The width of the channel a little way above
  the floor decides that: a V is near nothing there and a U is already a third of
  the mouth across, and the rate at which the width grows drops from that to the
  flank's lean as you rise. This is what makes it a drawn channel instead of a
  folded one, and it is the difference the user sees.
* **The wall is measured where it is a wall.** Across the flanks it must be the
  thickness asked for; down in the arc it is legitimately thicker, and a test
  that averages the two tells you nothing about either.

Specimens, all still in the directory, all rejected by name above:
`scoopula_84214874.3mf` is a sealed void with the handle across the trough (21
failures), `scoopula_264ece36.3mf` is an open trough with a *flat* bottom and
2 mm walls (5), `scoopula_4b5580e3.3mf` is a round channel whose back is open
and whose floor landed 1 um off the handle's top face (6). The current
`scoopula_1fec0f3f.3mf` passes.

Usage: python3 verify_scoopula.py scoopula_xxxxxxxx.3mf [--angle 45]
"""

import math
import sys

import checks as K
import csglib as csg
import memcap
import printability as PR

# --- What was asked for, restated so the verifier does not share the builder's
# --- geometry -------------------------------------------------------------
IN_SPEC = dict(scoop_len=62.0, hnd_l=98.0, rim_w=14.0, depth=8.0, wall=1.4,
               floor=1.6, lean=0.36, back_wall=3.0, hnd_w=7.0, hnd_w_tip=6.0,
               hnd_t=1.2, h_lip=1.6, grid=0.4, nozzle=0.4)
TOL = 0.1           # what a measurement may be out by and still pass
Y_STEP = 0.30       # mm, for the independent volume by column integration
PROBE = (0.0137, 0.0173)   # stand every probe off the axis and off whole mm


def hits(verts, tris, p0, p1):
    """t values where the segment p0->p1 crosses the surface, sorted.

    Every face here is planar, so this is an exact segment/triangle test rather
    than a raycast through a soup: it does not care how a polygon was
    triangulated, and a 1e-7 mm^2 sliver cannot flip an answer the way it can a
    parity point test.
    """
    dx, dy, dz = (p1[k] - p0[k] for k in range(3))
    out = []
    for a, b, c in tris:
        v0, v1, v2 = verts[a], verts[b], verts[c]
        e1 = (v1[0] - v0[0], v1[1] - v0[1], v1[2] - v0[2])
        e2 = (v2[0] - v0[0], v2[1] - v0[1], v2[2] - v0[2])
        px = dy * e2[2] - dz * e2[1]
        py = dz * e2[0] - dx * e2[2]
        pz = dx * e2[1] - dy * e2[0]
        det = e1[0] * px + e1[1] * py + e1[2] * pz
        if abs(det) < 1e-12:                      # parallel: no clean crossing
            continue
        inv = 1.0 / det
        tx, ty, tz = p0[0] - v0[0], p0[1] - v0[1], p0[2] - v0[2]
        u = (tx * px + ty * py + tz * pz) * inv
        if u < 0.0 or u > 1.0:
            continue
        qx = ty * e1[2] - tz * e1[1]
        qy = tz * e1[0] - tx * e1[2]
        qz = tx * e1[1] - ty * e1[0]
        vv = (dx * qx + dy * qy + dz * qz) * inv
        if vv < 0.0 or u + vv > 1.0:
            continue
        tt = (e2[0] * qx + e2[1] * qy + e2[2] * qz) * inv
        if 1e-9 < tt < 1.0 - 1e-9:
            out.append(tt)
    out.sort()
    merged = []
    for tt in out:
        if not merged or tt - merged[-1] > 1e-7:
            merged.append(tt)
    return merged


def cuts(verts, tris, p0, p1, ax):
    """Material intervals along axis `ax` as ([(lo, hi), ...], [crossings]).

    A closed surface is crossed an even number of times and the crossings
    alternate empty/material from outside, so consecutive pairs are the material.
    An odd count means the line grazed a face or the surface is not closed, and
    comes back as None with the crossings it did find.
    """
    coords = sorted(p0[ax] + tt * (p1[ax] - p0[ax])
                    for tt in hits(verts, tris, p0, p1))
    if len(coords) % 2:
        return None, coords
    return [(coords[i], coords[i + 1]) for i in
            range(0, len(coords) - 1, 2)], coords


def column(verts, tris, x, y, z0, z1):
    return cuts(verts, tris, (x, y, z0), (x, y, z1), 2)


def row(verts, tris, x, z, y0, y1):
    return cuts(verts, tris, (x, y0, z), (x, y1, z), 1)


def main():
    memcap.cap()
    path, angle = K.cli(sys.argv[1:] or [None], None)
    if path is None:
        print(__doc__)
        return 2
    angle = 45.0 if angle is None else angle

    rep = K.Report("scoopula")
    objects, _, _ = K.load_3mf(path, rep)
    verts, tris = objects[0]["verts"], objects[0]["tris"]
    shells = K.shell_triangles(verts, tris)
    meta = K.read_params(path)

    wrong = {k: (meta.get(k), v) for k, v in IN_SPEC.items()
             if meta.get(k) is None or abs(float(meta[k]) - v) > 1e-9}
    rep.check("the file's own metadata is the specification", not wrong,
              f"stated {sorted(wrong.items())}" if wrong else
              f"{len(IN_SPEC)} spec params read back out of the part")
    P = {k: float(meta.get(k, v)) for k, v in IN_SPEC.items()}
    L, HL = P["scoop_len"], P["hnd_l"]
    WALL, LEAN, RIM_W = P["wall"], P["lean"], P["rim_w"]
    DEPTH, HT, GRID = P["depth"], P["hnd_t"], P["grid"]
    NOZ, LIP = P["nozzle"], P["h_lip"]
    FLOOR, BACK = P["floor"], P["back_wall"]
    inset = WALL * math.hypot(1.0, LEAN)           # the cut's horizontal offset
    cos_lean = 1.0 / math.hypot(1.0, LEAN)

    xs = [v[0] for v in verts]
    ys = [v[1] for v in verts]
    zs = [v[2] for v in verts]
    x_lo, x_hi = min(xs), max(xs)
    y_lo, y_hi = min(ys), max(ys)
    z_bot, z_top = min(zs), max(zs)

    # --- one solid, and sound ---------------------------------------------
    rep.check("the handle and the trough are one solid", len(shells) == 1,
              f"{len(shells)} shells in the mesh")
    gaps = csg.real_gaps(verts, tris)
    rep.check("no hole in the surface", not gaps,
              f"{len(gaps)} edges used once with nothing along them" if gaps else
              f"{csg.crack_length(verts, tris):.2e} mm of edge used once, all of "
              f"it in T-junctions")
    vol = K.signed_volume(verts, tris)
    rep.check("the winding is outward", vol > 0.0, f"{vol:.1f} mm^3")

    # --- the overall shape, off the mesh ----------------------------------
    rep.near("the back is the flat of the part, on the bed", z_bot, 0.0, 1e-6)
    rep.near("the part is as long as the trough plus the handle",
             x_hi - x_lo, HL + L, TOL)
    rep.near("the part is as wide as the trough's rim", y_hi - y_lo,
             RIM_W + 2.0 * inset, TOL)
    rep.near("the part is as tall as the floor plus the channel depth", z_top,
             FLOOR + DEPTH, 0.01 + GRID)
    bed = PR.bed_contact(verts, tris)
    rep.check("it stands on the bed", bed > 800.0,
              f"{bed:.0f} mm^2 of contact with the bed, on a flat back")

    # --- an independent volume: integrate the columns over the footprint ---
    n_odd, n_cols = 0, 0
    stations = [x_lo + (x_hi - x_lo) * i / 24.0 for i in range(25)]
    areas = []
    for x in stations:
        tot, n = 0.0, 0
        y = y_lo - 0.5
        while y <= y_hi + 0.5:
            sp, _ = column(verts, tris, x + PROBE[0], y + PROBE[1],
                           z_bot - 1.0, z_top + 1.0)
            n_cols += 1
            if sp is None:
                n_odd += 1
            else:
                tot += sum(b - a for a, b in sp)
                n += 1
            y += Y_STEP
        areas.append(tot / max(1, n) * (y_hi - y_lo + 1.0))
    col_vol = sum(0.5 * (areas[i] + areas[i + 1]) * (stations[i + 1] - stations[i])
                  for i in range(len(stations) - 1))
    rep.check("every test column crosses the surface an even number of times",
              n_odd == 0,
              f"{n_odd} of {n_cols} columns met an odd crossing count"
              if n_odd else f"{n_cols} columns crossed the surface evenly: no "
                            f"grazing, no open edge, no parity lie")
    rep.check("the volume stands up to a column integration",
              abs(col_vol - vol) < 0.04 * vol,
              f"mesh {vol:.1f} mm^3, columns {col_vol:.1f} mm^3 "
              f"({100.0 * (col_vol - vol) / vol:+.1f}%)")

    # --- where is the trough open? cut it across, above the floor ---------
    z_cut = FLOOR + 1.0                       # down in the round bottom
    z_wall = FLOOR + DEPTH - 0.8              # up in the flanks, where a wall is
    open_x, solid_x, walls, mouths, odd = [], [], [], [], []
    flanks = []
    x = 0.5
    while x < L:
        sp, cr = row(verts, tris, x + PROBE[0], z_cut, y_lo - 0.5, y_hi + 0.5)
        if sp is None:
            odd.append(x)
        elif len(sp) == 2:                       # two walls: a trough
            open_x.append(x)
            walls.append(min(b - a for a, b in sp) * cos_lean)
            mouths.append((x, sp[1][0] - sp[0][1]))
        elif len(sp) == 1:                       # one band: solid heel or toe
            solid_x.append(x)
        x += 0.5
    rep.check("nothing grazed a face while cutting the trough across",
              not odd, f"{len(odd)} sections read oddly at "
                       f"{[round(o, 1) for o in odd[:4]]}" if odd else
                       f"{int(L / 0.5)} sections of the trough cut at "
                       f"z={z_cut:g} mm")
    cav_from, cav_to = (min(open_x), max(open_x)) if open_x else (0.0, 0.0)
    cav_len = cav_to - cav_from
    rep.check("the scoop is longer than it is wide, along the handle's axis",
              cav_len > 2.0 * RIM_W,
              f"open for {cav_len:.1f} mm along x against a {RIM_W:g} mm mouth: "
              f"a scoopula's trough and handle share one axis")
    rep.check("the trough is closed at the heel",
              any(s < cav_from - 0.5 for s in solid_x),
              f"{cav_from:.1f} mm of solid heel block from x=0")
    rep.check("the trough is closed at the toe",
              any(s > cav_to + 0.5 for s in solid_x),
              f"{L - cav_to:.1f} mm of solid toe wedge")
    # The leaf shape belongs on the bed, where the trough's width is the band
    # it stands on: down in the arc the channel is nearly the same width all
    # the way along, and a test there measures the loft's wobble, not the shape.
    sil = []
    for xb in [0.5 + 0.5 * i for i in range(int(L / 0.5))]:
        sp, _ = row(verts, tris, xb + PROBE[0], 0.8, y_lo - 0.5, y_hi + 0.5)
        if sp and len(sp) == 1:
            sil.append((xb, sp[0][1] - sp[0][0]))
    w_max = max((w for _, w in sil), default=0.0)
    mid = [w for xb, w in sil if abs(xb - 0.5 * L) < 1.0]
    toe = [w for xb, w in sil if xb > 0.92 * L]
    heel = [w for xb, w in sil if xb < 0.08 * L]
    rep.check("the trough is widest through its middle and runs out at the "
              "toe: a leaf, not a rectangle",
              sil and mid and toe and w_max - max(mid) < 0.15 and
              max(toe) < 0.70 * w_max,
              f"bed band {w_max:.2f} mm at its widest, "
              f"{max(mid):.2f} mm through the middle, {min(heel or [0]):.2f} mm "
              f"at the heel where it merges into the handle and {max(toe):.2f} "
              f"mm at the toe it runs out to" if sil and mid else
              "no section stood on one band at bed level")
    rep.check("the channel is never cut to a slit the nozzle cannot make",
              mouths and min(m[1] for m in mouths) >= 2.0 * NOZ,
              f"narrowest open channel {min(m[1] for m in mouths):.2f} mm at "
              f"z={z_cut:g} against a {NOZ:g} mm nozzle" if mouths else
              "no open section read two walls")
    x = 0.5
    while x < L:
        sp, _ = row(verts, tris, x + PROBE[0], z_wall, y_lo - 0.5, y_hi + 0.5)
        if sp and len(sp) == 2:
            flanks.append(min(b - a for a, b in sp) * cos_lean)
        x += 0.5
    rep.check("the flanks are the wall asked for, across themselves",
              len(flanks) > 10 and max(abs(w - WALL) for w in flanks) < TOL,
              f"{len(flanks)} sections cut at z={z_wall:.1f} mm, "
              f"{min(flanks):.3f} to {max(flanks):.3f} mm across against "
              f"{WALL:g} mm asked" if flanks else "no flank section read two "
              "walls at any height")
    rep.check("the wall is never thinner than asked, round bottom included",
              walls and min(walls) >= WALL - 0.05,
              f"{len(walls)} sections cut at z={z_cut:.1f} mm, thinnest "
              f"{min(walls):.3f} mm across: thicker than the wall is asked, "
              f"which is what a tangent arc leaves at the corners"
              if walls and min(walls) >= WALL - 0.05 else
              f"thinnest {min(walls):.3f} mm against {WALL:g} mm asked")

    # --- the heel closed where the powder is, not just under the floor ----
    z_back = z_bot + 0.75 * (z_top - z_bot)
    leaky = []
    for xb in [0.5, 1.0, 1.5, 2.0, BACK - 0.2]:
        sp, _ = row(verts, tris, xb + PROBE[0], z_back, y_lo - 0.5, y_hi + 0.5)
        if not sp or len(sp) != 1:
            leaky.append((round(xb, 2), "nothing there" if sp == [] else
                         f"{len(sp)} bands" if sp else "odd crossings"))
    rep.check("the heel is solid at the height the trough carries powder at",
              not leaky and cav_from >= BACK - 0.6,
              f"open back found at {leaky}" if leaky else
              f"solid in one band at z={z_back:.1f} mm for {cav_from:.1f} mm of "
              f"heel, and the channel does not start until past the "
              f"{BACK:g} mm back wall")

    # --- the nose: the channel runs out to a thin edge --------------------
    sp_nose, _ = row(verts, tris, x_hi - 0.4 + PROBE[0], 0.8,
                     y_lo - 0.5, y_hi + 0.5)
    rep.check("the trough converges to a nose much narrower than its mouth",
              sp_nose and len(sp_nose) == 1 and
              (sp_nose[0][1] - sp_nose[0][0]) < 0.55 * RIM_W,
              f"{sp_nose[0][1] - sp_nose[0][0]:.2f} mm across at the toe against "
              f"a {RIM_W:g} mm mouth, at {LIP:g} mm of height" if sp_nose and
              len(sp_nose) == 1 else f"{sp_nose and len(sp_nose or [])} bands "
              "at the toe")

    # --- the cavity is OPEN: fire every column upward ---------------------
    blocked = []
    for x in [cav_from + i * cav_len / 11.0 for i in range(12)]:
        sp, _ = row(verts, tris, x + PROBE[0], z_cut, y_lo - 0.5, y_hi + 0.5)
        if not sp or len(sp) != 2:
            blocked.append((round(x, 1), "no mouth here"))
            continue
        gap = (sp[0][1] + sp[1][0]) / 2.0          # the mouth's centre
        for fr in (0.0, 0.35, -0.35, 0.65, -0.65):
            y = gap + fr * 0.5 * (sp[1][0] - sp[0][1])
            col, cr = column(verts, tris, x + PROBE[0], y + PROBE[1],
                             z_bot - 1.0, z_top + 1.0)
            if col is None or len(col) != 1:
                blocked.append((round(x, 1), round(y, 1),
                                f"{len(col) if col else len(cr)} bands"))
            elif col[0][1] > z_top - 0.5 * DEPTH:
                blocked.append((round(x, 1), round(y, 1),
                                f"material up to {col[0][1]:.2f} mm, where the "
                                f"cavity should stand empty to {z_top:.1f}"))
    rep.check("nothing stands over the cavity, at 60 columns", not blocked,
              f"roof found at {blocked[:3]}" if blocked else
              "every column inside the mouth meets the bed, then the floor, then "
              "air: an open trough, not a sealed void with a roof")

    # --- the mouth: as wide as asked, and wider the higher you cut --------
    z_rim = z_top - 0.05
    profile = []
    for zz in (FLOOR + 0.4, FLOOR + 1.0, FLOOR + 2.0, z_wall, z_rim):
        sp, _ = row(verts, tris, 0.5 * L + PROBE[0], zz, y_lo - 0.5,
                     y_hi + 0.5)
        if sp and len(sp) == 2:
            profile.append((zz, sp[1][0] - sp[0][1], sp[0][1] - sp[0][0]))
    if len(profile) >= 4:
        m_lo, m_hi = profile[0][1], profile[-1][1]
        rep.near("the opening is the width asked for, at the rim", m_hi,
                 RIM_W, TOL)
        rep.check("the channel opens upward at every height, so no face of it "
                  "needs support",
                  all(b[1] > a[1] for a, b in zip(profile, profile[1:])),
                  " -> ".join(f"{w:.2f} mm @ z={z:.1f}" for z, w, _ in profile))
        rate_low = (profile[1][1] - profile[0][1]) / (profile[1][0] -
                                                     profile[0][0]) / 2.0
        rate_hi = (profile[-1][1] - profile[-2][1]) / (profile[-1][0] -
                                                      profile[-2][0]) / 2.0
        rep.check("the bottom is a round channel and not a V",
                  m_lo > 0.30 * m_hi and rate_low > 3.0 * LEAN,
                  f"{m_lo:.2f} mm across just {profile[0][0] - FLOOR:.1f} mm "
                  f"above the floor, a third of the {m_hi:.2f} mm mouth: a V "
                  f"would be {(profile[-1][1] * 0.4 / DEPTH):.2f} mm there, and "
                  f"the width grows {rate_low:.2f} mm per mm of rise down there "
                  f"against {rate_hi:.2f} up in the flanks")
        rep.near("the flanks lean outward at the rate asked for", rate_hi,
                 LEAN, 0.06)
        rep.near("the rim stands as a thin land, the wall's own secant",
                 profile[-1][2], inset, TOL)
    else:
        rep.check("the mouth could be profiled at five heights", False,
                  f"{len(profile)} of 5 heights read two walls")

    # --- the floor: one band, the thickness asked for, from the bed up ----
    floors = []
    for x in [cav_from + i * cav_len / 5.0 for i in range(6)]:
        col, _ = column(verts, tris, x + PROBE[0], PROBE[1], z_bot - 1.0,
                        z_top + 1.0)
        if col and len(col) == 1:
            floors.append(col[0])
    rep.check("the floor is a single band, on the bed, under the whole cavity",
              len(floors) == 6 and all(a < TOL for a, _ in floors),
              f"{len(floors)} columns read; lowest underside "
              f"{min((a for a, _ in floors), default=float('nan')):.3f} mm")
    rep.near("the floor is the thickness asked for, bed face to channel",
             max(b for _, b in floors[:4]), FLOOR, 0.01 + GRID)
    rep.check(f"the floor is a whole number of {GRID:g} mm layers",
              abs(FLOOR / GRID - round(FLOOR / GRID)) < 1e-9,
              f"{FLOOR:g} mm of floor is {FLOOR / GRID:g} layers")

    # --- the handle: in line, on the bed, the size asked for --------------
    rep.check("the handle is beyond the heel, on the far side from the toe",
              x_lo < -0.8 * HL and x_hi > 0.95 * L,
              f"tip at x={x_lo:.1f}, heel at x=0, toe at x={x_hi:.1f}")
    for x_at, want, name in ((x_lo + 1.0, P["hnd_w_tip"], "free end"),
                             (-2.0, P["hnd_w"], "joint")):
        sp, _ = row(verts, tris, x_at + PROBE[0], 0.8, y_lo - 0.5, y_hi + 0.5)
        if sp and len(sp) == 1:
            a, b = sp[0]
            rep.near(f"the handle is the width asked for at its {name}",
                     b - a, want, TOL)
            rep.near(f"the handle is centred on the trough's axis at its {name}",
                     0.5 * (a + b), 0.0, 0.05)
        else:
            rep.check(f"the handle is there at its {name}", False,
                      f"{len(sp) if sp else 'odd'} bands across at x={x_at:.1f}")
    col, _ = column(verts, tris, -0.5 * L + PROBE[0], PROBE[1], z_bot - 1.0,
                    z_top + 1.0)
    if col and len(col) == 1:
        rep.near("the handle is the thickness asked for", col[0][1] - col[0][0],
                 HT, TOL)
        rep.near("the handle lies on the bed with the trough", col[0][0], 0.0,
                 TOL)
        hnd_top = col[0][1]
        floor_top = max(b for _, b in floors[:4]) if floors else 0.0
        rep.check("the handle is thinner than the floor, so none of it stands "
                  "in the channel as a fin", HT < FLOOR,
                  f"{HT:g} mm of handle under a {FLOOR:g} mm floor")
        rep.check("and its top face stands a clear bead below the channel's "
                  "floor, so the joint has a step and not a coincidence",
                  floor_top - hnd_top >= NOZ,
                  f"handle top measured at z={hnd_top:.3f}, channel floor at "
                  f"z={floor_top:.3f}: {(floor_top - hnd_top) * 1000:.0f} um of "
                  f"step between them")
    else:
        rep.check("the handle is there mid-grip", False,
                  f"{col and len(col)} bands up at x={-0.5 * L:.1f}")
    lip, _ = column(verts, tris, x_hi - 0.4 + PROBE[0], PROBE[1], z_bot - 1.0,
                    z_top + 1.0)
    rep.near("the scraping lip stands the height asked for at the toe",
             lip[0][1] if lip and len(lip) == 1 else 0.0, LIP, 0.01 + GRID)
    rep.check("the lip's lowest edge is the bed plane, so it enters a pile at "
              "floor level", lip and lip[0][0] < TOL,
              f"toe column runs from z={lip[0][0]:.2f}" if lip else "no toe")

    # --- printing ----------------------------------------------------------
    flagged = PR.overhangs(verts, tris, shells=shells, exposed_only=True,
                           max_angle=angle)
    real = [f for f in flagged if f["area"] > 0.05]
    sliver = [f for f in flagged if f["area"] <= 0.05]
    # How much margin the part leaves, which is its own property and not the
    # machine's: the shallowest exposed face that looks down at all is the
    # printer limit this part stops clearing without support.
    every = PR.overhangs(verts, tris, shells=shells, exposed_only=True,
                         max_angle=89.0)
    dips = [90.0 - math.degrees(math.asin(min(1.0, -f["n"][2])))
            for f in every if f["area"] > 0.05 and f["n"][2] < 0.0]
    rep.check(f"no face a slicer has to carry is shallower than {angle:g} deg",
              not real,
              f"{len(real)} faces, {sum(f['area'] for f in real):.1f} mm^2, "
              f"lowest at z={real[0]['z']:.2f}" if real else
              f"the steepest exposed face looking down is {min(dips):.1f} deg "
              f"from horizontal, so support starts only on a printer shallower "
              f"than that; {len(sliver)} boolean slivers of "
              f"{sum(f['area'] for f in sliver):.1e} mm^2 are noise")
    off = []
    for (n, a, zmin, _zmax) in PR.triangle_normals(verts, tris):
        if a > 1e-6 and abs(n[2]) > 0.999:
            k = zmin / GRID
            if abs(k - round(k)) * GRID > 0.01:
                off.append((zmin, a))
    rep.check("every horizontal feature face lies on the layer grid", not off,
              f"{len(off)} horizontal faces off the {GRID:g} mm grid, nearest "
              f"z={off[0][0]:.3f}" if off else
              "bed, floor, rim land, handle top and cutter exits all land on "
              "whole layers (and on 0.2 mm ones too); vertices off the grid sit "
              "only on sloped faces, where a T-junction split has no height of "
              "its own to keep")
    rep.check("nothing is under two beads",
              all(v >= 2.0 * NOZ for v in (WALL, FLOOR, HT, LIP)),
              f"wall {WALL:g} mm = {WALL / NOZ:.1f} beads, floor {FLOOR:g} mm = "
              f"{FLOOR / NOZ:.1f}, handle {HT:g} mm = {HT / NOZ:.1f}, lip "
              f"{LIP:g} mm = {LIP / NOZ:.1f}, rim land {inset:.2f} mm = "
              f"{inset / NOZ:.1f}")

    print(f"  part size         : {x_hi - x_lo:.1f} x {y_hi - y_lo:.1f} x "
          f"{z_top - z_bot:.1f} mm")
    print(f"  open cavity       : x {cav_from:.1f} to {cav_to:.1f} mm, "
          f"{cav_len:.1f} mm of mouth, floor {FLOOR:g} mm above the bed")
    print(f"  overhang angle    : {angle:g} deg from horizontal "
          f"(the printer's limit, not the part's)")
    print(f"  vertices / tris   : {len(verts)} / {len(tris)}")
    print(f"  volume            : {vol:.1f} mm^3")

    # renders for human eyes, not assertions
    K.render(verts, tris, "scoopula_top.png", view="xy", su=6.0)
    K.render(verts, tris, "scoopula_side.png", view="xz", su=6.0, flip_v=True)
    print("  renders           : scoopula_top.png, scoopula_side.png")
    return 0 if rep.finish() else 1


if __name__ == "__main__":
    sys.exit(main())

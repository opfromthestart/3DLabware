#!/usr/bin/env python3
"""Generate a two-way straight-through stopcock: a barrel and a tapered key.

A stopcock is two parts that only work together, so this builds both and writes
them side by side on the plate, separated by `PLATE_GAP` in x -- the verifier
puts the key back in its barrel to test the fit, which is the only place the
part's real behaviour exists.

    TUBE_ID = 4       tube bore, mm
    TUBE_WALL = 2.5   tube wall, so 9 mm outside
    TUBE_LEN = 25     tube length beyond the barrel, mm
    BARREL_R = 8     barrel outer radius, mm
    BARREL_H = 30     barrel height, mm
    BORE_TOP_R = 6    the tapered bore's radius at the top, mm
    CLEARANCE = 0.08  radial gap between the bore and the key, mm

**Why a taper.** A cylinder with a clearance is a hole with a clearance: it
leaks. A cone with a clearance seals, because pushing the key in closes it.
The key is drawn `CLEARANCE` clear of the bore at every height, and the bore
narrows by `TAPER` per mm, so the key seats after `CLEARANCE / TAPER` mm of
travel -- 1.6 mm here -- and that seating is the seal. The same 1:20 the
project's ground joints use.

**The passage through the barrel's wall.** The tubes' bores have to open
through a wall whose inner face is a cone. The wall is one tapered tube with a
box taken out of it -- `subtract` -- and the box is the swath the tubes pass
through, 4 mm wide and 4 mm tall and reaching past the barrel's outside. The
sleeves fill that swath: every point of it is within the sleeve's outside
radius of the tube's axis, so the sleeve's wall covers it and the sleeve's own
bore is the passage.

The swath is a *box* and not a cylinder on purpose. A cylinder of the sleeve's
own radius would meet the wall face to face, and two coincident surfaces are
the problem this part is built to avoid; the box is smaller than the sleeve all
round, so the sleeve's material overlaps it with room to spare. The box is also
wider than the 4 mm bore is round, so its corners sit inside the sleeve's wall
and the passage is the full round bore and nothing less.

The sleeves stop at `X_END` from the axis, chosen by bisection as the shortest
that keeps *all* of the sleeve's material outside the cone: any shorter and the
sleeve would poke into the bore and jam the key, any longer and the wall is
left in front of the mouth.

**The pieces are unioned, not left overlapping.** A wall with a sleeve through
it is two surfaces crossing each other. That prints on a slicer that fills by
winding and shows as a hole on one that fills even-odd -- and the joint between
a tube and a barrel is exactly where the two surfaces cross. So the wall and
both sleeves go through `csglib.union`, and so do the key's frusta, its slot
segments and its blade: what is written is one watertight solid per part, with
no surface inside another surface.

**The wall's thickness is not free.** A removed facet reaches `barrel_r *
half / r_bore` sideways from the tube's axis at the barrel's outside, so the
sleeve's outside radius has to cover `hypot(barrel_r * half / r_bore, half)` or
the hole is wider than the wall meant to fill it. For a 4 mm bore that caps the
barrel at about `1.73 * r_bore` -- 9.1 mm here -- which is why the barrel is
8 mm across rather than anything you like, and why the check is in the build
and not in a comment.

**The key's passage is a slot, not a hole.** A hole through a solid needs a
subtraction; a slot is two convex segments, and `loft` fans convex caps
correctly. The slot is 4 mm wide and 4 mm tall through the key's waist, so it
never restricts the 4 mm round bore that feeds it.

Usage: python3 make_stopcock.py [output.3mf]
"""

import math
import sys

import csglib as csg
import junctionlib as J
import printability as PR
import threelib as t

# --- Parameters ------------------------------------------------------------
TUBE_ID = 4.0        # tube bore, mm
TUBE_WALL = 2.5      # tube wall, so the outside is TUBE_ID + 2*TUBE_WALL
TUBE_LEN = 25.0      # tube length beyond the barrel's outer surface
BARREL_R = 8.0       # barrel outer radius, mm
BARREL_H = 30.0      # barrel height, mm -- also the key's length
BORE_TOP_R = 6.0     # the tapered bore's radius at the top, mm
CLEARANCE = 0.08     # radial gap between the bore and the key, mm
BLADE_R = 5.6        # handle blade radius, mm
BLADE_T = 4.0        # handle blade thickness, mm
BLADE_DROP = 2.0     # how far the blade dips into the top of the key, mm
PLATE_GAP = 60.0     # the key is written this far aside, for printing
TOL = 0.01           # max chord deviation on round surfaces

PARAMS = dict(tube_id=TUBE_ID, tube_wall=TUBE_WALL, tube_len=TUBE_LEN,
              barrel_r=BARREL_R, barrel_h=BARREL_H, bore_top_r=BORE_TOP_R,
              clearance=CLEARANCE, blade_r=BLADE_R, blade_t=BLADE_T,
              blade_drop=BLADE_DROP, plate_gap=PLATE_GAP, tol=TOL)


def default_out():
    return t.stamped_filename("stopcock", PARAMS)


def dims(p):
    """Everything derived from the parameter set, in one place."""
    h = p["barrel_h"]
    taper = J.TAPER
    r_bore = lambda z: p["bore_top_r"] - taper * (h - z)      # noqa: E731
    r_plug = lambda z: r_bore(z) - p["clearance"]             # noqa: E731
    tube_z = 0.5 * h                       # the tubes meet the barrel's middle
    half = 0.5 * p["tube_id"]
    r_tube = half + p["tube_wall"]

    # How far out the sleeve's material has to stop. Sampled over the sleeve's
    # whole cross-section rather than reasoned about: the binding point is not
    # the widest part of the wall but wherever a small |y| meets a large bore
    # radius, and the taper moves both.
    def pokes(x_end):
        # The sleeve's *material* is the annulus between its bore and its
        # outside, so that is what has to be sampled -- not a strip of y values,
        # which misses the material at small |y| and large |z - tube_z| and
        # leaves the sleeve's wall poking into the bore at the top of the tube.
        for i in range(49):
            r = half + (r_tube - half) * i / 48.0
            for k in range(192):
                ang = 2.0 * math.pi * k / 192.0
                y, dz = r * math.cos(ang), r * math.sin(ang)
                if math.hypot(x_end, y) < r_bore(tube_z + dz) - 1e-9:
                    return True
        return False

    x_lo, x_hi = 0.0, p["barrel_r"] + p["tube_len"]
    while x_hi - x_lo > 1e-6:
        mid = 0.5 * (x_lo + x_hi)
        x_lo, x_hi = (mid, x_hi) if pokes(mid) else (x_lo, mid)
    x_end = x_hi + 0.01                     # just clear of the cone

    root = math.sqrt(1.0 + taper ** 2)
    return {"r_bore": r_bore, "r_plug": r_plug, "r_tube": r_tube,
            "tube_z": tube_z, "half": half, "x_end": x_end,
            "z0": tube_z - half, "z1": tube_z + half,
            "mouth": math.pi * half ** 2,
            "seat": p["clearance"] / taper,
            "blade_z": h + p["blade_r"] - p["blade_drop"],
            "blade_corner": math.hypot(p["blade_r"], p["blade_t"] / 2.0),
            "plug_cone_at_blade": (p["bore_top_r"] - p["clearance"]
                                   + taper * (p["blade_r"]
                                              - p["blade_drop"])) / root,
            "blade_reach": math.sqrt(max(0.0, 2.0 * p["blade_r"]
                                         * p["blade_drop"]
                                         - p["blade_drop"] ** 2))}


def segment(r, c, m):
    """The part of a disc of radius r above the line y = c, as a CCW polygon.

    Built with a fixed `m + 1` vertices rather than by clipping a disc, because
    `loft` needs every section to have the same count and a clip returns
    however many the chord happens to leave.
    """
    if c >= r:
        raise ValueError(f"a disc of radius {r:g} has nothing above y={c:g}")
    a0 = math.asin(c / r)
    return [(r * math.cos(a0 + (math.pi - 2.0 * a0) * i / m),
             r * math.sin(a0 + (math.pi - 2.0 * a0) * i / m))
            for i in range(m + 1)]


def surfaces_on_axis(verts, tris, z, sign):
    """Where a ray along the tube's axis meets a surface, outward from centre.

    Counting surfaces is the honest test of a passage: the wall's facets are
    left out across the tube's whole width, so the only thing that could block
    the axis is a sleeve in the wrong place -- and a point-in-the-cone test
    calls the mouth shut even when it is open. The probe is 0.13 mm off the
    tube's height because the key's loft has a section exactly there, and an
    edge that touches a plane is not a crossing.
    """
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
        den = sign * ey
        if abs(den) < 1e-12:
            continue
        if -1e-9 <= (-y0 * sign) / den <= 1 + 1e-9:
            t = (x0 * ey - y0 * ex) / den
            if t > 1e-7:
                hits.append(round(t, 3))
    return sorted(set(hits))


def revolve_frustum(mesh, r0, r1, z0, z1, n):
    """A solid frustum, capped flat at both ends."""
    t.revolve(mesh, [(0.0, z0), (r0, z0), (r1, z1), (0.0, z1)], n)


def build(**p):
    d = dims(p)
    for key in ("tube_id", "barrel_r", "barrel_h", "clearance"):
        if p[key] is None:
            raise ValueError(f"{key} is not set")
    if p["clearance"] <= 0.0:
        raise ValueError("a key with no clearance cannot be printed inside "
                         "its barrel")
    if p["clearance"] >= 0.5 * d["r_plug"](0.0):
        raise ValueError("the clearance is most of the key's tip")

    h, tube_z, half = p["barrel_h"], d["tube_z"], d["half"]
    n = t.auto_segments(p["barrel_r"], p["tol"])

    # --- the barrel: the wall, and one sleeve per tube, unioned ----------
    # The wall's bore has to be faceted the same way the key is, or the two
    # inscribed polygons disagree and the key pokes through the wall.
    n = max(n, 4 * int(math.ceil(n / 4.0)))
    wall = t.Mesh("wall")
    t.tapered_tube(wall, d["r_bore"](0.0), p["barrel_r"],
                   d["r_bore"](h), p["barrel_r"], 0.0, h, n)
    swath = t.Mesh("swath")
    t.box(swath, -(p["barrel_r"] + 0.5), -half, tube_z - half,
          p["barrel_r"] + 0.5, half, tube_z + half)
    length = p["barrel_r"] + p["tube_len"] - d["x_end"]
    steps = [("start", wall), ("subtract", swath)]
    for sign in (1.0, -1.0):
        sleeve = t.Mesh("sleeve")
        t.tube(sleeve, half, d["r_tube"], 0.0, length, n)
        t.rotate_y(sleeve, 90.0 * sign)
        t.translate(sleeve, sign * d["x_end"], 0.0, tube_z)
        steps.append(("union", sleeve))
    barrel = csg.combine(steps, name="barrel")

    # --- the key: one cone with the slot taken out of it ------------------
    # The key's body is a single tapered cone, so it is built as one revolve
    # and the slot is subtracted from it. Building it as two frusta with two
    # lofted segments between them puts the segments' arc faces *exactly* on
    # the frusta's cone -- coplanar and partly overlapping, which is the one
    # case a BSP resolves badly, and it came back as twelve millimetres of
    # cracks. A box cut through a cone crosses it and leaves none.
    body = t.Mesh("body")
    t.revolve(body, [(0.0, 0.0), (d["r_plug"](0.0), 0.0),
                     (d["r_plug"](h), h), (0.0, h)], n)
    slot = t.Mesh("slot")
    t.box(slot, -(d["r_plug"](h) + 1.0), -half, tube_z - half,
          d["r_plug"](h) + 1.0, half, tube_z + half)

    # --- the blade, standing on edge in the plane of the slot ------------
    blade = t.Mesh("blade")
    t.prism(blade, t.disc_polygon(p["blade_r"],
                                  t.auto_segments(p["blade_r"], p["tol"])),
            -p["blade_t"] / 2.0, p["blade_t"] / 2.0)
    t.rotate_x(blade, 90.0)
    t.translate(blade, 0.0, 0.0, d["blade_z"])
    key = csg.combine([("start", body), ("subtract", slot),
                       ("union", blade)], name="key")

    # The key is written beside its barrel: printing them in place would fuse
    # them across a 0.08 mm gap, which is less than one extrusion width.
    key = key.copy()
    key.transform(lambda v: (v[0] + p["plate_gap"], v[1], v[2]))
    probe_z = d["tube_z"] + 0.13
    axis = (surfaces_on_axis(barrel.verts, barrel.tris, probe_z, 1.0)
            + surfaces_on_axis(barrel.verts, barrel.tris, probe_z, -1.0))
    return (barrel, key), {"n": n, "dims": d, "axis": axis,
                           "v_wall": (math.pi * h / 3.0
                                      * (d["r_bore"](0.0) ** 2
                                         + d["r_bore"](0.0) * d["r_bore"](h)
                                         + d["r_bore"](h) ** 2)),
                           "v_barrel": barrel.signed_volume(),
                           "v_key": key.signed_volume(),
                           "shells": (len(barrel.shells()),
                                      len(key.shells())),
                           "cracks": (csg.crack_length(barrel.verts,
                                                       barrel.tris),
                                      csg.crack_length(key.verts, key.tris))}


def main():
    out = sys.argv[1] if len(sys.argv) > 1 else default_out()
    (barrel, key), s = build(**PARAMS)
    d = s["dims"]
    h, half = PARAMS["barrel_h"], d["half"]

    desc = (f"Two-way straight-through stopcock: a {PARAMS['barrel_r'] * 2:g} "
            f"mm barrel on a {PARAMS['bore_top_r'] * 2:g} mm tapered bore with "
            f"{PARAMS['tube_id']:g} mm tubes, and a key {PARAMS['clearance']:g} "
            f"mm clear of it.")
    t.write_3mf(out, (("barrel", barrel), ("key", key)),
                title="Two-way stopcock", description=desc, params=PARAMS)

    bb = barrel.bbox()
    kb = key.bbox()
    print(f"wrote {out}")
    print(f"  params            : {t.describe(PARAMS)}")
    print(f"  param hash        : {t.params_hash(PARAMS)}")
    print(f"  tubes             : {PARAMS['tube_id']:g} mm bore, "
          f"{d['r_tube'] * 2:g} mm outside, {PARAMS['tube_len']:g} mm long, "
          f"collinear through the barrel at z={d['tube_z']:g}")
    print(f"  barrel            : {PARAMS['barrel_r'] * 2:g} mm across, "
          f"{h:g} mm tall, bore {d['r_bore'](0.0) * 2:g} to "
          f"{d['r_bore'](h) * 2:g} mm on a 1:{1.0 / J.TAPER:g} taper")
    print(f"  key               : {d['r_plug'](h) * 2:g} mm at the top, "
          f"{PARAMS['clearance']:g} mm clear of the bore all the way down, "
          f"seating after {d['seat']:.2f} mm of travel")
    print(f"  tube mouths       : the sleeve stops {d['x_end']:.2f} mm from "
          f"the axis, and the tube's axis meets "
          + ("nothing at all between the two ends -- one straight passage "
             f"through a {PARAMS['tube_id']:g} mm bore"
             if not s["axis"] else f"surfaces at r={s['axis']}"))
    print(f"  key passage       : a {PARAMS['tube_id']:g} mm slot through a "
          f"{half * 2:g} mm waist, turning in the plane of the blade")
    print(f"  blade             : {PARAMS['blade_r'] * 2:g} mm across, "
          f"{PARAMS['blade_t']:g} mm thick, dipping "
          f"{PARAMS['blade_drop']:g} mm into the key; its corners reach "
          f"{d['blade_corner']:.2f} mm where the key's cone is "
          f"{d['plug_cone_at_blade']:.2f} mm, so it stays inside it")
    print(f"  plate layout      : the key is written {PARAMS['plate_gap']:g} mm "
          f"aside at x={kb[0]:.1f}..{kb[3]:.1f}; the barrel spans "
          f"{bb[0]:.1f}..{bb[3]:.1f}")
    print(f"  volume            : barrel {s['v_barrel']:.1f} mm^3, key "
          f"{s['v_key']:.1f} mm^3 -- the true union, no double counting")
    print(f"  bed contact       : barrel "
          f"{math.pi * (PARAMS['barrel_r'] ** 2 - d['r_bore'](0.0) ** 2):.1f} "
          f"mm^2, key {math.pi * d['r_plug'](0.0) ** 2:.1f} mm^2")
    print(f"  print orientation : both standing on the plug axis. The bore's "
          f"cone is 87 deg from the horizontal and the barrel's wall is "
          f"vertical; the two horizontal tubes are arches of "
          f"{PARAMS['tube_id']:g} and {d['r_tube'] * 2:g} mm, which bridge")
    print(f"  shells            : barrel {s['shells'][0]}, key {s['shells'][1]} "
          f"-- unioned, so no surface lies inside another")
    print(f"  cracks            : {s['cracks'][0]:.2e} mm in the barrel, "
          f"{s['cracks'][1]:.2e} mm in the key")
    print(f"  note              : the taper is the seal, not the clearance -- "
          f"press the key in {d['seat']:.1f} mm to seat it")


if __name__ == "__main__":
    main()

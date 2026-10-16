#!/usr/bin/env python3
"""Verify a water aspirator 3MF by measuring the file that was written.

Nothing here asks the build script what it thinks it made. The specification is
restated below as its own literals, the file is parsed out of its XML, and every
number comes off the mesh with rays. Where the two disagree, the part is wrong
-- or this file is, which amounts to the same thing.

The three structural demands, and how each is measured:

  * the throat starts before the gas inlet and ends after it: the narrow span
    is found from the measured bore, then compared against the inlet's own
    footprint;
  * the wall reads as a Z along its length: the outside leans out from the
    inlet tube to the chamber, ahead of the inside, which draws in from the
    bore to the throat, and across that join the section is one continuous
    band of plastic with no void between the two walls;
  * around the gas inlet the section is two rings, not one: a ray leaving the
    axis at the inlet's height must cross four surfaces -- a nozzle wall and a
    shell wall, with air between them.

Rays are the only honest way to measure a bore. A facetted circle puts a vertex
every ninth degree, and a ray running down a vertex column crosses each surface
twice and reads as empty space, so every probe here stands half a millimetre
off the axis and off the whole millimetres of the station it claims to measure.
`checks.inside`, which answers a single point, fails the same way: it was used
once here to confirm the air inlet was drilled, reported the hole as solid, and
was wrong -- a point sitting on the hole's own axis grazes every one of its
facets. It is used below for the section drawing and nothing else, and only a
quarter millimetre off a surface, where it agrees with the rays.

Run: python3 verify_aspirator.py [file.3mf] [--angle 45]
"""

import math
import sys
import zipfile
import xml.etree.ElementTree as ET

import csglib as csg
import checks as K
import memcap
import printability as PR
import threelib as t
import make_aspirator as gen        # for the default filename, and nothing else

# --- The specification, restated here. Radii in mm; z is the flow axis, which
# is vertical in the print, and the gas inlet points along +x. ---------------
INLET_BORE, INLET_OUT = 2.6, 4.6    # the water supply tube
INLET_TOP = 28.0                    # where its wall stops
THROAT, NOZ_WALL = 1.3, 1.2         # the jet and the tube it runs through
JET_RELEASE = 0.0                   # the channel gains nothing between the
                                    # ...throat and the tip: the water only ever
                                    # ...narrows, and the nozzle is a converging
                                    # ...tube and nothing else
NOZ0, NOZ_TOP = 24.0, 49.0          # the nozzle's foot, and its tip
FLOOR_BORE = 3.2                    # the chamber floor's hole, which grips the
                                    # ...nozzle rather than passing it
FOOT_BOTTOM = 9.8                   # where the shoulder starts leaning
FLOOR_TOP = 30.0                    # and where it stops, being the chamber's floor
CHAMBER_R, CHAMBER_OUT = 14.0, 16.0
CHAMBER_LO, CHAMBER_HI = 26.0, 56.0
NK0, NK2 = 54.0, 66.0               # the neck: starts in the shell, ends in the
                                    # ...outlet's wall
OUTLET_BORE, OUTLET_OUT = 2.6, 3.6
OUTLET_LO, OUTLET_HI = 60.0, 130.0  # and it starts 32 mm above the inlet's end
AIR_BORE, AIR_BOSS, AIR_OUT = 2.6, 3.6, 20.0
AIR_Z = 46.0                        # the gas inlet's station on the flow axis
NARROW_LO, NARROW_HI = 38.0, 49.0   # where the bore holds at the throat -- the
                                    # ...hold ends at the tip, not before it
PROBE = 0.5                         # how far the rays stand off the axis
OFF = 0.17                          # ...and a little off the whole millimetre,
                                    # ...so no probe lands on a station's plane
IN_SPEC = dict(inlet_bore=INLET_BORE, inlet_out=INLET_OUT, throat=THROAT,
               nozzle_wall=NOZ_WALL, chamber_r=CHAMBER_R,
               chamber_out=CHAMBER_OUT, outlet_bore=OUTLET_BORE,
               outlet_out=OUTLET_OUT, air_bore=AIR_BORE, air_boss=AIR_BOSS,
               air_out=AIR_OUT, z_air=AIR_Z, z_cone0=INLET_TOP,
               z_cone1=NARROW_LO, z_throat1=NARROW_HI, z_noz1=NOZ_TOP,
               z_ch0=CHAMBER_LO, z_ch1=CHAMBER_HI, z_out1=OUTLET_HI,
               jet_release=JET_RELEASE)

def profile(pts, z):
    """Straight lines through (station, value) pairs, held outside the ends."""
    if z <= pts[0][0]:
        return pts[0][1]
    for (za, ra), (zb, rb) in zip(pts, pts[1:]):
        if z <= zb:
            return ra + (z - za) / (zb - za) * (rb - ra)
    return pts[-1][1]





# What the part is made of, as washers: the same pieces the build script
# describes, written out again here so a slip in one is not repeated in the
# other. Each is (bottom, top, bore, skin), and they overlap on purpose -- the
# union keeps the *union* of these, which is what a band merge gives.
FOOT_SKIN = ((FOOT_BOTTOM, INLET_OUT - 0.3), (CHAMBER_LO, CHAMBER_OUT - 0.3),
             (30.0, CHAMBER_OUT - 0.3))
PIECES = [(0.0, INLET_TOP, INLET_BORE, INLET_OUT),           # inlet tube
          (FOOT_BOTTOM, 30.0, FLOOR_BORE, FOOT_SKIN),        # shoulder + floor
          (CHAMBER_LO, CHAMBER_HI, CHAMBER_R, CHAMBER_OUT),  # the shell
          (OUTLET_LO, OUTLET_HI, OUTLET_BORE, OUTLET_OUT)]   # outlet tube
# The nozzle holds at the throat right out to its tip: no final station, because
# JET_RELEASE is 0 and the channel is never asked to widen.
NOZZLE = ((NOZ0, INLET_BORE + 0.3), (INLET_TOP, INLET_BORE), (NARROW_LO, THROAT))
NECK = ((NK0, CHAMBER_R + 0.3, CHAMBER_OUT - 0.3), (CHAMBER_HI, CHAMBER_R + 0.3,
        CHAMBER_OUT - 0.3), (OUTLET_LO, OUTLET_BORE, CHAMBER_OUT - 4.0),
        (NK2, OUTLET_BORE + 0.3, OUTLET_OUT))               # closes the chamber on


def design_bands(z):
    """The radial bands of material the design leaves standing at station z."""
    spans = []
    for z0, z1, ri, ro in PIECES:
        if z0 <= z <= z1:
            # The shoulder's skin is a profile, not a constant.
            spans.append((ri, profile(ro, z) if isinstance(ro, tuple) else ro))
    if NOZ0 <= z <= NOZ_TOP:
        b = profile(NOZZLE, z)
        spans.append((b, b + NOZ_WALL))
    if NK0 <= z <= NK2:
        spans.append((profile([(a, r) for a, r, _ in NECK], z),
                      profile([(a, r) for a, _, r in NECK], z)))
    spans.sort()
    out = []
    for lo, hi in spans:
        if out and lo <= out[-1][1] + 1e-9:
            out[-1][1] = max(out[-1][1], hi)
        else:
            out.append([lo, hi])
    return [(a, b) for a, b in out]


def closed_form_volume(step=0.05):
    """Washers over `design_bands`, plus the boss, less the drill.

    Coarse where pieces bury one another -- this exists to catch a bore that was
    never drilled or a wall that was never built, which is thousands of cubic
    millimetres, not a rounding slip.
    """
    total = 0.0
    z = 0.0
    while z < OUTLET_HI:
        for lo, hi in design_bands(z):
            total += math.pi * (hi * hi - lo * lo) * step
        z += step
    base = CHAMBER_R + 0.4
    total += math.pi * AIR_BOSS ** 2 * (AIR_OUT - base)          # the boss
    y = -AIR_BOSS
    while y < AIR_BOSS:                       # less its foot, inside the wall
        h = math.sqrt(max(0.0, AIR_BOSS ** 2 - y * y))
        lo = max(base, math.sqrt(max(0.0, CHAMBER_R ** 2 - y * y)))
        hi = min(AIR_OUT, math.sqrt(max(0.0, CHAMBER_OUT ** 2 - y * y)))
        if hi > lo:
            total -= 2 * h * (hi - lo) * 0.05
        y += 0.05
    y = -AIR_BORE
    while y < AIR_BORE:                       # less the hole itself
        h = math.sqrt(max(0.0, AIR_BORE ** 2 - y * y))
        lo = min(base, math.sqrt(max(0.0, CHAMBER_R ** 2 - y * y)))
        total -= 2 * h * (AIR_OUT - lo) * 0.05
        y += 0.05
    return total


def expected_inner(z):
    """The radius the water's passage should have at station z: the inner edge
    of the innermost band the model leaves standing, so the chamber reads as the
    28 mm void it is rather than as the jet that runs through it."""
    bands = design_bands(z)
    return bands[0][0] if bands else None


# --- Rays ------------------------------------------------------------------
def prep(verts, tris):
    return [(verts[a],
             (verts[b][0] - verts[a][0], verts[b][1] - verts[a][1],
              verts[b][2] - verts[a][2]),
             (verts[c][0] - verts[a][0], verts[c][1] - verts[a][1],
              verts[c][2] - verts[a][2]))
            for a, b, c in tris]


def cast(data, o, d):
    """Every distance along o + t*d, t > 0, where the ray crosses a triangle."""
    ox, oy, oz = o
    dx, dy, dz = d
    hits = []
    for a, e1, e2 in data:
        px = dy * e2[2] - dz * e2[1]
        py = dz * e2[0] - dx * e2[2]
        pz = dx * e2[1] - dy * e2[0]
        det = e1[0] * px + e1[1] * py + e1[2] * pz
        if -1e-12 < det < 1e-12:
            continue
        inv = 1.0 / det
        tx, ty, tz = ox - a[0], oy - a[1], oz - a[2]
        u = inv * (tx * px + ty * py + tz * pz)
        if u < -1e-9 or u > 1 + 1e-9:
            continue
        qx = ty * e1[2] - tz * e1[1]
        qy = tz * e1[0] - tx * e1[2]
        qz = tx * e1[1] - ty * e1[0]
        v = inv * (dx * qx + dy * qy + dz * qz)
        if v < -1e-9 or u + v > 1 + 1e-9:
            continue
        tt = inv * (e2[0] * qx + e2[1] * qy + e2[2] * qz)
        if tt > 1e-9:
            hits.append(tt)
    hits.sort()
    return hits


def material(hits, start):
    """Material intervals along a ray whose start end is outside the part."""
    return [(start + hits[i], start + hits[i + 1])
            for i in range(0, len(hits) - 1, 2)]


def radial_bands(data, z, probe=PROBE):
    """The material's radial bands at height z, read outward from the axis.

    None if the number of crossings is odd -- the ray hit a surface dead-on --
    and [] if the station turns out to be hollow all the way to the sky.
    """
    hits = cast(data, (probe, -80.0, z), (0.0, 1.0, 0.0))
    if len(hits) % 2:
        return None
    rs = [math.hypot(probe, -80.0 + h) for h in hits if -80.0 + h > 0.0]
    return [(rs[i], rs[i + 1]) for i in range(0, len(rs) - 1, 2)]


def ledge(data, c):
    """How wide a downward-facing ledge is, measured from under it.

    Stand just below the face, in the void it bounds, and feel sideways towards
    and away from the axis: the nearer wall is the one holding the ledge up, and
    twice that distance is the ledge's width. A ray along the wall is safe here
    because it starts off the face's own plane.
    """
    x, y, z = c
    r = math.hypot(x, y)
    if r < 1e-6:
        return 0.0
    ux, uy = x / r, y / r
    o = (x + 0.0 * ux, y, z - 0.08)
    out = []
    for d in ((ux, uy, 0.0), (-ux, -uy, 0.0)):
        hits = cast(data, o, d)
        out.append(hits[0] if hits else 99.0)
    return 2.0 * min(out)


def edge_census(V, T):
    """Edges whose triangle multiplicity is not two, grouped by multiplicity.

    A boolean's seam can leave an edge with four triangles around it -- two
    surfaces meeting along one line, watertight, a crease -- or with one, which
    is a hole. `real_gaps` cannot tell those apart on its own, because a crease
    has neighbours lying along it like a T-junction does, so count them here.
    """
    count = {}
    for tri in T:
        for a, b in ((tri[0], tri[1]), (tri[1], tri[2]), (tri[2], tri[0])):
            k = (a, b) if a < b else (b, a)
            count[k] = count.get(k, 0) + 1
    by = {}
    for (a, b), n in count.items():
        if n == 2:
            continue
        d = by.setdefault(n, [0, 0.0, 0.0])
        length = math.dist(V[a], V[b])
        d[0] += 1
        d[1] += length
        d[2] = max(d[2], length)
    return by


def read_params(path):
    """The Parameter.* values the file carries about itself."""
    NS = {"c": "http://schemas.microsoft.com/3dmanufacturing/core/2015/02"}
    with zipfile.ZipFile(path) as z:
        root = ET.fromstring(z.read("3D/3dmodel.model"))
    out = {}
    for m in root.findall("c:metadata", NS):
        name = m.get("name") or ""
        if name.startswith("Parameter."):
            try:
                out[name.split(".", 1)[1]] = float(m.text)
            except (TypeError, ValueError):
                pass
    return out


def bbox(V):
    xs = [v[0] for v in V]
    ys = [v[1] for v in V]
    zs = [v[2] for v in V]
    return (min(xs), min(ys), min(zs), max(xs), max(ys), max(zs))


def tri_area(V, tri):
    ax, ay, az = V[tri[0]]
    ux, uy, uz = V[tri[1]][0] - ax, V[tri[1]][1] - ay, V[tri[1]][2] - az
    vx, vy, vz = V[tri[2]][0] - ax, V[tri[2]][1] - ay, V[tri[2]][2] - az
    nx, ny, nz = uy * vz - uz * vy, uz * vx - ux * vz, ux * vy - uy * vx
    return 0.5 * math.sqrt(nx * nx + ny * ny + nz * nz)


# --- The report ------------------------------------------------------------
def main():
    memcap.cap()
    path, angle = K.cli(sys.argv[1:], gen.default_out())
    R = K.Report("water aspirator")
    print(f"  file              : {path}")
    objs, _, _unit = K.load_3mf(path, R)
    R.check("one object in the file", len(objs) == 1,
            f"{len(objs)}: {[o['name'] for o in objs]}")
    if len(objs) != 1:
        return R.finish()
    V, T = objs[0]["verts"], objs[0]["tris"]
    data = prep(V, T)

    # The chord bound: a facetted circle is inscribed in the circle it stands
    # for, so a measured radius can only read small, by the sagitta at most.
    meta = read_params(path)
    tol = meta.get("tol", 0.05)
    nseg = t.auto_segments(CHAMBER_OUT, tol)
    SAG = CHAMBER_OUT * (1.0 - math.cos(math.pi / nseg))
    # A hole has its own sagitta, off its own radius and its own facet count, and
    # the wall it is drilled through has that one too: measured across a drilled
    # circle, both bite, so neither bound alone covers the reading.
    n_hole = t.auto_segments(AIR_BORE, tol)
    HOLE_SAG = 2 * AIR_BORE * (1.0 - math.cos(math.pi / n_hole))

    def within(got, want):
        return want - SAG - 1e-9 <= got <= want + 1e-6

    wrong = {k: (meta.get(k), v) for k, v in IN_SPEC.items()
             if meta.get(k) is None or abs(meta[k] - v) > 1e-9}
    R.check("the file's own metadata is the specification", not wrong,
            f"disagrees: {wrong}" if wrong else
            f"{len(IN_SPEC)} params, {nseg} segments, chord under {SAG:.3f} mm")

    R.check("triangle indices in range",
            all(0 <= i < len(V) for tri in T for i in tri))
    R.check("every coordinate finite",
            all(all(math.isfinite(c) for c in v) for v in V))
    R.check("no triangle repeats a vertex", all(len(set(tri)) == 3 for tri in T))
    shells = K.shell_triangles(V, T)
    gaps = csg.real_gaps(V, T)
    sliv = [tri for tri in T if tri_area(V, tri) < 1e-9]
    census = edge_census(V, T)
    R.check("one closed shell", len(shells) == 1, f"{len(shells)} shells")
    R.check("the surface has no real gaps", not gaps,
            f"{len(sliv)} zero-area slivers left on the boolean seams, "
            f"{sum(d[0] for d in census.values())} edges not used twice "
            f"(see the next line)"
            if not gaps else f"{len(gaps)} open edges, e.g. {gaps[:3]}")
    R.check("nothing stands open: no edge is used once", 1 not in census,
            "the edges that are not used twice are "
            + ", ".join(f"{d[0]} with {n} triangles ({d[1]:.1f} mm, longest "
                        f"{d[2]:.1f} mm)" for n, d in sorted(census.items()))
            + " -- creases where two curved surfaces met along one line, "
              f"crack_length {csg.crack_length(V, T):.1f} mm")

    bb = bbox(V)
    want = (-CHAMBER_OUT, -CHAMBER_OUT, 0.0, AIR_OUT, CHAMBER_OUT, OUTLET_HI)
    R.check("bounding box", all(abs(a - b) < SAG + 0.02 for a, b in zip(bb, want)),
            f"{bb[3]-bb[0]:.1f} x {bb[4]-bb[1]:.1f} x {bb[5]-bb[2]:.1f} mm; the "
            f"gas inlet reaches {bb[3]:.1f} mm out from the axis")
    vol = K.signed_volume(V, T)
    cf = closed_form_volume()
    R.near("volume against the closed form", vol, cf, 0.025 * cf, fmt="{:.0f}")
    R.check("it stands along its own axis, tall way up",
            (bb[5] - bb[2]) > max(bb[3] - bb[0], bb[4] - bb[1]) + 50.0,
            f"{bb[5]-bb[2]:.0f} mm tall, {max(bb[3]-bb[0], bb[4]-bb[1]):.0f} mm "
            f"across")

    # -- the channel, every half millimetre of it
    bands_at, shut, off = {}, [], []
    z = 0.5
    while z <= OUTLET_HI - 0.5:
        bands = radial_bands(data, z + OFF)
        if not bands:
            shut.append(round(z, 1))
        else:
            bands_at[round(z, 1)] = bands
            # Where the wall is conical the ray reads it at the height it
            # actually crossed, which is OFF above the station's own name.
            want = expected_inner(z + OFF)
            if not within(bands[0][0], want):
                off.append((round(z, 1), bands[0][0], want))
        z += 0.5
    R.check("the water passage is open the whole way through", not shut,
            f"{len(bands_at)} of {len(bands_at) + len(shut)} stations open"
            + (f", shut at {shut[:6]}" if shut else ""))
    R.check("the bore follows the promised profile", not off,
            f"{len(off)} stations out by more than the chord: "
            + ", ".join(f"z={a}: {b:.2f} vs {c:.2f}" for a, b, c in off[:5])
            if off else f"{len(bands_at)} stations, each within {SAG:.3f} mm of "
            f"{2 * INLET_BORE:g} mm narrowing to {2 * THROAT:g} mm")

    stations = sorted(bands_at)
    bores = {z: bands_at[z][0][0] for z in stations}

    # -- demand 1: the narrow part must bracket the gas inlet
    narrow = [z for z in stations if bores[z] <= THROAT + 0.02]
    inlet_lo, inlet_hi = AIR_Z - AIR_BORE, AIR_Z + AIR_BORE
    R.near("the throat is as narrow as asked", min(bores.values()), THROAT, SAG)
    # A ray crossing the jet reads the channel, not the tip's end face, so the
    # sampled narrow span stops about half a millimetre short of it. Where
    # JET_RELEASE holds the throat to the tip, the tip is the throat's end.
    narrow_end = max(narrow) if JET_RELEASE > 0.0 else NOZ_TOP
    R.check("the narrow part starts before the gas inlet and ends after it",
            min(narrow) < inlet_lo and narrow_end > inlet_hi,
            f"narrow over z={min(narrow):.1f}..{narrow_end:g}, and the "
            f"{2 * AIR_BORE:g} mm inlet spans z={inlet_lo:g}..{inlet_hi:g}, so "
            f"the throat clears it by {inlet_lo - min(narrow):.1f} mm below and "
            f"{narrow_end - inlet_hi:.2f} mm above"
            + (f" (sampled to {max(narrow):.1f})" if narrow_end != max(narrow)
               else ""))
    R.check("it holds the throat longer than the inlet is wide",
            max(narrow) - min(narrow) >= 2 * AIR_BORE,
            f"{max(narrow) - min(narrow):.1f} mm of throat round a "
            f"{2 * AIR_BORE:g} mm inlet")
    cone = [(z, bores[z]) for z in stations if INLET_TOP <= z <= NARROW_LO]
    R.check("the water is squeezed steadily, with no step in the cone",
            all(b[1] <= a[1] + 1e-9 for a, b in zip(cone, cone[1:]))
            and within(cone[0][1], expected_inner(cone[0][0] + OFF))
            and within(cone[-1][1], expected_inner(cone[-1][0] + OFF)),
            f"{2 * cone[0][1]:.1f} mm at z={cone[0][0]:g} down to "
            f"{2 * cone[-1][1]:.1f} mm at z={cone[-1][0]:g}")
    # JET_RELEASE is the one knob that decides whether the water widens again
    # before it leaves the nozzle. At 0.0 the channel must not grow by anything:
    # the jet is still at its throat speed when it comes out, and it slows in the
    # open chamber, past the inlet's mouth, which is what carries the air away.
    tail = [(z, bores[z]) for z in stations
            if NARROW_LO <= z <= NOZ_TOP and bores[z] < CHAMBER_R / 2.0]
    R.check("and the water only ever narrows -- it never widens again",
            tail and all(b[1] <= a[1] + 1e-9 for a, b in zip(tail, tail[1:]))
            and max(b for _, b in tail) <= THROAT + SAG + 1e-9,
            f"from z={NARROW_LO:g} to the tip at z={NOZ_TOP:g} the channel "
            f"stays {2 * THROAT:g} mm across, the widest reading in it being "
            f"{2 * max(b for _, b in tail):.2f} mm")
    R.check("so the jet is still a jet past the mouth of the air inlet",
            NOZ_TOP > AIR_Z + AIR_BORE and NOZ_TOP < CHAMBER_HI,
            f"the tip is at z={NOZ_TOP:g}, {NOZ_TOP - (AIR_Z + AIR_BORE):.2f} mm "
            f"above the inlet's far edge and {CHAMBER_HI - NOZ_TOP:g} mm below "
            f"the chamber's top")

    # -- demand 3: two rings at the gas inlet
    for zz in (inlet_lo + 1.0, AIR_Z, inlet_hi - 1.0):
        bands = radial_bands(data, zz + OFF)
        ok = (bands and len(bands) == 2
              and within(bands[0][0], THROAT)
              and within(bands[0][1], THROAT + NOZ_WALL)
              and within(bands[1][0], CHAMBER_R)
              and within(bands[1][1], CHAMBER_OUT))
        R.check(f"two rings across the section at z={zz:g}", ok,
                ("bands " + " ".join(f"[{a:.2f},{b:.2f}]" for a, b in bands)
                 if bands else "this section is closed"))
    at_inlet = radial_bands(data, AIR_Z + OFF)
    R.near("the nozzle's wall is the thickness asked for",
           at_inlet[0][1] - at_inlet[0][0], NOZ_WALL, 2 * SAG)
    annulus = at_inlet[1][0] - at_inlet[0][1]
    R.check("the air has an annulus wider than the hole that feeds it",
            annulus > 2 * AIR_BORE, f"{annulus:.2f} mm between the two rings")
    R.check("the chamber is a chamber, not a slot",
            at_inlet[1][0] > THROAT + NOZ_WALL + 2 * AIR_BORE,
            f"{2 * at_inlet[1][0]:.1f} mm of void round a "
            f"{2 * (THROAT + NOZ_WALL):.1f} mm nozzle")

    # -- demand 2: read lengthwise, the wall is a Z. The outside has to change
    # its radius at a different station from the inside, and the two diagonals
    # have to be joined by wall that is solid all the way across.
    outer = {z: bands_at[z][-1][1] for z in stations}
    # The shoulder's own flank: skin radius between the inlet tube and the
    # chamber's shell, below the shell's bottom edge.
    flank = [(z, outer[z]) for z in stations
             if INLET_OUT + 0.4 < outer[z] < CHAMBER_OUT - 0.4 and z < CHAMBER_LO]
    lean = 90.0
    if len(flank) > 2:
        (z0, r0), (z1, r1) = flank[0], flank[-1]
        lean = math.degrees(math.atan2(z1 - z0, max(1e-9, r1 - r0)))
    step_out = flank[0][0] if flank else None
    step_in = min((z for z in stations if bores[z] < INLET_BORE - 0.4),
                  default=None)
    R.check("the outside leans out from the inlet tube to the chamber",
            len(flank) > 2,
            f"the skin runs {outer[flank[0][0]]:.1f} to "
            f"{outer[flank[-1][0]]:.1f} mm over z={flank[0][0]:g}..{flank[-1][0]:g}, "
            f"{lean:.0f} degrees off the horizontal" if len(flank) > 2 else
            "the skin never leaves the inlet tube")
    R.check("the inside draws in from the bore towards the throat",
            step_in is not None,
            f"the bore is a full {2 * INLET_BORE:g} mm down to z={step_in - 0.5:g} "
            f"and narrowing by z={step_in:g}" if step_in else "the bore never narrows")
    R.check("and the two diagonals sit at different stations, which is the Z",
            step_out is not None and step_in is not None and step_out < step_in,
            f"outside out at z={flank[0][0]:g}..{flank[-1][0]:g}, inside in at "
            f"z={step_in:g}..{NARROW_LO:g}" if step_in else "no inner step")
    R.check("the shoulder is steep enough that the printer carries it alone",
            lean > 45.0, f"{lean:.0f} degrees off the horizontal over "
            f"{len(flank)} stations")
    # The join is where the shoulder's skin and the nozzle's cone are the same
    # piece of plastic. Above the floor's top (Z_FOOT1) they legitimately part
    # into two rings, and that gap is the air's passage, so it is not a fault.
    # Strictly below the floor's top: exactly on a piece's end plane a ray
    # crosses the face, not the wall, and the band that is there one layer down
    # simply is not readable.
    join = [z for z in stations
            if (step_out or 0) <= z <= (step_in or 0) and z < FLOOR_TOP]
    apart = [z for z in join if len(bands_at[z]) != 1]
    one = [z for z in join if len(bands_at[z]) == 1]
    mid = bands_at[one[len(one) // 2]][0] if one else (0.0, 0.0)
    R.check("across the join it is one continuous band, not two walls apart",
            join and not apart and len(join) >= 4,
            f"{len(join)} stations from z={join[0]:g} to {join[-1]:g}: material "
            f"from {mid[0]:.1f} mm to {mid[1]:.1f} mm, with nothing hollow in "
            f"between" if not apart else
            f"{len(apart)} of {len(join)} read two walls, first at z="
            f"{apart[0]:g}: " + " ".join(f"[{a:.2f},{b:.2f}]"
                                         for a, b in bands_at[apart[0]]))
    upper = [z for z in stations if NOZ_TOP + 1.0 < z < CHAMBER_HI - 0.5]
    R.check("past the nozzle's tip the shell is still one plain ring",
            upper and all(len(bands_at[z]) == 1
                          and within(bands_at[z][0][0], CHAMBER_R)
                          and within(bands_at[z][0][1], CHAMBER_OUT)
                          for z in upper),
            f"{len(upper)} stations above z={NOZ_TOP:g} read "
            f"{2 * bands_at[upper[-1]][0][0]:.1f} mm void, "
            f"{2 * bands_at[upper[-1]][0][1]:.1f} mm outside")

    # -- the gas inlet: drilled, open, and short of the jet
    for yy, zz in ((1.0, 46.6), (-1.2, 47.4), (0.5, 46.0), (2.2, 44.9)):
        hits = cast(data, (8.0, yy, zz), (1.0, 0.0, 0.0))
        R.check(f"the annulus opens straight out of the side at y={yy:g}, "
                f"z={zz:g}", not hits,
                f"{len(hits)} surfaces in the way at x=" +
                " ".join(f"{8.0 + h:.1f}" for h in hits[:3])
                if hits else "nothing between the jet's skin and the air")
    col = cast(data, (CHAMBER_OUT - 1.0, PROBE, -80.0), (0.0, 0.0, 1.0))
    mats = [(a, b) for a, b in material(col, -80.0)
            if b > CHAMBER_LO - 3 and a < CHAMBER_HI + 3]
    R.check("the hole is drilled through the chamber's own wall",
            len(mats) == 2,
            " ".join(f"[{a:.2f},{b:.2f}]" for a, b in mats)
            + f" -- {len(mats)} lands of material up that column")
    if len(mats) == 2:
        gap = mats[1][0] - mats[0][1]
        R.near("round and the diameter asked for", gap, 2 * AIR_BORE,
               HOLE_SAG + 2 * SAG)
        R.near("centred on the throat's station",
               (mats[1][0] + mats[0][1]) / 2.0, AIR_Z, 2 * SAG)
        R.check("with wall left below it and above it",
                mats[0][1] - CHAMBER_LO > 4 and CHAMBER_HI - mats[1][0] > 4,
                f"{mats[0][1] - CHAMBER_LO:.1f} mm of wall under the hole, "
                f"{CHAMBER_HI - mats[1][0]:.1f} mm over it, and the shell's own "
                f"ends are {CHAMBER_LO:g} and {CHAMBER_HI:g}")
    boss = material(cast(data, (0.0, AIR_BOSS - 0.5, AIR_Z + 0.13),
                         (1.0, 0.0, 0.0)), 0.0)
    R.check("the boss stands proud of the chamber and is fused into it",
            len(boss) == 1 and abs(boss[0][1] - AIR_OUT) < SAG + 0.02
            and boss[0][0] < CHAMBER_OUT + 1.5,
            f"one land of material from x={boss[0][0]:.2f} to "
            f"{boss[0][1]:.2f} mm" if boss else "nothing out that way")
    skin_above = radial_bands(data, AIR_Z + AIR_BORE - 0.6 + OFF)
    # And on the face the hole comes from: run a ray along the axis line at the
    # inlet's height. It crosses the jet's wall twice, either side of the axis,
    # and the hole's floor is beyond that -- so if both crossings land where the
    # nozzle's own surfaces are, the drill never touched the tube it runs in.
    # This ray runs along +x at y = PROBE, so it reads positions on a chord, not
    # radii: convert back with the offset it was cast at, or a jet that is
    # exactly THROAT reads a tenth of a millimetre too small and looks wrong.
    facing = sorted((math.hypot(min(abs(a), abs(b)), PROBE),
                     math.hypot(max(abs(a), abs(b)), PROBE)) for a, b in
                    material(cast(data, (-80.0, PROBE, AIR_Z + OFF),
                                   (1.0, 0.0, 0.0)), -80.0)
                    if max(abs(a), abs(b)) < CHAMBER_R - 1.0)
    R.check("the drill left the jet's own tube unbroken",
            all(within(b[0][1], THROAT + NOZ_WALL)
                for b in (at_inlet, skin_above)) and
            within(skin_above[0][0], THROAT) and len(facing) == 2 and
            all(within(x[1], THROAT + NOZ_WALL) and within(x[0], THROAT)
                for x in facing),
            f"the nozzle still reads a {2 * at_inlet[0][0]:.1f} mm jet inside a "
            f"{2 * at_inlet[0][1]:.1f} mm tube above the hole, and its wall "
            f"crosses this ray at {facing[0][0]:.2f} and {facing[0][1]:.2f} mm "
            f"on either side of the axis, where it should"
            if len(facing) == 2 else f"{len(facing)} walls on the air's line")

    # -- printability, standing on the inlet end as the file lies
    faces = [f for f in PR.overhangs(V, T, max_angle=angle, exposed_only=True,
                                     shells=shells) if f["area"] > 0.05]
    groups = {}
    for f in faces:
        x, y, z = f["c"]
        from_axis = math.hypot(y, z - AIR_Z)
        if x > CHAMBER_R - 2 and from_axis <= AIR_BORE + 0.4 and z > AIR_Z:
            k = "air roof"          # the hole's ceiling, inside the chamber
        elif x > CHAMBER_R - 2 and from_axis <= AIR_BOSS + 0.4:
            k = "boss"              # the inlet's own pad, outside the part
        elif math.hypot(x, y) < CHAMBER_OUT:
            # Anything else shallow and facing down is inside the part's own
            # silhouette. How wide is the ledge? Stand just under its ceiling
            # and feel sideways for the wall that holds it up: a lip too thin
            # to matter, or a shelf the printer has to carry.
            k = "lip" if ledge(data, (x, y, z)) < 1.2 else "ceiling"
        else:
            k = "other"
        g = groups.setdefault(k, {"n": 0, "area": 0.0, "shallow": 999.0,
                                  "where": None, "widest": 0.0})
        g["n"] += 1
        g["area"] += f["area"]
        if k == "lip":
            g["widest"] = max(g["widest"], ledge(data, (x, y, z)))
        if f["angle"] < g["shallow"]:
            g["shallow"] = f["angle"]
            g["where"] = (round(x, 1), round(y, 1), round(z, 1))
    label = {"air roof": "the air inlet's roof, inside the chamber",
             "boss": "the boss's underside, outside the part",
             "lip": "joint lips inside the silhouette, too thin to matter",
             "ceiling": "INTERNAL SHELF, would need trapped support",
             "other": "unclassified"}
    for k in sorted(groups, key=lambda g: -groups[g]["area"]):
        g = groups[k]
        print(f"  {g['area']:7.1f} mm^2 over {g['n']:4d} faces, shallowest "
              f"{min(g['shallow'], 180.0):3.0f} deg at {g['where']} -- {label[k]}")
    roof = groups.get("air roof", {"area": 0.0, "n": 0})
    boss_g = groups.get("boss", {"area": 0.0, "n": 0})
    lip = groups.get("lip", {"area": 0.0, "n": 0, "widest": 0.0})
    R.check("nothing inside needs a support that could not be got at",
            not groups.get("ceiling") and not groups.get("other"),
            f"the only ceiling inside is the air inlet's own roof, "
            f"{roof['area']:.0f} mm^2 over {roof['n']} faces, and it bridges a "
            f"{2 * AIR_BORE:g} mm through-hole; the rest of what faces down "
            f"inside the silhouette is {lip['n']} lips {lip['widest']:.2f} mm "
            f"across at the joints, thinner than the nozzle that lays them")
    R.check("the outside support is only under the boss",
            boss_g["area"] < 250.0,
            f"{boss_g['area']:.0f} mm^2 under the inlet's pad, removable")
    bc = PR.bed_contact(V, T)
    R.check("the first layer sits on the inlet's end face", bc > 20.0,
            f"{bc:.0f} mm^2 of {(INLET_OUT - INLET_BORE):.1f} mm-wide ring, "
            f"{2 * INLET_OUT:.1f} mm across -- thin for a part this tall, so "
            f"splay it with a brim")

    findings, summ = PR.assess(V, T, overhang_deg=angle)
    for f in findings:
        print(f"  [{f['level']:5s}] {f['check']:9s} {f['msg']}")
    for k in sorted(summ):
        print(f"  {k:21s}: {summ[k]}")

    section(data)
    return R.finish()


def section(data):
    """The axial section, from the same rays as the measurements: 1 mm per
    column, 3 mm per row, read at y = PROBE so nothing runs down a vertex
    column. The gas inlet's hole shows as the gap in the wall on the right."""
    x0, x1 = -18.0, 28.0
    cols = int(x1 - x0)
    print(f"\n  section through the axis, read at y={PROBE} mm, z upward")
    print("   z  " + "".join("+" if x % 10 == 0 else
                             "|" if x % 5 == 0 else " "
                             for x in range(int(x0), int(x1))))
    z = OUTLET_HI - 1.0
    while z > 0:
        row = ["."] * cols
        for a, b in material(cast(data, (x0, PROBE, z), (1.0, 0.0, 0.0)), x0):
            for i in range(max(0, int(a - x0)), min(cols, int(b - x0) + 1)):
                row[i] = "#"
        mark = ""
        if abs(z - AIR_Z) < 1.6:
            mark = "   <- the air inlet"
        print(f"{z:4.0f}  " + "".join(row) + mark)
        z -= 3.0


if __name__ == "__main__":
    if not main():
        sys.exit(1)

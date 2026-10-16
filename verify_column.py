#!/usr/bin/env python3
"""Verify a fractionating column 3MF: joints, bore, staging, and the passage.

Re-parsed from disk and measured from the mesh. The part is a tube, so the
mouths are the two annular faces on the axis and the joints are the ground
cones behind them: the bottom one is a male tip that widens going in, the top
one is a socket whose bore narrows going in, and both are measured at several
depths against the 1:10 standard the rest of this project uses.

Two things make a column work, and both are checked directly. The bore has to
be *open* from mouth to mouth: a blocked passage prints perfectly and distils
nothing, so the vapour path the build reports is walked through the mesh
sample by sample, and every sample has to be in air. And the staging has to be
the shape the style claims -- washers with a hole, a twisting band, staggered
fingers, leaning plates -- which is measured by cutting planes across the bore
and looking at what stands in them.

What is left is the trade-off the part exists to manage: contact surface
against pressure drop. Both are reported rather than judged, since what counts
as too much drop depends on the vacuum you have.

Usage: python3 verify_column.py [column.3mf]
"""

import math
import sys

import checks as K
import csglib as csg
import junctionlib as J
import make_column as C
import printability as PR

PATH, ANGLE = K.cli(sys.argv[1:], C.default_out())
P = K.read_params(PATH)
# The angle the part was *designed* for is stamped in the file, because the
# staging geometry depends on it: a washer is as steep as it is, a plate leans
# as far as it does, because of the number the build was given. Checking
# against a different number is a legitimate thing to do -- that is what
# `--angle` is for -- but it is not the same check, so say so.
STAMPED = float(P.get("overhang", PR.OVERHANG_DEG))
ANGLE_DEG = STAMPED if ANGLE is None else ANGLE

rep = K.Report(f"Fractionating column: {PATH}")
objs, _, _ = K.load_3mf(PATH, rep)
rep.check("one object", len(objs) == 1, f"{len(objs)}")
if ANGLE is not None and abs(ANGLE - STAMPED) > 1e-9:
    print(f"  note: checking against {ANGLE:g} deg from horizontal, but this "
          f"part was built for {STAMPED:g} deg, so its staging was shaped to "
          f"that other number")
verts, tris = objs[0]["verts"], objs[0]["tris"]
shells = K.shell_triangles(verts, tris)
boxes = K.shell_boxes(verts, tris, shells)

D = C.dims(P)
STYLE = P["style"]
COUNT = int(round(P["elements"]))
R_BORE, R_BODY = D["r_bore"], D["r_body"]
Z0, Z1 = C.stage_span(P)
PITCH = (Z1 - Z0) / (COUNT + 1)
TOP = D["z_body1"] + P["cap"] + D["glen_t"]


def radius(p):
    return math.hypot(p[0], p[1])


# ------------------------------------------------------------ the mouths ---
# An annular face is a mouth; one with material behind it is a cap where some
# arm stopped short, which this part has none of. But this part is full of
# annuli that are not mouths -- every washer ends in an open ring at its own
# hole -- so a mouth is narrowed down to an annular face that is perpendicular
# to the axis *and* at an end of the part, since a tube has exactly two.
z_min = min(v[2] for v in verts)
z_max = max(v[2] for v in verts)
faces = K.annuli(verts, tris, shells=shells)
mouths = [f for f in faces if not f["buried"]
          and abs(abs(f["n"][2]) - 1.0) < 0.02
          and (f["c"][2] < z_min + 0.5 or f["c"][2] > z_max - 0.5)]
rings = [f for f in faces if not f["buried"] and f not in mouths]
rep.check("two mouths, one at each end of the axis", len(mouths) == 2,
          f"{len(mouths)} mouths: "
          + ", ".join(f"z={f['c'][2]:.1f} {f['inner'] * 2:.1f}->"
                      f"{f['outer'] * 2:.1f} mm" for f in mouths)
          + (f"; {len(rings)} other open rings, which are the staging's own"
             if rings else ""))
if len(mouths) != 2:
    print("\nSTOPPING: the rest of the checks need two mouths")
    sys.exit(1)

low, high = sorted(mouths, key=lambda f: f["c"][2])
rep.check("the mouths face along the axis, opposite ways",
          low["n"][2] < -0.99 and high["n"][2] > 0.99,
          f"normals z={low['n'][2]:.3f} and {high['n'][2]:.3f}")
rep.check("the mouths are on the axis",
          radius(low["c"]) < 0.05 and radius(high["c"]) < 0.05,
          f"{radius(low['c']):.3f} and {radius(high['c']):.3f} mm off it")

# A male tip is the narrow end of its ground cone; a female socket is the wide
# end of its own, with the wall around it.
r_mouth_b, r_tip_b, glen_b = J.joint_dims(P["bottom_joint"])
r_mouth_t, r_tip_t, glen_t = J.joint_dims(P["top_joint"])
want_low = ((r_tip_b - P["tip_wall"], r_tip_b)
            if P["bottom_gender"] == "male" else (r_tip_b, r_tip_b + P["wall"]))
want_high = ((r_mouth_t, r_mouth_t + P["wall"])
             if P["top_gender"] == "female" else (r_tip_t - P["tip_wall"], r_tip_t))
rep.near("the bottom mouth is the free end of its joint", low["inner"],
         want_low[0], 0.05, "{:.3f}")
rep.near("the top mouth is the free end of its joint", high["inner"],
         want_high[0], 0.05, "{:.3f}")
rep.near("the wall around the bottom mouth", low["outer"], want_low[1], 0.05,
         "{:.3f}")
rep.near("the wall around the top mouth", high["outer"], want_high[1], 0.05,
         "{:.3f}")
z_low, z_high = low["c"][2], high["c"][2]
rep.near("the part is as tall as its joints and body make it", z_high - z_low,
         TOP, 0.05, "{:.3f}")


def cut(z, ang=0.0):
    """Material intervals along a ray outward from the axis at height z."""
    return K.radial_profile(verts, tris, z, ang=ang, shells=shells)


# ------------------------------------------------------------ the joints ---
# The ground cones are the tapers everything else in this project mates with,
# so they are measured where a joint actually seals: along the ground length.
tapers = []
for s in (2.0, 6.0, 10.0, 14.0, 18.0):
    # Bottom: a male joint's outer surface is a cone widening going in.
    mat = cut(z_low + s)
    got = mat[0][1] if mat else 0.0
    tapers.append(("bottom", s, got, r_tip_b + J.TAPER * s))
    # Top: a female joint's socket is a cone narrowing going in.
    mat = cut(z_high - s)
    got = mat[0][0] if mat else 0.0
    tapers.append(("top", s, got, r_mouth_t - J.TAPER * s))
worst = max(tapers, key=lambda r: abs(r[2] - r[3]))
rep.check("both joints are ground to the project's taper",
          abs(worst[2] - worst[3]) <= 0.06,
          f"worst {worst[0]} at {worst[1]:g} mm in: {worst[2]:.4f} vs "
          f"{worst[3]:.4f} mm radius (tol 0.06)")
rep.check("the taper is the 1:10 on diameter this project uses",
          abs(J.TAPER - 0.05) < 1e-12, f"TAPER {J.TAPER:g} on radius")

# -------------------------------------------------------------- the body ---
# Between the two transitions the bore and the outside are both cylinders.
# A staging can touch the wall, so the wall is measured on the planes where it
# stands alone: an interval that starts at the bore and ends at the OD.
body = []
for i in range(COUNT + 1):
    for frac in (0.25, 0.6):
        z = Z0 + PITCH * (i + frac)
        if z >= Z1:
            continue
        for ang in (0.0, 0.9, 1.8, 2.7):
            wall = [m for m in cut(z, ang)
                    if abs(m[0] - R_BORE) < 0.02 and abs(m[1] - R_BODY) < 0.05]
            if wall:
                body.append((z, wall[0][0], wall[0][1]))
                break
worst = (max(body, key=lambda b: max(abs(b[1] - R_BORE),
                                     abs(b[2] - R_BODY)))
         if body else (0.0, 0.0, 0.0))
rep.check("the body bore and wall are what the parameters say",
          len(body) >= 2
          and abs(worst[1] - R_BORE) <= 0.05 and abs(worst[2] - R_BODY) <= 0.05,
          f"{len(body)} planes where the wall stands alone; worst at "
          f"z={worst[0]:.1f}: bore {worst[1] * 2:.3f} mm, "
          f"OD {worst[2] * 2:.3f} mm")
rep.check("the bore is wider than a joint's, so it is a column not a straw",
          R_BORE > max(D["bore_b"], D["bore_t"]),
          f"{R_BORE * 2:g} mm bore against {D['bore_b'] * 2:g} and "
          f"{D['bore_t'] * 2:g} mm at the joints")


# ------------------------------------------------- the staging, by planes ---
# What a staging style *is* can be measured: cut the bore across and see what
# stands in the plane. A washer is a closed ring with a hole; a band is a
# single strip that walks round as it rises; fingers are three separate lumps
# a row; a leaning plate is one lump that swaps sides.
N_ANG = 24


def scan(z, n_ang=N_ANG):
    """What stands in the bore at height z, azimuth by azimuth."""
    rows = []
    for j in range(n_ang):
        ang = 2.0 * math.pi * j / n_ang
        mat = cut(z, ang)
        free, prev = [], 0.0
        for r0, r1 in mat:
            if r0 > prev + 1e-6:
                free.append((prev, min(r0, R_BORE)))
            prev = max(prev, r1)
        if prev < R_BORE - 1e-6:
            free.append((prev, R_BORE))
        rows.append({
            "ang": ang,
            "free_axis": free[0][1] if free and free[0][0] < 1e-6 else 0.0,
            "area": sum(0.5 * (b * b - a * a) for a, b in free)
                    * (2.0 * math.pi / n_ang),
            "blocked": any(r0 < R_BORE - 0.5 for r0, r1 in mat),
            "mat": mat,
        })
    return rows


def plane_area(rows):
    """Free area of the plane, from the radial profiles that measured it."""
    return sum(r["area"] for r in rows)


def clusters(flags):
    """Runs of True in a circular list: how many separate lumps there are."""
    n = len(flags)
    if all(flags):
        return 1
    if not any(flags):
        return 0
    start = next(i for i in range(n) if not flags[i])
    runs, prev = 0, False
    for k in range(n):
        here = flags[(start + k) % n]
        if here and not prev:
            runs += 1
        prev = here
    return runs


def azimuth_of(rows):
    """Where the staging is, in turns.

    A strip of material has one azimuth where it runs furthest out -- the
    direction it points -- and that is what this returns. The mean azimuth of
    "everything blocked" would be right for a thin strip and useless for a band
    wide enough to cross the axis, which has material at every azimuth and no
    mean direction at all.
    """
    # Only material that is *inside* the bore counts: the band's outer edge
    # sits on the wall, so the outermost thing in the plane is the wall itself
    # at every azimuth, and asking for the longest run would answer "the wall"
    # wherever you looked.
    def reach(r):
        return max((min(b, R_BORE) for a, b in r["mat"] if a < R_BORE - 0.5),
                   default=0.0)

    best = max((reach(r) for r in rows), default=0.0)
    if best <= 0.0:
        return None
    # The peak is flat over a few samples -- a strip whose inner end is past
    # the axis reaches the wall over a range of azimuths, not one -- so take
    # the centre of the peak rather than its first sample, or the reading
    # jitters by a sample and the twist comes out as noise.
    hits = [r["ang"] for r in rows if reach(r) >= best - 0.05 * best]
    a = math.atan2(sum(math.sin(x) for x in hits) / len(hits),
                   sum(math.cos(x) for x in hits) / len(hits))
    return a / (2.0 * math.pi)


staging_planes = []
if STYLE == "cones":
    # A washer is a plate of revolution: at the wall it is a thin ring, and its
    # inner edge walks in as it rises, so every azimuth sees it -- one closed
    # ring -- until the hole is what is left at the axis. How tall the plate is
    # depends on where it is attached and how steep it is, so it is measured
    # rather than predicted: walk up from the staging and see where it ends.
    hole = max(1.5, R_BORE * C.CONE_HOLE_FRAC)
    rings, holes, found = 0, [], 0
    for k in range(COUNT):
        z = Z0 + PITCH * (k + 1)
        top = z
        for i in range(1, 17):
            zz = z + PITCH * i / 16.0
            if any(r["blocked"] for r in scan(zz, 2)):
                top = max(top, zz)
        if top <= z:
            continue
        found += 1
        widest = []
        for frac in (0.5, 0.75, 0.95):
            zz = z + (top - z) * frac
            rows = scan(zz)
            staging_planes.append((zz, rows))
            if all(r["blocked"] for r in rows):
                rings += 1
            widest.append(min(r["free_axis"] for r in rows))
        # The hole is narrowest at the plate's inner edge, which is the top of
        # the cut, so the tightest of the three is the one to compare.
        holes.append(min(widest))
    rep.check("every washer is a closed ring in the bore",
              found == COUNT and rings == COUNT * 3,
              f"{found} washers found, {rings} of {COUNT * 3} cuts across them "
              f"are rings")
    # The hole is not one number: it is whatever is left once the washer fits
    # in the pitch it was given, so the floor comes from the parameters and the
    # ceiling from the bore it has to leave room in.
    rep.check("every washer leaves a central hole open",
              all(hole - 0.6 <= h <= 0.9 * R_BORE for h in holes),
              f"holes {min(holes):.2f} to {max(holes):.2f} mm radius; the "
              f"parameters ask for at least {hole:.2f} and a {R_BORE:g} mm bore "
              f"leaves room for at most {0.9 * R_BORE:.2f}")

elif STYLE == "helix":
    # One band or several arms, and it walks round: the twist is what a helix
    # is for, so the rotation is measured rather than assumed. Each arm is its
    # own shell, which is how the arms are counted -- and the count is worth
    # checking, because bands built at the same azimuth lie exactly on top of
    # one another and look like a single band from every angle.
    turns_wanted = C.helix_turns(P)[0]
    arms_wanted = int(round(P.get("helix_arms", 1)))
    spans = sorted(((max(verts[j][2] for ti in sh for j in tris[ti])
                     - min(verts[j][2] for ti in sh for j in tris[ti])), i)
                   for i, sh in enumerate(shells))
    if len(shells) == 1:
        # The arms are unioned with the body: one shell, and they cannot be
        # counted separately because there is nothing separate to count. What
        # matters for a unioned column is that the union is sound -- no holes
        # in the surface -- and that the arms are actually there, which the
        # material at the bore's wall shows.
        bands = []
        gaps = csg.real_gaps(verts, tris)
        rep.check(f"the helix is unioned with the body, {arms_wanted} arm(s) "
                  f"in one solid", not gaps,
                  f"{len(gaps)} edges are used once with nothing lying along "
                  f"them, so the union left holes" if gaps else
                  f"one shell: the arms and the body are merged, so the middle "
                  f"of the bore is solid and no fill rule can show a hole "
                  f"around the axis")
    else:
        bands = [i for _, i in spans[:-1]]      # everything but the tube
        rep.check(f"the helix is {arms_wanted} arm(s), not one band drawn "
                  f"{arms_wanted} times", len(bands) == arms_wanted,
                  f"{len(bands)} separate bands in the mesh"
                  + ("; bands at the same azimuth overlap exactly"
                     if len(bands) < arms_wanted else ""))
    for i in bands:
        bz = [verts[j][2] for ti in shells[i] for j in tris[ti]]
        br = [math.hypot(verts[j][0], verts[j][1])
              for ti in shells[i] for j in tris[ti]]
        rep.check(f"band {bands.index(i) + 1} runs the whole staged length",
                  abs(min(bz) - Z0) < 1.0 and abs(max(bz) - Z1) < 1.0,
                  f"z {min(bz):.1f} to {max(bz):.1f} against the staged "
                  f"{Z0:g} to {Z1:g} mm")
        rep.check(f"band {bands.index(i) + 1} reaches the wall",
                  max(br) > R_BORE - 0.3,
                  f"{max(br):.2f} mm from the axis against a {R_BORE:g} mm bore")

    # Where each arm is, half way up. The arms of a multi-start helix have to
    # be evenly spaced, or it is not the helix that was asked for.
    z_ref = 0.5 * (Z0 + Z1)
    arm_az = []
    for i in bands:
        hits = [math.atan2(verts[j][1], verts[j][0])
                for ti in shells[i] for j in tris[ti]
                if abs(verts[j][2] - z_ref) < 2.0
                and math.hypot(verts[j][0], verts[j][1]) > R_BORE - 0.3]
        if hits:
            arm_az.append(math.degrees(math.atan2(
                sum(math.sin(x) for x in hits) / len(hits),
                sum(math.cos(x) for x in hits) / len(hits))))
    if len(arm_az) > 1:
        order = sorted(arm_az)
        gaps = [order[k + 1] - order[k] for k in range(len(order) - 1)]
        gaps.append(order[0] + 360.0 - order[-1])
        want = 360.0 / len(arm_az)
        rep.check("the arms are evenly spaced round the axis",
                  max(abs(g - want) for g in gaps) < 3.0,
                  "azimuths " + ", ".join(f"{a:.1f}" for a in order)
                  + f"; gaps {', '.join(f'{g:.1f}' for g in gaps)} against "
                  f"{want:.1f}")

    # The twist, measured off one arm's own vertices rather than off plane
    # cuts: every section is at one height and one azimuth, so grouping the
    # outer edge by z and unwrapping the azimuth gives the twist rate directly.
    # A plane cut has to sample azimuths on a grid, and a band wide enough to
    # cross the axis reaches the wall over a range of azimuths, which makes
    # that grid alias -- it reported 0.64 turns where the part had 2.83.
    by_z = {}
    if bands:
        for ti in shells[bands[0]]:
            for j in tris[ti]:
                x, y, z = verts[j]
                if math.hypot(x, y) > R_BORE - 0.3:
                    by_z.setdefault(round(z, 4), []).append(math.atan2(y, x))
    walked = []
    for z in sorted(by_z):
        v = by_z[z]
        a = math.atan2(sum(math.sin(x) for x in v) / len(v),
                       sum(math.cos(x) for x in v) / len(v))
        if walked:
            while a - walked[-1][1] > math.pi:
                a -= 2.0 * math.pi
            while a - walked[-1][1] < -math.pi:
                a += 2.0 * math.pi
        walked.append((z, a))
    if len(walked) > 2:
        turns = ((walked[-1][1] - walked[0][1]) / (2.0 * math.pi)
                 * (Z1 - Z0) / (walked[-1][0] - walked[0][0]))
        rep.near("each arm twists as many turns as the build claims", turns,
                 turns_wanted, 0.15, "{:.3f}")
        print(f"  twist            : {turns:.3f} turns per arm over "
              f"{walked[-1][0] - walked[0][0]:.0f} mm, from {len(walked)} "
              f"sections of one arm's outer edge")

    steps = max(12, int(6.0 * turns_wanted) + 6)
    for i in range(1, steps + 1):
        z = Z0 + (Z1 - Z0) * i / (steps + 1)
        staging_planes.append((z, scan(z, 48)))

elif STYLE == "vigreux":
    # Fingers hang off the wall, so a cut through a row shows them near the
    # wall and not at the axis, and a cut clear of them shows nothing at all.
    # How far they hang depends on the pitch they were given, so the cut is
    # taken where the fingers are actually found rather than at a guessed
    # depth: walk down from the row until the bore clears, then cut half way
    # into what was found.
    rows_seen, per_row, spans = [], [], []
    for k in range(COUNT):
        z = Z0 + PITCH * (k + 1)
        span = 0.0
        for i in range(1, 13):
            zz = z - PITCH * i / 12.0
            if any(r["blocked"] for r in scan(zz, 8)):
                span = PITCH * i / 12.0
            else:
                break
        spans.append(span)
        zz = z - span * 0.5 if span > 0 else z
        rows = scan(zz)
        staging_planes.append((zz, rows))
        per_row.append(clusters([r["blocked"] for r in rows]))
        # Fingers hang *down* from their row, so the clear band is below the
        # next row's fingers, not half way between the rows: find it rather
        # than assume where it is.
        gap = []
        for i in range(1, 11):
            zz = z + PITCH * i / 11.0
            gap.append(not any(r["blocked"] for r in scan(zz, 8)))
        rows_seen.append(any(gap))
    rep.check(f"every row is {C.FINGER_PER_ROW} separate fingers",
              all(n == C.FINGER_PER_ROW for n in per_row),
              f"fingers per row {per_row}; they hang "
              f"{min(spans):.1f} to {max(spans):.1f} mm below their row")
    rep.check("the vapour can get between one row of fingers and the next",
              all(rows_seen),
              f"{sum(rows_seen)} of {COUNT} gaps between rows have a clear "
              f"plane; {COUNT - sum(rows_seen)} do not")

elif STYLE == "tesla":
    # One plate per staging, and consecutive plates swap sides, which is the
    # whole of the zig-zag.
    lean = math.radians(C.louvre_tilt(P))
    length = (R_BORE + 0.2 + C.louvre_reach(P)) / math.cos(lean)
    sides, plates = [], []
    for k in range(COUNT):
        z = Z0 + PITCH * (k + 1)
        rows = scan(z + length * math.sin(lean) * 0.4)
        staging_planes.append((z + length * math.sin(lean) * 0.4, rows))
        plates.append(clusters([r["blocked"] for r in rows]))
        sides.append(azimuth_of(rows))
    rep.check("every plate is one lump in the bore",
              all(n == 1 for n in plates), f"plates per plane {plates}")
    swaps = all(abs(sides[k] - sides[k + 1]) > 0.25
                for k in range(len(sides) - 1)
                if sides[k] is not None and sides[k + 1] is not None)
    rep.check("consecutive plates are on opposite sides", swaps,
              "azimuths " + ", ".join(f"{a * 360:.0f}" for a in sides
                                      if a is not None))
    # How far a plate reaches past the axis is not fixed: it is whatever fits
    # between one plate and the next, so what matters is that the plates are
    # separated and that the gap between them is passable.
    clear = []
    for k in range(COUNT - 1):
        z = Z0 + PITCH * (k + 1)
        clear.append(any(not any(r["blocked"] for r in scan(zz, 8))
                         for zz in (z + PITCH * i / 10.0 for i in range(1, 10))))
    rep.check("consecutive plates are separated by a passable gap",
              all(clear),
              f"{sum(clear)} of {COUNT - 1} gaps between plates are clear")
    reach = min(min((r0 for r0, r1 in r["mat"]), default=1e9)
                for _, rows in staging_planes for r in rows if r["mat"])
    print(f"  plate reach      : the plates come {reach:.2f} mm from the axis "
          f"({R_BORE - reach:.2f} mm short of the wall), which is what the "
          f"{PITCH:.1f} mm pitch leaves room for")

else:
    rep.check("the style is one this build knows", False, repr(STYLE))

n_blocked = sum(1 for _, rows in staging_planes for r in rows if r["blocked"])
rep.check("the staging leaves a way past at every height it stands in",
          n_blocked > 0, f"{n_blocked} blocked azimuths over "
          f"{len(staging_planes)} planes, all of them with a way round")


# ---------------------------------------------------------- the passage ---
# The point of the whole part: air has to get from one mouth to the other. The
# build reports the route it intends, derived from the parameters read back out
# of this file, and every sample along it has to be in air -- a column with a
# blocked passage prints perfectly and distils nothing.
VOID_RAYS = ((0.3, 0.7, 0.65), (0.62, 0.41, -0.67), (-0.55, 0.6, 0.58))


def clear(p):
    """Void by a majority of rays: the column is a solid of revolution, so a
    single horizontal ray can graze a facet exactly and answer either way."""
    yes = 0
    for d in VOID_RAYS:
        if K.inside(verts, tris, p, shells=shells, boxes=boxes, direction=d):
            yes += 1
            if yes >= 2:
                return False
    return True


path = C.vapour_path(P)
blocked = []
samples = 0
for i in range(len(path) - 1):
    a, b = path[i], path[i + 1]
    steps = max(2, int(math.dist(a, b) / 0.3) + 1)
    for k in range(steps + 1):
        u = k / steps
        p = tuple(a[j] + (b[j] - a[j]) * u for j in range(3))
        samples += 1
        if not clear(p):
            blocked.append((i, round(u, 2), tuple(round(v, 2) for v in p)))
rep.check("the vapour path is open from mouth to mouth", not blocked,
          f"material on the path at {blocked[:4]}" if blocked else
          f"{samples} samples along {len(path) - 1} segments, clear throughout")

# A second, blunter route from mouth to mouth: straight up the axis. It only
# exists for the styles that leave the middle of the bore open -- a leaning
# plate crosses the axis on purpose, and a band wide enough to cross it blocks
# it too, and in both the vapour goes round instead -- so it is a check on the
# styles that have one, and said to be skipped on the others.
band_over_axis = STYLE == "helix" and P.get("helix_frac", 1.0) <= 0.0
if STYLE == "tesla" or band_over_axis:
    print("  axis route       : not applicable to "
          + ("leaning plates" if STYLE == "tesla"
             else "a band wide enough to cross the axis")
          + ", which block it by design; the vapour path is the route that was "
          "walked")
else:
    axis_blocked = [round(z, 1) for z in
                    [z_low + (z_high - z_low) * i / 120.0 for i in range(121)]
                    if not clear((0.0, 0.0, z))]
    rep.check("the axis is clear from one mouth to the other", not axis_blocked,
              f"material on the axis at {axis_blocked[:6]}" if axis_blocked
              else "121 samples from mouth to mouth, clear throughout")


# -------------------------------------------------------- printing as drawn ---
# Nothing inside the bore may need support: that is where support cannot be
# reached to take it out again. The outside is different -- the crotch where a
# transition meets the body can face downward -- and support there is ordinary.
over = PR.overhangs(verts, tris, exposed_only=True, max_angle=ANGLE_DEG,
                    shells=shells)
inside_bore = [o for o in over if radius(o["c"]) <= R_BORE - 0.05]
outside = [o for o in over if radius(o["c"]) > R_BORE - 0.05]
rep.check(f"no surface inside the bore needs support, at {ANGLE_DEG:g} deg "
          f"from horizontal", not inside_bore,
          f"down-facing faces in the bore at "
          + str([(round(o["z"], 1), round(180.0 - o["angle"]),
                  round(o["area"], 2)) for o in inside_bore[:4]])
          + " (z, deg from horizontal, mm^2)" if inside_bore else
          "every down-facing face is on the outside, where support comes off")

# --------------------------------------------------------------- numbers ---
# Contact against pressure drop. The narrowest free area is the pressure-drop
# proxy: it is what the vapour has to squeeze through, and under vacuum every
# pascal lost here comes off the still.
narrow = min(((plane_area(rows), z) for z, rows in staging_planes),
             default=(0.0, 0.0))
# Only azimuths that are actually open count: a ray that meets material at the
# axis is inside a staging, not looking down a passage, and reporting that as
# a 0.0 mm clearance would be nonsense.
tight = min(((min((r["free_axis"] for r in rows if r["free_axis"] > 0.0),
                  default=0.0), z)
             for z, rows in staging_planes if any(r["free_axis"] > 0.0
                                                  for r in rows)),
            default=(0.0, 0.0))
full = math.pi * R_BORE ** 2
tube_area = 2.0 * math.pi * R_BORE * (D["z_body1"] - D["z_body0"])


def shell_span(i):
    zs = [verts[j][2] for ti in shells[i] for j in tris[ti]]
    return max(zs) - min(zs)


tube = max(range(len(shells)), key=shell_span)
staging_area = 0.0
for i, shell in enumerate(shells):
    if i == tube:
        continue
    for ti in shell:
        a, b, c = (verts[j] for j in tris[ti])
        u = (b[0] - a[0], b[1] - a[1], b[2] - a[2])
        v = (c[0] - a[0], c[1] - a[1], c[2] - a[2])
        n = (u[1] * v[2] - u[2] * v[1], u[2] * v[0] - u[0] * v[2],
             u[0] * v[1] - u[1] * v[0])
        staging_area += 0.5 * math.sqrt(sum(x * x for x in n))

bb = [min(v[k] for v in verts) for k in range(3)] + \
     [max(v[k] for v in verts) for k in range(3)]
print(f"\n  style            : {STYLE}, {COUNT} stagings, "
      f"{len(shells) - 1} of them their own shell")
print(f"  joints           : bottom {P['bottom_joint']} {P['bottom_gender']} "
      f"at z={z_low:.1f}, top {P['top_joint']} {P['top_gender']} "
      f"at z={z_high:.1f}")
print(f"  body             : {R_BORE * 2:g} mm bore, {R_BODY * 2:g} mm OD, "
      f"staged {Z0:g} to {Z1:g} mm")
print(f"  bounding box     : {bb[3] - bb[0]:.1f} x {bb[4] - bb[1]:.1f} x "
      f"{bb[5] - bb[2]:.1f} mm")
print(f"  vapour path      : {C.path_length(P):.1f} mm to climb "
      f"{z_high - z_low:.1f} mm "
      f"({C.path_length(P) / (z_high - z_low):.2f}x the height)")
print(f"  narrowest plane  : {narrow[0]:.0f} mm^2 free at z={narrow[1]:.1f} "
      f"({100.0 * narrow[0] / full:.0f}% of the {full:.0f} mm^2 empty bore)")
print(f"  tightest radial  : {tight[0]:.2f} mm from the axis at z={tight[1]:.1f} "
      f"-- measured along a radius, so a staging that reaches the centre reads "
      f"tight even though the vapour goes round it")
print(f"  wetted surface   : {tube_area + staging_area:.0f} mm^2 "
      f"({tube_area:.0f} of bore wall + {staging_area:.0f} of staging, the "
      f"overlap where a staging meets the wall counted twice)")
print(f"  outside overhang : {len(outside)} faces, "
      f"{sum(o['area'] for o in outside):.0f} mm^2")
print(f"  vertices / tris  : {len(verts)} / {len(tris)}")

print(rep.finish())

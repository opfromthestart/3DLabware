#!/usr/bin/env python3
"""Standard taper joints, and the branching junction they meet at.

Shared by `make_wye.py` and `make_collector.py`; both are several arms leaving
one junction, and the hard part is the same in each.

THE JUNCTION. Several tubes cannot be glued into a branched passage at one
point: an arm's wall runs straight through the other arms' bores and plugs
them, and what comes out is a sealed node with, at best, a crescent of a gap
left over. What a branch *is*, and what `junction` builds, is

    (the arms' outer solids u a blend ball) - (the arms' bores)

with the bores taken out *after* the union, so nothing is left standing inside
them. Every one of those primitives is convex and contains the junction, so a
ray leaving the junction crosses each of them in one interval; the material is
therefore star-shaped about the junction -- one interval rho(u)..R(u) per
direction, rho being the furthest a bore reaches along u and R the furthest an
outer solid or the ball reaches. That is what makes the junction meshable
radially, with no boolean kernel behind it.

Each arm's tube stops short of the junction by its `trim`: the distance at
which its wall would start to lie inside another arm's bore. Out at that plane
the tube adds nothing the junction does not already carry, and any closer in it
would plug a passage.

The ball, of the tube's own radius, comes with the outer solids: the arms are
tangent to it, and in an arrangement where every arm leans the same way it is
the ball that closes the back of the junction. Arms spanning more than a
half-space -- an opposite pair -- already cover it with their own outer solids.

The taper joints are 1:10 on the diameter, the same rule `make_elbow.py`
states and every part here shares; a real ground joint is 1:20, so these mate
with each other rather than with glassware.
"""

import math

import threelib as t

# Standard taper joints: mouth diameter / ground length, in mm.
JOINTS = {
    "5/13": (5.0, 13.0), "7/16": (7.0, 16.0), "10/18": (10.0, 18.0),
    "12/21": (12.0, 21.0), "14/20": (14.0, 20.0), "19/22": (19.0, 22.0),
    "24/40": (24.0, 40.0), "29/42": (29.0, 42.0), "34/45": (34.0, 45.0),
    "40/50": (40.0, 50.0), "45/50": (45.0, 50.0), "50/50": (50.0, 50.0),
}
TAPER = 1.0 / 20.0

# Construction margins, not things to choose a part by.
TRIM_MARGIN = 0.05   # how far past the plugging distance each wall stops
CLIP_MARGIN = 0.4    # how far the junction reaches past the arms' trimmed ends

# How far the junction's solid is inflated, on both surfaces, so that a boolean
# has something to resolve. The blend follows each arm's profile, so its
# surfaces land *exactly* on the arm's -- two surfaces in the same place, which
# a BSP cannot resolve, because it needs them to cross.
#
# The direction matters and was got wrong first. *Shrinking* the cavity pushes
# it into the bores: near each arm the junction's cavity is that arm's bore, so
# the inlet came out 5.80 mm against 6.00, its taper read 1:8.3 instead of 1:10,
# and the centre where four bores meet filled in solid. *Inflating* the solid --
# cavity out as well as skin out -- leaves the bores their own width and puts
# the junction's inner surface inside each arm's wall, where the union removes
# it and the arm's bore survives as the boundary.
JUNCTION_BITE = 0.25


def parse_joint(spec):
    if isinstance(spec, str):
        dia, length = spec.split("/")
        return float(dia), float(length)
    return float(spec[0]), float(spec[1])


def joint_dims(spec):
    """(mouth radius, tip radius, ground length)."""
    dia, length = parse_joint(spec)
    r_mouth = dia / 2.0
    return r_mouth, r_mouth - length * TAPER, length


def end_bore(spec, gender, tip_wall):
    """Bore radius at that arm's mouth."""
    _, r_tip, _ = joint_dims(spec)
    return r_tip if gender == "female" else r_tip - tip_wall


def arm(tag, spec, gender, length, direction, wall, tip_wall):
    """One arm of a junction: a tapered joint on a tube aimed at the centre."""
    r_mouth, r_tip, glen = joint_dims(spec)
    return {"tag": tag, "spec": spec, "gender": gender, "len": length,
            "dir": direction, "r_mouth": r_mouth, "r_tip": r_tip, "glen": glen,
            "bore": end_bore(spec, gender, tip_wall)}


def inner_end_outer(a, wall):
    """Outside radius where that arm's joint ends and its tube would begin."""
    return a["r_tip"] + wall if a["gender"] == "female" else a["r_mouth"]


def tube_radius(arms, wall, bore_j):
    """The tube's outside radius: as fat as the widest joint allows."""
    return max(max(inner_end_outer(a, wall) for a in arms), bore_j + wall)


def arm_profile(a, s, wall, flare, tube_r, bore_j):
    """(bore radius, outer radius) of one arm at arc length `s` from its mouth.

    Piecewise linear in `s`, breaking where a part's swept stations break: the
    end of the ground zone, the end of the flare, and the junction. The
    junction piece measures its rays against this same profile, so the two
    agree about where each surface is.
    """
    fem = a["gender"] == "female"
    glen, length = a["glen"], a["len"]
    if s <= glen:
        # The bore is the socket's cone (female) or straight through (male);
        # the outer surface is the joint's own cone either way.
        bore = (a["r_mouth"] - TAPER * s) if fem else a["bore"]
        cone = (a["r_mouth"] - TAPER * s + wall) if fem else (a["r_tip"] + TAPER * s)
        return bore, cone
    # Past the ground zone the bore runs to the junction's, so a female socket
    # ends at its tip radius with no shoulder behind it.
    start = a["r_tip"] if fem else a["bore"]
    bore = start + (bore_j - start) * (s - glen) / (length - glen)
    cone_end = (a["r_mouth"] - TAPER * glen + wall) if fem else (a["r_tip"] + TAPER * glen)
    span = flare if flare > 0.0 else (length - glen)
    if s < glen + span:
        # A straight ramp from the boss down to the tube.
        return bore, cone_end + (tube_r - cone_end) * (s - glen) / span
    return bore, tube_r


def arm_waypoints(a, wall, flare, tube_r, bore_j):
    """The arc lengths along one arm where its profile changes slope."""
    glen = a["glen"]
    span = flare if flare > 0.0 else (a["len"] - glen)
    return [0.0, glen, min(a["len"], glen + span), a["len"]]


def reach(a, u, bore, limit, wall, flare, tube_r, bore_j):
    """How far along `u` that arm's bore, or its outer solid, reaches.

    Measured from the junction, which every one of those solids contains. Round
    the junction the profiles only widen going out, so a ray from there crosses
    each solid exactly once and a bisection finds where.

    The closed end planes matter: a ray lying in the junction plane is still
    inside an arm's bore at the radius of its end disc, and reporting zero
    there would drag the junction mesh right across the passage.
    """
    c = sum(u[k] * a["dir"][k] for k in range(3))
    sp = math.sqrt(max(0.0, 1.0 - c * c))
    hi = limit if c <= 1e-12 else min(limit, a["len"] / c)
    lo = 0.0
    for _ in range(26):
        mid = 0.5 * (lo + hi)
        s = a["len"] - mid * c
        ri, ro = arm_profile(a, s, wall, flare, tube_r, bore_j)
        if -1e-9 <= s <= a["len"] + 1e-9 and mid * sp <= (ri if bore else ro) + 1e-9:
            lo = mid
        else:
            hi = mid
    # The bisection assumes the ray crosses the solid once, so check that it
    # did. A profile that narrows going out and then widens again -- which is
    # just what a female socket's outer cone does, along with every joint's
    # shaped end -- can cross a ray three times, and the answer would then
    # depend on where the bisection happened to start. Whether that can
    # actually happen depends on the ray as well as the profile: the crossing
    # is unique while the ray spreads sideways faster than the profile grows.
    held = True
    for k in range(49):
        r = hi * k / 48.0
        s = a["len"] - r * c
        ri, ro = arm_profile(a, s, wall, flare, tube_r, bore_j)
        inside = (-1e-9 <= s <= a["len"] + 1e-9 and
                  r * sp <= (ri if bore else ro) + 1e-9)
        if inside and not held:
            raise ValueError(
                f"arm {a['tag']}'s {'bore' if bore else 'outer'} surface is "
                f"crossed more than once by a ray from the junction, so the "
                f"junction cannot be swept radially over it; steady the "
                f"profile of its {a['spec']} joint with a longer arm or a "
                f"shorter ground zone")
        held = inside
    return lo


def junction_radii(arms, u, r_cap, tube_r, wall, flare, bore_j,
                   bite=JUNCTION_BITE):
    """(inner, outer) radius of the junction solid along the direction `u`.

    `bite` inflates the junction -- the skin out and the cavity out -- so that
    the arms' surfaces, which the blend follows exactly, no longer coincide
    with it. See `JUNCTION_BITE` for why the direction matters.
    """
    limit = r_cap + 0.5
    inner = 0.0
    outer = tube_r                       # the blend ball reaches its own radius
    for a in arms:
        inner = max(inner, reach(a, u, True, limit, wall, flare, tube_r, bore_j))
        outer = max(outer, reach(a, u, False, limit, wall, flare, tube_r, bore_j))
    if inner > 0.0:
        inner = inner + bite
    return inner, min(r_cap + bite, outer + bite)


def junction(mesh, arms, seg, rings, r_cap, tube_r, wall, flare, bore_j):
    """Mesh the junction: (the outer solids u the blend ball) - (the bores).

    Swept radially, direction by direction: a vertex out on each ray where a
    bore ends and another where the outer skin ends. Where no bore reaches the
    ray, the inner radius is zero and the material starts at the junction
    itself, which is how the ball closes the back of the part.

    Directions where the bores cover the ray's whole length -- the openings --
    carry no surface, so cells there are dropped and the rim between the
    passage wall and the outer skin is fenced off along the boundary. That
    fence sits inside the material, right under an arm's own wall, so it never
    shows in a passage.
    """
    U, MAT, RHO, RHI = [], [], [], []
    for i in range(rings + 1):
        th = math.pi * i / rings
        for j in range(seg):
            ph = 2.0 * math.pi * j / seg
            u = (math.sin(th) * math.cos(ph), math.sin(th) * math.sin(ph),
                 math.cos(th))
            inner, outer = junction_radii(arms, u, r_cap, tube_r, wall, flare,
                                          bore_j)
            U.append(u)
            RHO.append(inner)
            RHI.append(outer)
            MAT.append(inner < outer - 1e-9)
    IDX = [[i * seg + j for j in range(seg)] for i in range(rings + 1)]

    verts, faces = [], []
    inner_id, outer_id = {}, {}

    def vertex(kind, i, j):
        cache = inner_id if kind == "i" else outer_id
        pole = i == 0 or i == rings
        key = (kind, i, 0 if pole else j)
        if key in cache:
            return cache[key]
        k = IDX[i][j]
        r = RHO[k] if kind == "i" else RHI[k]
        cache[key] = len(verts)
        verts.append(tuple(x * r for x in U[k]))
        return len(verts) - 1

    def cell_material(i, j):
        j2 = (j + 1) % seg
        return (MAT[IDX[i][j]] and MAT[IDX[i + 1][j]] and
                MAT[IDX[i + 1][j2]] and MAT[IDX[i][j2]])

    def add(ids, ref):
        """One face, oriented so its normal points along `ref`."""
        ids = [v for n, v in enumerate(ids) if n == 0 or v != ids[n - 1]]
        if len(ids) > 2 and ids[0] == ids[-1]:
            ids.pop()
        if len(ids) < 3:
            return
        a, b, c = (verts[v] for v in ids[:3])
        e1 = [b[k] - a[k] for k in range(3)]
        e2 = [c[k] - a[k] for k in range(3)]
        nrm = (e1[1] * e2[2] - e1[2] * e2[1],
               e1[2] * e2[0] - e1[0] * e2[2],
               e1[0] * e2[1] - e1[1] * e2[0])
        cen = [sum(verts[v][k] for v in ids) / len(ids) for k in range(3)]
        if sum(nrm[k] * ref[k] for k in range(3)) < 0.0:
            ids = list(reversed(ids))
        faces.append(tuple(ids))

    for i in range(rings):
        for j in range(seg):
            if not cell_material(i, j):
                continue
            j2 = (j + 1) % seg
            corners = [(i, j), (i + 1, j), (i + 1, j2), (i, j2)]
            out = [sum(U[IDX[a][b]][k] for a, b in corners) / 4.0 for k in range(3)]
            add([vertex("o", a, b) for a, b in corners], out)       # outer skin
            add([vertex("i", a, b) for a, b in corners],
                [-x for x in out])                                  # passage wall
            # Fence the rim: every cell edge whose far side carries no surface.
            for (p, q), (ni, nj) in (((0, 1), (i, j - 1)), ((1, 2), (i + 1, j)),
                                     ((2, 3), (i, (j + 1) % seg)),
                                     ((3, 0), (i - 1, j))):
                if not (0 <= ni < rings) or not (0 <= nj < seg):
                    continue
                if cell_material(ni, nj):
                    continue
                a, b = corners[p], corners[q]
                ia, ib = vertex("i", *a), vertex("i", *b)
                oa, ob = vertex("o", *a), vertex("o", *b)
                if ia == ib and oa == ob:
                    continue
                nj2 = (nj + 1) % seg
                hole = [U[IDX[ni][nj]][k] + U[IDX[ni + 1][nj]][k] +
                        U[IDX[ni + 1][nj2]][k] + U[IDX[ni][nj2]][k]
                        for k in range(3)]
                add([ia, ib, ob, oa], hole)
    mesh.add_solid(verts, faces)
    return len(verts), len(faces)


def trim_of(a, others, wall, flare, tube_r, bore_j):
    """How far that arm's wall must stop short of the junction.

    The wall reaches into another arm's bore wherever the two are closer than
    that bore's radius; walk out from the junction and find the last station
    that still does, over every other arm. Out past that plane the wall is
    clean, so it may run to the mouth; any closer in it would be a ridge across
    a passage.
    """
    need = 0.0
    for b in others:
        c = sum(a["dir"][k] * b["dir"][k] for k in range(3))
        c = max(-1.0, min(1.0, c))
        sin_t = math.sqrt(max(0.0, 1.0 - c * c))
        if sin_t < 1e-9:
            continue                      # collinear: one shared cylinder
        steps = int(a["len"] / 0.01) + 1
        for m in range(steps, -1, -1):
            axial = min(a["len"], m * 0.01)
            ri, ro = arm_profile(a, a["len"] - axial, wall, flare, tube_r, bore_j)
            lo, hi = ri * c, ro * c
            if min(lo, hi) < axial * sin_t + bore_j and max(lo, hi) > axial * sin_t - bore_j:
                need = max(need, axial)
                break
    return need


def plan(arms, wall, flare, bore_j, tol):
    """Trim lengths, junction reach and facet counts for a set of arms.

    Everything that depends on the arms as a group rather than on one of them:
    how far each wall pulls back, how far the junction must then reach to meet
    them, and how finely to facet the round surfaces.
    """
    cuts = [trim_of(a, [b for b in arms if b is not a], wall, flare,
                    tube_radius(arms, wall, bore_j), bore_j) + TRIM_MARGIN
            for a in arms]
    tube_r = tube_radius(arms, wall, bore_j)
    r_cap = 0.0
    for a, cut in zip(arms, cuts):
        r_cap = max(r_cap, math.hypot(cut, arm_profile(
            a, a["len"] - cut, wall, flare, tube_r, bore_j)[1]))
    r_cap += CLIP_MARGIN
    return cuts, r_cap, tube_r, t.auto_segments(tube_r, tol)


def check_arms(arms, cuts, r_cap, wall, flare, tube_r, bore_j):
    """Reject the arrangements the radial junction cannot describe.

    What has to hold for the radial sweep is checked where it is used: `reach`
    confirms, ray by ray, that the ray crosses each arm's surfaces once. What
    is left here is that an arm still has a tube to sweep after its wall has
    been pulled back off the junction.
    """
    for a, cut in zip(arms, cuts):
        if a["len"] - cut <= 0.0:
            raise ValueError(f"arm {a['tag']} is too short for the junction")


def sweep_arm(mesh, a, cut, n, wall, flare, tube_r, bore_j):
    """Sweep one arm's tube from its mouth in to the junction, stopping short.

    The stations are the profile's break points plus the trimmed end; the
    profile runs straight between breaks, so landing the last station on that
    same line keeps the arm's surface exactly where the junction piece looks
    for it.
    """
    length = a["len"]
    pts = [s for s in arm_waypoints(a, wall, flare, tube_r, bore_j)
           if s < length - cut - 1e-9]
    pts.append(length - cut)
    # The arm runs from its mouth *inwards* to the junction, so arc length s is
    # measured from the mouth and the station sits (length - s) out along the
    # arm's direction. The tangent is therefore -dir: the loft advances in
    # +tangent, and pointing it the other way would turn every face inside out.
    back = tuple(-c for c in a["dir"])
    stations = []
    for s in pts:
        ri, ro = arm_profile(a, s, wall, flare, tube_r, bore_j)
        pos = tuple(a["dir"][k] * (length - s) for k in range(3))
        stations.append((pos, back, ri, ro))
    t.sweep(mesh, stations, n, ref=(0.0, 1.0, 0.0))

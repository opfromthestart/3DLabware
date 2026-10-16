#!/usr/bin/env python3
"""Verify a fraction collector 3MF: joints, the rotation axis, the passage.

Re-parsed from disk and measured from the mesh. The four arms are found without
being told where they are: each ends in a mouth annulus, and a mouth's plane
normal is that arm's axis. The female socket is told from the male outlets by
its radii, and the three outlets have to be evenly spaced about the inlet's
axis, because that is what makes one index position the next.

The two things that make this part work are checked directly. Every pair of
mouths has to be joined through the junction -- an adapter whose bores stop
short of each other is not a fraction collector -- and the whole thing has to
print as it stands, with no support inside the passage.

Usage: python3 verify_collector.py [collector.3mf]
"""

import math
import sys

import checks as K
import junctionlib as J
import make_collector as C
import printability as PR

PATH, ANGLE = K.cli(sys.argv[1:], C.default_out())
ANGLE_DEG = PR.OVERHANG_DEG if ANGLE is None else ANGLE
P = K.read_params(PATH)

rep = K.Report(f"Fraction collector: {PATH}")
objs, _, _ = K.load_3mf(PATH, rep)
rep.check("one object", len(objs) == 1, f"{len(objs)}")
verts, tris = objs[0]["verts"], objs[0]["tris"]
shells = K.shell_triangles(verts, tris)


def sub(a, b):
    return (a[0] - b[0], a[1] - b[1], a[2] - b[2])


def dot(a, b):
    return a[0] * b[0] + a[1] * b[1] + a[2] * b[2]


def cross(a, b):
    return (a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2],
            a[0] * b[1] - a[1] * b[0])


def unit(a):
    m = math.sqrt(dot(a, a))
    return (a[0] / m, a[1] / m, a[2] / m)


# ---------------------------------------------------------------- mouths ---
# A mouth is an annular face at the end of an arm; a cap where an arm's wall
# stops short of the junction is one too. They are told apart by whether
# material lies behind the face: step off a mouth and you are in open air.
faces = K.annuli(verts, tris, shells=shells)
outside = [f for f in faces if not f["buried"]]
outside.sort(key=lambda f: (-f["n_tri"], -f["area"]))
rep.check("one mouth annulus per arm", len(outside) == 1 + P["count"],
          f"{len(faces)} annular faces, {len(outside)} of them on the outside")


def end_bore(spec, gender, tip_wall):
    _, r_tip, _ = J.joint_dims(spec)
    return r_tip if gender == "female" else r_tip - tip_wall


def signature(spec, gender, tip_wall):
    """The (outer, inner) radii a mouth of that joint should show."""
    r_mouth, r_tip, _ = J.joint_dims(spec)
    if gender == "male":
        return r_tip, end_bore(spec, gender, tip_wall)
    return r_mouth + P["wall"], r_mouth


def mismatch(face, spec, gender):
    want = signature(spec, gender, P["tip_wall"])
    return abs(face["outer"] - want[0]) + abs(face["inner"] - want[1])


inlet = min(outside, key=lambda f: mismatch(f, P["inlet_joint"], "female"))
if mismatch(inlet, P["inlet_joint"], "female") > 0.2:
    inlet = None
    rep.check("the inlet is a female socket", False, "no matching mouth")
else:
    rep.check("the inlet is a female socket", True,
              f"{P['inlet_joint']} socket, mouth {inlet['outer'] * 2:.2f} mm "
              f"outside {inlet['inner'] * 2:.2f} mm inside")
outlets = [f for f in outside if f is not inlet]
costs = [mismatch(f, P["out_joint"], "male") for f in outlets]
rep.check("the other mouths are the outlet joints",
          len(outlets) == P["count"] and max(costs) < 0.2,
          f"total mismatch {sum(costs):.4f} mm across {len(outlets)} outlets")

inlet_axis = inlet["n"]                    # out of the socket, along +z
outlet_axes = [f["n"] for f in outlets]    # out of each outlet, pointing down

# ------------------------------------------------------- the rotation axis ---
# The adapter turns on the inlet, so the inlet's axis *is* the axis of
# rotation, and it is the part's +z by construction.
rep.near("the inlet stands on the +z axis", inlet_axis[2], 1.0, 0.002,
         "{:.4f}")
tilt = math.degrees(math.acos(max(-1.0, min(1.0, inlet_axis[2]))))
rep.check("the rotation axis is upright", tilt < 0.2, f"{tilt:.4f} deg off +z")

# The two angles the part is asked for, measured where they name something:
# inlet to outlet, and outlet to its neighbour.
angles = [math.degrees(math.acos(max(-1.0, min(1.0, dot(u, inlet_axis)))))
          for u in outlet_axes]
rep.check("every outlet hangs at the angle asked for off the inlet",
          max(abs(a - P["in_out_angle"]) for a in angles) < 0.1,
          "/".join(f"{a:.3f}" for a in angles)
          + f" vs {P['in_out_angle']:g} deg")

# Which outlet neighbours which is a matter of the order they were written in,
# and the axes are only found by their mouth planes, so take the smallest
# angle in the ring rather than trusting the order.
n_out = len(outlet_axes)
ring_pairs = [(i, (i + 1) % n_out, math.degrees(math.acos(max(-1.0, min(
    1.0, dot(outlet_axes[i], outlet_axes[(i + 1) % n_out]))))))
    for i in range(n_out)]
# Sort the axes around the axis and re-pair, so "the next one round" means
# what it says even if the outlet list came back shuffled.
perp = []
for u in outlet_axes:
    v = tuple(u[k] - dot(u, inlet_axis) * inlet_axis[k] for k in range(3))
    perp.append(math.atan2(dot(v, (0.0, 1.0, 0.0)), dot(v, (1.0, 0.0, 0.0)))
                % (2.0 * math.pi))
order = sorted(range(n_out), key=lambda i: perp[i])
steps = []
for i in range(n_out):
    a, b = order[i], order[(i + 1) % n_out]
    steps.append(math.degrees((perp[b] - perp[a]) % (2.0 * math.pi)))
    ring_pairs[i] = (a, b, math.degrees(math.acos(max(-1.0, min(
        1.0, dot(outlet_axes[a], outlet_axes[b]))))))
neighbour = [p[2] for p in ring_pairs]
pitch = neighbour[:-1]                 # the fan's own steps, in order
wrap = neighbour[-1]                   # last back round to first
rep.check("consecutive outlets stand apart by the angle asked for",
          max(abs(v - P["out_out_angle"]) for v in pitch) < 0.1,
          "/".join(f"{v:.3f}" for v in pitch)
          + f" vs {P['out_out_angle']:g} deg")

# Around the axis is not a parameter either: the two angles imply it, and what
# the implication buys you is whether the outlets close a ring or stand in a
# fan. Only the fan's own steps can match; the wrap-round one is whatever is
# left of the circle, and saying so is the point of the check.
want_step = C.azimuth_step(P["in_out_angle"], P["out_out_angle"], P["count"])
even = 360.0 / int(round(P["count"]))
want_wrap = 360.0 - (len(pitch)) * want_step
ok = (max(abs(s - want_step) for s in steps[:-1]) < 0.1 and
      abs(steps[-1] - want_wrap) < 0.1)
closes = abs(want_step - even) < 1e-6
rep.check("the outlets sit around the axis where those angles put them", ok,
          "/".join(f"{s:.3f}" for s in steps)
          + f" vs {want_step:.3f} steps and {want_wrap:.3f} closing the circle"
          + ("; the ring closes, evenly spaced" if closes else
             f"; a fan, {abs(want_wrap - want_step):.3f} deg off even"))
rep.near("the wrap-round pair is the angle that leaves", wrap,
         math.degrees(math.acos(max(-1.0, min(1.0,
             math.cos(math.radians(P["in_out_angle"])) ** 2 +
             math.sin(math.radians(P["in_out_angle"])) ** 2 *
             math.cos(math.radians(abs(steps[-1]))))))), 0.1, "{:.3f} deg")

radii = []
for f in outlets:
    p = tuple(f["c"][k] - dot(f["c"], inlet_axis) * inlet_axis[k]
              for k in range(3))
    radii.append(math.sqrt(dot(p, p)))
rep.near("the outlet mouths stand the same distance out", max(radii) - min(radii),
         0.0, 0.05, "{:.4f} mm")
rep.near("neighbouring mouths stand this far apart", C.spacing(
    P["out_arm"], P["in_out_angle"], P["out_out_angle"], P["count"]),
    2.0 * (sum(radii) / len(radii)) * math.sin(math.radians(want_step) / 2.0),
    0.3, "{:.2f} mm")

# ------------------------------------------------------------ the joints ---
# The arms as the parameters describe them: needed both to check the joint
# tapers (which uses the tube radius as a fence) and to check that the
# junction left the ground zones alone.
# The arms as the parameters describe them: needed both to check the joint
# tapers (which uses the tube radius as a fence) and to check that the
# junction left the ground zones alone.
# The ground zone has to survive the junction: the junction pulls each arm's
# wall back off the centre, and if it bit into the joint there would be no
# cone left to seat on. Recompute the trims from the parameters on file.
arms = [J.arm("in", P["inlet_joint"], "female", P["inlet_arm"], (0, 0, 1),
              P["wall"], P["tip_wall"])]
for k, u in enumerate(C.outlet_dirs(P["in_out_angle"],
                                    P["out_out_angle"],
                                    P["count"])):
    arms.append(J.arm(f"out{k+1}", P["out_joint"], "male", P["out_arm"], u,
                      P["wall"], P["tip_wall"]))
bore_j = max(a["bore"] for a in arms)
cuts, r_cap, tube_r, _ = J.plan(arms, P["wall"], P["flare"], bore_j, P["tol"])
clear = [a["len"] - cut - a["glen"] for a, cut in zip(arms, cuts)]
rep.check("the junction does not eat into any ground zone",
          min(clear) > 0.0,
          "tube left past the ground zones: "
          + "/".join(f"{c:.2f}" for c in clear) + " mm")


# Which arm a mouth belongs to has to be settled by its axis, not by the order
# the mouths came back in: the outlets are identical, so any labelling of them
# is as good as any other, and the file cannot know which one you call out1.
tube_r = J.plan(arms, P["wall"], P["flare"], max(a["bore"] for a in arms),
                P["tol"])[2]


def arm_for(face):
    axis = face["n"]                       # outward, along that arm
    return max(arms, key=lambda a: dot(a["dir"], axis))


for tag, face, spec, gender in (
        [("in", inlet, P["inlet_joint"], "female")] +
        [(arm_for(f)["tag"], f, P["out_joint"], "male") for f in outlets]):
    mine = arm_for(face)
    r_mouth, r_tip, glen = J.joint_dims(spec)
    axis = tuple(-n for n in face["n"])       # into the arm, from its mouth

    def on_taper(d):
        """That arm's own surface at depth `d`.

        A plane across one outlet also cuts its neighbours, so the radii have
        to be fenced. How close a neighbour can come is geometry, not
        guesswork: the arms all meet at the junction, so the distance from
        this arm's axis point to the neighbour's axis is a dot product, and
        that arm's wall can be no nearer than its own outside radius at its
        own station. A neighbour whose solid does not even reach this plane
        -- the inlet, well above the junction by here -- drops out.

        Each ring is a polygon, not a circle, because the cone is faceted, so
        the radius is taken as the mean of the ring. An extreme of it is
        biased in or out by the facet angle -- enough to tip a 1:10 check over
        its tolerance.
        """
        at = tuple(face["c"][k] + d * axis[k] for k in range(3))
        limit = float("inf")
        for other in arms:
            if other["tag"] == tag:
                continue
            u = other["dir"]
            along = dot(at, u)
            s = other["len"] - along
            if s < 0.0 or s > other["len"]:
                continue
            foot = tuple(along * u[k] for k in range(3))
            far = math.sqrt(dot(sub(at, foot), sub(at, foot)))
            limit = min(limit, far - J.arm_profile(
                other, s, P["wall"], P["flare"], tube_r, bore_j)[1])
        rs = sorted(r for r in K.plane_radii(verts, tris, face["c"], axis, d)
                    if r < limit - 0.5)
        if not rs:
            raise AssertionError(
                f"{tag}: no surface of its own found at depth {d:g} mm; "
                f"neighbours are within {limit:.2f} mm")
        inner_ring = [r for r in rs if r < rs[0] + 0.5]
        rest = [r for r in rs if r >= rs[0] + 0.5]
        if not rest:
            return sum(inner_ring) / len(inner_ring)
        outer_ring = [r for r in rest if r < rest[0] + 0.5]
        inner = sum(inner_ring) / len(inner_ring)
        outer = sum(outer_ring) / len(outer_ring)
        return inner if gender == "female" else outer

    fracs = (0.05, 0.3, 0.5, 0.7, 0.95)
    taps = [on_taper(f * glen) for f in fracs]
    ds = [f * glen for f in fracs]
    slope = (taps[-1] - taps[0]) / (ds[-1] - ds[0])
    icept = taps[0] - slope * ds[0]
    resid = max(abs(v - (icept + slope * d)) for v, d in zip(taps, ds))
    end0, end1 = ((r_tip, r_mouth) if gender == "male" else (r_mouth, r_tip))
    rep.check(f"{tag} joint taper is straight over {glen:g} mm", resid < 0.02,
              " ".join(f"{v:.3f}" for v in taps) + f"  (residual {resid:.4f})")
    rep.near(f"{tag} free-end radius is {end0:g} mm", icept, end0, 0.05,
             "{:.3f} mm")
    rep.near(f"{tag} inboard radius is {end1:g} mm", icept + slope * glen, end1,
             0.05, "{:.3f} mm")
    rep.near(f"{tag} taper is 1:10 on the diameter",
             abs(glen / ((icept - (icept + slope * glen)) * 2)), 10.0, 0.05)

# --------------------------------------------------------- the passage ---
# The inlet has to reach all three outlets with nothing standing in the way.
# Walk each arm's axis into the junction and out along every other arm: that
# path *is* the passage, so every sample of it should be void.
VOID_RAYS = ((0.3, 0.7, 0.65), (-0.4, 0.9, 0.2), (0.11, 0.21, -0.97))


def clear_at(p):
    """Void by a majority of rays: this part is symmetric about the model's
    own z, and the checker's default ray is horizontal, so single-ray answers
    on the axis planes are a coin toss."""
    yes = 0
    for d in VOID_RAYS:
        if K.inside(verts, tris, p, shells=shells, direction=d, boxes=boxes):
            yes += 1
            if yes >= 2:
                return False
    return True


boxes = K.shell_boxes(verts, tris, shells)

axes = {}
mouths = {"in": inlet}
for i, f in enumerate(outlets):
    mouths[f"out{i + 1}"] = f
for tag, f in mouths.items():
    axes[tag] = tuple(-n for n in f["n"])

junction = (0.0, 0.0, 0.0)
rep.check("the junction is open, not filled in", clear_at(junction),
          f"centre ({junction[0]:.3f}, {junction[1]:.3f}, {junction[2]:.3f})")

lengths = {"in": P["inlet_arm"]}
for i in range(int(P["count"])):
    lengths[f"out{i + 1}"] = P["out_arm"]

blocked = []
tags = list(mouths)
for a in tags:
    for b in tags:
        if a == b:
            continue
        for k in range(41):
            f = k / 40.0
            if f <= 0.5:
                d = min(2.0 * f * lengths[a], lengths[a] - 0.05)
                p = tuple(mouths[a]["c"][i] + d * axes[a][i] for i in range(3))
            else:
                d = max(2.0 * (1.0 - f) * lengths[b], 0.05)
                p = tuple(mouths[b]["c"][i] + d * axes[b][i] for i in range(3))
            if not clear_at(p):
                blocked.append((a, b, round(f, 2)))
rep.check("every mouth is joined to every other through the junction",
          not blocked,
          f"material on the path at {blocked[:4]}" if blocked else
          f"all {len(tags) * (len(tags) - 1)} mouth-to-mouth paths walked, "
          f"clear throughout")

# ------------------------------------------------------- printing as drawn ---
# Nothing inside a passage may need support: that is where support cannot be
# reached to take it out again. The junction's own underside does have
# down-facing surfaces -- the crotch under the outlets, where a sphere-like
# blend always does -- but that is outside the part, where support is ordinary
# and comes off. So the check is on the passages, and the outside is reported.
over = PR.overhangs(verts, tris, exposed_only=True, max_angle=ANGLE_DEG,
                    shells=K.shell_triangles(verts, tris))
inside_bore = []
for o in over:
    # What matters is the air the face holds up, so test the point just off
    # it, not the face itself. The bore widens from the mouth towards the
    # junction, so the radius to compare against is the profile's, not the
    # joint's: using the mouth radius misses the overhangs nearest the
    # junction, which are exactly the ones a shallow outlet angle creates.
    air = tuple(o["c"][k] + 0.15 * o["n"][k] for k in range(3))
    for a in arms:
        u = a["dir"]
        along = dot(air, u)
        if not (0.05 <= along <= a["len"] - 0.05):
            continue
        perp = math.sqrt(max(0.0, dot(air, air) - along * along))
        r_bore = J.arm_profile(a, a["len"] - along, P["wall"], P["flare"],
                               tube_r, bore_j)[0]
        if perp <= r_bore - 0.05:
            inside_bore.append((a["tag"], round(perp, 2), round(r_bore, 2)))
            break
rep.check(f"no passage surface needs support, at {ANGLE_DEG:g} deg from "
          f"horizontal", not inside_bore,
          f"down-facing surfaces inside a bore at {inside_bore[:4]}"
          if inside_bore else
          "every down-facing face is on the outside, where support comes off")
outside_area = sum(o["area"] for o in over)
print(f"  note: {len(over)} down-facing faces on the outside of the junction, "
      f"{outside_area:.0f} mm^2; support there is reachable, so it is not a "
      f"defect, but this part does not print support-free")

print(rep.finish())

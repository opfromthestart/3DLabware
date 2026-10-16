#!/usr/bin/env python3
"""Verify a three-way adapter 3MF: joints, angles, arm lengths, connectivity.

Re-parsed from disk and measured from the mesh. The arm directions are found
without being told where they are: each arm ends in a mouth annulus, and a
mouth's plane normal is that arm's axis. The angles between those normals are
the pairwise angles of the part.

Connectivity is checked two ways. The wall has to be there from the junction
out to each mouth: the part is a union of separate closed shells -- three arms
and the junction piece -- and nothing about that construction guarantees they
actually overlap; a gap between them would still form a mesh, still have a
plausible volume and still not print. And the *passage* has to be clear, from
each mouth through the junction and out of every other mouth, because a
three-way adapter whose bores do not meet is not a three-way adapter.

Usage: python3 verify_wye.py [wye.3mf]
"""

import math
import sys
import xml.etree.ElementTree as ET
import zipfile

import checks as K
import make_wye as W
import threelib as t

PATH = sys.argv[1] if len(sys.argv) > 1 else W.default_out()
P = K.read_params(PATH)


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


rep = K.Report(f"Three-way adapter: {PATH}")
objs, _, _ = K.load_3mf(PATH, rep)
rep.check("one object", len(objs) == 1, f"{len(objs)}")
verts, tris = objs[0]["verts"], objs[0]["tris"]
shells = K.shell_triangles(verts, tris)
boxes = K.shell_boxes(verts, tris, shells)


def tri_area(tri):
    a, b, c = verts[tri[0]], verts[tri[1]], verts[tri[2]]
    n = cross(sub(b, a), sub(c, a))
    return 0.5 * math.sqrt(dot(n, n))


def section(origin, axis, d):
    """Radii where the mesh crosses the plane at depth d along `axis`."""
    rs = []
    for tri in tris:
        for i, j in ((0, 1), (1, 2), (2, 0)):
            a, b = verts[tri[i]], verts[tri[j]]
            da = dot(sub(a, origin), axis) - d
            db = dot(sub(b, origin), axis) - d
            if (da < 0.0) == (db < 0.0):
                continue
            u = da / (da - db)
            q = tuple(a[k] + u * (b[k] - a[k]) for k in range(3))
            radial = sub(sub(q, origin), tuple(d * axis[k] for k in range(3)))
            rs.append(math.sqrt(dot(radial, radial)))
    return rs


# ---------------------------------------------------------------- mouths ---
# A mouth is an annular face at the end of an arm. A cap where an arm's wall
# stops short of the junction is one too, and a bigger one -- so the two are
# told apart by whether material lies behind the face: step off a mouth and you
# are in open air, step off an internal cap and you are inside the junction.
faces = K.annuli(verts, tris)
faces.sort(key=lambda f: (f["buried"], -f["n_tri"], -f["area"]))
rep.check("one mouth annulus per arm",
          len([f for f in faces if not f["buried"]]) == 3,
          f"{len(faces)} annular faces, "
          f"{len([f for f in faces if not f['buried']])} of them on the outside")
mouths = [f for f in faces if not f["buried"]][:3]


def end_bore(spec, gender):
    _, r_tip, _ = W.joint_dims(spec)
    return r_tip if gender == "female" else r_tip - P["tip_wall"]


def signature(spec, gender, tip_wall):
    r_mouth, r_tip, _ = W.joint_dims(spec)
    if gender == "male":
        return r_tip, end_bore(spec, gender)
    return r_mouth + P["wall"], r_mouth


declared = {tag: (P["joint_" + tag], P["gender_" + tag]) for tag in "123"}
best, best_cost = None, None
for perm in ((0, 1, 2), (0, 2, 1), (1, 0, 2), (1, 2, 0), (2, 0, 1), (2, 1, 0)):
    cost = 0.0
    for tag, k in zip("123", perm):
        spec, gender = declared[tag]
        want = signature(spec, gender, P["tip_wall"])
        cost += abs(mouths[k]["outer"] - want[0]) + abs(mouths[k]["inner"] - want[1])
    if best_cost is None or cost < best_cost:
        best, best_cost = perm, cost
paired = {tag: mouths[best[i]] for i, tag in enumerate("123")}
rep.check("each arm's mouth matches its declared joint", best_cost < 0.2,
          f"total mismatch {best_cost:.4f} mm across the three arms")

# ---------------------------------------------------------------- angles ---
axes = {}
for tag in "123":
    fin = paired[tag]
    axes[tag] = (-fin["n"][0], -fin["n"][1], -fin["n"][2])

# With arms that differ, each mouth matches its own joint and the labels mean
# something. With three identical arms there is nothing to match on, the
# labelling is arbitrary, and demanding a particular one would fail whenever a
# coin flip went the other way. Compare the *set* of angles instead: that is
# what the part actually has to satisfy, and it holds either way.
sigs = [signature(declared[t][0], declared[t][1], P["tip_wall"]) for t in "123"]
distinct = len({(round(a, 3), round(b, 3)) for a, b in sigs}) == 3
measured = [math.degrees(math.acos(max(-1.0, min(1.0, dot(axes[a], axes[b])))))
            for a, b in (("1", "2"), ("1", "3"), ("2", "3"))]
wanted = [P["angle_12"], P["angle_13"], P["angle_23"]]
worst = max(abs(g - w) for g, w in zip(sorted(measured), sorted(wanted)))
rep.check("the three pairwise angles are the requested set", worst < 0.05,
          "measured " + "/".join(f"{v:.4f}" for v in sorted(measured))
          + "  requested " + "/".join(f"{v:g}" for v in sorted(wanted)))
if distinct:
    for a, b in (("1", "2"), ("1", "3"), ("2", "3")):
        got = math.degrees(math.acos(max(-1.0, min(1.0, dot(axes[a], axes[b])))))
        rep.near(f"angle {a}{b}", got, P["angle_" + a + b], 0.05, "{:.4f} deg")
else:
    print("  note: all three arms are identical, so which is 1 and which is 3 "
          "is arbitrary;")
    print("  note: only the set of pairwise angles is meaningful for this part")

# The three axes must be concurrent: that point is the junction.
a1, d1 = paired["1"]["c"], axes["1"]
a2, d2 = paired["2"]["c"], axes["2"]
nrm = cross(d1, d2)
m = math.sqrt(dot(nrm, nrm))
if m < 1e-9:
    rep.check("arm axes are not parallel", False, "arms 1 and 2 are collinear")
    junction = paired["1"]["c"]
else:
    w = sub(a2, a1)
    s = dot(cross(w, d2), nrm) / (m * m)
    junction = tuple(a1[k] + s * d1[k] for k in range(3))
    off3 = math.sqrt(dot(sub(paired["3"]["c"], junction),
                         sub(paired["3"]["c"], junction)))
    # arm 3's axis passes through the junction if its mouth lies on it: the
    # perpendicular distance is what matters, so measure that.
    v = sub(paired["3"]["c"], junction)
    along = dot(v, axes["3"])
    perp = math.sqrt(max(0.0, dot(v, v) - along * along))
    rep.check("the three arm axes meet at one point", perp < 0.05,
              f"arm 3 misses the junction by {perp:.4f} mm")

lengths = sorted(math.sqrt(dot(sub(paired[t]["c"], junction),
                                 sub(paired[t]["c"], junction))) for t in "123")
want_len = sorted(P["arm_" + t] for t in "123")
rep.check("arm lengths (mouth to junction) match the request",
          max(abs(g - w) for g, w in zip(lengths, want_len)) < 0.05,
          "/".join(f"{v:.3f}" for v in lengths) + "  vs  "
          + "/".join(f"{v:g}" for v in want_len))

# ------------------------------------------------------------ the joints ---
for tag in "123":
    spec, gender = declared[tag]
    fin = paired[tag]
    r_mouth, r_tip, glen = W.joint_dims(spec)
    axis = axes[tag]
    pick = max if gender == "male" else min
    fracs = (0.05, 0.3, 0.5, 0.7, 0.95)
    taps = [pick(section(fin["c"], axis, f * glen)) for f in fracs]
    ds = [f * glen for f in fracs]
    slope = (taps[-1] - taps[0]) / (ds[-1] - ds[0])
    icept = taps[0] - slope * ds[0]
    resid = max(abs(v - (icept + slope * d)) for v, d in zip(taps, ds))
    end0, end1 = ((r_tip, r_mouth) if gender == "male" else (r_mouth, r_tip))
    rep.check(f"arm {tag} joint taper is straight over {glen:g} mm", resid < 0.02,
              " ".join(f"{v:.3f}" for v in taps) + f"  (residual {resid:.4f})")
    rep.near(f"arm {tag} free-end radius is {end0:g} mm", icept, end0, 0.05,
             "{:.3f} mm")
    rep.near(f"arm {tag} inboard radius is {end1:g} mm", icept + slope * glen,
             end1, 0.05, "{:.3f} mm")
    measured = abs(glen / ((icept - (icept + slope * glen)) * 2))
    rep.near(f"arm {tag} taper is 1:10 on the diameter", measured, 10.0, 0.05)

    # Bore just past the ground zone, where a male tip keeps its wall.
    if gender == "male":
        cut = section(fin["c"], axis, 0.5)
        rep.near(f"arm {tag} tip keeps a {P['tip_wall']:g} mm wall",
                 max(cut) - min(cut), P["tip_wall"], 0.06, "{:.3f} mm")

# --------------------------------------------------------- connectivity ---
# Walk each arm from the junction out to its mouth. The part is separate
# shells that were built to overlap; if any two of them missed, the walk would
# step outside the material and this is where that shows up.
# Sample across the arm at each depth and require *some* material. Walking the
# axis itself would find nothing: the axis is the bore, which is meant to be
# empty. What has to be continuous is the wall.
holes = []
for tag in "123":
    fin = paired[tag]
    axis = axes[tag]
    perp = unit(cross(axis, (0.0, 0.0, 1.0) if abs(axis[2]) < 0.9
                      else (1.0, 0.0, 0.0)))
    for i in range(61):
        # Start just inside: at d=0 every sample would lie exactly in the
        # mouth's own plane, on the boundary, where the parity test is
        # undefined and arbitrarily reports either answer.
        d = 0.05 + (P["arm_" + tag] - 0.1) * i / 60.0
        # axis points inward, so depth is added to the mouth's centroid.
        base = tuple(fin["c"][k] + d * axis[k] for k in range(3))
        if not any(K.inside(verts, tris, tuple(
                base[k] + r * perp[k] for k in range(3)),
                shells=shells, boxes=boxes)
                for r in [0.1 * j for j in range(1, 140)]):
            holes.append((tag, round(d, 2)))
rep.check("each arm is solid material from the junction to its mouth", not holes,
          f"no material at {holes[:5]}" if holes else
          "3 arms walked from the junction to the mouth, wall present throughout")

# ------------------------------------------------------ the junction ---
# The junction is meant to be *open*: the three bores have to run into one
# another with nothing left standing in the way. Glue three tubes together at
# the junction point instead and each arm's wall plugs the other two bores --
# two of the three pairs of mouths come out shut and the third is pinched.
#
# The junction centre is tested along a ray that is not parallel to any arm.
# This model's arm axes lie in one plane, and the checker's default ray
# direction is horizontal, which is that same plane: there the parity test
# straddles the surface and returns either answer.
VOID_RAY = (0.3, 0.7, 0.65)
rep.check("the junction is open, not filled in",
          not K.inside(verts, tris, junction, direction=VOID_RAY,
                       shells=shells, boxes=boxes),
          f"centre ({junction[0]:.3f}, {junction[1]:.3f}, {junction[2]:.3f})")

# From each mouth, through the junction, out of every other mouth. Every
# sample of that path is meant to be void, because that path *is* the passage
# the part exists for.
def depth_point(tag, d):
    return tuple(paired[tag]["c"][k] + d * axes[tag][k] for k in range(3))


blocked = []
for a in "123":
    for b in "123":
        if a == b:
            continue
        for i in range(41):
            f = i / 40.0
            if f <= 0.5:
                d = 2.0 * f * P["arm_" + a]
                p = depth_point(a, min(d, P["arm_" + a] - 0.05))
            else:
                d = 2.0 * (1.0 - f) * P["arm_" + b]
                p = depth_point(b, max(d, 0.05))
            if K.inside(verts, tris, p, shells=shells, boxes=boxes):
                blocked.append((a, b, round(f, 2)))
rep.check("each pair of mouths is joined through the junction", not blocked,
          f"material on the path at {blocked[:4]}" if blocked else
          "all 6 mouth-to-mouth paths walked, clear throughout")

print(rep.finish())

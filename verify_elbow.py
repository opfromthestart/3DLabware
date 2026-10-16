#!/usr/bin/env python3
"""Verify an elbow 3MF: joint sizes, taper, angle and arm lengths.

The file is re-parsed from disk and everything is measured from the mesh, so a
build that silently produces the wrong taper or the wrong angle fails here.

The two ground joints are found without being told where they are: they are the
two largest *planar* faces on the part, because every other surface is either
part of a cylinder or part of the swept bend. Their planes give the arm axes
and so the included angle; their centres give the arm lengths.

Usage: python3 verify_elbow.py [elbow.3mf]
"""

import math
import sys
import xml.etree.ElementTree as ET
import zipfile

import checks as K
import make_elbow as E
import threelib as t

PATH = sys.argv[1] if len(sys.argv) > 1 else E.default_out()
CORE = "{http://schemas.microsoft.com/3dmanufacturing/core/2015/02}"


def read_params(path):
    """The parameters the writer stamped into the 3MF.

    Check against these, not against the module's globals: the point of the
    stamp is that a file carries its own account of itself, so verifying a
    75 degree elbow built from a module defaulting to 105 has to compare
    against 75. Every parameter is a string in the XML; numbers are converted
    back, joint designations like "14/20" are left alone.
    """
    root = ET.fromstring(zipfile.ZipFile(path).read("3D/3dmodel.model"))
    out = {}
    for m in root.findall(f"{CORE}metadata"):
        name = m.get("name", "")
        if not name.startswith("Parameter."):
            continue
        val = (m.text or "").strip()
        try:
            out[name[10:]] = float(val)
            continue
        except ValueError:
            pass
        # Non-numeric parameters are written as their repr, quotes and all.
        if len(val) >= 2 and val[0] == val[-1] and val[0] in "'\"":
            val = val[1:-1]
        out[name[10:]] = val
    return out


P = read_params(PATH)
R_MOUTH_A, R_TIP_A, LEN_A = E.joint_dims(P["joint_a"])
R_MOUTH_B, R_TIP_B, LEN_B = E.joint_dims(P["joint_b"])
# The bore is a profile, not a number: a female end's socket finishes at its
# tip radius and the pipe carries on at that size, a male end's tip is that
# same radius on the outside so its bore is smaller by the tip wall.
def end_bore(spec, gender):
    _, r_tip, _ = E.joint_dims(spec)
    return r_tip if gender == "female" else r_tip - P["tip_wall"]


BORE_A = end_bore(P["joint_a"], P["gender_a"])
BORE_B = end_bore(P["joint_b"], P["gender_b"])
DELTA_R = math.radians(180.0 - P["angle"])
ARC = P["bend_r"] * DELTA_R
TOTAL = P["arm_a"] + ARC + P["arm_b"]
FLARE = P["flare"]
TUBE_R = max(R_TIP_A + P["wall"] if P["gender_a"] == "female" else R_MOUTH_A,
             R_TIP_B + P["wall"] if P["gender_b"] == "female" else R_MOUTH_B,
             max(BORE_A, BORE_B) + P["wall"])


def bore_at(depth):
    """Bore radius at arc length `depth` from joint A's mouth."""
    if depth <= LEN_A:
        return (R_MOUTH_A - E.TAPER * depth
                if P["gender_a"] == "female" else BORE_A)
    if depth >= TOTAL - LEN_B:
        d = depth - (TOTAL - LEN_B)
        return (R_TIP_B + E.TAPER * d
                if P["gender_b"] == "female" else BORE_B)
    u = (depth - LEN_A) / ((TOTAL - LEN_B) - LEN_A)
    return BORE_A + (BORE_B - BORE_A) * u

rep = K.Report(f"Elbow geometry: {PATH}")
objs, _, _ = K.load_3mf(PATH, rep)
rep.check("one object", len(objs) == 1, f"{len(objs)}")
verts, tris = objs[0]["verts"], objs[0]["tris"]
sh = K.shells(verts, tris)
rep.check("part is a single connected solid", len(sh) == 1, f"{len(sh)} shells")


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


def planar_faces(verts, tris, tol=1e-3):
    """{(normal, offset): [tri, ...]} for every coplanar group of triangles."""
    groups = {}
    for i, tri in enumerate(tris):
        a, b, c = verts[tri[0]], verts[tri[1]], verts[tri[2]]
        n = cross(sub(b, a), sub(c, a))
        m = math.sqrt(dot(n, n))
        if m < 1e-12:
            continue
        n = (n[0] / m, n[1] / m, n[2] / m)
        key = (round(n[0], 3), round(n[1], 3), round(n[2], 3), round(dot(n, a), 3))
        groups.setdefault(key, []).append(i)
    return groups


def tri_area(verts, tri):
    a, b, c = verts[tri[0]], verts[tri[1]], verts[tri[2]]
    n = cross(sub(b, a), sub(c, a))
    return 0.5 * math.sqrt(dot(n, n))


def section(verts, tris, origin, axis, d):
    """Radii where the mesh crosses the plane at depth d along `axis`.

    Sampling vertices cannot measure a taper: the sweep puts stations only at
    the ends of each straight run, so a ground cone has vertices at its mouth
    and at its tip and nowhere between. Cutting the mesh finds the surface
    wherever it actually is.
    """
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


faces = planar_faces(verts, tris)


def face_info(key, idxs):
    """Describe a planar group, and say whether it is an annulus.

    Area alone does not identify a mouth: a straight tube's facet is a single
    long planar quad, 175 mm by half a millimetre, which beats the mouth
    annulus. What separates them is shape -- an annulus has all its vertices at
    exactly two radii from its centre, a facet does not.
    """
    n = key[:3]
    pts = sorted({v for i in idxs for v in tris[i]})
    ring = [verts[v] for v in pts]
    c = tuple(sum(p[k] for p in ring) / len(ring) for k in range(3))
    rad = sorted(math.sqrt(dot(sub(p, c), sub(p, c))) for p in ring)
    gaps = [(rad[i + 1] - rad[i], i) for i in range(len(rad) - 1)]
    gap, at = max(gaps) if gaps else (0.0, 0)
    annulus = (gap > 0.5 and at + 1 >= 8 and len(rad) - at - 1 >= 8)
    return {"n": n, "c": c, "outer": rad[-1], "inner": rad[0],
            "area": sum(tri_area(verts, tris[i]) for i in idxs),
            "n_tri": len(idxs), "annulus": annulus}


cands = [face_info(k, v) for k, v in faces.items() if len(v) >= 8]
cands = [c for c in cands if c["annulus"]]
# Rank by triangle count first: a mouth is a full ring of them, a stray flat
# face that happens to be annular has few.
cands.sort(key=lambda c: (-c["n_tri"], -c["area"]))
rep.check("two annular faces found", len(cands) >= 2, f"{len(cands)} annuli")

ends = cands[:2]

ang = math.degrees(math.acos(max(-1.0, min(1.0, dot(ends[0]["n"], ends[1]["n"])))))
rep.near("included angle between the arms", ang, P["angle"], 0.05, "{:.4f} deg")

# Which end is which. Do not assume the arms came out in parameter order:
# with two joints of the same size the only thing telling them apart is the
# shape of the mouth annulus, so match on that.
def signature(spec, gender, bore_r, wall):
    """(outer, inner) radii of a joint's mouth annulus.

    A male joint's mouth runs from the bore out to the joint radius -- the
    ground cone is its outside. A female's runs from the joint radius out to
    the joint radius plus the wall, because its mouth is a hole of the joint's
    own size.
    """
    r_mouth, r_tip, _ = E.joint_dims(spec)
    if gender == "male":
        # A male's free end is its narrow tip: the ground cone is its outside.
        return r_tip, bore_r
    # A female's free end is its wide mouth, and that mouth is a hole of the
    # joint's own size, so the wall shows around it.
    return r_mouth + wall, r_mouth


declared = {"A": (P["joint_a"], P["gender_a"]),
            "B": (P["joint_b"], P["gender_b"])}
best, best_cost = (0, 1), None
for perm in ((0, 1), (1, 0)):
    cost = 0.0
    for tag, k in zip(("A", "B"), perm):
        spec, gender = declared[tag]
        want = signature(spec, gender, end_bore(spec, gender), P["wall"])
        cost += abs(ends[k]["outer"] - want[0]) + abs(ends[k]["inner"] - want[1])
    if best_cost is None or cost < best_cost:
        best, best_cost = perm, cost
paired = {"A": ends[best[0]], "B": ends[best[1]]}

for tag in ("A", "B"):
    spec, gender = declared[tag]
    fin = paired[tag]
    r_mouth, r_tip, ln = E.joint_dims(spec)
    want_o, want_i = signature(spec, gender, end_bore(spec, gender), P["wall"])
    rep.check(f"joint {tag} is the {gender} one", best_cost is not None,
              f"mouth r {fin['outer']:.2f} outer / {fin['inner']:.2f} inner "
              f"matches {gender} {spec}")
    rep.near(f"joint {tag} mouth outer radius is {spec}", fin["outer"],
             want_o, 0.05, "{:.3f} mm")

    # Depth along the arm from the mouth's plane, radius from the arm's axis.
    # A male joint is the outer surface, a female the socket bore.
    axis = (-fin["n"][0], -fin["n"][1], -fin["n"][2])
    pick = max if gender == "male" else min
    # Keep the cuts off the mouth and shoulder planes: a cut coplanar with a
    # face has both edge endpoints at zero depth and the crossing test, which
    # is strict, skips every one of them. Fit the taper instead of sampling its
    # ends, which also tests that it is straight rather than just correct at
    # two points.
    fracs = (0.05, 0.3, 0.5, 0.7, 0.95)
    taps = [pick(section(verts, tris, fin["c"], axis, f * ln)) for f in fracs]
    ds = [f * ln for f in fracs]
    slope = (taps[-1] - taps[0]) / (ds[-1] - ds[0])
    icept = taps[0] - slope * ds[0]
    resid = max(abs(t - (icept + slope * d)) for t, d in zip(taps, ds))
    rep.check(f"joint {tag} taper is straight over its {ln:g} mm",
              resid < 0.02,
              " ".join(f"{v:.3f}" for v in taps) + f"  (residual {resid:.4f} mm)")
    end0, end1 = ((r_tip, r_mouth) if gender == "male" else (r_mouth, r_tip))
    rep.near(f"joint {tag} free-end radius is {end0:g} mm", icept, end0, 0.05,
             "{:.3f} mm")
    rep.near(f"joint {tag} inboard radius is {end1:g} mm",
             icept + slope * ln, end1, 0.05, "{:.3f} mm")
    # Magnitude only: a male's radius grows with depth and a female's shrinks,
    # so the signed ratio flips sign between them while the taper is the same.
    measured = abs(ln / ((icept - (icept + slope * ln)) * 2))
    rep.near(f"joint {tag} taper is 1:10 on the diameter", measured, 10.0, 0.05)


# Arm lengths: fit each arm's axis from the straight run, intersect the two
# axes to get the corner, and compare the mouth-to-corner distance with the
# geometry of a circular bend tangent to both arms.
delta = math.radians(180.0 - P["angle"])
for tag in ("A", "B"):
    fin, arm = paired[tag], E.PARAMS["arm_" + tag.lower()]
    axis = (-fin["n"][0], -fin["n"][1], -fin["n"][2])
    far = 0.4 * arm
    pts = []
    for v in verts:
        d = dot(sub(v, fin["c"]), axis)
        if 2.0 <= d <= far:
            pts.append(dot(v, axis))
    # Fit only needs the axis offset, which the mouth plane already fixes, so
    # measure the centroid drift of slabs perpendicular to the fitted axis.
    slabs = {}
    for v in verts:
        d = dot(sub(v, fin["c"]), axis)
        if 2.0 <= d <= far:
            slabs.setdefault(round(d), []).append(v)
    cents = []
    for d in sorted(slabs):
        g = slabs[d]
        cents.append(tuple(sum(p[k] for p in g) / len(g) for k in range(3)))
    if len(cents) >= 2:
        d0 = unit(sub(cents[-1], cents[0]))
        ends_axis = (fin["c"], d0)
    else:
        ends_axis = (fin["c"], axis)
    fin["axis"] = ends_axis

(a1, d1), (a2, d2) = paired["A"]["axis"], paired["B"]["axis"]
# Closest approach of the two axis lines is their intersection point (they are
# coplanar): solve a1 + s*d1 = a2 + u*d2 in the plane.
nrm = cross(d1, d2)
m = math.sqrt(dot(nrm, nrm))
if m < 1e-9:
    rep.check("arms are not parallel", False, "fitted axes are parallel")
else:
    w = sub(a2, a1)
    s = dot(cross(w, d2), nrm) / (m * m)
    corner = tuple(a1[k] + s * d1[k] for k in range(3))
    reach = P["bend_r"] * math.tan(delta / 2.0)
    for tag in ("A", "B"):
        fin, arm = paired[tag], E.PARAMS["arm_" + tag.lower()]
        got = math.sqrt(dot(sub(fin["c"], corner), sub(fin["c"], corner)))
        rep.near(f"arm {tag} length (mouth to turn)", got - reach, arm, 0.6,
                 "{:.3f} mm")

# Bore and wall through the plain tube. Measured by cutting the mesh, and
# compared against the bore profile at that arc length rather than one number,
# because the bore is meant to change along the elbow.
axis = paired["A"]["axis"][1]
c0 = paired["A"]["c"]
mid = (LEN_A + P["flare"] + P["arm_a"]) / 2.0
radii = section(verts, tris, c0, axis, mid)
inner, outer = min(radii), max(radii)
rep.near(f"bore diameter at {mid:g} mm along", inner * 2, bore_at(mid) * 2,
         0.06, "{:.3f} mm")
rep.near("outer diameter through the arm", outer * 2, TUBE_R * 2, 0.06,
         "{:.3f} mm")
# `wall` is the wall at the female end. It is meant to be a little thicker
# further along, because the bore shrinks there while the outside does not.
rep.check("wall is at least the specified thickness all the way",
          outer - inner >= P["wall"] - 0.01,
          f"{outer - inner:.3f} mm at {mid:g} mm along "
          f"(nominal {P['wall']:g}, bore {inner * 2:.3f})")

# The female socket must run straight into the pipe. A step here is a ledge for
# solids to settle on and something for an incoming joint to bottom out on.
if "female" in (P["gender_a"], P["gender_b"]):
    tag = "A" if P["gender_a"] == "female" else "B"
    fin = paired[tag]
    depth = LEN_A if tag == "A" else LEN_B
    r_axis = (-fin["n"][0], -fin["n"][1], -fin["n"][2])
    inside = section(verts, tris, fin["c"], r_axis, depth - 0.5)
    just_in = section(verts, tris, fin["c"], r_axis, depth + 0.5)
    step = abs(min(inside) - min(just_in))
    rep.check(f"female socket runs into the pipe with no shoulder", step < 0.06,
              f"inner r {min(inside):.3f} at {depth - 0.5:g} mm vs "
              f"{min(just_in):.3f} at {depth + 0.5:g} mm (step {step:.3f} mm)")

# And at a male end the tip keeps a real wall, which is the whole reason the
# bore is allowed to shrink in the first place.
for tag in ("A", "B"):
    gender = P["gender_" + tag.lower()]
    if gender != "male":
        continue
    fin = paired[tag]
    r_axis = (-fin["n"][0], -fin["n"][1], -fin["n"][2])
    ln = LEN_A if tag == "A" else LEN_B
    # At the tip itself, not part way along the cone: the wall thickens as the
    # ground cone widens, and it is the tip that has to survive.
    cut = section(verts, tris, fin["c"], r_axis, 0.5)
    got = end_bore(P["joint_" + tag.lower()], gender)
    rep.near(f"male joint {tag} tip keeps a {P['tip_wall']:g} mm wall",
             max(cut) - min(cut), P["tip_wall"], 0.06, "{:.3f} mm")
    rep.near(f"male joint {tag} tip bore is {got * 2:g} mm",
             min(cut) * 2, got * 2, 0.06, "{:.3f} mm")

# Every vertex must sit on the nominal centreline within the tube's radii. This
# is the check that ties arms, bend and radii together: a wrong bend radius, a
# wrong arm length or a wrong taper all show up as outliers here.
def centreline_pts(n=400):
    total = P["arm_a"] + P["bend_r"] * delta + P["arm_b"]
    return [E._centreline(total * k / n, P["arm_a"], P["bend_r"],
                          delta)[0] for k in range(n + 1)]


cl = centreline_pts()
worst_rad = 0.0
for v in verts:
    best = min(math.sqrt(dot(sub(v, p), sub(v, p))) for p in cl)
    worst_rad = max(worst_rad, best)
lo, hi = min(BORE_A, BORE_B), max(TUBE_R, ends[0]["outer"],
                                    ends[1]["outer"]) + P["wall"]
lo *= 0.98
rep.check("every vertex lies on the nominal centreline path",
          lo - 0.1 <= worst_rad <= hi + 0.1,
          f"farthest {worst_rad:.2f} mm from the centreline "
          f"(tube spans {lo:.2f}..{hi:.2f})")

print(rep.finish())

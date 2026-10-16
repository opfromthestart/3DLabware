"""checks -- verification helpers for scripted 3D models.

The rule these enforce: **verify the artifact, not the generator.** Every
loader here re-reads the file that was written, so a bug in the build
script cannot hide behind itself.
"""

from __future__ import annotations

import math
import xml.etree.ElementTree as ET
import zipfile

NS = {"c": "http://schemas.microsoft.com/3dmanufacturing/core/2015/02"}
OPC_PARTS = ("[Content_Types].xml", "_rels/.rels", "3D/3dmodel.model")


# ---------------------------------------------------------------------------
# Reporting
# ---------------------------------------------------------------------------

def cli(argv, default_path):
    """(path, angle) from a verifier's arguments.

    A verifier takes an optional 3MF path and an optional overhang angle. The
    angle is the shallowest face that still prints without support, which is a
    property of the printer and not of the part, so it is settable per run:
    `--angle 20` for a machine that copes with shallow overhangs, or the
    OVERHANG_ANGLE environment variable for a whole session. `None` means use
    the library's own default.
    """
    path, angle = default_path, None
    i = 0
    while i < len(argv):
        a = argv[i]
        if a in ("--angle", "-a"):
            i += 1
            if i >= len(argv):
                raise SystemExit("--angle needs a number of degrees")
            angle = float(argv[i])
        elif a.startswith("--angle="):
            angle = float(a.split("=", 1)[1])
        elif a.startswith("-"):
            raise SystemExit(f"unknown option {a}")
        else:
            path = a
        i += 1
    return path, angle


class Report:
    def __init__(self, title):
        self.title = title
        self.rows = []
        print(f"== {title} ==")

    def check(self, label, ok, detail=""):
        ok = bool(ok)
        self.rows.append((label, ok, detail))
        print(f"  [{'PASS' if ok else 'FAIL'}] {label}"
              + (f" -- {detail}" if detail else ""))
        return ok

    def near(self, label, got, want, tol=1e-6, fmt="{:.6g}"):
        return self.check(label, abs(got - want) <= tol,
                          f"{fmt.format(got)} vs {fmt.format(want)} (tol {tol:g})")

    def finish(self):
        bad = [r[0] for r in self.rows if not r[1]]
        print("\n" + ("ALL CHECKS PASSED" if not bad
                      else f"{len(bad)} FAILURE(S): {bad}"))
        return not bad


# ---------------------------------------------------------------------------
# Loading (independent of the writing code)
# ---------------------------------------------------------------------------

def load_3mf(path, report=None):
    """Parse a 3MF from scratch. Returns (objects, build_ids, unit).

    Each object is a dict: {name, id, verts, tris}.
    """
    zf = zipfile.ZipFile(path)
    if report:
        report.check("zip integrity", zf.testzip() is None)
        names = zf.namelist()
        for part in OPC_PARTS:
            report.check(f"contains {part}", part in names)
    for part in OPC_PARTS:
        try:
            ET.fromstring(zf.read(part))
            if report:
                report.check(f"{part} parses", True)
        except ET.ParseError as exc:
            if report:
                report.check(f"{part} parses", False, str(exc))

    root = ET.fromstring(zf.read("3D/3dmodel.model"))
    objects = []
    for obj in root.findall(".//c:resources/c:object", NS):
        mesh = obj.find("c:mesh", NS)
        verts = [(float(v.get("x")), float(v.get("y")), float(v.get("z")))
                 for v in mesh.findall("c:vertices/c:vertex", NS)]
        tris = [(int(t.get("v1")), int(t.get("v2")), int(t.get("v3")))
                for t in mesh.findall("c:triangles/c:triangle", NS)]
        objects.append({"name": obj.get("name"), "id": obj.get("id"),
                        "verts": verts, "tris": tris})
    build_ids = [i.get("objectid") for i in root.findall(".//c:build/c:item", NS)]
    if report:
        report.check("unit is millimetre", root.get("unit") == "millimeter",
                     root.get("unit"))
        report.check("build list references every object",
                     sorted(build_ids) == sorted(o["id"] for o in objects),
                     f"build={build_ids}")
    return objects, build_ids, root.get("unit")


# ---------------------------------------------------------------------------
# Mesh analysis
# ---------------------------------------------------------------------------

def signed_volume(verts, tris):
    total = 0.0
    for a, b, c in tris:
        ax, ay, az = verts[a]
        bx, by, bz = verts[b]
        cx, cy, cz = verts[c]
        total += (ax * (by * cz - bz * cy)
                  - ay * (bx * cz - bz * cx)
                  + az * (bx * cy - by * cx))
    return total / 6.0


def shell_triangles(verts, tris):
    """Triangle index lists grouped into connected shells.

    A model here is a *soup*: shapes are combined by overlapping them, so one
    physical solid is normally several closed shells sharing no vertices. That
    matters for any point-in-solid test, which is only meaningful per shell.
    """
    parent = list(range(len(verts)))

    def find(a):
        while parent[a] != a:
            parent[a] = parent[parent[a]]
            a = parent[a]
        return a

    for a, b, c in tris:
        for x, y in ((a, b), (b, c), (c, a)):
            rx, ry = find(x), find(y)
            if rx != ry:
                parent[rx] = ry

    groups = {}
    for i, tri in enumerate(tris):
        groups.setdefault(find(tri[0]), []).append(i)
    return list(groups.values())


def shell_boxes(verts, tris, shells=None):
    """Axis-aligned bounding box per shell, for pruning point tests.

    A part built by overlapping shapes is several shells, and a point inside
    one of them is still tested against all the others: without this, a dense
    walk through a staged column costs one ray test per staging per point.
    """
    if shells is None:
        shells = shell_triangles(verts, tris)
    out = []
    for indices in shells:
        lo = [float("inf")] * 3
        hi = [float("-inf")] * 3
        for i in indices:
            for j in tris[i]:
                v = verts[j]
                for k in range(3):
                    if v[k] < lo[k]:
                        lo[k] = v[k]
                    if v[k] > hi[k]:
                        hi[k] = v[k]
        out.append((lo, hi))
    return out


def _inside_one_shell(verts, tris, indices, p, direction):
    """Odd/even ray crossings against a single closed shell."""
    dx, dy, dz = direction
    n = math.sqrt(dx * dx + dy * dy + dz * dz)
    dx, dy, dz = dx / n, dy / n, dz / n
    count = 0
    for i in indices:
        tri = tris[i]
        v0, v1, v2 = verts[tri[0]], verts[tri[1]], verts[tri[2]]
        e1 = tuple(v1[k] - v0[k] for k in range(3))
        e2 = tuple(v2[k] - v0[k] for k in range(3))
        px = dy * e2[2] - dz * e2[1]
        py = dz * e2[0] - dx * e2[2]
        pz = dx * e2[1] - dy * e2[0]
        det = e1[0] * px + e1[1] * py + e1[2] * pz
        if abs(det) < 1e-12:
            continue
        inv = 1.0 / det
        tx, ty, tz = p[0] - v0[0], p[1] - v0[1], p[2] - v0[2]
        u = (tx * px + ty * py + tz * pz) * inv
        if u < 0.0 or u > 1.0:
            continue
        qx = ty * e1[2] - tz * e1[1]
        qy = tz * e1[0] - tx * e1[2]
        qz = tx * e1[1] - ty * e1[0]
        v = (dx * qx + dy * qy + dz * qz) * inv
        if v < 0.0 or u + v > 1.0:
            continue
        t = (e2[0] * qx + e2[1] * qy + e2[2] * qz) * inv
        if t > 1e-9:
            count += 1
    return count % 2 == 1


def inside(verts, tris, p, shells=None,
           direction=(1.0, 0.7548776662466927, 0.0), boxes=None):
    """Is point `p` inside the solid?

    Union semantics: true if inside *any* closed shell. Counting parity across
    the whole triangle soup instead computes the symmetric difference, so a
    point where two shapes overlap -- the exact region a union is meant to
    fill -- comes out as "outside". That silently reports a solid part as
    hollow wherever it was built by overlapping two shells, which is how every
    part here is built.

    The ray direction is deliberately irrational-looking. These solids are
    full of axis-aligned edges and vertices, and a ray along an axis strikes
    them exactly: such a crossing is counted twice by the two triangles
    sharing the edge and cancels. A skewed direction misses every edge it
    does not truly cross.

    Pass `shells` (from `shell_triangle`/`shell_triangles`) and `boxes` (from
    `shell_boxes`) when testing many points: deriving them per call costs a
    union-find over every vertex and a pass over every triangle, and a shell
    whose box misses `p` cannot contain it.
    """
    if shells is None:
        shells = shell_triangles(verts, tris)
    if boxes is None:
        boxes = shell_boxes(verts, tris, shells)
    for s, (lo, hi) in zip(shells, boxes):
        if (p[0] < lo[0] or p[0] > hi[0] or p[1] < lo[1] or p[1] > hi[1]
                or p[2] < lo[2] or p[2] > hi[2]):
            continue
        if _inside_one_shell(verts, tris, s, p, direction):
            return True
    return False


def radial_profile(verts, tris, z, ang=0.0, shells=None):
    """Material intervals along the ray outward from the axis at height z.

    Returns [(r0, r1), ...] sorted, in mm. Each crossing is decided by a real
    point-in-solid test at the interval midpoint rather than by parity, so
    coincident surfaces -- two parts touching on a shared face, or a bar
    flush with a wall -- cannot produce a wrong answer.

    Pass `shells` when taking many profiles of the same mesh; deriving it per
    call costs a union-find over every vertex, which dominates on a model with
    a few thousand bars.
    """
    if shells is None:
        shells = shell_triangles(verts, tris)
    dx, dy = math.cos(ang), math.sin(ang)
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
        den = dx * ey - dy * ex
        if abs(den) < 1e-12:
            continue
        t = (x0 * ey - y0 * ex) / den
        u = (x0 * dy - y0 * dx) / den
        if -1e-9 <= u <= 1 + 1e-9 and t > 1e-7:
            hits.append(t)
    hits.sort()
    merged = []
    for h in hits:
        if not merged or h - merged[-1] > 1e-6:
            merged.append(h)

    out, prev_r = [], 0.0
    for r in merged:
        rm = 0.5 * (prev_r + r)
        if r - prev_r > 1e-6 and inside(verts, tris, (rm * dx, rm * dy, z),
                                         shells=shells):
            out.append((prev_r, r))
        prev_r = r
    return out


def read_params(path):
    """The parameter set a part was written with, read back from its 3MF.

    Every part here stamps its parameters into the file, so a verifier can
    check the artifact against what was asked for without being handed the
    generator's own state -- and without trusting it.
    """
    root = ET.fromstring(zipfile.ZipFile(path).read("3D/3dmodel.model"))
    out = {}
    for m in root.findall("{%s}metadata" % NS["c"]):
        name = m.get("name", "")
        if not name.startswith("Parameter."):
            continue
        val = (m.text or "").strip()
        try:
            out[name[10:]] = float(val)
        except ValueError:
            if len(val) >= 2 and val[0] == val[-1] and val[0] in "'\"":
                val = val[1:-1]
            if val in ("True", "False"):
                val = val == "True"
            out[name[10:]] = val
    return out


def plane_radii(verts, tris, origin, u, d):
    """Radii, about the axis `u` through `origin`, where the mesh crosses the
    plane `d` along it. Measuring a taper means asking this at several depths
    and fitting a line: unlike a silhouette it sees the real surface, and it
    cannot be fooled by a curved profile.
    """
    rs = []
    for tri in tris:
        for i, j in ((0, 1), (1, 2), (2, 0)):
            a, b = verts[tri[i]], verts[tri[j]]
            da = sum((a[k] - origin[k]) * u[k] for k in range(3)) - d
            db = sum((b[k] - origin[k]) * u[k] for k in range(3)) - d
            if (da < 0.0) == (db < 0.0):
                continue
            t = da / (da - db)
            q = tuple(a[k] + t * (b[k] - a[k]) for k in range(3))
            radial = tuple(q[k] - origin[k] - d * u[k] for k in range(3))
            rs.append(math.sqrt(sum(v * v for v in radial)))
    return rs


def annuli(verts, tris, direction=None, min_triangles=8, min_gap=0.5,
           probe=3.0, shells=None):
    """Ring-shaped faces: coplanar faces with a hole in the middle.

    A tube's mouth is such a face, and so is a cap where an arm stops short of
    something. Grouping triangles by plane and normal finds them; a ring is a
    group with a clear radial gap, `min_triangles` or more triangles, and at
    least that many vertices on each side of the gap.

    Each one comes back with the point `on_face` -- on the ring itself, half
    way between its inner and outer radius, because the mean of an annulus's
    vertices is the hole -- and with `buried`: whether stepping off the face
    along its own normal lands in material. A mouth is not buried; a cap that
    runs into something is.

    `probe` has to be long enough to see past a passage. The junction's fence
    faces sit at the rim of a bore with their normals pointing along it, so a
    quarter of a millimetre step lands in the void and they read as mouths --
    which they are not, and the wye's verifier counted four where it wanted
    one. Three millimetres reaches the wall beyond.

    Coplanarity alone is not enough to call a group a ring: facets on the far
    side of a cylinder can land in the same plane by coincidence. A ring also
    has to go all the way round, so the group's vertices are checked for
    azimuthal coverage about the plane's centre.

    Pass `shells` from `shell_triangles` when calling this on a big mesh; the
    point-in-solid test derives it otherwise, which costs a union-find.
    """
    def sub(a, b):
        return (a[0] - b[0], a[1] - b[1], a[2] - b[2])

    def dot(a, b):
        return a[0] * b[0] + a[1] * b[1] + a[2] * b[2]

    def cross(a, b):
        return (a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2],
                a[0] * b[1] - a[1] * b[0])

    if shells is None:
        shells = shell_triangles(verts, tris)

    def tri_area(t):
        a, b, c = verts[t[0]], verts[t[1]], verts[t[2]]
        n = cross(sub(b, a), sub(c, a))
        return 0.5 * math.sqrt(dot(n, n))

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

    out = []
    for key, idxs in groups.items():
        if len(idxs) < min_triangles:
            continue
        pts = sorted({v for i in idxs for v in tris[i]})
        ring = [verts[v] for v in pts]
        # The centre of a ring is its *area* centroid, not the mean of its
        # vertices. A boolean splits a face unevenly, and the vertex mean then
        # walks off the axis -- the wye's mouth came back 1.25 mm off centre and
        # the radial gap test failed on a face that was perfectly sound. Area
        # weighting does not depend on how the face was triangulated.
        total = 0.0
        acc = [0.0, 0.0, 0.0]
        for i in idxs:
            t0, t1, t2 = tris[i]
            a, b, d = verts[t0], verts[t1], verts[t2]
            area = tri_area(tris[i])
            total += area
            for k in range(3):
                acc[k] += area * (a[k] + b[k] + d[k]) / 3.0
        c = (tuple(x / total for x in acc) if total > 0.0
             else tuple(sum(p[k] for p in ring) / len(ring) for k in range(3)))
        rad = sorted(math.sqrt(dot(sub(p, c), sub(p, c))) for p in ring)
        if rad[-1] <= min_gap:
            continue
        mid = 0.5 * (rad[0] + rad[-1])

        def radius(p):
            return math.sqrt(dot(sub(p, c), sub(p, c)))

        # Go all the way round? A ring's vertices cover the azimuth evenly; a
        # coincidental cluster of coplanar facets leaves gaps.
        n = key[:3]
        seed = (1.0, 0.0, 0.0) if abs(n[0]) < 0.9 else (0.0, 1.0, 0.0)
        e1 = tuple(seed[k] - dot(seed, n) * n[k] for k in range(3))
        m1 = math.sqrt(dot(e1, e1))
        e1 = tuple(v / m1 for v in e1)
        e2 = cross(n, e1)
        ax = [math.atan2(dot(sub(p, c), e2), dot(sub(p, c), e1))
              for p in ring if radius(p) > 1e-9]
        if len(ax) < min_triangles:
            continue
        # Is the material a ring? Scan the face's own plane on a polar grid
        # about the centre and ask which cells the triangles cover. A ring
        # covers every azimuth, contiguously in radius, with a hole all the way
        # round; a cluster of coplanar facets does not. This replaces "the
        # largest gap between consecutive vertex radii", which depended on how
        # the face was triangulated -- a boolean splits a mouth into hundreds of
        # coplanar triangles whose new vertices sit *inside* the ring's width,
        # and the wye's three mouths stopped being recognised at all.
        flat = {}
        for i in idxs:
            for v in tris[i]:
                if v not in flat:
                    q = sub(verts[v], c)
                    flat[v] = (dot(q, e1), dot(q, e2))
        outer = max(math.hypot(x, y) for x, y in flat.values())

        def covered(x, y):
            for i in idxs:
                t0, t1, t2 = tris[i]
                (x0, y0), (x1, y1), (x2, y2) = (flat[t0], flat[t1], flat[t2])
                d = (y1 - y2) * (x0 - x2) + (x2 - x1) * (y0 - y2)
                if abs(d) < 1e-15:
                    continue
                a = ((y1 - y2) * (x - x2) + (x2 - x1) * (y - y2)) / d
                b = ((y2 - y0) * (x - x2) + (x0 - x2) * (y - y2)) / d
                if a >= -1e-9 and b >= -1e-9 and a + b <= 1.0 + 1e-9:
                    return True
            return False

        spokes, steps = 72, 40
        inner, ok = None, True
        holes = []
        for k in range(spokes):
            ang = 2.0 * math.pi * k / spokes
            seen = [covered(outer * m / steps * math.cos(ang),
                            outer * m / steps * math.sin(ang))
                    for m in range(1, steps + 1)]
            if not any(seen):
                ok = False                     # this azimuth is uncovered
                break
            first = seen.index(True)
            last = len(seen) - 1 - seen[::-1].index(True)
            if not all(seen[first:last + 1]):
                ok = False                     # material, a hole, then material
                break
            here = outer * (first + 1) / steps
            holes.append(here)
            inner = here if inner is None else min(inner, here)
        if not ok or inner is None or inner <= min_gap:
            continue
        # A ring's hole is a circle, so the inner radius cannot wander from one
        # azimuth to the next. Without this a patch of coplanar facets -- the
        # junction's fence faces, for instance -- passes the scan and comes back
        # as a mouth: the wye grew a spurious annulus 7.18 mm out and its
        # verifier read the arm's length as 7.18 mm instead of 45.
        if max(holes) > 1.5 * min(holes):
            continue
        anchor = min(ring, key=lambda p: abs(radius(p) - mid))
        ra = radius(anchor)
        on_face = tuple(c[k] + (mid / ra) * (anchor[k] - c[k]) for k in range(3))
        # `inside` has its own default ray; only override it when asked.
        ray = {} if direction is None else {"direction": direction}
        # Several depths, because just behind a face that stops short of
        # something there can be a sliver of void -- the bore next door, where
        # nothing belongs -- with solid material beyond it. One step of
        # `probe` would call such a face a mouth.
        buried = False
        for depth in (probe, 2.0 * probe, 4.0 * probe, 8.0 * probe):
            step = tuple(on_face[k] + depth * key[k] for k in range(3))
            if inside(verts, tris, step, shells=shells, **ray):
                buried = True
                break
        out.append({"n": key[:3], "c": c, "inner": rad[0], "outer": rad[-1],
                    "on_face": on_face, "n_tri": len(idxs),
                    "area": sum(tri_area(tris[i]) for i in idxs),
                    "buried": buried})
    return out


def shells(verts, tris):
    """Split into connected shells: {indices, bbox, volume, max_radius}."""
    parent = list(range(len(verts)))

    def find(a):
        while parent[a] != a:
            parent[a] = parent[parent[a]]
            a = parent[a]
        return a

    for a, b, c in tris:
        for u, v in ((a, b), (b, c)):
            ru, rv = find(u), find(v)
            if ru != rv:
                parent[ru] = rv

    members = {}
    for i in range(len(verts)):
        members.setdefault(find(i), []).append(i)

    out = []
    for root_idx, idxs in members.items():
        pts = [verts[i] for i in idxs]
        xs, ys, zs = zip(*pts)
        mset = set(idxs)
        vol = 0.0
        for t in tris:
            if t[0] in mset:
                a, b, c = (verts[i] for i in t)
                vol += (a[0] * (b[1] * c[2] - b[2] * c[1])
                        - a[1] * (b[0] * c[2] - b[2] * c[0])
                        + a[2] * (b[0] * c[1] - b[1] * c[0]))
        out.append({
            "indices": idxs, "n": len(idxs),
            "bbox": (min(xs), min(ys), min(zs), max(xs), max(ys), max(zs)),
            "volume": vol / 6.0,
            "max_radius": max(math.hypot(x, y) for x, y, _ in pts),
        })
    return out


def aabb_overlap_pairs(shell_list, tol=1e-9):
    """Pairs of shells whose AABBs overlap by more than `tol` in every axis.

    Valid only when every shell is an axis-aligned box -- a single annulus
    shell spans the whole disc, so its AABB tells you nothing useful.
    """
    hits = []
    n = len(shell_list)
    for i in range(n):
        ai = shell_list[i]["bbox"]
        for j in range(i + 1, n):
            aj = shell_list[j]["bbox"]
            if all(min(ai[k + 3], aj[k + 3]) - max(ai[k], aj[k]) > tol
                   for k in range(3)):
                hits.append((i, j))
    return hits


def render(verts, tris, path, view="xy", su=10.0, sv=None, window=None,
           flip_v=False):
    """Rasterise the mesh to a PNG. `view` picks which axes map to (u, v).

    xy -> top view; yz -> side elevation. `su`/`sv` are pixels per mm and may
    differ, which is how you exaggerate a thin feature for inspection.
    """
    from PIL import Image, ImageDraw

    sv = su if sv is None else sv
    ax = {"xy": (0, 1), "yz": (1, 2), "xz": (0, 2)}[view]
    iu, iv = ax
    pts = [(v[iu], v[iv]) for v in verts]
    if window:
        u0, u1, v0, v1 = window
        sel = [p for p in pts if u0 <= p[0] <= u1 and v0 <= p[1] <= v1]
        if sel:
            pts = sel
    us = [p[0] for p in pts]
    vs = [p[1] for p in pts]
    umin, umax, vmin, vmax = min(us), max(us), min(vs), max(vs)
    umin, vmin = min(umin, 0.0), min(vmin, 0.0)
    umax, vmax = max(umax, 0.0), max(vmax, 0.0)
    w = max(1, int((umax - umin) * su) + 1)
    h = max(1, int((vmax - vmin) * sv) + 1)
    if w * h > 40_000_000:
        su = sv = min(su, sv) * math.sqrt(40_000_000 / (w * h))
        w = max(1, int((umax - umin) * su) + 1)
        h = max(1, int((vmax - vmin) * sv) + 1)

    def to_px(idx):
        u, v = verts[idx][iu], verts[idx][iv]
        px = (u - umin) * su
        py = (v - vmin) * sv
        if flip_v:
            py = (vmax - v) * sv
        return (px, py)

    img = Image.new("L", (w, h), 0)
    draw = ImageDraw.Draw(img)
    for t in tris:
        draw.polygon([to_px(i) for i in t], fill=255)
    img.save(path)
    return img, (w, h)


def raster_area_mm2(img, su, sv):
    """Filled area in mm^2, given the per-axis scales."""
    filled = sum(1 for p in img.getdata() if p > 127)
    return filled / (su * sv)


def scanline_runs(img, world_row, v_origin, sv, axis="u"):
    """Count bright runs along a scanline -- robust where area is not.

    Total pixel area is badly biased by pixel quantisation when features are
    only a few pixels wide, but the *number* of features along a line is not.
    """
    px = img.load()
    w, h = img.size
    row = int(round((world_row - v_origin) * sv))
    row = max(0, min(h - 1, row))
    vals = [px[x, row] > 127 for x in range(w)]
    runs, prev = 0, False
    for v in vals:
        if v and not prev:
            runs += 1
        prev = v
    return runs

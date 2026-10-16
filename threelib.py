"""threelib -- small, dependency-free helpers for generating 3D meshes.

Written for scripted/CAD-style modelling where you build geometry from
parametric primitives and export to 3MF or STL. Everything is a triangle
mesh; there is no B-rep and no exact boolean engine. Shapes are combined by
*overlap* and unified by the slicer at slice time (see "Combining shapes"
in GUIDE.md).

Conventions
-----------
* Right-handed, Z-up. Units are millimetres unless you say otherwise.
* Triangles wound counter-clockwise seen from outside => outward normals.
  `Mesh.signed_volume()` is positive only if every shell is wound outward;
  use it as a cheap global winding check.
* Polygons are given counter-clockwise in XY and must be convex when passed
  to `clip_convex`.
* Every vertex is (x, y, z) with floats; every triangle is a 3-tuple of
  vertex indices.

Only the standard library is required.
"""

from __future__ import annotations

import math
import struct
import zipfile
import zipfile as _zipfile  # noqa: F401  (documented dependency)
import zlib  # noqa: F401
import hashlib
from xml.sax.saxutils import escape

__all__ = [
    "Mesh", "Scene",
    "box", "prism", "cylinder", "tube", "tapered_tube", "sweep", "cone", "sphere", "torus", "wedge",
    "ring", "disc_polygon", "rect_polygon", "auto_segments",
    "clip_convex", "clip_halfplane", "clean_polygon", "polygon_area",
    "ensure_ccw", "rotate_z", "rotate_x", "rotate_y", "translate",
    "union_find_shells", "params_hash", "stamped_filename", "describe",
    "write_3mf", "write_stl",
]


# ---------------------------------------------------------------------------
# Geometry helpers
# ---------------------------------------------------------------------------

def _cross(ax, ay, bx, by):
    return ax * by - ay * bx


def _sub(a, b):
    return (a[0] - b[0], a[1] - b[1])


def _add(a, b):
    return (a[0] + b[0], a[1] + b[1])


def _scale(a, k):
    return (a[0] * k, a[1] * k)


def _dot(a, b):
    return a[0] * b[0] + a[1] * b[1]


def _dot3(a, b):
    return a[0] * b[0] + a[1] * b[1] + a[2] * b[2]


def _cross3(t, u):
    """3D cross product. Note the name: `_cross` above is the 2D scalar one."""
    return (t[1] * u[2] - t[2] * u[1],
            t[2] * u[0] - t[0] * u[2],
            t[0] * u[1] - t[1] * u[0])


def polygon_area(poly):
    """Signed shoelace area; positive when the polygon is CCW."""
    n = len(poly)
    if n < 3:
        return 0.0
    s = 0.0
    for i in range(n):
        x0, y0 = poly[i]
        x1, y1 = poly[(i + 1) % n]
        s += x0 * y1 - x1 * y0
    return s / 2.0


def ensure_ccw(poly):
    """Return the polygon counter-clockwise."""
    return list(poly) if polygon_area(poly) >= 0 else list(reversed(poly))


def clean_polygon(poly, dup_eps=1e-9, collinear_eps=1e-12):
    """Drop consecutive duplicates and collinear vertices.

    Clipping can leave near-coincident points; turning those into triangles
    creates degenerate geometry that slicers dislike.
    """
    out = []
    for p in poly:
        if not out or _dist(p, out[-1]) > dup_eps:
            out.append(p)
    if len(out) > 1 and _dist(out[0], out[-1]) <= dup_eps:
        out.pop()
    changed = True
    while changed and len(out) > 3:
        changed = False
        keep = []
        n = len(out)
        for i in range(n):
            a, b, c = out[i - 1], out[i], out[(i + 1) % n]
            ab, bc = _sub(b, a), _sub(c, b)
            cr = _cross(ab[0], ab[1], bc[0], bc[1])
            # Only drop b if it lies between a and c (true for convex output).
            if abs(cr) <= collinear_eps and _dot(ab, bc) >= 0:
                changed = True
                continue
            keep.append(b)
        out = keep
    return out


def _dist(a, b):
    return math.hypot(a[0] - b[0], a[1] - b[1])


# ---------------------------------------------------------------------------
# Polygon clipping
# ---------------------------------------------------------------------------

def clip_halfplane(poly, a, b):
    """Clip `poly` (list of 2D points) to the left of the directed line a->b.

    "Left" means cross(b-a, p-a) >= 0, i.e. the CCW interior of a polygon
    that traverses a->b as one of its edges.
    """
    if not poly:
        return []
    dx, dy = b[0] - a[0], b[1] - a[1]
    out = []
    n = len(poly)
    for i in range(n):
        p = poly[i]
        q = poly[(i + 1) % n]
        pin = _cross(dx, dy, p[0] - a[0], p[1] - a[1]) >= 0.0
        qin = _cross(dx, dy, q[0] - a[0], q[1] - a[1]) >= 0.0
        if pin:
            out.append(p)
        if pin != qin:
            denom = _cross(dx, dy, q[0] - p[0], q[1] - p[1])
            if denom != 0.0:
                t = _cross(dx, dy, a[0] - p[0], a[1] - p[1]) / denom
                out.append((p[0] + t * (q[0] - p[0]), p[1] + t * (q[1] - p[1])))
    return out


def clip_convex(subject, clip):
    """Clip a polygon against a *convex* CCW polygon (Sutherland-Hodgman).

    Both inputs should be simple polygons; because the clip region is convex
    the result is convex, so it can be triangulated with a simple fan.

    To cut a shape to a circle, pass a *faceted* disc polygon as `clip` and
    use the very same polygon to build the round solid -- then the clipped
    edge coincides exactly with the curved surface instead of poking through
    it. `disc_polygon` and `ring` are designed for exactly that.
    """
    out = list(subject)
    n = len(clip)
    for i in range(n):
        if not out:
            return []
        out = clip_halfplane(out, clip[i], clip[(i + 1) % n])
    return out


# ---------------------------------------------------------------------------
# Tessellation
# ---------------------------------------------------------------------------

def auto_segments(radius, tol=0.01, multiple=4, lo=12, hi=4096):
    """Segments for a full circle so the chord sagitta stays under `tol`.

    sagitta = r*(1 - cos(pi/n)); solve for the smallest n, then round up to
    a multiple of `multiple` so that flats land on the axes.
    """
    if radius <= 0 or tol <= 0:
        return lo
    ratio = max(-1.0, min(1.0, 1.0 - tol / radius))
    n = math.ceil(math.pi / math.acos(ratio))
    n = max(lo, min(hi, n))
    if multiple > 1:
        n = int(math.ceil(n / multiple) * multiple)
    return max(multiple, n)


def ring(radius, n, phase=0.0, cx=0.0, cy=0.0):
    """n CCW points on a circle. Same vertex set as `disc_polygon`."""
    return [
        (cx + radius * math.cos(phase + 2 * math.pi * i / n),
         cy + radius * math.sin(phase + 2 * math.pi * i / n))
        for i in range(n)
    ]


def disc_polygon(radius, n, phase=0.0, cx=0.0, cy=0.0):
    """Alias of `ring`, for use as a convex clip region."""
    return ring(radius, n, phase, cx, cy)


def rect_polygon(x0, y0, x1, y1):
    """CCW axis-aligned rectangle."""
    return [(x0, y0), (x1, y0), (x1, y1), (x0, y1)]


# ---------------------------------------------------------------------------
# Parameter fingerprints
# ---------------------------------------------------------------------------

def params_hash(params, length=8):
    """Stable short hash of a parameter dict.

    Scripts that default to a fixed output filename will silently overwrite
    an earlier model when a parameter is edited. Naming outputs
    `prefix_<hash>.3mf` means only *identical* parameters collide, so every
    variant you generate stays on disk.

    Stability comes from sorting keys and using repr() for values, so
    0.5 and 0.50 hash identically and dict insertion order does not matter.
    """

    def norm(v):
        # Round to 9 significant digits first, so that hair-level float
        # noise in a derived parameter does not fork the filename while a
        # deliberate change still does.
        if isinstance(v, float):
            return repr(float(f"{v:.9g}"))
        return repr(v)

    blob = "\n".join(f"{k}={norm(params[k])}" for k in sorted(params))
    return hashlib.sha256(blob.encode("utf-8")).hexdigest()[:length]


def stamped_filename(prefix, params, ext=".3mf", length=8):
    """`prefix_<8-char hash>.3mf` -- unique per parameter set."""
    return f"{prefix}_{params_hash(params, length)}{ext}"


def describe(params):
    """One-line summary of a parameter dict, for printouts."""
    return ", ".join(f"{k}={params[k]!r}" for k in sorted(params))


# ---------------------------------------------------------------------------
# Mesh container
# ---------------------------------------------------------------------------

class Mesh:
    """A triangle soup with optional per-shell bookkeeping.

    There is deliberately no boolean union: build shapes that overlap and
    let the slicer merge them. `shells()` recovers the individual solids
    when you need to reason about them separately.
    """

    def __init__(self, name="mesh"):
        self.name = name
        self.verts = []
        self.tris = []

    # -- construction --------------------------------------------------
    def add_verts(self, verts):
        base = len(self.verts)
        self.verts.extend((float(v[0]), float(v[1]), float(v[2])) for v in verts)
        return base

    def add_tris(self, tris, base=0):
        for t in tris:
            self.tris.append((base + t[0], base + t[1], base + t[2]))

    def add_solid(self, verts, faces):
        """Append an indexed polyhedron. `verts` is 3D, `faces` index loops."""
        base = self.add_verts(verts)
        for f in faces:
            for k in range(1, len(f) - 1):  # fan-triangulate
                self.tris.append((base + f[0], base + f[k], base + f[k + 1]))

    def merge(self, other, offset=(0.0, 0.0, 0.0)):
        ox, oy, oz = offset
        base = len(self.verts)
        self.verts.extend((v[0] + ox, v[1] + oy, v[2] + oz) for v in other.verts)
        self.tris.extend((a + base, b + base, c + base) for a, b, c in other.tris)
        return self

    def transform(self, fn):
        self.verts = [fn(v) for v in self.verts]
        return self

    def copy(self):
        m = Mesh(self.name)
        m.verts = list(self.verts)
        m.tris = list(self.tris)
        return m

    def __iadd__(self, other):
        self.merge(other)
        return self

    # -- measurement ---------------------------------------------------
    def bbox(self):
        if not self.verts:
            return None
        xs = [v[0] for v in self.verts]
        ys = [v[1] for v in self.verts]
        zs = [v[2] for v in self.verts]
        return (min(xs), min(ys), min(zs), max(xs), max(ys), max(zs))

    def signed_volume(self):
        """Sum of signed tetrahedron volumes. Positive iff all shells face out."""
        total = 0.0
        for a, b, c in self.tris:
            ax, ay, az = self.verts[a]
            bx, by, bz = self.verts[b]
            cx, cy, cz = self.verts[c]
            total += (ax * (by * cz - bz * cy)
                      - ay * (bx * cz - bz * cx)
                      + az * (bx * cy - by * cx))
        return total / 6.0

    def area_xy(self):
        """Area of the silhouette projected on XY (union is NOT computed)."""
        total = 0.0
        for a, b, c in self.tris:
            ax, ay, _ = self.verts[a]
            bx, by, _ = self.verts[b]
            cx, cy, _ = self.verts[c]
            total += abs((bx - ax) * (cy - ay) - (cx - ax) * (by - ay)) / 2.0
        return total

    def face_areas_3d(self):
        out = []
        for a, b, c in self.tris:
            ax, ay, az = self.verts[a]
            bx, by, bz = self.verts[b]
            cx, cy, cz = self.verts[c]
            ux, uy, uz = bx - ax, by - ay, bz - az
            vx, vy, vz = cx - ax, cy - ay, cz - az
            nx, ny, nz = uy * vz - uz * vy, uz * vx - ux * vz, ux * vy - uy * vx
            out.append(0.5 * math.sqrt(nx * nx + ny * ny + nz * nz))
        return out

    def max_radius_xy(self, cx=0.0, cy=0.0):
        return max(math.hypot(v[0] - cx, v[1] - cy) for v in self.verts) if self.verts else 0.0

    def check(self):
        """Basic integrity report as a list of (ok, message)."""
        checks = []
        nv = len(self.verts)
        checks.append(("no NaN/inf verts",
                       all(all(math.isfinite(c) for c in v) for v in self.verts)))
        checks.append(("indices in range",
                       all(0 <= i < nv for t in self.tris for i in t)))
        checks.append(("no degenerate triangles",
                       all(len(set(t)) == 3 for t in self.tris)))
        checks.append(("outward winding",
                       self.signed_volume() > 0 or not self.tris))
        return checks

    def shells(self):
        """Group vertex indices into connected shells -> list of (root, indices)."""
        parent = list(range(len(self.verts)))

        def find(a):
            while parent[a] != a:
                parent[a] = parent[parent[a]]
                a = parent[a]
            return a

        for a, b, c in self.tris:
            for u, v in ((a, b), (b, c)):
                ru, rv = find(u), find(v)
                if ru != rv:
                    parent[ru] = rv
        groups = {}
        for i in range(len(self.verts)):
            groups.setdefault(find(i), []).append(i)
        return list(groups.items())


def union_find_shells(mesh):
    """Same as Mesh.shells() but returns shells as dicts of extent + volume."""
    out = []
    tri_of = {}
    for t in mesh.tris:
        tri_of.setdefault(t[0], []).append(t)
    for root, members in mesh.shells():
        mset = set(members)
        xs = [mesh.verts[i][0] for i in members]
        ys = [mesh.verts[i][1] for i in members]
        zs = [mesh.verts[i][2] for i in members]
        vol = 0.0
        for t in mesh.tris:
            if t[0] in mset:
                a, b, c = (mesh.verts[i] for i in t)
                vol += (a[0] * (b[1] * c[2] - b[2] * c[1])
                        - a[1] * (b[0] * c[2] - b[2] * c[0])
                        + a[2] * (b[0] * c[1] - b[1] * c[0]))
        out.append({
            "root": root,
            "indices": members,
            "n_verts": len(members),
            "bbox": (min(xs), min(ys), min(zs), max(xs), max(ys), max(zs)),
            "volume": vol / 6.0,
        })
    return out


# ---------------------------------------------------------------------------
# Transforms (applied to every vertex of a Mesh)
# ---------------------------------------------------------------------------

def translate(mesh, dx=0.0, dy=0.0, dz=0.0):
    return mesh.transform(lambda v: (v[0] + dx, v[1] + dy, v[2] + dz))

def rotate_z(mesh, deg):
    a = math.radians(deg)
    ca, sa = math.cos(a), math.sin(a)
    return mesh.transform(lambda v: (v[0] * ca - v[1] * sa,
                                    v[0] * sa + v[1] * ca, v[2]))


def rotate_x(mesh, deg):
    a = math.radians(deg)
    ca, sa = math.cos(a), math.sin(a)
    return mesh.transform(lambda v: (v[0], v[1] * ca - v[2] * sa,
                                    v[1] * sa + v[2] * ca))


def rotate_y(mesh, deg):
    a = math.radians(deg)
    ca, sa = math.cos(a), math.sin(a)
    return mesh.transform(lambda v: (v[0] * ca + v[2] * sa, v[1],
                                    -v[0] * sa + v[2] * ca))


# ---------------------------------------------------------------------------
# Primitives
# ---------------------------------------------------------------------------

def box(mesh, x0, y0, z0, x1, y1, z1):
    """Axis-aligned box, wound outward."""
    v = [(x0, y0, z0), (x1, y0, z0), (x1, y1, z0), (x0, y1, z0),
         (x0, y0, z1), (x1, y0, z1), (x1, y1, z1), (x0, y1, z1)]
    faces = [
        (0, 3, 2, 1),  # bottom (-Z)
        (4, 5, 6, 7),  # top (+Z)
        (0, 1, 5, 4),  # -Y
        (1, 2, 6, 5),  # +X
        (2, 3, 7, 6),  # +Y
        (3, 0, 4, 7),  # -X
    ]
    mesh.add_solid(v, faces)
    return mesh


def prism(mesh, poly, z0, z1):
    """Extrude a CCW convex polygon between z0 and z1."""
    poly = ensure_ccw(poly)
    n = len(poly)
    verts = [(x, y, z0) for x, y in poly] + [(x, y, z1) for x, y in poly]
    faces = [tuple(range(n - 1, -1, -1)),          # bottom, reversed -> -Z
             tuple(range(n, 2 * n))]               # top -> +Z
    for i in range(n):
        j = (i + 1) % n
        faces.append((i, j, n + j, n + i))          # side, CCW -> outward
    mesh.add_solid(verts, faces)
    return mesh


def cylinder(mesh, r, z0, z1, n=64, cx=0.0, cy=0.0, phase=0.0):
    return prism(mesh, ring(r, n, phase, cx, cy), z0, z1)


def tapered_tube(mesh, r_in0, r_out0, r_in1, r_out1, z0, z1, n=64,
                 cx=0.0, cy=0.0, phase=0.0):
    """Hollow tapered tube: inner and outer radii given separately at z0 and z1.

    Unlike `tube`, both the bore and the outside may change radius with height.
    Keep the bore vertical (r_in0 == r_in1) when printing a part with this: a
    bore that narrows as the part rises leaves its ceiling unsupported.
    """
    inner_b = [(x, y, z0) for x, y in ring(r_in0, n, phase, cx, cy)]
    inner_t = [(x, y, z1) for x, y in ring(r_in1, n, phase, cx, cy)]
    outer_b = [(x, y, z0) for x, y in ring(r_out0, n, phase, cx, cy)]
    outer_t = [(x, y, z1) for x, y in ring(r_out1, n, phase, cx, cy)]
    verts = inner_b + inner_t + outer_b + outer_t
    IB, IT, OB, OT = 0, n, 2 * n, 3 * n
    faces = []
    for i in range(n):
        j = (i + 1) % n
        faces.append((OB + i, OB + j, OT + j, OT + i))   # outer, outward
        faces.append((IB + i, IT + i, IT + j, IB + j))   # bore, toward axis
        faces.append((OT + i, OT + j, IT + j, IT + i))   # top annulus, +Z
        faces.append((IB + i, IB + j, OB + j, OB + i))   # bottom annulus, -Z
    mesh.add_solid(verts, faces)
    return mesh


def revolve(mesh, profile, n=64, phase=0.0):
    """Revolve a profile about the z axis into a solid of revolution.

    `profile` is the outline of the solid's cross-section, run from the axis
    outward and back again: a list of (radius, z) in order of increasing z,
    starting and ending on the axis for a solid with two poles. A point at
    radius 0 becomes a single vertex on the axis rather than a ring of
    coincident ones -- which is the difference between this and sweeping an
    annulus with an inner radius of zero, where the collapsed inner surface
    contributes three degenerate triangles for every real one.

    A profile that starts or ends at radius > 0 gets a flat disc there, so a
    capped cylinder is [(0, z0), (r, z0), (r, z1), (0, z1)] and a cone on a
    disc is [(0, z0), (r, z0), (r1, z1)].

    This is the primitive for anything solid and round -- a stopper, a knob, a
    bulb, a flask -- where `sweep` only does tubes and `cone` only does one
    straight run.
    """
    if len(profile) < 2:
        raise ValueError("a revolve needs at least two profile points")
    for k in range(len(profile) - 1):
        if profile[k + 1][1] < profile[k][1]:
            raise ValueError("the profile must run in order of increasing z")
        if profile[k][0] < 0.0:
            raise ValueError("a radius cannot be negative")
    verts, rings = [], []
    for r, z in profile:
        if r <= 0.0:
            rings.append((True, len(verts)))
            verts.append((0.0, 0.0, z))
        else:
            first = len(verts)
            for i in range(n):
                a = phase + 2.0 * math.pi * i / n
                verts.append((r * math.cos(a), r * math.sin(a), z))
            rings.append((False, first))
    faces = []
    for k in range(len(profile) - 1):
        a_pole, a = rings[k]
        b_pole, b = rings[k + 1]
        for i in range(n):
            j = (i + 1) % n
            if a_pole and b_pole:
                raise ValueError("a profile cannot run along the axis")
            elif a_pole:
                faces.append((a, b + j, b + i))
            elif b_pole:
                faces.append((a + i, a + j, b))
            else:
                faces.append((a + i, a + j, b + j, b + i))
    # A flat disc where the profile stops short of the axis. Seen from +z the
    # bottom one is wound the other way, so both faces point out of the solid.
    if not rings[0][0]:
        c = len(verts)
        verts.append((0.0, 0.0, profile[0][1]))
        for i in range(n):
            j = (i + 1) % n
            faces.append((c, rings[0][1] + j, rings[0][1] + i))
    if not rings[-1][0]:
        c = len(verts)
        verts.append((0.0, 0.0, profile[-1][1]))
        for i in range(n):
            j = (i + 1) % n
            faces.append((c, rings[-1][1] + i, rings[-1][1] + j))
    mesh.add_solid(verts, faces)
    return mesh


def sweep(mesh, stations, n=64, phase=0.0, ref=(0.0, 1.0, 0.0)):
    """Loft an annular section along a planar centreline.

    `stations` is a list of ((x, y, z), (tx, ty, tz), r_in, r_out). Consecutive
    stations are joined by a ring of quads, so the centreline may bend; the
    two ends are capped with annuli.

    Tangents are given explicitly rather than derived from the neighbouring
    points. A station that only changes radius -- the flat floor of a socket,
    where the bore steps in -- has no direction of travel of its own, and
    deriving one there produces a ring skew to both its neighbours.

    The centreline must lie in the xz plane. Rings are built on a fixed
    out-of-plane reference, which is what stops a bend from rolling: a frame
    carried along a curve accumulates twist, and across a quarter turn the far
    end arrives rotated.

    `r_in` is the bore and `r_out` the outside, both perpendicular to the
    centreline. A zero `r_in` degenerates the inner surface to the axis, which
    is what a solid rod wants.

    `ref` is the out-of-plane reference the rings are built on, and it is
    projected perpendicular to the tangent before use. Do not skip that step:
    for a centreline in the xz plane the default reference already *is*
    perpendicular so it looks like a no-op, but against an oblique axis it is
    not, the rings come out as ellipses, and the solid is silently the wrong
    size -- a diagonal rod came out 18% light.
    """
    rings = []
    for p, t, r_in, r_out in stations:
        # Project the reference perpendicular to the tangent, then complete
        # the right-handed frame with the cross product.
        proj = _dot3(ref, t)
        e1 = (ref[0] - proj * t[0], ref[1] - proj * t[1], ref[2] - proj * t[2])
        m1 = math.sqrt(_dot3(e1, e1))
        if m1 < 1e-9:
            alt = (1.0, 0.0, 0.0) if abs(t[0]) < 0.9 else (0.0, 0.0, 1.0)
            proj = _dot3(alt, t)
            e1 = (alt[0] - proj * t[0], alt[1] - proj * t[1], alt[2] - proj * t[2])
            m1 = math.sqrt(_dot3(e1, e1))
        e1 = (e1[0] / m1, e1[1] / m1, e1[2] / m1)
        e2 = _cross3(t, e1)
        m = math.sqrt(_dot3(e2, e2))
        if m < 1e-12:
            raise ValueError("sweep tangent is parallel to the reference axis")
        e2 = (e2[0] / m, e2[1] / m, e2[2] / m)

        def at(radius, k):
            a = 2.0 * math.pi * k / n + phase
            ca, sa = math.cos(a), math.sin(a)
            return tuple(p[i] + radius * (ca * e1[i] + sa * e2[i])
                         for i in range(3))

        rings.append(([at(r_in, k) for k in range(n)],
                      [at(r_out, k) for k in range(n)]))

    verts, faces = [], []
    for r_in, r_out in rings:
        verts.extend(r_in)
        verts.extend(r_out)
    # Station k: bore at 2kn, outside at 2kn + n.
    for k in range(len(rings) - 1):
        ib, ob = 2 * k * n, 2 * k * n + n
        it, ot = ib + 2 * n, ob + 2 * n
        for i in range(n):
            j = (i + 1) % n
            faces.append((ob + i, ob + j, ot + j, ot + i))   # outside
            faces.append((ib + i, it + i, it + j, ib + j))   # bore
    ib, ob = 0, n
    for i in range(n):                                       # start cap
        j = (i + 1) % n
        faces.append((ib + i, ib + j, ob + j, ob + i))
    last = 2 * (len(rings) - 1) * n
    ib, ob = last, last + n
    for i in range(n):                                       # end cap
        j = (i + 1) % n
        faces.append((ob + i, ob + j, ib + j, ib + i))
    mesh.add_solid(verts, faces)
    return mesh


def loft(mesh, sections, close=True):
    """Loft a closed polygon through a list of sections.

    `sections` is a list of point lists, all the same length, each going round
    the same way. Consecutive sections are joined ring by ring and the ends are
    capped with the polygon itself, so what comes out is one closed solid.

    Nothing here assumes the sections are alike -- only that they have the same
    number of vertices and go round in the same direction. That is what makes
    this the primitive for shapes whose cross-section *turns* as it rises: a
    twisted band, a helical fin, a swept channel wall. `sweep` cannot do those,
    because its rings all face the same way by design (a frame carried along a
    curve accumulates twist, and for a bent tube that is exactly what you do
    not want).

    Section k must be wound counter-clockwise seen from the direction of
    travel, or the solid comes out inside-out; check `signed_volume()`.
    """
    if len(sections) < 2:
        raise ValueError("a loft needs at least two sections")
    counts = {len(s) for s in sections}
    if len(counts) != 1:
        raise ValueError(f"sections have different vertex counts: {sorted(counts)}")
    n = counts.pop()
    if n < 3:
        raise ValueError("a section needs at least three points")
    verts = [tuple(p) for s in sections for p in s]
    faces = []
    for k in range(len(sections) - 1):
        a, b = k * n, (k + 1) * n
        for i in range(n):
            j = (i + 1) % n
            faces.append((a + i, a + j, b + j, b + i))
    if close:
        faces.append(tuple(range(n - 1, -1, -1)))
        faces.append(tuple(range(len(verts) - n, len(verts))))
    mesh.add_solid(verts, faces)
    return mesh


def cone(mesh, r_bot, r_top, z0, z1, n=64, cx=0.0, cy=0.0, phase=0.0):
    """Truncated cone / frustum. r_top=0 gives a cone, r_bot=r_top a cylinder."""
    bottom = ring(r_bot, n, phase, cx, cy) if r_bot > 0 else [(cx, cy)] * n
    top = ring(r_top, n, phase, cx, cy) if r_top > 0 else [(cx, cy)] * n
    verts = [(x, y, z0) for x, y in bottom] + [(x, y, z1) for x, y in top]
    faces = [tuple(range(n - 1, -1, -1)), tuple(range(n, 2 * n))]
    for i in range(n):
        j = (i + 1) % n
        faces.append((i, j, n + j, n + i))
    mesh.add_solid(verts, faces)
    return mesh


def tube(mesh, r_in, r_out, z0, z1, n=64, cx=0.0, cy=0.0, phase=0.0):
    """Annular prism (a pipe wall). Open bore, capped top and bottom."""
    inner_b = ring(r_in, n, phase, cx, cy)
    inner_t = [(x, y, z1) for x, y in inner_b]
    inner_b = [(x, y, z0) for x, y in inner_b]
    outer_b = [(x, y, z0) for x, y in ring(r_out, n, phase, cx, cy)]
    outer_t = [(x, y, z1) for x, y in ring(r_out, n, phase, cx, cy)]
    verts = inner_b + inner_t + outer_b + outer_t
    IB, IT, OB, OT = 0, n, 2 * n, 3 * n
    faces = []
    for i in range(n):
        j = (i + 1) % n
        faces.append((OB + i, OB + j, OT + j, OT + i))   # outer, outward
        faces.append((IB + i, IT + i, IT + j, IB + j))   # inner, toward axis
        faces.append((OT + i, OT + j, IT + j, IT + i))   # top annulus, +Z
        faces.append((IB + i, IB + j, OB + j, OB + i))   # bottom annulus, -Z
    mesh.add_solid(verts, faces)
    return mesh


def sphere(mesh, r, cz=0.0, seg=32, rings=16, cx=0.0, cy=0.0):
    verts = [(cx, cy, cz + r)]
    for i in range(1, rings):
        phi = math.pi * i / rings
        z = cz + r * math.cos(phi)
        rad = r * math.sin(phi)
        for j in range(seg):
            th = 2 * math.pi * j / seg
            verts.append((cx + rad * math.cos(th), cy + rad * math.sin(th), z))
    verts.append((cx, cy, cz - r))
    faces = []
    top, bottom = 0, len(verts) - 1
    for j in range(seg):
        faces.append((top, 1 + j, 1 + (j + 1) % seg))
    for i in range(rings - 2):
        base0, base1 = 1 + i * seg, 1 + (i + 1) * seg
        for j in range(seg):
            k = (j + 1) % seg
            faces.append((base0 + j, base1 + j, base1 + k, base0 + k))
    base = 1 + (rings - 2) * seg
    for j in range(seg):
        faces.append((bottom, base + (j + 1) % seg, base + j))
    mesh.add_solid(verts, faces)
    return mesh


def torus(mesh, r_major, r_minor, seg=48, sides=16, cz=0.0):
    verts, faces = [], []
    for i in range(seg):
        a = 2 * math.pi * i / seg
        ca, sa = math.cos(a), math.sin(a)
        for j in range(sides):
            b = 2 * math.pi * j / sides
            rr = r_major + r_minor * math.cos(b)
            verts.append((rr * ca, rr * sa, cz + r_minor * math.sin(b)))
    for i in range(seg):
        i2 = (i + 1) % seg
        for j in range(sides):
            j2 = (j + 1) % sides
            faces.append((i * sides + j, i2 * sides + j,
                          i2 * sides + j2, i * sides + j2))
    mesh.add_solid(verts, faces)
    return mesh


def _perm_parity(p):
    """0 if p is an even permutation of (0,1,2), else 1."""
    perm = list(p)
    swaps = 0
    for i in range(3):
        while perm[i] != i:
            j = perm.index(i)
            perm[i], perm[j] = perm[j], perm[i]
            swaps += 1
    return swaps % 2


def wedge(mesh, x0, y0, z0, x1, y1, z1, axis="x"):
    """A ramp: full height at one end, tapering to zero thickness at the other.

    `axis` is the direction the taper runs along. The result is a triangular
    prism with half the volume of the bounding box.

    Built in a canonical orientation (taper along X) and then permuted. An
    odd permutation mirrors the space and would invert every normal, so the
    face winding is reversed to compensate.
    """
    # Canonical: taper along X, extruded along Y, rising in Z.
    v = [
        (x0, y0, z0),   # 0
        (x1, y0, z0),   # 1
        (x1, y0, z1),   # 2
        (x0, y1, z0),   # 3
        (x1, y1, z0),   # 4
        (x1, y1, z1),   # 5
    ]
    faces = [
        (0, 3, 4, 1),   # bottom, -Z
        (0, 1, 2),      # front,  -Y
        (3, 5, 4),      # back,   +Y
        (1, 4, 5, 2),   # right,  +X
        (0, 2, 5, 3),   # sloped face, (-X, +Z)
    ]
    perms = {"x": (0, 1, 2), "y": (1, 0, 2), "z": (2, 1, 0)}
    if axis not in perms:
        raise ValueError(f"axis must be 'x', 'y' or 'z' (got {axis!r})")
    p = perms[axis]
    verts = [tuple(c[k] for k in p) for c in v]
    if _perm_parity(p):
        faces = [tuple(reversed(f)) for f in faces]
    mesh.add_solid(verts, faces)
    return mesh


# ---------------------------------------------------------------------------
# Scene / export
# ---------------------------------------------------------------------------

class Scene:
    """Named meshes exported as separate objects in one 3MF."""

    def __init__(self):
        self.objects = []  # list of (name, Mesh)

    def add(self, mesh, name=None):
        self.objects.append((name or mesh.name, mesh))
        return mesh

    def merged(self):
        m = Mesh("merged")
        for _, mesh in self.objects:
            m.merge(mesh)
        return m


_CONTENT_TYPES = """<?xml version="1.0" encoding="UTF-8"?>
<Types xmlns="http://schemas.openxmlformats.org/package/2006/content-types">
  <Default Extension="rels"
           ContentType="application/vnd.openxmlformats-package.relationships+xml"/>
  <Default Extension="model"
           ContentType="application/vnd.ms-package.3dmanufacturing-3dmodel+xml"/>
</Types>
"""

_RELS = """<?xml version="1.0" encoding="UTF-8"?>
<Relationships
    xmlns="http://schemas.openxmlformats.org/package/2006/relationships">
  <Relationship Id="rel0"
      Target="/3D/3dmodel.model"
      Type="http://schemas.microsoft.com/3dmanufacturing/2013/01/3dmodel"/>
</Relationships>
"""

_ZIP_DATE = (1980, 1, 1, 0, 0, 0)  # fixed -> reproducible archives


def _model_xml(meshes, title="model", description="", params=None):
    parts = ['<?xml version="1.0" encoding="UTF-8"?>',
             '<model unit="millimeter" xml:lang="en-US"',
             '       xmlns="http://schemas.microsoft.com/3dmanufacturing/core/2015/02">',
             f'  <metadata name="Title">{escape(title)}</metadata>']
    if description:
        parts.append(f'  <metadata name="Description">{escape(description)}</metadata>')
    if params:
        # Record the exact parameter set that produced this file, so a model
        # found on disk can be traced back to the script that made it.
        for k in sorted(params):
            parts.append(f'  <metadata name="Parameter.{k}">'
                         f'{escape(repr(params[k]))}</metadata>')
        parts.append(f'  <metadata name="ParameterHash">'
                     f'{params_hash(params)}</metadata>')
    parts.append('  <resources>')
    for oid, (name, mesh) in enumerate(meshes, start=1):
        # repr() gives the shortest string that round-trips exactly. Fixed
        # precision (e.g. %.6g) silently perturbs every coordinate -- enough
        # to break exact-fit geometry and to make radii differ per-vertex.
        v = "\n".join(f'        <vertex x="{x!r}" y="{y!r}" z="{z!r}"/>'
                      for x, y, z in mesh.verts)
        t = "\n".join(f'        <triangle v1="{a}" v2="{b}" v3="{c}"/>'
                      for a, b, c in mesh.tris)
        parts.append(f'    <object id="{oid}" name="{escape(name)}" type="model">')
        parts.append('      <mesh>')
        parts.append(f'        <vertices>\n{v}\n        </vertices>')
        parts.append(f'        <triangles>\n{t}\n        </triangles>')
        parts.append('      </mesh>')
        parts.append('    </object>')
    parts.append('  </resources>')
    parts.append('  <build>')
    for oid, _ in enumerate(meshes, start=1):
        parts.append(f'    <item objectid="{oid}"/>')
    parts.append('  </build>')
    parts.append('</model>')
    return "\n".join(parts) + "\n"


def write_3mf(path, meshes, title="model", description="", params=None):
    """Write one or more (name, Mesh) pairs to a 3MF file.

    `meshes` may be a single (name, Mesh) pair, a bare Mesh, or a list of
    either. Objects land in the build list in order. Pass the generating
    `params` dict to embed them in the file's metadata.
    """
    if isinstance(meshes, Mesh):
        meshes = [(meshes.name, meshes)]
    elif (isinstance(meshes, tuple) and len(meshes) == 2
          and isinstance(meshes[1], Mesh)):
        meshes = [meshes]
    meshes = [(n, m) if isinstance(m, Mesh) else (n, Mesh(n)) for n, m in meshes]
    xml = _model_xml(meshes, title, description, params)
    with zipfile.ZipFile(path, "w", zipfile.ZIP_DEFLATED) as z:
        for name, data in (("[Content_Types].xml", _CONTENT_TYPES),
                           ("_rels/.rels", _RELS),
                           ("3D/3dmodel.model", xml)):
            info = zipfile.ZipInfo(name, date_time=_ZIP_DATE)
            info.compress_type = zipfile.ZIP_DEFLATED
            info.external_attr = 0o644 << 16
            z.writestr(info, data)
    return path


def write_stl(path, mesh, binary=True):
    """Export a Mesh to STL (binary by default)."""
    tris = mesh.tris
    if binary:
        with open(path, "wb") as fh:
            fh.write(b"threelib STL".ljust(80, b"\0"))
            fh.write(struct.pack("<I", len(tris)))
            for a, b, c in tris:
                ax, ay, az = mesh.verts[a]
                bx, by, bz = mesh.verts[b]
                cx, cy, cz = mesh.verts[c]
                ux, uy, uz = bx - ax, by - ay, bz - az
                vx, vy, vz = cx - ax, cy - ay, cz - az
                nx, ny, nz = uy * vz - uz * vy, uz * vx - ux * vz, ux * vy - uy * vx
                ln = math.sqrt(nx * nx + ny * ny + nz * nz) or 1.0
                fh.write(struct.pack("<3f", nx / ln, ny / ln, nz / ln))
                fh.write(struct.pack("<9f", ax, ay, az, bx, by, bz, cx, cy, cz))
                fh.write(struct.pack("<H", 0))
    else:
        with open(path, "w") as fh:
            fh.write("solid mesh\n")
            for a, b, c in tris:
                pts = [mesh.verts[i] for i in (a, b, c)]
                ux = [pts[1][k] - pts[0][k] for k in range(3)]
                vx = [pts[2][k] - pts[0][k] for k in range(3)]
                nx = [ux[1] * vx[2] - ux[2] * vx[1],
                      ux[2] * vx[0] - ux[0] * vx[2],
                      ux[0] * vx[1] - ux[1] * vx[0]]
                ln = math.sqrt(sum(n * n for n in nx)) or 1.0
                fh.write("  facet normal %.6g %.6g %.6g\n    outer loop\n"
                         % tuple(n / ln for n in nx))
                for p in pts:
                    fh.write("      vertex %.6g %.6g %.6g\n" % p)
                fh.write("    endloop\n  endfacet\n")
            fh.write("endsolid mesh\n")
    return path

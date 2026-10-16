#!/usr/bin/env python3
"""Boolean operations on closed meshes, by BSP tree.

Why this exists. Every part in this repo is built by overlapping shells: a
sleeve's wall passes through a barrel's wall, a disc dips into a cone, a bar
lies across a junction. `checks.inside` unions those shells when it answers a
point query, but the *mesh* is still two surfaces crossing each other. A slicer
that fills by winding unions it; a slicer or viewer that fills even-odd counts
the overlap twice and shows the joint as a hole. A hole in the model is a hole
in the print, so the union has to be made, not assumed.

The algorithm is the classic BSP one -- Naylor, Amanatides and Thibault, in the
shape `csg.js` popularised. Each solid becomes a tree of planes; `clip_to`
cuts one tree's polygons against another's volume; the union is A's polygons
outside B plus B's polygons outside A.

Every face in this repo is planar -- curved surfaces are faceted -- so the
result is exact up to the faceting rather than approximate. The stopcock's
clearance is 0.08 mm and a voxel method would not hold it.

What to expect:

  * the output is welded, so the pieces really are one shell and
    `Mesh.shells()` says so; without welding a union of two shells comes back
    as two shells that happen to touch, which is the bug this module exists to
    fix;
  * coplanar faces are resolved by normal direction: a face shared by two
    solids is dropped, which is what makes two boxes side by side one box;
  * the input must be closed and wound outward (`signed_volume() > 0`), and the
    faces of one solid must not be coincident-and-opposite with the other's,
    which is a zero-thickness gap rather than a solid;
  * it is O(n log n) on a good day and O(n^2) on a bad one. A part of a few
    thousand triangles is seconds; ten thousand is minutes.

Usage:

    import csglib as csg
    solid = csg.union(shell_a, shell_b, shell_c)
    bored = csg.subtract(solid, hole)
"""

import math
import sys

import checks
import threelib as t

EPS = 1e-9
COPLANAR, FRONT, BACK, SPANNING = 0, 1, 2, 3

# A BSP built by inserting polygons one at a time is as deep as the polygon
# count in the worst case -- a linked list of planes, tens of thousands deep for
# a part this size. Every tree walk below is iterative for that reason: the
# recursive form blows the interpreter's stack, and raising the limit just
# moves the crash into C.


def _sub(a, b):
    return (a[0] - b[0], a[1] - b[1], a[2] - b[2])


def _cross(a, b):
    return (a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2],
            a[0] * b[1] - a[1] * b[0])


def _dot(a, b):
    return a[0] * b[0] + a[1] * b[1] + a[2] * b[2]


def _newell(pts):
    """A polygon's area vector, by Newell's method.

    The cross product of the first three vertices is the cheap way and the
    wrong one: a polygon that came out of a split can easily start with three
    collinear points, and then the cross product is zero and the polygon has no
    plane. Newell sums over every edge, so collinear runs anywhere in the
    polygon cannot take the normal away.
    """
    n = [0.0, 0.0, 0.0]
    for i in range(len(pts)):
        a, b = pts[i], pts[(i + 1) % len(pts)]
        n[0] += (a[1] - b[1]) * (a[2] + b[2])
        n[1] += (a[2] - b[2]) * (a[0] + b[0])
        n[2] += (a[0] - b[0]) * (a[1] + b[1])
    return (n[0], n[1], n[2])


def _unit(v):
    n = math.sqrt(_dot(v, v))
    if n < 1e-15:
        raise ValueError("a degenerate polygon has no plane")
    return (v[0] / n, v[1] / n, v[2] / n)


class Poly:
    """A planar convex polygon, with the plane it lies in."""

    __slots__ = ("pts", "normal", "w")

    def __init__(self, pts, normal=None, w=None):
        self.pts = [tuple(p) for p in pts]
        if normal is None:
            normal = _unit(_newell(self.pts))
            w = _dot(normal, self.pts[0])
        self.normal, self.w = normal, w

    def flipped(self):
        n = self.normal
        return Poly(list(reversed(self.pts)), (-n[0], -n[1], -n[2]), -self.w)

    def flip(self):
        self.pts.reverse()
        n = self.normal
        self.normal, self.w = (-n[0], -n[1], -n[2]), -self.w

    def split(self, plane, cof, cob, front, back):
        """Sort this polygon into the plane's four buckets, splitting if needed."""
        pn, pw = plane
        types, kind = [], COPLANAR
        for v in self.pts:
            d = _dot(pn, v) - pw
            ty = BACK if d < -EPS else (FRONT if d > EPS else COPLANAR)
            kind |= ty
            types.append(ty)
        if kind == COPLANAR:
            (cof if _dot(pn, self.normal) > 0.0 else cob).append(self)
        elif kind == FRONT:
            front.append(self)
        elif kind == BACK:
            back.append(self)
        else:
            f, b = [], []
            n = len(self.pts)
            for i in range(n):
                j = (i + 1) % n
                ti, tj = types[i], types[j]
                vi, vj = self.pts[i], self.pts[j]
                if ti != BACK:
                    f.append(vi)
                if ti != FRONT:
                    b.append(vi)
                if (ti | tj) == SPANNING:
                    di = _dot(pn, vi) - pw
                    dj = _dot(pn, vj) - pw
                    if abs(di - dj) > 1e-15:
                        u = di / (di - dj)
                        x = (vi[0] + u * (vj[0] - vi[0]),
                             vi[1] + u * (vj[1] - vi[1]),
                             vi[2] + u * (vj[2] - vi[2]))
                        f.append(x)
                        b.append(x)
            # A piece of a polygon lies in the same plane as the polygon, so
            # the plane is inherited rather than recomputed -- recomputing is
            # how a cut face ends up with no plane at all.
            if len(f) >= 3:
                front.append(Poly(f, self.normal, self.w))
            if len(b) >= 3:
                back.append(Poly(b, self.normal, self.w))


def _choose_plane(polys, samples=8):
    """The plane that splits the fewest polygons, and balances best.

    csg.js takes the first polygon's plane. On a tube that is arbitrary, and
    the tree comes out as a list: each facet's plane cuts most of the others,
    so the polygon count multiplies at every level and a boolean on a part of
    a few thousand facets stops finishing. Sampling a handful of candidates and
    scoring each costs one pass over the polygons per candidate, and it is the
    difference between a build and a hang.
    """
    # Only polygons that lie in their own plane. Two things stop that being
    # true. A degenerate polygon -- zero area, from a plane grazing a vertex --
    # has a meaningless normal and a meaningless w. And a *non-planar* polygon
    # has a perfectly good normal but vertices that are not all at w, which is
    # what a cut face comes out as when its corners are rounded differently.
    # Either way its own plane does not absorb it: it is pushed to a child, that
    # child picks the same plane, and the build never terminates. Measured: two
    # polygons cycling forever at 10.8 million nodes, and a body built from a
    # cone at 21.7 million. The tolerance has to be the one `split` classifies
    # with, or a polygon can pass here and still fail to be absorbed.
    def usable(p):
        n = p.normal
        if _dot(n, n) <= 0.5:
            return False
        for v in p.pts:
            if abs(_dot(n, v) - p.w) > EPS:
                return False
        return True

    good = [p for p in polys if usable(p)]
    if not good:
        return None
    if len(good) <= 2:
        return (good[0].normal, good[0].w)
    step = max(1, len(good) // samples)
    best, best_score = (good[0].normal, good[0].w), None
    for i in range(0, len(good), step):
        n, w = good[i].normal, good[i].w
        front = back = span = 0
        for q in polys:
            kind = COPLANAR
            for v in q.pts:
                d = _dot(n, v) - w
                kind |= BACK if d < -EPS else (FRONT if d > EPS else COPLANAR)
                if kind == SPANNING:
                    break
            if kind == SPANNING:
                span += 1
            elif kind == FRONT:
                front += 1
            elif kind == BACK:
                back += 1
        score = span * 8 + abs(front - back)
        if best_score is None or score < best_score:
            best, best_score = (n, w), score
    return best


class Node:
    """A BSP node: one plane, the polygons in it, and the two half-spaces."""

    __slots__ = ("plane", "front", "back", "polys")

    def __init__(self, polys=None):
        self.plane = None
        self.front = None
        self.back = None
        self.polys = []
        if polys:
            self.build(polys)

    def invert(self):
        stack = [self]
        while stack:
            node = stack.pop()
            for p in node.polys:
                p.flip()
            if node.plane is not None:
                n, w = node.plane
                node.plane = ((-n[0], -n[1], -n[2]), -w)
            node.front, node.back = node.back, node.front
            if node.front is not None:
                stack.append(node.front)
            if node.back is not None:
                stack.append(node.back)

    def clip_polygons(self, polys):
        out = []
        stack = [(self, polys)]
        while stack:
            node, items = stack.pop()
            if not items:
                continue
            if node.plane is None:
                out.extend(items)
                continue
            front, back = [], []
            for p in items:
                p.split(node.plane, front, back, front, back)
            if node.front is not None:
                stack.append((node.front, front))
            else:
                out.extend(front)
            if node.back is not None:
                stack.append((node.back, back))
            # No back child means the back half-space is inside the solid, and
            # what is inside it is not part of the surface.
        return out

    def clip_to(self, other):
        stack = [self]
        while stack:
            node = stack.pop()
            node.polys = other.clip_polygons(node.polys)
            if node.front is not None:
                stack.append(node.front)
            if node.back is not None:
                stack.append(node.back)

    def all_polygons(self):
        out, stack = [], [self]
        while stack:
            node = stack.pop()
            out.extend(node.polys)
            if node.front is not None:
                stack.append(node.front)
            if node.back is not None:
                stack.append(node.back)
        return out

    def build(self, polys):
        stack = [(self, polys)]
        while stack:
            node, items = stack.pop()
            if not items:
                continue
            if node.plane is None:
                node.plane = _choose_plane(items)
            if node.plane is None:
                # Nothing here has a usable plane -- the rest is degenerate
                # slivers. They belong to this node; recursing on them is what
                # does not terminate.
                node.polys.extend(items)
                continue
            front, back = [], []
            for p in items:
                p.split(node.plane, node.polys, node.polys, front, back)
            if front:
                if node.front is None:
                    node.front = Node()
                stack.append((node.front, front))
            if back:
                if node.back is None:
                    node.back = Node()
                stack.append((node.back, back))


# --------------------------------------------------------------- conversions ---
def _shell_triangles(mesh):
    """Connected components of the mesh, without importing checks."""
    parent = list(range(len(mesh.verts)))

    def find(a):
        while parent[a] != a:
            parent[a] = parent[parent[a]]
            a = parent[a]
        return a

    for tri in mesh.tris:
        r = find(tri[0])
        for j in tri[1:]:
            parent[find(j)] = r
    groups = {}
    for k, tri in enumerate(mesh.tris):
        groups.setdefault(find(tri[0]), []).append(k)
    return list(groups.values())


def _shells_touch(mesh, a, b, samples=2000):
    """Do two shells actually interpenetrate? A few vertices decide.

    Overlapping bounding boxes do not: the column's three helix bands are
    radial strips 120 degrees apart, so their boxes overlap near the axis while
    the bands themselves never meet. What does decide is a vertex of one,
    nudged along the inward normal of a triangle it belongs to, landing inside
    the other. A handful of samples per shell is enough to catch a real
    intersection and cheap enough to run before every boolean.
    """
    normals = {}
    for k in a:
        tri = mesh.tris[k]
        p0, p1, p2 = (mesh.verts[j] for j in tri)
        u = _sub(p1, p0)
        w = _sub(p2, p0)
        n = _cross(u, w)
        length = math.sqrt(_dot(n, n))
        if length < 1e-12:
            continue
        n = (n[0] / length, n[1] / length, n[2] / length)
        for j in tri:
            acc = normals.setdefault(j, [0.0, 0.0, 0.0])
            for d in range(3):
                acc[d] += n[d]
    for j, n in normals.items():
        if samples <= 0:
            break
        samples -= 1
        length = math.sqrt(sum(x * x for x in n))
        if length < 1e-9:
            continue
        v = mesh.verts[j]
        probe = tuple(v[d] - 0.02 * n[d] / length for d in range(3))
        if checks.inside(mesh.verts, mesh.tris, probe, shells=[b],
                         boxes=[_box_of(mesh, b)]):
            return True
    return False


def _box_of(mesh, shell):
    pts = [mesh.verts[j] for k in shell for j in mesh.tris[k]]
    return ([min(p[k] for p in pts) for k in range(3)],
            [max(p[k] for p in pts) for k in range(3)])


def open_edges(tris):
    """Edges used by an odd number of triangles -- the mesh is not closed.

    A boolean needs a closed solid, and `signed_volume() > 0` does not prove
    one: `t.cylinder` builds a tube with no caps, and an open cylinder still
    integrates to the right volume because the caps' contributions cancel.
    Feeding that to a BSP gives a confidently wrong answer, so it is checked.
    """
    count = {}
    for a, b, c in tris:
        for i, j in ((a, b), (b, c), (c, a)):
            key = (i, j) if i < j else (j, i)
            count[key] = count.get(key, 0) + 1
    return [e for e, n in count.items() if n != 2]


def to_polygons(mesh, skip_degenerate=True, require_closed=True):
    """Every triangle as a polygon, with its plane. The mesh must be closed.

    Zero-area triangles are dropped rather than fatal. They carry no volume and
    no surface, and the primitives that sweep an annulus with an inner radius
    of zero -- `t.tube(mesh, 0, r, ...)` -- emit hundreds of them at the axis.
    Dropping them cannot change the solid, because a triangle with no area
    bounds nothing.
    """
    if mesh.signed_volume() <= 0.0:
        raise ValueError("the mesh is inside-out or not closed: signed volume "
                         "is not positive")
    if require_closed:
        # Each mesh must be one solid that does not intersect itself. A mesh of
        # several shells is fine when they are disjoint -- that is a union
        # already -- but when they interpenetrate, a BSP built over the whole
        # soup answers inside/outside wrongly, and the boolean returns a
        # confidently wrong solid. The collector did exactly that: unioning a
        # mesh of four overlapping arms with its junction returned the arms and
        # silently deleted the junction.
        shells = _shell_triangles(mesh)
        if len(shells) > 1:
            boxes = [_box_of(mesh, s_) for s_ in shells]
            for i in range(len(shells)):
                for j in range(i + 1, len(shells)):
                    if all(boxes[i][1][k] >= boxes[j][0][k]
                           and boxes[j][1][k] >= boxes[i][0][k]
                           for k in range(3)) \
                            and _shells_touch(mesh, shells[i], shells[j]):
                        raise ValueError(
                            f"the mesh has {len(shells)} shells and shells "
                            f"{i} and {j} overlap: a boolean needs one solid "
                            f"that does not intersect itself, so union the "
                            f"pieces one at a time instead. Feeding a "
                            f"self-intersecting mesh in returns a wrong solid "
                            f"without complaining")
        gaps = real_gaps(mesh.verts, mesh.tris)
        if gaps:
            raise ValueError(
                f"the mesh is not closed: {len(gaps)} edges are used once and "
                f"have no vertex anywhere along them, so there is a gap in the "
                f"surface, e.g. {gaps[:3]}. A boolean needs a closed solid -- "
                f"`t.tube(mesh, 0, r, ...)` with an inner radius of zero is one "
                f"of the ways to make a hole in one")
    out = []
    for tri in mesh.tris:
        pts = [mesh.verts[i] for i in tri]
        n = _cross(_sub(pts[1], pts[0]), _sub(pts[2], pts[0]))
        if _dot(n, n) < 1e-18:
            if not skip_degenerate:
                raise ValueError("a degenerate polygon has no plane")
            continue
        out.append(Poly(pts))
    return out


def weld(verts, tris, tol=1e-6):
    """Snap vertices together and drop the triangles that collapse."""
    index, remap = {}, []
    for v in verts:
        key = (round(v[0] / tol), round(v[1] / tol), round(v[2] / tol))
        if key not in index:
            index[key] = len(index)
        remap.append(index[key])
    out_verts = [None] * len(index)
    for key, i in index.items():
        out_verts[i] = (key[0] * tol, key[1] * tol, key[2] * tol)
    out_tris = []
    for a, b, c in tris:
        a, b, c = remap[a], remap[b], remap[c]
        if a != b and b != c and a != c:
            out_tris.append((a, b, c))
    return out_verts, out_tris


def repair_tjunctions(verts, tris, tol=1e-6, cell=1.0, snap=1e-5):
    """Put the missing vertices into the long edges that pass through them.

    A boolean's output is geometrically watertight but topologically cracked.
    Where one big face meets two smaller ones, the big face's edge runs from P
    to R while the small ones use P-Q and Q-R: every one of those three edges
    is used by exactly one triangle, so an edge count says the surface is open
    even though the surface has no gap in it. A slicer sees a hairline crack
    and may fill it as a hole, which is the whole point of doing the union.

    Splitting the long edge at Q closes it. The split polygon stays convex --
    the inserted points are collinear with the edge they lie on -- so the fan
    triangulation is still exact.

    `snap` is how far off an edge a vertex may be and still count as a point on
    it. Raising it to the faceting scale -- 0.02 mm -- does close the gaps a
    boolean leaves where a plane cuts a curved surface, and it was measured
    doing so on the scoopula. It also *bends* the surface by up to that much,
    which the closed-form tests catch: a drilled block's volume came out
    0.0012 mm^3 wrong against a tolerance of 0.00004. Exactness is worth more
    than four cosmetic gaps, so it stays at the weld tolerance. Pass a larger
    one only when a mesh is known to be a display artefact rather than a part.

    Some of the triangles this makes have no area: where a face has points
    inserted on two of its edges, the fan from any corner is flat along those
    two edges. They are kept, and they have to be. A flat triangle is the
    counterpart of the neighbour's split edge, and dropping it puts the crack
    straight back -- which is exactly what happened the first time this was
    written, and it showed up as a union that refused its own output for being
    open.
    """
    grid = {}
    for i, v in enumerate(verts):
        key = (int(math.floor(v[0] / cell)), int(math.floor(v[1] / cell)),
               int(math.floor(v[2] / cell)))
        grid.setdefault(key, []).append(i)

    def near_segment(a, b):
        va, vb = verts[a], verts[b]
        d = _sub(vb, va)
        length = math.sqrt(_dot(d, d))
        if length < tol:
            return []
        steps = max(1, int(length / (0.5 * cell)) + 1)
        found = {}
        for s in range(steps + 1):
            t = s / steps
            mid = (va[0] + d[0] * t, va[1] + d[1] * t, va[2] + d[2] * t)
            key = (int(math.floor(mid[0] / cell)),
                   int(math.floor(mid[1] / cell)),
                   int(math.floor(mid[2] / cell)))
            for dx in (-1, 0, 1):
                for dy in (-1, 0, 1):
                    for dz in (-1, 0, 1):
                        for i in grid.get((key[0] + dx, key[1] + dy,
                                           key[2] + dz), ()):
                            if i == a or i == b or i in found:
                                continue
                            vp = verts[i]
                            u = (_dot(_sub(vp, va), d)) / (length * length)
                            if not (tol < u < 1.0 - tol):
                                continue
                            closest = (va[0] + d[0] * u, va[1] + d[1] * u,
                                       va[2] + d[2] * u)
                            off = _sub(vp, closest)
                            if math.sqrt(_dot(off, off)) <= snap:
                                found[i] = u
        return [i for i, _ in sorted(found.items(), key=lambda kv: kv[1])]

    # Only the edges that are used once can be T-junctions; walking all three
    # edges of every triangle costs three times what it needs to on a mesh
    # where the great majority of edges are already paired.
    unpaired = set()
    count = {}
    for tri in tris:
        for i, j in ((tri[0], tri[1]), (tri[1], tri[2]), (tri[2], tri[0])):
            key = (i, j) if i < j else (j, i)
            count[key] = count.get(key, 0) + 1
    for key, n in count.items():
        if n == 1:
            unpaired.add(key)

    out = []
    for tri in tris:
        poly, extra = [], 0
        for k in range(3):
            a, b = tri[k], tri[(k + 1) % 3]
            poly.append(a)
            if ((a, b) if a < b else (b, a)) in unpaired:
                for m in near_segment(a, b):
                    poly.append(m)
                    extra += 1
        if not extra:
            out.append(tri)
            continue
        # Fan from a corner that gives no zero-area triangle. The inserted
        # points are collinear with the edge they sit on, so a fan from a
        # corner adjacent to an insertion produces slivers -- and a sliver is
        # not harmless: it is a face with no area that a slicer still has to
        # make sense of.
        best = None
        for apex in range(len(poly)):
            order = poly[apex:] + poly[:apex]
            fan, ok = [], True
            for k in range(1, len(order) - 1):
                a, b, c = order[0], order[k], order[k + 1]
                if a == b or b == c or a == c:
                    ok = False
                    break
                n = _cross(_sub(verts[b], verts[a]), _sub(verts[c], verts[a]))
                if _dot(n, n) < 1e-16:
                    ok = False
                    break
                fan.append((a, b, c))
            if ok:
                best = fan
                break
            if best is None:
                best = fan
        out.extend(best)
    return verts, out


def real_gaps(verts, tris, tol=1e-6, cell=1.0, snap=1e-5):
    """Edges used once with nothing lying along them: a gap, not a T-junction.

    A boolean's output is full of T-junctions -- an edge used once whose length
    is covered by the edges of two or three neighbours. That is not a hole, and
    refusing it would refuse every boolean this module produces. A gap is an
    edge used once with no vertex anywhere on it, and that is worth refusing.

    `snap` is how far off the edge a vertex may be and still count as lying on
    it, and it has to match what `repair_tjunctions` used: a gap that the repair
    would have closed is not a gap, and judging it with a tighter tolerance than
    the repair would report holes in a surface that has none.
    """
    verts, tris = weld(verts, tris)
    grid = {}
    for i, v in enumerate(verts):
        grid.setdefault((int(math.floor(v[0] / cell)),
                         int(math.floor(v[1] / cell)),
                         int(math.floor(v[2] / cell))), []).append(i)
    count = {}
    for tri in tris:
        for i, j in ((tri[0], tri[1]), (tri[1], tri[2]), (tri[2], tri[0])):
            key = (i, j) if i < j else (j, i)
            count[key] = count.get(key, 0) + 1
    gaps = []
    for (a, b), n in count.items():
        if n != 1:
            continue
        va, vb = verts[a], verts[b]
        d = _sub(vb, va)
        length2 = _dot(d, d)
        if length2 < 1e-18:
            continue
        covered = False
        steps = max(1, int(math.sqrt(length2) / (0.5 * cell)) + 1)
        for s in range(steps + 1):
            t = s / steps
            mid = (va[0] + d[0] * t, va[1] + d[1] * t, va[2] + d[2] * t)
            key = (int(math.floor(mid[0] / cell)),
                   int(math.floor(mid[1] / cell)),
                   int(math.floor(mid[2] / cell)))
            for i in grid.get(key, ()):
                if i == a or i == b:
                    continue
                u = _dot(_sub(verts[i], va), d) / length2
                if not (1e-9 < u < 1.0 - 1e-9):
                    continue
                closest = (va[0] + d[0] * u, va[1] + d[1] * u,
                           va[2] + d[2] * u)
                off = _sub(verts[i], closest)
                if _dot(off, off) <= snap * snap:
                    covered = True
                    break
            if covered:
                break
        if not covered:
            gaps.append((a, b))
    return gaps


def to_mesh(polys, name="union", repair=True):
    """Triangulate by fan -- the polygons are convex, so a fan is exact."""
    verts, tris = [], []
    for p in polys:
        base = len(verts)
        verts.extend(p.pts)
        for k in range(1, len(p.pts) - 1):
            tris.append((base, base + k, base + k + 1))
    verts, tris = weld(verts, tris)
    tris = dedupe_faces(verts, tris)
    if repair:
        # Repair first, then merge. Merging first and repairing after breaks
        # the connectivity: the repair's fans introduce vertices the merged
        # faces do not use, and a box union that was one shell came back as
        # three. Merging last leaves the surface connected -- the merged faces
        # keep every boundary vertex -- and the volume is unchanged to 1e-6.
        verts, tris = repair_tjunctions(verts, tris)
        verts, tris = merge_coplanar(verts, tris)
        # And again: merging removes interior vertices, and one of those can be
        # sitting on a *neighbouring* patch's edge. That leaves an edge used
        # once -- 1.1 mm of it on the drilled-block test -- which a second
        # repair closes.
        verts, tris = repair_tjunctions(verts, tris)
    mesh = t.Mesh(name)
    mesh.verts = verts
    mesh.tris = tris
    return mesh


def _combine_polys(pa, pb, how):
    """The three operations, straight out of csg.js, on polygon lists."""
    na, nb = Node(pa), Node(pb)
    if how == "union":
        na.clip_to(nb)
        nb.clip_to(na)
        nb.invert()
        nb.clip_to(na)
        nb.invert()
        na.build(nb.all_polygons())
        return na.all_polygons()
    if how == "subtract":
        na.invert()
        na.clip_to(nb)
        nb.clip_to(na)
        nb.invert()
        nb.clip_to(na)
        nb.invert()
        na.build(nb.all_polygons())
        na.invert()
        return na.all_polygons()
    if how == "intersect":
        na.invert()
        nb.clip_to(na)
        nb.invert()
        na.clip_to(nb)
        nb.clip_to(na)
        na.build(nb.all_polygons())
        na.invert()
        return na.all_polygons()
    raise ValueError(f"unknown operation {how!r}")


def combine(operations, name="solid"):
    """A chain of booleans that repairs once, at the end.

    `operations` is a list of ("union" | "subtract" | "intersect", mesh), the
    first of which is the solid to start from. Every intermediate stays as
    polygons, so the mesh is welded and its T-junctions closed exactly once --
    on the answer. Doing it between steps multiplies the polygon count that the
    next boolean has to chew through, and the stopcock's barrel is the part
    where that stopped being an optimisation and became the difference between
    a build and a hang.
    """
    polys = None
    for op, mesh in operations:
        if polys is None:
            if op != "start":
                raise ValueError("the first operation must be 'start'")
            polys = to_polygons(mesh)
            continue
        polys = _combine_polys(polys, to_polygons(mesh), op)
    if polys is None:
        raise ValueError("a combine needs at least one mesh")
    out = to_mesh(polys, name)
    out.name = name
    return out


def union(*meshes, name="union"):
    """One closed solid covering exactly the material of all of them.

    The intermediates stay as polygons. Building a mesh between steps would
    weld and repair T-junctions on a result that is about to be cut up again,
    which is both wasted work and a much worse input: the subtract of a tube
    from a block goes from 512 triangles to 6260 after a repair, and the next
    boolean then has twelve times as much to do.
    """
    if not meshes:
        raise ValueError("a union needs at least one mesh")
    polys = to_polygons(meshes[0])
    for m in meshes[1:]:
        polys = _combine_polys(polys, to_polygons(m), "union")
    out = to_mesh(polys, name)
    out.name = name
    return out


def subtract(a, b, name="subtract"):
    """A with B's material taken out of it."""
    return to_mesh(_combine_polys(to_polygons(a), to_polygons(b), "subtract"),
                   name)


def intersect(a, b, name="intersect"):
    """Only the material both of them have."""
    return to_mesh(_combine_polys(to_polygons(a), to_polygons(b), "intersect"),
                   name)


def _plane_of(verts, tri):
    a, b, c = (verts[j] for j in tri)
    n = _cross(_sub(b, a), _sub(c, a))
    length = math.sqrt(_dot(n, n))
    if length < 1e-12:
        return None
    n = (n[0] / length, n[1] / length, n[2] / length)
    return n, _dot(n, a)


def _boundary_loop(comp, verts, tris):
    """The single loop round a patch of coplanar triangles, or None.

    Edges used once inside the patch are its boundary; a simple patch's
    boundary is one closed loop, and walking it gives the vertices in the
    triangles' own winding. Anything else -- a hole, two loops, a pinch -- is
    not something this can describe, and the caller keeps the triangles.
    """
    used = {}
    for k in comp:
        tri = tris[k]
        for i in range(3):
            a, b = tri[i], tri[(i + 1) % 3]
            key = (a, b) if a < b else (b, a)
            used[key] = used.get(key, 0) + 1
    nxt = {}
    for k in comp:
        tri = tris[k]
        for i in range(3):
            a, b = tri[i], tri[(i + 1) % 3]
            key = (a, b) if a < b else (b, a)
            if used[key] != 1:
                continue
            if a in nxt:
                return None                 # two ways out of one vertex
            nxt[a] = b
    if len(nxt) < 3:
        return None
    start = next(iter(nxt))
    loop, cur = [start], nxt[start]
    while cur != start:
        if cur not in nxt or len(loop) > len(nxt):
            return None
        loop.append(cur)
        cur = nxt[cur]
    if len(loop) != len(nxt):
        return None
    return loop


def _triangulate(loop, verts, normal):
    """Fan a convex loop, ear-clip anything else. Indices out."""
    seed = (1.0, 0.0, 0.0) if abs(normal[0]) < 0.9 else (0.0, 1.0, 0.0)
    e1 = _sub(seed, tuple(_dot(seed, normal) * x for x in normal))
    m = math.sqrt(_dot(e1, e1))
    e1 = tuple(x / m for x in e1)
    e2 = _cross(normal, e1)
    pts = [(_dot(verts[j], e1), _dot(verts[j], e2)) for j in loop]
    n = len(loop)

    def cross2(o, a, b):
        return ((pts[a][0] - pts[o][0]) * (pts[b][1] - pts[o][1])
                - (pts[a][1] - pts[o][1]) * (pts[b][0] - pts[o][0]))

    convex = all(cross2(i, (i + 1) % n, (i + 2) % n) >= -1e-12
                 for i in range(n))
    if convex:
        return [(loop[0], loop[k], loop[k + 1]) for k in range(1, n - 1)]

    def inside(p, a, b, c):
        d1 = cross2(a, b, p)
        d2 = cross2(b, c, p)
        d3 = cross2(c, a, p)
        return (d1 >= -1e-12 and d2 >= -1e-12 and d3 >= -1e-12)

    idx = list(range(n))
    out = []
    guard = 0
    while len(idx) > 3 and guard < 4 * n:
        guard += 1
        for i in range(len(idx)):
            a = idx[i - 1]
            b = idx[i]
            c = idx[(i + 1) % len(idx)]
            if cross2(a, b, c) <= 1e-12:
                continue                    # reflex or flat
            if any(inside(j, a, b, c) for j in idx
                   if j not in (a, b, c)):
                continue
            out.append((loop[a], loop[b], loop[c]))
            idx.pop(i)
            break
        else:
            return None                     # no ear found: leave it alone
    if len(idx) == 3:
        out.append((loop[idx[0]], loop[idx[1]], loop[idx[2]]))
    return out or None


def dedupe_faces(verts, tris):
    """Drop triangles that another triangle covers exactly.

    When two solids share a boundary plane -- a handle and a scoop both standing
    on the bed, a plate and a cradle both starting at z=0 -- their faces in that
    plane overlap, and a union keeps *both*: the coplanar rule sorts them by
    normal direction, and two down-facing faces sort the same way. What comes
    out is a surface with every edge used twice in the same place, which the
    T-junction repair cannot close and `real_gaps` reads as holes.

    A face is a duplicate if the same three vertices bound another triangle.
    The first one is kept.
    """
    seen = set()
    out = []
    for tri in tris:
        key = tuple(sorted(tri))
        if key in seen:
            continue
        seen.add(key)
        out.append(tri)
    return out


def merge_coplanar(verts, tris, tol=1e-6, min_gain=2):
    """Merge coplanar triangles that share an edge into single faces.

    A boolean splits every face it crosses, so a tube's facets come back as
    tens of thousands of triangles lying in the same few planes: the surface is
    the same and the mesh is thirty times bigger, which is what a run of the
    column's union does to memory. Merging each coplanar patch back into one
    face and re-triangulating it returns the surface with the interior vertices
    gone.

    A patch whose boundary is not a single simple loop -- a ring, a patch with a
    hole -- is left as it was. `min_gain` is how many triangles a merge has to
    save before it is worth doing.
    """
    # No weld here: the caller has already welded, and welding again drops
    # triangles whose vertices coincide -- which the repair's fans produce, and
    # which are what closes a T-junction. That alone turned a one-shell box
    # union into nine.
    groups = {}
    flat = []
    for k, tri in enumerate(tris):
        pl = _plane_of(verts, tri)
        if pl is None:
            # A triangle with no area belongs to no plane, and it is *kept*:
            # dropping it puts a hole in the surface. The repair's fans are
            # full of them -- flat triangles are what closes a T-junction --
            # and dropping them turned a one-shell box union into three, and
            # then nine.
            flat.append(k)
            continue
        n, w = pl
        key = (round(n[0], 4), round(n[1], 4), round(n[2], 4), round(w, 4))
        groups.setdefault(key, []).append(k)

    out = [tris[k] for k in flat]
    for key, idxs in groups.items():
        if len(idxs) < 2:
            out.extend(tris[k] for k in idxs)
            continue
        # Connected patches within the plane, joined by shared edges.
        parent = {k: k for k in idxs}

        def find(a):
            while parent[a] != a:
                parent[a] = parent[parent[a]]
                a = parent[a]
            return a

        edge_owner = {}
        for k in idxs:
            tri = tris[k]
            for i in range(3):
                a, b = tri[i], tri[(i + 1) % 3]
                key2 = (a, b) if a < b else (b, a)
                if key2 in edge_owner:
                    parent[find(k)] = find(edge_owner[key2])
                else:
                    edge_owner[key2] = k
        comps = {}
        for k in idxs:
            comps.setdefault(find(k), []).append(k)
        for comp in comps.values():
            if len(comp) < min_gain + 1:
                out.extend(tris[k] for k in comp)
                continue
            loop = _boundary_loop(comp, verts, tris)
            fan = _triangulate(loop, verts, key[:3]) if loop else None
            if not fan or len(fan) >= len(comp):
                out.extend(tris[k] for k in comp)
            else:
                out.extend(fan)
    # Drop the vertices the merge left behind. Merging removes interior
    # vertices, and `Mesh.shells()` groups *vertices*: an unreferenced one
    # becomes a shell of its own, which turned a sound box union into nine
    # shells of one vertex each.
    used = {j for tri in out for j in tri}
    if len(used) != len(verts):
        remap = {}
        keep = []
        for j in range(len(verts)):
            if j in used:
                remap[j] = len(keep)
                keep.append(verts[j])
        out = [(remap[a], remap[b], remap[c]) for a, b, c in out]
        verts = keep
    return verts, out


def gap_stats(verts, tris):
    """Count, median, worst and total length of the real gaps in a surface.

    A count says nothing about whether a surface is printable. A boolean
    between faceted surfaces leaves gaps where their facet planes disagree, and
    those are the faceting error -- a hundredth of a millimetre, far below what
    a nozzle resolves. A gap of two millimetres is a hole. Reporting the sizes
    is what lets a caller tell one from the other instead of counting them
    equally.
    """
    lengths = []
    for a, b in real_gaps(verts, tris):
        d = _sub(verts[a], verts[b])
        lengths.append(math.sqrt(_dot(d, d)))
    if not lengths:
        return {"count": 0, "median": 0.0, "worst": 0.0, "total": 0.0}
    lengths.sort()
    return {"count": len(lengths), "median": lengths[len(lengths) // 2],
            "worst": lengths[-1], "total": sum(lengths)}


def crack_length(verts, tris):
    """Total length of the edges used by only one triangle.

    A boolean between two *curved* surfaces cannot close every edge: the
    intersection curve is a polyline on each surface and the two polylines only
    agree where their planes happened to coincide. The gaps are the faceting
    error -- micrometres at the resolutions this repo uses -- and what matters
    is their length, not their existence. Planar geometry should come out at
    zero.
    """
    total = 0.0
    for i, j in open_edges(tris):
        d = _sub(verts[i], verts[j])
        total += math.sqrt(_dot(d, d))
    return total


def volume(mesh):
    return mesh.signed_volume()


def overlap_volume(a, b):
    """How much material the two have in common, exactly."""
    return intersect(a, b).signed_volume()

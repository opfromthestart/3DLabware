# Scripted 3D design — a guide for LLMs (and humans) generating meshes

This is the working reference for producing printable 3D models from
parametric code. It is built around `threelib.py` in this directory, but the
reasoning generalises to any triangle-mesh toolchain.

Everything here was learned the hard way on a real task. The failure modes
in [§8](#8-pitfalls-that-will-bite-you) are all things that pass silently
and produce a file that *looks* fine and slices wrong.

## 0. Models in this repo

| Script | Part | Verified figures |
|---|---|---|
| `make_grid.py` | 150×150 mm two-layer grid | 300 disjoint bars, 4530.000 mm³, 22,650 mm² bed |
| `make_grid_shell.py` | grid floor in a 110 mm bore | 147 shells, disc trimmed flush to the bore |
| `make_filter.py` | Buchner filter cartridge | 70 mm disc, 30 staggered layers, 10 mm headroom, hollow funnel |
| `make_holder.py` | holder / funnel cradle | 90 mm sheet, 18 mm hole, 10.4 mm cradle, 1 mm seat |
| `make_elbow.py` | distillation elbow, standard taper joints | 14/20 male+female, 200+200 mm arms, 105° and 75° |
| `make_wye.py` | three-way adapter, standard taper joints | 14/20 all male, 45 mm arms, 75/180/105°, all three bores joined |
| `make_collector.py` | fraction collector (rotating "cow") | 14/20 female inlet on the rotation axis, 3 male outlets 155° off the inlet and 43° apart, all four bores joined |
| `make_column.py` | fractionating column, four staging styles | 16 mm bore × 150 mm staged, ground joints both ends, `cones`/`helix`/`vigreux`/`tesla`, vapour path proved open |
| `make_stopcock.py` | two-way straight-through stopcock | 4 mm tubes collinear, 9 to 12 mm bore on 1:20, key 0.08 mm clear and seating after 1.6 mm, passage proved open at 0° and blocked by the key at 90° |
| `make_pfunnel_body.py` / `_cap.py` / `_screw.py` | powder-addition funnel, three parts | 50/35 mouth + 29/42 stem, internal 3.5 mm equalizer bore proved open tip-to-headspace, 14 mm auger tube over a 7 mm throat, two-start 42 mm-pitch flight 0.45 mm clear, cap seats on a common 1:20 generatrix over 24.9 mm of ground |
| `make_stopper.py` | ground stopper, solid, disc handle | 14/25 male solid cone, 13.2 mm disc on edge inside the taper's implied cone, fused 2.5 mm deep over 18 mm², prints with no support |
| `make_aspirator.py` | water aspirator (Venturi vacuum pump) | 5.2 mm supply bore squeezed to a 2.6 mm throat that starts before, and ends after, the 5.2 mm air inlet; the section at the inlet is *two rings* — a 2.6 mm jet in a 5.0 mm nozzle, 11.5 mm of air annulus, inside a 28 mm chamber 30 mm tall; the 9.2→32 mm shoulder leans out at 55° *ahead* of the inner cone, and the join is one continuous band; `JET_RELEASE = 0` so the water only ever narrows, and the outlet starts 32 mm above where the inlet ends; 36×32×130 mm, 17.8 cm³, 1 shell, 0 gaps, stands on a 45 mm² ring (brim it), only the boss wants removable support |

| `make_scoopula.py` | scoopula: round channel **in line with** a flat handle | 160 x 17.0 x 9.6 mm; 62 mm channel, 14 mm round mouth 8 mm deep over a 5.9 mm-radius tangent arc, 1.4 mm wall measured *across itself* at every height, 1.6 mm floor, open for 46 mm of its length behind a 3.1 mm back wall and running out to a 6.2 x 1.6 mm nose, no roof over the channel at 60 upward columns, heel solid at carrying height, 1201 mm^2 of bed, steepest exposed down-facing face 70.2 deg from horizontal, 3.49 cm^3, 1 shell, 0 gaps |
Each has a matching verifier (`verify_grid.py`, `verify_shell.py`,
`verify_parts.py`, `verify_elbow.py`, `verify_wye.py`, `verify_collector.py`,
`verify_column.py`, `verify_stopper.py`, `verify_stopcock.py`,
`verify_aspirator.py`, `verify_weighboat.py`,
`verify_scoopula.py`) that re-parses the written 3MF from disk. The filter and holder are verified *together*,
including the fit between them — see [§7](#7-verification-discipline).

**Measured glass.** The printed joints above mate with each other at the
nominal table in `junctionlib.JOINTS`, but real ware in the cupboard runs
loose against it: the user's 14/20 glass fits best at the setting
**"14.5/22"** (mouth 14.5 mm, ground 22 mm) rather than "14/20". Any
joint spec parses as a `"<mouth>/<ground length>"` string
(`junctionlib.parse_joint`), so a part aimed at that glass should stamp
`joint="14.5/22"` in its PARAMS; a fit verifier comparing against real
dimensions should use those, not the nameplate.

---

## 1. What the file formats actually are

You do not need a CAD kernel or `lib3mf`. Both common formats are trivial to
write with the standard library.

**3MF is a ZIP of XML.** Three parts, and all three matter:

| Part | Purpose |
|---|---|
| `[Content_Types].xml` | declares MIME types per extension |
| `_rels/.rels` | points at the model part |
| `3D/3dmodel.model` | the geometry |

The model part is a single `<mesh>` of `<vertex>` and `<triangle>` lists:

```xml
<model unit="millimeter" xml:lang="en-US"
       xmlns="http://schemas.microsoft.com/3dmanufacturing/core/2015/02">
  <resources>
    <object id="1" type="model">
      <mesh>
        <vertices>
          <vertex x="0" y="0" z="0"/>
          ...
        </vertices>
        <triangles>
          <triangle v1="0" v2="1" v3="2"/>
          ...
        </triangles>
      </mesh>
    </object>
  </resources>
  <build>
    <item objectid="1"/>
  </build>
</model>
```

Traps that produce files that look valid but will not open:

- `unit` must be `millimeter` (or `millimeter` with a scale). If you write
  coordinates in mm and declare the unit as `inch`, the part is 25.4× wrong.
- `objectid` in `<build>` must match an `<object id>`.
- Element order is fixed: `resources` then `build`; `vertices` then
  `triangles`.
- The namespace URI is not optional.

**STL** is simpler still (a triangle list with optional normals) but carries
no units, no colours, and no assembly structure. Use 3MF by default.

---

## 2. Mental model

- A model is a **triangle soup**. No shared edges are required.
- A **shell** is a set of triangles connected through shared vertices. A box
  is one shell of 8 vertices and 12 triangles. A cylinder wall is one shell.
- **Winding** decides which way a face points. Counter-clockwise seen from
  outside = outward normal. Getting this backwards makes every slicer treat
  your solid as inside-out.
- Signed volume is your global sanity check: for correctly wound meshes it
  is positive and equals the real volume.
- There is **no boolean engine** here. You combine shapes by letting them
  overlap and trusting the slicer to union them.

That last point is the single biggest decision in this workflow, so treat it
deliberately — see [§6](#6-combining-shapes-overlap-or-not).

---

## 3. Designing for the printer

The mesh is not the part. The part is whatever a nozzle of a given diameter
can build, one layer at a time. Most bad prints are not mesh bugs — they are
perfectly valid meshes that asked the machine for something it cannot
represent.

`python3 check_printability.py <file.3mf> [nozzle_mm] [layer_mm]` runs the
checks below and reports them for a specific process.

**Design in the orientation it will be printed.** Build Z-up with the part
flat on `z = 0`. Rotating a mesh afterwards does not rotate the layer lines,
and re-orienting can turn a self-supporting part into one full of overhangs.
If a part needs a different orientation, redesign it — don't rotate.

**Respect the layer grid.** Every `z` coordinate should be a multiple of the
layer height. A feature straddling a layer boundary gets sliced across, so its
real dimension depends on whichever layer height you happened to pick. The
grid here sits at `z = 0.2 / 0.4`, which is exactly 1 and 2 layers at 0.2 mm
— but lands **off-grid at 0.3 mm layers**, 0.1 mm out. If a part is meant to
print at several layer heights, quantise its `z` values to the coarsest one
you intend to support.

**Minimum feature sizes** (FDM rules of thumb, not laws):

| Feature | Rule of thumb |
|---|---|
| lateral wall — screen, lattice, light fit | >= 1x nozzle (one bead is normal) |
| lateral wall — load-bearing | >= 2x nozzle |
| holes | >= 3x nozzle (they print undersize) |
| any wall in Z | >= 2 layers |

The grid in this repo is 0.5 mm bars on a 0.4 mm nozzle — **1.25 beads**, i.e.
one extrusion per bar. That is exactly right for a screen or lattice and
prints cleanly; do not treat it as a defect. Two caveats are worth stating
rather than panicking about: the bars won't carry meaningful load, and their
real width follows your extrusion calibration, so don't rely on 0.5 mm as a
tolerance. The same model *does* become unprintable on a 0.6 mm nozzle, where
0.5 mm is narrower than a single bead.

**Don't confuse "below the nozzle" with "below two nozzles."** These are
completely different verdicts and conflating them produces false alarms on
perfectly good models. A sub-nozzle feature cannot be extruded at all; a
single-bead feature is normal. `check_printability.py` reports them at
different severities for exactly this reason.

**Holes come out smaller than you draw them.** Material shrinks as it cools,
and the nozzle squishes into the opening. Design holes 0.1–0.2 mm oversize
for a nominal fit. The grid's 0.5 mm holes will probably measure 0.4 mm or
less — fine for a screen or filter, wrong if you need a precise aperture.

**Overhangs.** A downward-facing surface shallower than the printer's
overhang angle needs support. Vertical walls and steep arches are free, and
the angle is a property of the printer, not of the part: the usual 45° rule is
conservative, and a well-tuned machine holds far shallower. It is
`printability.OVERHANG_DEG`, settable with `OVERHANG_ANGLE=20` or per call,
and every verifier that looks at overhangs takes `--angle 20` and says in the
check which value it used — a check that quietly assumed 45 would pass a part
that will not print, and a build that quietly assumed 45 would refuse one that
will.

**The 45° rule says nothing about horizontal faces.** A bar bridging a
0.5 mm gap and a 70 mm ceiling have the *same* normal, `-Z`, so any
normal-based test flags both. They are not remotely equivalent. So
`printability.overhangs` tags every face:

- `sloped` — a real overhang, reported as **WARN**
- `flat` — a horizontal ceiling, reported as **INFO** with the caveat that a
  lattice bridges it but a wide ceiling does not

Apply the same scepticism to how you *measure* an overhang. A face is
"resting on the bed" only if the whole face is in the first layer; test its
**highest** point, not its lowest. Testing the lowest point silently hides
every overhang on a part that simply touches the bed at its base — a
10.4 mm conical wall has all 912 of its downward faces at `zmin = 0`, and a
`zmin`-based test reports a clean part that in fact needs a conical skirt of
support. That single comparison is the whole check.

**Print orientation is a design decision, not a rendering preference.** Two
choices in this repo show the trade-off:

- The **filter** is modelled cone-down but printed grid-side down. Printed
  as drawn, its cone would narrow *upward* and leave its bore ceiling
  unsupported; flipped, every layer is smaller than the one below, and the
  bore is a straight cylinder. The same part is unprintable in one
  orientation and support-free in the other.
- The **holder's** cradle is an upward-opening funnel. Its wall necessarily
  flares outward as z rises, so it reports **304 sloped overhang faces** and
  genuinely needs support. That is inherent to the shape, not a modelling
  error — make the outer wall a vertical cylinder and it becomes
  support-free, at the cost of a thick, heavy base. Report the requirement;
  do not quietly redesign around it.

**Bridges and membranes.** A thin floor spanning holes prints fine as long as
the pattern is *continuous* — each new layer is supported by the one below.
Isolated islands of material bridging a gap are the risky case. The grid's
0.4 mm floor is a continuous lattice so it bridges naturally, but a 0.4 mm
membrane spanning a 110 mm bore is mechanically fragile in *use*, however
well it prints.

**Stagger stacked layers.** Bars in a layer directly above the same bars
below bridge a whole pitch. Shift each same-direction layer by half a pitch
and every bar lands over the *gaps* instead, halving the span — for this
filter, 1.0 mm down to 0.25 mm. `add_bar_layer(..., shift=...)` does it:

```python
stagger = STAGGER * ((i // 2) % 2) * pitch    # 0, 0.5p, 0, 0.5p, ...
```

Test it from the mesh, not from the generator's own arithmetic: for each
pair of same-direction layers, assert every bar centre in one equals a bar
centre in the other plus the shift. Skip the outermost bars — a staggered
layer legitimately gains one bar at an edge with no partner.

Staggering creates a second problem at a round boundary. Bars now hang half
off the disc, and the offcut is a sliver 0.25 mm wide — narrower than a 0.4 mm
nozzle, so it prints as a blob or not at all. `min_width=` rejects a clipped
piece that is narrower than a full bar *along the bar's narrow axis only*, so
a bar clipped short along its length (the normal case) is kept.

**Quantise every z to the layer grid.** A collar top at 7.5 mm is 37.5 layers
at 0.2 mm: it gets sliced across, so its real height becomes whatever the
layer height happens to be. Snap it up — `ceil(h / layer) * layer` — and let
the realised value be what goes into the parameters, not the one you asked
for. Here that turns a requested 1.5 mm rise into an actual 1.6 mm.

And test "is this a whole number of layers" against the **nearest integer**,
not with `% 1.0`: `7.6 / 0.2` is `37.99999999999999`, so the modulo form
reports `1.0` for a perfectly good layer boundary.

**A follower must have the same profile shape as what it follows, measured
from the same datum.** The holder's cradle is meant to sit a constant 1.5 mm
off the filter's cone. It did not, three times over, in three different ways:

1. It restated the cone instead of asking for it, and measured the taper from
   the sheet rather than from where the taper starts. The cradle then crossed
   the filter's surface -- interpenetrating low down, standing 10 mm clear at
   the rim, touching nowhere.
2. Once it followed the real surface, it was still *one straight taper* while
   the filter's profile has a **knee** where a straight spout meets the cone.
   A taper that matches at both ends bulges away from a piecewise profile in
   between -- 13 mm of gap along the middle.
3. With the shape right, the datum was still wrong: the cradle measured its
   bore from the sheet's underside while the filter begins at the sheet's top
   face, so it accumulated a whole extra taper of gap (1.5 mm grew to 2.5 mm)
   before the filter started.

Each fix was invisible in the numbers the part printed about itself, because
those numbers came from the same formula. Measure a fit from the **meshes**,
at every height, and assert the spread -- a "constant clearance" claim is only
true if the spread is zero.

**An overrun must stay inside the wall that hides it.** Shapes are overlapped
rather than boolean-unioned, so the funnel is extended up into the cup to
avoid coincident faces. The funnel flares going down, so that overrun moves
outward as well as up: harmless through a 4 mm cup wall, but through a 1.5 mm
one it pushed 1.3 mm of funnel out through the side of the part. Make the
overrun conditional on the room available, and when there is none, butt the
two identical cross-sections instead.

**Let each part be its own shape, and let the verifier do the coupling.** A
holder derived from the filter by import looks elegant and fails badly: change
one number in the filter and the holder silently changes too, so a mismatch
can never be detected -- the holder is always "correct" by construction. Two
standalone scripts with their own parameters, plus a verifier that imports
both and checks the fit, means a change to either one shows up as a failed
check instead of a part that does not fit.

**A lofted part has vertices only where you put stations.** Measuring a
swept solid by sampling its vertices is measuring your own station list, not
the shape. Two verifier checks were wrong for this reason before either was
right: the ground cone of a joint has vertices at its mouth and its tip and
nowhere between, so an intermediate-depth sample found nothing at all and
reported `nan`; and reading the taper off the endpoints cannot tell a straight
taper from a curved one. Cut the mesh with a plane instead
(`verify_elbow.section`) and fit the result -- the fit then also tests that the
taper is *straight*, which endpoint sampling never could.

Keep the cuts off the faces you are measuring the ends of. A cut coplanar with
the mouth annulus has both endpoints of every edge at zero depth, and a strict
crossing test counts none of them, so the measurement at depth 0 is empty.

**Identify a feature by its shape, not its size.** The mouth annuli were found
by taking the largest planar faces. That is wrong: a straight tube's facet is a
single long planar quad -- 175 mm by half a millimetre -- which beats the mouth
annulus on area, so the check picked facets on opposite sides of the tube and
reported an "included angle" of exactly 180 degrees. An annulus is recognisable
because all its vertices sit at exactly two radii from its centre; rank on that.

**A verifier should read the part's parameters out of the part.** Comparing a
file against the module's globals only works for the default variant. The
second elbow, built at 75 degrees from a module defaulting to 105, failed every
angle and length check while being perfectly correct. The writer already stamps
every parameter into the 3MF, so verify against the stamped values -- and
expect strings there to carry their `repr` quotes.

**Bed contact.** The first layer needs area to stick. Both models are fine
(22,650 mm² and 10,568 mm²). A part touching the bed only on thin ribs will
curl or release mid-print.

**Slice before you commit.** `printability.py` catches the gross classes. A
real slicer preview catches what it cannot: travel moves, support placement,
infill interaction, and thin walls merging with each other.

---

## 4. The toolkit

### `threelib.py` — primitives, clipping, IO

```python
import threelib as t

mesh = t.Mesh("part")
t.tube(mesh, r_in=55, r_out=58, z0=0, z1=10, n=172)   # annulus / pipe wall
t.box(mesh, 0, 0, 0, 10, 5, 2)                       # axis-aligned box
t.cylinder(mesh, 10, 0, 20, n=128)
t.cone(mesh, 10, 0, 0, 30, n=128)                     # r_top=0 -> point
t.sphere(mesh, 10, cz=5)
t.torus(mesh, 20, 5)
t.wedge(mesh, 0, 0, 0, 10, 10, 10, axis="x")         # ramp / triangular prism
t.prism(mesh, polygon, z0, z1)                        # extrude a convex CCW polygon
t.revolve(mesh, profile, n)                           # solid of revolution from (r, z)

t.write_3mf("part.3mf", (mesh.name, mesh), title="Part")
t.write_stl("part.stl", mesh)                         # binary by default
```

`revolve` is the one for solids rather than tubes: `profile` is the outline from
the axis outward and back, `[(0, z0), (r, z0), (r1, z1), (0, z1)]`, and a point
at radius 0 becomes a single vertex on the axis. Sweeping an annulus with an
inner radius of zero instead emits three degenerate triangles for every real
one, and its signed volume is right, which is exactly the sort of thing that
passes a check and clogs a slicer.

Clipping and tessellation — the tools for curved boundaries:

```python
n = t.auto_segments(radius, tol=0.01)   # segments so chord error <= tol
disc = t.disc_polygon(radius, n)         # CCW regular n-gon, reusable
kept = t.clip_convex(subject_poly, disc) # Sutherland-Hodgman, convex clip
kept = t.clean_polygon(kept)             # drop dupes/collinear points
```

Measurement and structure:

```python
mesh.bbox()              # (minx, miny, minz, maxx, maxy, maxz)
mesh.signed_volume()     # positive iff all shells face outward
mesh.max_radius_xy()     # for radial parts
mesh.shells()            # connected components -> one per solid piece
mesh.area_xy()           # projected silhouette area
```

### `csglib.py` — booleans, so a joint can be a joint

```python
import csglib as csg

solid = csg.union(a, b, c)              # one closed solid covering all three
bored = csg.subtract(solid, hole)       # a with the hole's material gone
both  = csg.intersect(a, b)             # only what they share
solid = csg.combine([("start", base), ("subtract", box), ("union", cap)])

csg.overlap_volume(a, b)                # how much they have in common
csg.open_edges(tris)                    # edges used once -- T-junctions or gaps
csg.real_gaps(verts, tris)              # the subset that are gaps, not T-junctions
csg.crack_length(verts, tris)           # total length of those
```

A BSP boolean (Naylor/Amanatides/Thibault, the shape `csg.js` popularised). Every
face in this repo is planar -- curves are faceted -- so it is exact rather than
approximate, which matters when a fit is 0.08 mm. `combine` exists so a chain of
operations repairs its T-junctions once, at the end: repairing between steps
turned a 512-triangle tube into 6260 before the next boolean had a go at it.

Two things it will not do. It will not touch a mesh that is not closed, and it
will not close a crack between two *curved* surfaces where their facet planes
disagree -- that gap is the faceting error, reported by `crack_length` rather
than hidden.

### `check_overlaps.py` — which parts still need one

```bash
python3 check_overlaps.py part.3mf [part2.3mf ...]
```

Reports shell pairs that interpenetrate: a vertex of one, nudged into its own
solid along the inward normal, landing inside the other. A slicer that fills by
winding unions those shells; one that fills even-odd shows the joint as a hole,
which is a hole in the print.

### `gridlib.py` — the layered grid pattern

```python
import gridlib as g

spans = g.bar_spans(150.0, width=0.5, gap=0.5)
# -> 150 flush, symmetric intervals covering exactly 0..150

stats = g.add_grid(mesh, size=150, width=0.5, gap=0.5, height=0.2,
                   z0=0.0, center=(0, 0), clip_poly=disc)
# vertical bars   z0        .. z0+height
# horizontal bars z0+offset .. z0+offset+height   (offset defaults to height)
```

The two-layer layout is the important part. Bars in one direction sit at
`z=0..0.2`, bars in the other at `z=0.2..0.4`, so **nothing interpenetrates**.

### `checks.py` — verification helpers

```python
import checks as K

rep = K.Report("Part")
objs, build_ids, unit = K.load_3mf("part.3mf", rep)   # re-parses the file
shells = K.shells(objs[0]["verts"], objs[0]["tris"])
K.aabb_overlap_pairs(shells)      # pairwise overlap, for box-like solids
K.render(v, t, "top.png", view="xy", su=10)
K.scanline_runs(img, y_world, v_origin, sv)
```

### `printability.py` — process limits

```python
import printability as P

findings, summary = P.assess(verts, tris, shells,
                             nozzle=0.4, layer_height=0.2)
```

Checks overhangs against the printer's angle, laterally thin shells against
the nozzle, `z` values against the layer grid, and first-layer bed contact.
These are heuristics and they are documented as such in the module — but they
catch the failures that produce a file which verifies perfectly and still will
not print. See [§3](#3-designing-for-the-printer).

The overhang angle is `OVERHANG_DEG`, default 45 and settable with
`OVERHANG_ANGLE=20` in the environment or `max_angle=` per call. Verifiers
that look at overhangs take `--angle 20` and name the value they used in the
check; parts whose geometry depends on it stamp it into the file as a
parameter (`make_column.py`'s `overhang`), so a build can be re-checked
against the printer it was actually shaped for.

---

## 5. Workflow

1. **Restate the spec as numbers.** Write down every dimension with its unit
   and convert to a derived quantity: inner radius, pitch, total height.
2. **Ask about anything that changes topology.** Ambiguities that alter the
   shape are worth one question; see [§9](#9-when-to-ask-the-user).
3. **Build with closed forms in mind.** Before writing geometry, work out the
   answer analytically — total bar count, material area, expected volume.
   You cannot check a model against a number you never derived.
4. **Generate, then verify the artifact** by re-reading the file.
5. **Render and inspect** — but see the caveat about rasterised checks.
6. **Report the numbers**, not just "done": dimensions, bar counts, volume,
   and any assumption you made.

### 5.1 Output names must encode the parameters

A script with a hardcoded `grid_150mm.3mf` will silently overwrite an earlier
model the moment you edit a parameter, and the only evidence that anything
changed is the missing file. Worse, a name that lists *some* parameters
(`filter_70mm.3mf`) still collides as soon as an unlisted one moves.

Hash the whole parameter set instead:

```python
PARAMS = dict(grid_dia=70.0, layers=30, bar_w=0.5, bar_gap=0.5, ...)

def default_out():
    return t.stamped_filename("filter", PARAMS)   # filter_a77bf917.3mf
```

Three properties matter, and all three are tested:

| Property | Why |
|---|---|
| keys sorted | dict insertion order must not change the name |
| values via `repr` | `0.5` and `0.50` are the same parameter |
| floats rounded to 9 sig figs | float noise must not fork the filename, but a deliberate change still does |

Every generator exposes `default_out()`, and every verifier calls it rather
than repeating the filename. A verifier that hardcodes `grid_150mm.3mf` will
happily check a stale file forever; deriving the path from the same `PARAMS`
dict makes that failure impossible. The parameters are also written into the
3MF `<metadata>` block, so a file found on disk can be traced back to the
script state that produced it.

### 5.2 Derive mating parts from each other

When one part has to fit another, do not restate the shared dimension in both
scripts. Have the second import the first and compute:

```python
# make_holder.py -- the cradle is the filter's cone, offset outwards
def void_r(z, clearance=CLEARANCE):
    slope = (F.SHOULDER_R - F.STEM_R) / F.CONE_H
    return F.STEM_R + clearance + slope * (z - plate_h)
```

Restating it means a change to `SHOULDER_R` updates the filter and silently
leaves the holder behind — and the two then print as parts that no longer fit,
which is the failure mode least likely to be noticed. Include the dependency
in the holder's own parameters (`filter_hash=t.params_hash(F.PARAMS)`) so the
holder's filename changes too.

---

## 6. Combining shapes: overlap or not

Two options, and they are not equivalent.

**Overlap** — emit separate closed shells that intersect. Every slicer
(PrusaSlicer, Cura, Orca) unions them on load, so this just works. Cost:
your file reports the *sum* of volumes, not the union, and coincident faces
can occasionally confuse a slicer.

**Disjoint** — arrange so nothing interpenetrates. Then signed volume is
the true volume with no double counting, and there are no coincident faces.
Strongly preferred where you can achieve it.

**Subtraction is the third option, and overlap-land cannot do it.** Overlapping
shells only ever *add* material, so a passage that has to run *through* one of
them cannot be cut afterwards. `make_wye.py` is the worked case; the trap it
falls into and how it is escaped are in [§8](#8-pitfalls-that-will-bite-you).

The two-layer grid is the canonical example of choosing disjoint: because the
two bar directions occupy different `z` bands, all 300 bars are disjoint.
Compare the reported volume:

```
overlapping single-layer grid: 4530.000 mm^3  (= sum of boxes)
                                 (true union is 3389.950 mm^3 -- 25% lower)
two-layer disjoint grid:       4530.000 mm^3  (true union, exactly)
```

Same number, completely different meaning. Always say which one you mean.

---

## 7. Verification discipline

**Verify the artifact, not the generator.** Re-parse the written file with
code that shares nothing with the builder. If your builder and your checker
share a helper, a bug in that helper produces two agreeing wrong answers.

The layered checks that earned their keep, roughly in order of value:

1. **Package** — zip integrity, all three OPC parts present, XML parses,
   `unit` is `millimeter`, build list references real objects.
2. **Mesh sanity** — indices in range, no degenerate triangles, all
   coordinates finite, signed volume positive.
3. **Shell count and structure** — `len(mesh.shells())` equals the number of
   solids you intended; each shell has the vertex count you expect.
4. **Overlap** — `aabb_overlap_pairs` for box-like solids. Non-empty result
   means unintended interpenetration. (Meaningless for a single annulus
   shell, whose AABB spans the whole disc — use a radial containment test
   there instead.)
5. **Dimensions** — bounding box against the spec, plus per-layer/per-bar
   intervals. This is where parity bugs surface.
6. **Volume** — against the closed form you derived in step 3.
7. **Independent method** — re-derive one number by a *different* algorithm.
   For the trimmed grid, Simpson's rule on `2*sqrt(R²-y²) dy` validates the
   polygon clipping; for facets, compare against the inscribed-polygon area
   ratio.

### Verifying that two parts fit

Each part can pass every dimensional check on its own and still not fit its
partner. That class of bug survives to the printer, so make it a test.

`verify_parts.py` transforms the filter into its seated pose and checks it
against the holder's solid:

1. **Plane coincidence** — the filter's shoulder must land exactly on the
   rim plane (`SEAT_C - GRID_H == Z_SEAT`, to 1e-9, not "close").
2. **Radial overlap** — the rim must be *narrower* than the filter's disc or
   the disc falls through the cradle. Assert the seat is a usable width, not
   merely non-negative.
3. **Gap, measured on the facets** — not on the ideal cones. Interpolate both
   faceted polygons at a shared height and take the minimum radial
   difference over all angles. Two inscribed polygons with different segment
   counts give a gap that is *not* the nominal one (1.4940 mm measured
   against 1.5 mm nominal), because the faceting differs.
4. **Containment** — the decisive one. Transform every vertex of part A into
   part B's frame and test it for containment in B's solid. Non-empty result
   means the parts clash.

Two more subtleties, both of which produced confident wrong answers here.

**Contact is an area overlap, not a vertex hit.** Two annular faces meet over
the intersection of their radius ranges, which says nothing about where
vertices happen to land. With the holder wall at 4 mm the rim was
`r 31.5..35.5` and 356 filter vertices sat inside it, so a vertex test looked
fine. At 0.8 mm the rim is `r 31.5..32.3` and **zero** vertices land in it,
while the surfaces still touch over 0.8 mm. A vertex test reports contact when
it happens to exist and misses it when it matters. Compare intervals.

**A part built by overlapping shapes needs union semantics when you test a
point.** Ray-casting parity across the whole triangle soup computes the
*symmetric difference*: two triangles' worth of shell around a point that is
inside both shapes cancels to "outside". Every part here is a soup of
overlapping shells, so the region where a thick cup wall swallows a funnel's
thinner wall reported as void -- a solid part described as hollow exactly
where the union was meant to fill it. Test each shell separately and OR the
results (`checks.inside`). Coincident surfaces need the same care: a ray along
an axis hits shared edges exactly and counts them twice, so use a skewed
direction.

**A hollow part can pass every measurement and still be solid.** The filter's
outlet was a cone with a bore drilled through it: correct bounding box,
plausible volume, valid mesh, and it would not filter anything, because the
only void inside it was the bore and the whole 70 mm disc drained onto a solid
ceiling. Volume, bounding box and shell count cannot see this. Check the
*interior*: walk a radial profile at several heights and assert the passage is
open where it should be, and measure what fraction of the media actually has
somewhere to drain to -- here it was 6%.

**Where two features meet, check for a gap, not just for overlap.** The
funnel's outer wall stopped at r=30 while the cup's bore was r=35, leaving a
5 mm annular void under the media's edge. Anything passing through the disc
outboard of r=30 fell into it and ran down the outside of the part. No
measurement of either feature alone shows this; it is a property of the
junction. Sample the band between them and require it to be solid.

**Compare polygons as sets, not as lists.** The collar bore is generated under
one angle convention and the reference under another, so a positional
comparison reports a **70 mm** "offset" for two identical 136-gons — that is
a diameter, i.e. the wrap point, not an error. Measure each vertex against
its *nearest* counterpart instead.

And the original trap still stands: **contact is not interference.** The
filter's shoulder lies *on* the rim, so 356 of its 22,416 vertices sit
exactly on the holder's surface. A naive containment test calls that a clash.
Test strictly inside the material and treat the shared plane as contact — but
also assert the contact exists, so you have not merely suppressed the error
to make it go away.

The same trap, reversed, hides real errors: a face that merely *touches* the
bed at its lowest point is not supported by the bed. Test the highest point.

### Renders are for humans, not assertions

Pixel area is a bad assertion target. PIL fills polygons pixel-inclusively,
so every feature edge can claim up to half a pixel. On the 150 mm grid that
was a **1998 mm²** discrepancy on a 16950 mm² part — 12% — entirely from
rasterisation bias, with a per-pixel fill budget around 2220 mm². It looked
exactly like a geometry bug and was not one.

When you must measure a render:

- budget for it explicitly (`perimeter_px / 2`), and
- prefer **counting features** to measuring area. A scanline across a gap row
  counting 150 bright runs is robust where total area is not, because run
  *counts* are insensitive to sub-pixel edge position.

---

## 8. Pitfalls that will bite you

**Parity — the division that doesn't come out even.**
150 mm at 1.0 mm pitch is 150.5 bars. There is no exact all-0.5 mm layout.
`bar_spans` resolves it by widening the two end bars to `width + gap/2`,
which absorbs the remainder symmetrically. Alternatives: half-width edge
bars (too thin to slice), or letting the grid run long and clipping (which
`add_grid(clip_poly=...)` does). Pick one deliberately and say which.

**Faceted curves are inscribed, not exact.**
Sutherland-Hodgman cannot clip to a true circle, only to a polygon.
An `n`-gon has area `(n/2)·R²·sin(2π/n)`, which is *less* than `πR²`. A
trimmed grid measured against the ideal-circle integral is legitimately
~0.02% low; against the wrong reference it looks broken. Use
`auto_segments(r, tol)` and state the tolerance.

**Share the facet polygon between clip region and solid.**
If the grid is trimmed with a 172-gon but the tube is built from a 168-gon,
the trimmed bars either poke through the bore or float off it. Build the
bore polygon once and use it for both `tube()` and `clip_poly=`.

**Coordinate precision.**
`%.6g` on a coordinate like `54.9999711` silently rounds it. This made every
wall vertex report a slightly different radius, so "exactly two distinct
radii" became 43 spurious values and radius comparisons failed by 5e-5 mm.
Use `repr(float)` — shortest round-trip, exact, and no size penalty worth
mentioning.

**Slivers.**
A bar clipped near the tangent to a circle can produce a tiny sliver whose
triangles are numerically meaningless. Drop pieces below a minimum area
(`add_grid(min_area=...)`) and report how many were dropped.

**Centring.**
`bar_spans` returns intervals over `0..size`. To centre them you must offset
by `-size/2` *as well as* offsetting the perpendicular axis. Getting this
half-right produces a 225×225 mm grid from a 150 mm request, with the two
layers offset from each other — and a bounding-box assertion is what catches
it. **Assert your bounding box every single time.**

**Union area is not material.**
For a grid, `footprint − holes` is the *projection* of all bars. If the bars
split across two `z` layers, each layer holds only one direction, so material
per layer is `bar_sum × length` and the 3D volume is
`2 × bar_sum × length × thickness` — which is **not**
`(footprint − holes) × thickness`. On this grid that projection-based figure
is 6779.9 mm³ against a true 4530 mm³: 33% too large.

**A union of shells cannot branch.** Overlapping shells only ever *add*
material, so a passage that has to run *through* one of them cannot be cut
afterwards — it has to be missing from the start. The failure is quiet, and
easy to write by accident: glue three tubes together at one point for a tee and
each tube's wall runs straight through the other two bores and plugs them.
That is exactly what the first wye did — two of the three pairs of mouths came
out shut and the third was pinched to a 3 mm crescent, in a mesh that sliced,
had a plausible volume, and was not a tee. A part can be solid, watertight and
useless, which is why `verify_wye.py` now walks mouth-to-mouth *through* the
junction and not just around the walls. With no boolean kernel to hand, the way
out is that the solid wanted is `(outer solids u blend ball) − (bores)`, and if
every one of those primitives is convex and contains the junction, a ray out of
the junction crosses each in one interval: the material is star-shaped about
it, and can be swept radially as `rho(u)..R(u)` per direction with the bores
simply absent. `make_wye.py` does that, and pulls each arm's wall back by the
distance at which it would start to lie inside another arm's bore.

**An angle's reference direction is not obvious, and getting it wrong is
invisible in a render.** A Vigreux finger was meant to hang down and inward;
its drop was written as `-cos(drop)` in `z`, which makes `drop` the angle from
the *vertical*, so the fingers that were supposed to point 55° down actually
lay 35° off level, with undersides 2° from horizontal. The part looked like a
Vigreux column from every angle and was a stack of nearly flat shelves. The
only thing that caught it was a checker measuring face slopes; state the
reference in the parameter's name or in one line above it, and check the
resulting angle rather than the number you passed in.

**A staging parameter that fights the geometry has to be told so.** Asking
for Vigreux fingers that reach the middle of a 16 mm bore is asking for long
fingers, and a finger is long when its drop is shallow: at 25° below level it
reaches 91% of its length inwards, at 55° only 57%. But the two ways to make
one reach are not interchangeable — a shallow finger is a downward face the
printer has to bridge (which is what the overhang angle is for), and three
fingers a row at 120° come no closer to each other than `sqrt(3) × tip`, so a
tip near the axis caps the finger's thickness at about `0.58 × tip`. A first
cut of `make_column.py` *solved* for the drop instead of using it, and the
result was worse than wrong: a thick finger over a short length implies a
large taper, a large taper forces a steeper drop, a steeper drop shortens the
finger, and the loop converged on 3 mm stubs pointing 88° down — fingers that
reached 0.1 mm into the bore while the parameters said 2.6 mm thick and
dropping 25°. Now the drop is the parameter, the tip radius is the design
intent, and a request the geometry cannot honour is refused with the number
that would work: *"a 1.5 mm tip needs FINGER_R at most 0.80 mm, or FINGER_R
2.6 needs a tip of at least 3.58 mm"*.

**Staging has to fit the pitch it is given.** Elements per column is a
parameter people turn up, and every internal shape has a height that follows
from its own geometry: a washer from its slope and radial span, a leaning plate
from its reach, a finger from its length and drop. Crowd them and they meet
each other, and the bore that was open at nine stagings is sealed at
twenty-seven — a part that still slices, still has a plausible volume, and
cannot distil anything. Either scale the shape to the pitch (`make_column.py`
does: the washer's hole opens up, the plate's reach shrinks, the finger
shortens) or refuse the count with the number that would fit, and let
`verify_column.py` walk the vapour path to prove the result is still open.

**A blank visual check is not a test.** I shipped a grid that looked correct
in every render while 300 solids interpenetrated at every crossing. The
defect was only visible in `shells()` and an overlap test.

**"Correct geometry" is not "printable".** Every check in
[§7](#7-verification-discipline) can pass while the model puts a feature at
`z = 0.2` when the layer height is 0.3, or asks a 0.6 mm nozzle for a 0.5 mm
wall. Geometry verification and process verification are different questions;
run both. But the reverse mistake is just as common: a checker that cries
"unprintable" at a perfectly ordinary single-bead wall trains you to ignore
it. Calibrate the thresholds to the part, and make severity levels that let
"fine, with a caveat" be distinct from "will not print".

**A hollow profile, revolved, is not a hollow solid.**
`revolve` closes *every* end of the profile it is given: an outline that stops
short of the axis gets a flat disc, so a profile drawn as an annulus — bore and
skin, the way a wall looks on paper — comes out solid plastic with a plausible
bounding box, plausible shells, plausible everything. `tube`, `tapered_tube` and
`sweep` take the bore and the outside separately and leave the bore open, which
is what a wall needs. The tell is a volume measured against a washer integral,
plus one probe sitting in the supposed bore reading solid.

**A piece's stations have to ascend.** `sweep` and `tapered_tube` join the
stations they are handed, in the order they arrive: one that sits lower than the
one before it folds the shell back through itself, the piece comes out with
negative volume, and the boolean reports only that "the mesh is inside-out or not
closed" -- nothing about the parameter that did it. Two hand-edited numbers were
enough: a nozzle's tip pulled below its own throat, and two neck stations given
the same height. `make_aspirator.py` now walks every piece's station list in
`validate` and names the pair that does not rise; equal heights are allowed,
since that is a segment of no length, and the profile that owns it drops it.

**A ring whose bore is outside its own skin eats its neighbour.** An annulus with
`r_in > r_out` is not a thin wall, it is an inside-out solid, and the union then
swallows whatever it overlaps -- here ten millimetres of the inlet tube simply
vanished from the result, with a shorter part, one shell and no error. Shrink one
radius in a stack of them and the bore that used to sit comfortably inside a skin
does not stay there. Check every profile at every station, `bore < skin`, not just
the one you edited.

**Coincident curved surfaces leave a boolean nothing to cut with.**
Two cylinders of the same radius running alongside each other never intersect,
so there is no curve to split them: the result is a doubled wall with a hairline
gap inside, a chain that refuses to close ("edges are used once"), or a surface
that survives where it should have vanished. Where one piece buries another,
offset the buried surface — a bore made *wider* than the governing one, a skin
made *narrower* — by a fraction of a millimetre. The surface the design wants is
then the one that wins, and the other disappears inside plastic. This repo's
aspirator carries that offset as a named parameter (`BURIED = 0.3`) precisely so
it is a design decision and not a coincidence.

**A step in diameter is a ceiling; a lean is not.**
Joining a 21 mm tube to a 32 mm chamber with a flat shoulder made 456 mm² of
horizontal down-facing annulus — an outside surface, so support would have been
removable, but 5 mm wide and printed 25 mm up in mid-air. Slanting the same
material into a cone at 55° off the horizontal costs nothing in strength, reads
better as the intended Z-shaped wall, and leaves the printer nothing to hold up.
Wherever a part changes diameter, the joint's underside is a design choice.

**Rays and point tests both lie at a vertex column.**
An `n`-gon facetted circle puts a vertex every 360°/n — at `n = 40`, every 9°.
A probe that lands on such a column, or a ray that runs down it, crosses every
surface *twice*, and parity flips: an open bore reads as solid, and a solid wall
reads as empty. `checks.inside` is no exception. Stand every probe a fraction of
a millimetre off the axis and off the whole millimetres, and note that a chord
measures *small* — by at most the sagitta `r(1 − cos(π/n))` — so a tolerance
against a nominal radius must be one-sided, not `±`.

**A snapped height is a float divide, and floats eat halves.**
`round(z / grid) * grid` looks like quantising and sometimes *is*: 1.4 / 0.4 is 3.4999999999999996, so it rounds to 3 and the feature lands 0.2 mm low, silently. On the scoopula the floor was asked at 1.4 mm and came out 1.2 mm -- which happened to be exactly the handle's thickness, so the joint that every other paragraph here goes to trouble to keep *off* a coincidence got one. Two fixes, both of which belong in the generator: snap with `floor(z / grid + 0.5 + 1e-9)`, and refuse a feature height that is not already a whole number of layers, naming the two grid heights it sits between (`on_grid` in `make_scoopula.py`). A measurement tolerance in the verifier will not save you: 1.2 mm of floor passed a check on 1.4 mm that allowed for a grid of 0.4.

**Off-grid `z` from a boolean is not an off-grid feature.**
A lofted part that has been through a boolean reports hundreds of vertex `z` values off the layer grid, and `printability.layer_alignment` calls every one of them a warn. Most are T-junction splits: a vertex dropped onto a *sloped* face where another surface crossed it, which has no height of its own to keep and gets sliced away without changing the part. What has to be on the grid is the **feature heights** -- the planes that horizontal faces lie in: bed, floor, rim land, handle top. Test those (`verify_scoopula.py` walks the triangles and asserts every horizontal face above 1e-6 mm^2 sits on a whole layer) and ignore the vertex tally, or you will spend the afternoon chasing noise.

**Self-test the primitives.** `wedge` shipped producing zero volume and a
non-manifold mesh while every example that used it "worked". One table of
primitive-vs-closed-form tests caught it immediately.

---

## 9. When to ask the user

Ask when the ambiguity changes **topology or a dimension**, and the options
are qualitatively different — not when it merely changes a number you could
default sensibly. In practice that means:

- whether two co-linear features overlap or are stacked (changes the mesh's
  fundamental structure),
- whether a spec's numbers are compatible (0.5 mm lines at 0.5 mm pitch),
- whether a feature is trimmed at an inner or outer boundary,
- whether a floor is solid or open.

Do *not* ask about tolerance conventions, facet counts, or file format —
those are defaults you should choose and state.

---

## 10. Pre-flight checklist

Before calling a model done:

- [ ] Every spec dimension appears in the report with its unit.
- [ ] Bounding box matches the spec on all three axes.
- [ ] Signed volume is positive (outward winding everywhere).
- [ ] Volume matches the closed form derived by hand.
- [ ] Shell count matches the number of intended solids.
- [ ] No unintended overlap (or stated intentional overlap).
- [ ] Index range, degenerate triangles, finite coordinates.
- [ ] Zip/XML integrity, correct unit, build list resolves.
- [ ] Renders produced; area assertions budgeted, or replaced by run counts.
- [ ] `check_printability.py` run for the intended nozzle and layer height.
- [ ] No `warn`-level findings; `info`-level ones read and understood, not
      dismissed.
- [ ] No overhangs past the printer's angle (`--angle`), or supports planned
      and stated.
- [ ] All `z` values on the layer grid for the target layer height.
- [ ] Enough bed contact for first-layer adhesion.
- [ ] Any fit, hole or clearance allowance stated (holes print undersize).
- [ ] Assumptions and defaults stated explicitly.

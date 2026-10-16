#!/usr/bin/env python3
"""Generate a fractionating column: a tube of staged internals, one piece.

A fractionating column does not separate anything by itself. It gives the
vapour and the reflux a long, wet, shared path, so that the mixture climbing
it re-equilibrates again and again; what comes off the top is then richer in
the more volatile component than one evaporation would give. Two things matter
and they pull against each other:

* **contact** -- surface area, and a path that keeps the vapour moving over the
  wet surfaces rather than straight up the middle;
* **pressure drop** -- under vacuum, every pascal lost across the column is a
  pascal lost from the still, and it shifts everything in the pot. A column
  packed tight enough for good separation can make the distillation impossible.

The third constraint is the one that shapes this part: it is printed *as one
piece, standing up*, so nothing inside may need support -- support that goes
into a 16 mm bore cannot be got out again.

    BORE = 16        the open bore of the body, mm
    BODY = 150       the staged length between the joints, mm
    STYLE = "cones"  what the staging is made of

FOUR STAGING STYLES, all printable standing up for the same reason: every
surface inside makes at least 45 degrees with the horizontal. A cone whose
half-angle is at most 45 degrees has no surface steeper than that -- its
normals are never more than `half-angle` off vertical -- which is why cones are
the obvious shape here, and why a horizontal finger or plate never is (the
underside of anything horizontal is a ceiling).

    cones    conical washers, each attached to the wall, each with a central
             hole; the vapour climbs through the holes, the reflux runs back
             down the same cones. Sieve-plate behaviour with no flat plate.
    helix    one or more bands, standing on edge, twisting as they rise,
             touching the wall at the outer edge. Contact all the way up, a
             long path, and no staging to flood. HELIX_ARMS > 1 makes a
             multi-start helix -- a triangular one at three -- by phase-shifting
             the bands by 360/arms, which is the difference between three arms
             and one arm built three times.
    vigreux  fingers pointing down and inward, staggered, the classic pattern:
             the reflux clings to and re-evaporates off them.
    tesla    louvred plates, tilted no more than 45 degrees off vertical and
             alternating, so the vapour zig-zags. This borrows the *only*
             useful part of the Tesla valve idea here -- that its passages can
             be planar, and planar walls print vertical -- because the diode
             behaviour itself is nothing to do with fractionation.

The build prints the geometry and the vapour path; `verify_column.py` then
walks that path through the mesh to prove the bore really is open end to end,
and measures what the path is worth -- the surface it touches, and the
narrowest free cross-section, which is the pressure-drop proxy. A column with
a blocked passage prints perfectly and distils nothing.

Joints are per-end and independent, a male 14/20 at the bottom into the flask
and a female 14/20 at the top for the still head, on the same 1:10 ground
joint standard as everything else here.

Usage: python3 make_column.py [output.3mf] [style] [--elements N] [--angle N]

The style may be given instead of or as well as the output name, so all four
can be built and compared, and the two numbers that decide the staging can be
overridden without editing the file:

    python3 make_column.py column_helix.3mf helix --elements 27 --angle 20
    python3 make_column.py column.3mf vigreux --finger-r 0.8 --finger-tip 1.5

`--angle` is the shallowest face this printer holds without support, in degrees
from the horizontal; it is also `OVERHANG_ANGLE` in the environment, and it is
stamped into the file, because the staging is shaped around it.
"""

import math
import sys

import checks as K
import csglib as csg
import junctionlib as J
import printability as PR
import threelib as t

# --- Parameters ------------------------------------------------------------
BORE = 16.0          # bore of the body, mm
BODY = 150.0         # staged length between the joints, mm
STYLE = "helix"      # cones | helix | vigreux | tesla
ELEMENTS = 18         # stagings: plates, turns, finger rows or louvres
BOTTOM_JOINT, BOTTOM_GENDER = "14/25", "male"
TOP_JOINT, TOP_GENDER = "14/25", "female"
WALL = 2.5           # tube wall thickness
TIP_WALL = 1.0       # wall at a male tip
CAP = 6.0            # length over which the bore meets the body, mm
# The shallowest face that still prints without support, in degrees from the
# horizontal. It is a printer property, not a part property, but it decides
# real geometry here -- how steep a washer is, how far a plate leans, how far a
# finger drops -- so it is a parameter and goes into the file: a column has to
# say which printer it was designed for, or the numbers in it mean nothing.
OVERHANG = PR.OVERHANG_DEG

# The helix band's inner edge, as a fraction of the bore radius. Zero puts the
# inner edge on the axis and a negative value puts it past the centre, so the
# band spans the whole bore and the vapour has to go round it rather than up
# the middle -- which is a different column, not a broken one.
HELIX_FRAC = -0.1

# How many bands the helix has: one start, or several arms evenly spaced round
# the axis. Three arms is a triangular helix, and the arms are phase-shifted by
# 360/arms so they are three arms and not one band built three times.
HELIX_ARMS = 3

# Vigreux fingers. `FINGER_TIP` is the design intent -- how close to the middle
# of the bore the tips come -- and everything else follows from it: the length
# from the drop, the length again from the pitch it has to hang in. `FINGER_R`
# and `FINGER_TAPER` shape it: a taper of 0 is a rod, and any taper eats into
# the drop, because the underside of a tapered finger is `drop - taper` below
# level and that is the number the printer has to be able to bridge.
#
# A finger is long when the drop is shallow. At 25 degrees below level a finger
# reaches 91% of its length inwards, so 7.4 mm of finger spans the 6.7 mm from
# the wall to a tip 1.5 mm off the axis; at 55 degrees the same 7.4 mm reaches
# only 4.2 mm and stops a third of the way out. The catch is the underside: a
# finger that shallow is only printable on a printer that holds shallow
# overhangs, which is what `--angle` and `OVERHANG_ANGLE` are for.
FINGER_PER_ROW, FINGER_DROP, FINGER_R = 3, 25.0, 2.6
FINGER_TIP, FINGER_TAPER = 1.5, 0.0
FINGER_CLEARANCE = 0.8   # mm of air to leave between two fingers a row

TOL = 0.01           # max chord deviation on round surfaces

PARAMS = dict(bore=BORE, body=BODY, style=STYLE, elements=ELEMENTS,
              overhang=OVERHANG, helix_frac=HELIX_FRAC, helix_arms=HELIX_ARMS,
              finger_tip=FINGER_TIP, finger_drop=FINGER_DROP,
              finger_r=FINGER_R, finger_taper=FINGER_TAPER,
              bottom_joint=BOTTOM_JOINT, bottom_gender=BOTTOM_GENDER,
              top_joint=TOP_JOINT, top_gender=TOP_GENDER, wall=WALL,
              tip_wall=TIP_WALL, cap=CAP, tol=TOL)

STYLES = ("cones", "helix", "vigreux", "tesla")

# Staging detail. These are construction constants rather than parameters, but
# they live here because the vapour path has to be worked out from the same
# numbers the staging was built with -- if the two drift apart the path stops
# being a proof that the passage is open.
CONE_SLOPE, CONE_THICK, CONE_HOLE_FRAC = 50.0, 1.6, 0.1
BAND_THICK = 1.2

# How far a staging piece reaches into the wall. Every piece is attached to the
# wall, and stopping exactly at the bore's surface is a *coincidence*: the two
# surfaces occupy the same place, and a boolean union of two surfaces that
# coincide rather than cross leaves cracks along them -- the wye lost 195 mm of
# its surface that way. Reaching a little way in makes the joint a crossing,
# which the union resolves exactly, and the extra material is inside the wall
# where it does not show.
BITE = 0.5

# How far each helix arm's ends are staggered in z. Every arm is a twisted
# slab, and their end caps are flat: three arms starting at the same height put
# three coplanar quads in the same plane, partly overlapping -- the one case a
# BSP resolves badly. Measured on the default column: unioning two arms left 3
# real gaps of 0.477 mm and the third union was refused; staggering the ends by
# 0.05 mm gives 0 gaps. It is invisible and each arm still reaches the staged
# span.
BAND_STAGGER = 0.05

# How many staging pieces the build will union in one call. The helix has three
# and takes twelve seconds; the cones have eighteen and the vigreux thirty-six,
# and a single call over that many solids is not something this boolean can do.
UNION_MAX_PIECES = 8
LOUVRE_TILT, LOUVRE_THICK, LOUVRE_REACH = 40.0, 2.0, 3.0


def finger_geometry(p, pitch):
    """(length, drop, radius, tip) for a Vigreux finger, and why each is what it is.

    The drop is the parameter, used as given: it is the angle the finger hangs
    at, it is what sets how far inwards a given length reaches, and a build
    that quietly changed it would be building a different column from the one
    asked for. The length then follows from how close to the middle the tips
    are meant to come, and the only thing that overrides it is the pitch: a
    finger taller than the gap between one row and the next would grow into the
    row above, so past that point the tips fall short of the middle and the
    build says so rather than pretending.

    The taper is checked, not solved for: a finger's underside sits
    `drop - taper` below level, and if that is shallower than the printer can
    hold, no length fixes it -- the message names both ways out.
    """
    d = dims(p)
    r_base = d["r_bore"] + 0.2
    drop = p["finger_drop"]
    tip = p["finger_tip"]
    # Three fingers a row at 120 degrees: their axes come closest to each other
    # where they are closest to the axis, and that distance is `sqrt(3) * tip`.
    # Two fingers of radius r touch once that is 2r, which is a hard ceiling on
    # how near the middle a tip can come for a given thickness -- and worth
    # saying out loud, because "make the fingers reach the middle" and "make
    # the fingers thick" are the same request pulling opposite ways.
    r_finger = p["finger_r"]
    closest = math.sqrt(3.0) * tip
    if closest < 2.0 * r_finger + FINGER_CLEARANCE - 1e-9:
        raise ValueError(
            f"{FINGER_PER_ROW} fingers a row come no closer to each other than "
            f"{closest:.2f} mm, and {r_finger:g} mm fingers are "
            f"{2.0 * r_finger:g} mm thick: a {tip:g} mm tip needs FINGER_R at "
            f"most {(closest - FINGER_CLEARANCE) / 2.0:.2f} mm, or FINGER_R "
            f"{r_finger:g} needs a tip of at least "
            f"{(2.0 * r_finger + FINGER_CLEARANCE) / math.sqrt(3.0):.2f} mm")
    length = (r_base - tip) / math.cos(math.radians(drop))
    tallest = 0.85 * pitch / math.sin(math.radians(drop))
    clipped = length > tallest
    if clipped:
        length = tallest
    tip_got = r_base - length * math.cos(math.radians(drop))
    underside = drop - p["finger_taper"]
    if underside <= p["overhang"]:
        raise ValueError(
            f"a finger dropping {drop:g} deg with a {p['finger_taper']:g} deg "
            f"taper has its underside {underside:g} deg below the horizontal, "
            f"and this build is told the printer holds {p['overhang']:g} deg. "
            f"A finger that shallow needs a printer that bridges shallower "
            f"overhangs than that (--angle {underside - 1.0:g}, or "
            f"OVERHANG_ANGLE={underside - 1.0:g}) or a steeper FINGER_DROP")
    return length, drop, p["finger_r"], (tip_got, tip, clipped)


def louvre_reach(p):
    """How far past the axis a leaning plate reaches, given the pitch it has.

    A plate's height follows from how far it reaches, and it has to fit between
    one plate and the next: past a certain pitch the plates stop reaching the
    middle at all and become baffles along the wall, which is the honest thing
    to build and the honest thing to report. The vapour path is worked out from
    this same number, so the two cannot disagree about where the plates are.
    """
    d = dims(p)
    z0, z1 = stage_span(p)
    pitch = (z1 - z0) / (int(round(p["elements"])) + 1)
    lean = math.radians(louvre_tilt(p))
    return min(LOUVRE_REACH,
               (0.85 * pitch - LOUVRE_THICK * math.cos(lean)) / math.tan(lean)
               - d["r_bore"] - 0.2)


def louvre_tilt(p):
    """How far a leaning plate leans, from the printer's overhang angle.

    `rotate_y(tilt)` tips the plate until its long axis is `tilt` below the
    horizontal, and its two large faces contain that axis -- so the faces lean
    `tilt` off the horizontal, not off the vertical. The underside is the face
    that matters and it has to be steeper than the printer can bridge.
    """
    return max(LOUVRE_TILT, p["overhang"] + 5.0)


def default_out():
    return t.stamped_filename("column", PARAMS)


def dims(p):
    """Everything derived from the parameter set, in one place."""
    r_bore = p["bore"] / 2.0
    r_body = r_bore + p["wall"]
    r_mouth_b, r_tip_b, glen_b = J.joint_dims(p["bottom_joint"])
    r_mouth_t, r_tip_t, glen_t = J.joint_dims(p["top_joint"])
    bore_b = r_tip_b - p["tip_wall"] if p["bottom_gender"] == "male" else r_tip_b
    bore_t = r_tip_t - p["tip_wall"] if p["top_gender"] == "male" else r_tip_t
    return {"r_bore": r_bore, "r_body": r_body, "glen_b": glen_b,
            "glen_t": glen_t, "r_mouth_b": r_mouth_b, "r_tip_b": r_tip_b,
            "r_mouth_t": r_mouth_t, "r_tip_t": r_tip_t,
            "bore_b": bore_b, "bore_t": bore_t,
            # z of the body's first and last station
            "z_body0": glen_b + p["cap"], "z_body1": glen_b + p["cap"] + p["body"]}


def end_profile(p, which):
    """The two radii at each end of a ground joint: (glen, free, inward).

    A male joint is its cone as an *outer* surface, narrowest at its free end,
    with a bore of constant radius; a female joint is the same cone as a
    *socket*, so its bore is widest at the mouth and narrows going in, while
    its outside is the cone plus the wall. `s` is measured from the free end,
    exactly as `junctionlib` measures it, and here the free end is the bottom
    of the column at z = 0 or the top of it.

    Each pair is `(r_bore, r_outer)` at s = 0 and at s = glen.
    """
    r_mouth, r_tip, glen = J.joint_dims(p[which + "_joint"])
    gender, wall = p[which + "_gender"], p["wall"]
    if gender == "male":
        bore = r_tip - p["tip_wall"]
        return glen, (bore, r_tip), (bore, r_mouth)
    return glen, (r_mouth, r_mouth + wall), (r_tip, r_tip + wall)


def body_stations(p):
    """The tube itself, free end to free end, as swept stations.

    Every transition is a cone whose slope is at most 45 degrees from vertical,
    which is what `cap` is checked against: a shoulder would be a flat annular
    ceiling inside the bore, and the whole point of this part is that there is
    nothing inside it that needs support.
    """
    d = dims(p)
    r_bore, r_body, cap = d["r_bore"], d["r_body"], p["cap"]
    glen_b, free_b, in_b = end_profile(p, "bottom")
    glen_t, free_t, in_t = end_profile(p, "top")
    return [
        # Bottom: the free end, the end of the ground zone, then out to the body.
        (0.0, free_b[0], free_b[1]),
        (glen_b, in_b[0], in_b[1]),
        (glen_b + cap, r_bore, r_body),
        # Body.
        (d["z_body1"], r_bore, r_body),
        # Top: out to the end of the ground zone, then the free end.
        (d["z_body1"] + cap, in_t[0], in_t[1]),
        (d["z_body1"] + cap + glen_t, free_t[0], free_t[1]),
    ]


def sweep_body(mesh, p, n):
    """Sweep the body along +z from the bottom mouth to the top mouth."""
    stations = body_stations(p)
    out = []
    for z, r_in, r_out in stations:
        out.append(((0.0, 0.0, z), (0.0, 0.0, 1.0), r_in, r_out))
    t.sweep(mesh, out, n, ref=(1.0, 0.0, 0.0))
    return stations


# ---------------------------------------------------------------------------
# Staging
# ---------------------------------------------------------------------------

def plaque(mesh, z_outer, r_outer, r_inner, slope_deg, thickness, n,
           overhang, up=True):
    """One conical washer: a plate attached at the wall, with a central hole.

    A washer is an annulus swept along the axis with its radii sloping, so the
    rings stay horizontal -- which leaves `thickness` a horizontal wall, and the
    wall a little thinner than that where the plate is steep. `slope_deg` is
    measured from the horizontal and has to beat the printer's overhang angle,
    because the plate's underside is a downward face: a plate flat enough to
    puddle is a plate flat enough to need support.
    """
    if slope_deg <= overhang:
        raise ValueError(
            f"a washer at {slope_deg:g} deg from horizontal is shallower than "
            f"the {overhang:g} deg this build is told the printer holds, "
            f"so its underside would need support inside the bore")
    dz = abs(r_outer - r_inner) * math.tan(math.radians(slope_deg))
    z_inner = z_outer + dz if up else z_outer - dz
    if z_inner < z_outer:
        z_outer, z_inner = z_inner, z_outer
    stations = [((0.0, 0.0, z_outer), (0.0, 0.0, 1.0),
                 r_outer - thickness, r_outer),
                ((0.0, 0.0, z_inner), (0.0, 0.0, 1.0),
                 r_inner, r_inner + thickness)]
    t.sweep(mesh, stations, n, ref=(1.0, 0.0, 0.0))
    return z_inner


def finger(mesh, z, azimuth_deg, r_base, length, drop_deg, radius, tip, n):
    """A Vigreux finger: a rod off the wall, pointing down and inward.

    Its base sits inside the wall so the union is one solid. The axis drops
    `drop_deg` below horizontal, which is what keeps it printable: the
    underside of a horizontal finger is a ceiling, and a finger hanging at 25
    degrees has its underside at 25 degrees. `radius` and `tip` are the two
    ends; a tip of 0 makes a point, which is what a tapered finger wants,
    because a flat tip facing down is an overhang in its own right.
    """
    a = math.radians(azimuth_deg)
    drop = math.radians(drop_deg)
    # `drop_deg` is below the horizontal: at 90 the finger points straight
    # down, at 0 it lies flat. Getting this backwards is not a small error --
    # a finger meant to hang down at 70 degrees ends up 20 degrees off level,
    # with an underside only a couple of degrees from horizontal, which is a
    # ceiling the printer has to bridge.
    axis = (-math.cos(a) * math.cos(drop), -math.sin(a) * math.cos(drop),
            -math.sin(drop))
    base = (math.cos(a) * r_base, math.sin(a) * r_base, z)
    bump = t.Mesh("finger")
    t.cone(bump, radius, tip, 0.0, length, n=n)
    phi = math.degrees(math.acos(max(-1.0, min(1.0, axis[2]))))
    theta = math.degrees(math.atan2(axis[1], axis[0]))
    t.rotate_y(bump, phi)
    t.rotate_z(bump, theta)
    t.translate(bump, *base)
    mesh.merge(bump)
    return base, tuple(base[k] + length * axis[k] for k in range(3))


def louvre(mesh, z, tilt_deg, r_bore, depth, thickness, reach, side=1.0):
    """A leaning baffle: buried in the wall on one side, reaching past the axis.

    Tilting is the only way to put a wall across a bore without a ceiling, and
    `tilt_deg` is measured off the *horizontal*: the plate's long axis tips
    that far below level, and its two large faces contain that axis, so they
    lean by the same angle. The underside is one of them, which is why the
    tilt has to beat the printer's overhang angle rather than the other way
    round -- `louvre_tilt` is where that is worked out.

    `reach` is how far past the axis the free end goes, which is what sets the
    gap the vapour gets through on that side. A plate that reached *both* walls
    would seal the bore, which is the trap this style walks into: the useful
    part of the Tesla idea here is only that its walls can be planar and lean,
    since planar walls print and the diode behaviour is nothing to do with
    fractionation.
    """
    lean = math.radians(tilt_deg)
    length = (r_bore + 0.2 + reach) / math.cos(lean)
    plate = t.Mesh("louvre")
    t.box(plate, -length / 2.0, -depth / 2.0, -thickness / 2.0,
          length / 2.0, depth / 2.0, thickness / 2.0)
    t.rotate_y(plate, tilt_deg)   # the +x end drops, so that is the buried end
    t.translate(plate, r_bore + 0.2 - (length / 2.0) * math.cos(lean), 0.0,
                z + (length / 2.0) * math.sin(lean))
    if side < 0.0:
        t.rotate_z(plate, 180.0)
    mesh.merge(plate)
    return length


def stage_span(p):
    """The staged stretch of the body: where the internals live."""
    d = dims(p)
    return d["z_body0"] + 2.0, d["z_body1"] - 2.0


def helix_turns(p):
    """(turns used, turns asked for, most that print) for the helix style.

    The band's outer edge leans `atan(r * dtheta/dz)` off vertical, and a
    surface that leans far enough becomes a downward face the printer has to
    bridge: the ceiling is `90 - printability.OVERHANG_DEG` degrees of lean, so
    the twist rate has a hard ceiling at any radius. It is reported rather than
    quietly obeyed -- a column that dropped half its staging would still look
    right, and would fractionate nothing like as well -- and it moves with the
    printer, since a machine that holds shallow overhangs can take a finer
    twist and get a longer path out of the same tube.
    """
    d = dims(p)
    z0, z1 = stage_span(p)
    r_out = d["r_bore"] + 0.2
    want = int(round(p["elements"])) / 3.0
    ceiling = math.tan(math.radians(90.0 - p["overhang"]))
    most = (z1 - z0) * ceiling / (2.0 * math.pi * r_out)
    return min(want, most), want, most


def twisted_band(mesh, z0, z1, r_in, r_out, turns, thickness, overhang,
                 offset_a=0.0, sections_per_turn=48, lean_r=None):
    """A band standing on edge, twisting as it rises.

    A band's surfaces are near vertical, so the twist is free to be fine -- and
    the finer it is the longer the path the vapour has to follow. It is not free
    to be *unlimited*: every section is a set of distinct planes, and a boolean
    between two of these bands splits every facet by every plane of the other.
    Cutting it to 16 sections per turn was tried and is worse, not better:
    the arms' surfaces then differ by more than the boolean can reconcile and
    the union came back as three shells and 597 gaps, where the fine ones leave
    three. The cost has to be paid somewhere else. The limit is
    the surface angle at the outer edge: the surface leans by
    `atan(r * d(theta)/dz)`, so the twist rate is chosen to keep the band's
    faces steeper than the printer's overhang angle at `r_out`.

    `offset_a` starts the band at a different azimuth, which is how a
    multi-start helix is made: three bands at 0, 120 and 240 degrees are three
    arms of one helix, and three bands at 0, 0 and 0 are one band built three
    times. Pass it by name.
    """
    height = z1 - z0
    rate = 2.0 * math.pi * turns / height
    # The lean is checked at the radius where the band is *exposed*, which is
    # not necessarily where it ends: the band reaches BITE into the wall so the
    # union has a crossing, and those last tenths of a millimetre are buried.
    # Checking them would refuse a band whose visible faces are all fine.
    lean = math.degrees(math.atan((r_out if lean_r is None else lean_r) * rate))
    limit = 90.0 - overhang
    if lean > limit + 1e-9:   # inclusive: helix_turns lands exactly on it
        raise ValueError(
            f"{turns:g} turns in {height:g} mm leans the band's outer edge "
            f"{lean:.1f} deg off vertical; {limit:g} is where its faces start "
            f"needing support, so use fewer turns or a smaller bore -- or tell "
            f"the build the printer holds shallower overhangs than "
            f"{overhang:g} deg")
    steps = max(8, int(abs(turns) * sections_per_turn))
    half = thickness/2
    sections = []
    for k in range(steps + 1):
        u = k / steps
        a = rate * height * u + offset_a
        z = z0 + height * u
        ca, sa = math.cos(a), math.sin(a)
        pts = []
        # Counter-clockwise seen from +z, which is the direction of travel up
        # the loft. Wound the other way the band is inside-out: its signed
        # volume comes out negative, and while that was invisible with the band
        # merged into the tube's own positive volume, a slicer filling by
        # winding reads a reversed shell as material to *remove*.
        for x, y in ((r_in, -half), (r_out, -half), (r_out, half), (r_in, half)):
            pts.append((x * ca - y * sa, x * sa + y * ca, z))
        sections.append(pts)
    t.loft(mesh, sections)
    return lean


def mesh_of(mesh, shell):
    """One shell of a mesh, as a mesh of its own."""
    one = t.Mesh(mesh.name)
    one.verts = mesh.verts
    one.tris = [mesh.tris[k] for k in shell]
    return one


def stage_column(mesh, p, n):
    """Add the staging named by STYLE; return a description of it."""
    d = dims(p)
    r_bore = d["r_bore"]
    z0, z1 = stage_span(p)
    count = int(round(p["elements"]))
    style = p["style"]
    if style not in STYLES:
        raise ValueError(f"unknown style {style!r}; try one of {', '.join(STYLES)}")
    if count < 1:
        raise ValueError("a column needs at least one staging")

    info = {"style": style, "elements": count}
    pitch = (z1 - z0) / (count + 1)

    if style == "cones":
        # Conical washers: each one attached to the wall, each with a central
        # hole, all opening upward so the reflux runs out to the wall and the
        # vapour comes up the middle.
        r_attach = r_bore + 1.5
        hole = max(1.5, r_bore * CONE_HOLE_FRAC)
        # A washer's height follows from its slope and its radial span, so
        # crowding the stagings together has to be paid for at the hole: a
        # plate that does not fit in the pitch would run into the next one.
        fit = 0.8 * pitch / math.tan(math.radians(CONE_SLOPE))
        if r_attach - hole > fit:
            hole = r_attach - fit
        if hole > r_bore * 0.9:
            raise ValueError(
                f"{count} washers over {z1 - z0:g} mm leaves {pitch:.1f} mm of "
                f"pitch each, and a {CONE_SLOPE:g} deg washer that fits in that "
                f"has a hole of {hole * 2:.1f} mm in a {r_bore * 2:g} mm bore: "
                f"at most {int((z1 - z0) * math.tan(math.radians(CONE_SLOPE)) / (0.8 * (r_attach - 1.5)) ) - 1} "
                f"of these fit, so use fewer elements or a longer body")
        for k in range(count):
            z = z0 + pitch * (k + 1)
            plaque(mesh, z, r_attach, hole, CONE_SLOPE, CONE_THICK, n,
                   p["overhang"], up=True)
        info["hole"] = hole
        info["r_attach"] = r_attach

    elif style == "helix":
        turns, want, most = helix_turns(p)
        arms = int(round(p["helix_arms"]))
        if arms < 1:
            raise ValueError("a helix needs at least one arm")
        for i in range(arms):
            lean = twisted_band(mesh, z0 + i * BAND_STAGGER,
                                z1 - (arms - 1 - i) * BAND_STAGGER,
                                r_bore * p["helix_frac"], r_bore + BITE,
                                turns, BAND_THICK, p["overhang"],
                                offset_a=i * 2.0 * math.pi / arms,
                                lean_r=r_bore)
        info["arms"] = arms
        if arms > 1 and r_bore * p["helix_frac"] < 0.0:
            print(f"  note: {arms} arms that each reach past the axis cross one "
                  f"another near it, so the middle of the bore is solid and "
                  f"the vapour goes round the outside")
        info["turns"] = turns
        info["turns_wanted"] = want
        info["turns_max"] = most
        info["lean"] = lean

    elif style == "vigreux":
        rows = count
        length, drop, radius, (tip, tip_want, clipped) = finger_geometry(p, pitch)
        # A tapered finger comes to a point; a rod keeps its radius, and then
        # its flat tip is a disk facing down at `90 - drop`, which for a
        # shallow finger is steep enough not to matter.
        tip_r = (0.0 if p["finger_taper"] > 0.0
                 else max(0.0, radius - length * math.tan(
                     math.radians(p["finger_taper"]))))
        for row in range(rows):
            z = z0 + pitch * (row + 1)
            for j in range(FINGER_PER_ROW):
                az = 360.0 * j / FINGER_PER_ROW + (60.0 if row % 2 else 0.0)
                finger(mesh, z, az, r_bore + BITE, length, drop, radius,
                       tip_r, n)
        info["rows"] = rows
        info["per_row"] = FINGER_PER_ROW
        info["tip"] = tip
        info["tip_wanted"] = tip_want
        info["clipped"] = clipped
        info["length"] = length
        info["drop"] = drop
        info["radius"] = radius

    elif style == "tesla":
        # Leaning baffles across the bore, buried in the wall on one side and
        # reaching past the axis on the other, alternating: the vapour has to
        # cross over between them, and every face is within 45 degrees of
        # vertical, which is what lets the plates print at all.
        # A leaning plate's height is `length * sin(tilt)`, and its length
        # follows from how far past the axis it reaches -- so a tight pitch
        # makes the plates shorter, and past a point they stop reaching the
        # middle at all and become baffles along the wall.
        tilt = louvre_tilt(p)
        reach = louvre_reach(p)
        length = 0.0
        for k in range(count):
            z = z0 + pitch * (k + 1)
            length = louvre(mesh, z, tilt, r_bore + BITE, min(r_bore, 8.0),
                            LOUVRE_THICK, reach=reach,
                            side=1.0 if k % 2 == 0 else -1.0)
        info["tilt"] = tilt
        info["reach"] = reach
        info["reach_wanted"] = LOUVRE_REACH
        info["length"] = length
        info["gap"] = r_bore - reach
    return info


# ---------------------------------------------------------------------------
# The vapour path
# ---------------------------------------------------------------------------

def vapour_path(p):
    """A polyline that the open passage must contain, bottom mouth to top.

    Not the *only* way through -- the verifier is checking the passage is open,
    not that the vapour obeys this line -- but every staging style here has an
    intended route, and if the mesh does not leave it open then either the
    staging was built wrong or it choked the column.
    """
    d = dims(p)
    r_bore = d["r_bore"]
    count = int(round(p["elements"]))
    z0, z1 = stage_span(p)
    pitch = (z1 - z0) / (count + 1)
    pts = [(0.0, 0.0, 0.0), (0.0, 0.0, d["glen_b"] * 0.5),
           (0.0, 0.0, d["z_body0"])]
    style = p["style"]
    if style == "cones":
        # Up the middle, through each hole.
        for k in range(count):
            pts.append((0.0, 0.0, z0 + pitch * (k + 1)))
    elif style == "helix":
        # Following the channel between turns, at the band's mid radius. The
        # band is a radial strip at azimuth a(z), so a point at a(z) + 180 deg
        # is always on the far side of the axis from it, and clear of it.
        turns = helix_turns(p)[0]
        r_mid = r_bore * 0.6
        # Enter from the far side of the axis, below where the band starts: a
        # band whose inner edge is past the centre blocks the axis itself, and
        # the way in is round it.
        pts.append((-r_mid, 0.0, z0 - 1.0))
        for k in range(1, 25):
            u = k / 24.0
            a = 2.0 * math.pi * turns * u + math.pi
            pts.append((r_mid * math.cos(a), r_mid * math.sin(a),
                        z0 + (z1 - z0) * u))
    elif style == "vigreux":
        # Between the fingers, one sector at a time. Each row's fingers are 120
        # degrees apart and the next row's are staggered 60, so the gap the
        # vapour climbs walks round by 60 degrees per row, and the crossing from
        # one sector to the next happens in the clear band above a row: the
        # fingers hang *down*, so between the top of one row and the tips of
        # the next there is a band with nothing in it at any azimuth, and that
        # is the only place in the pitch where a crossing fits.
        length, drop, radius, (tip, _, _) = finger_geometry(p, pitch)
        r_mid = 0.5 * (r_bore + 0.2 + tip)
        band = pitch - length * math.sin(math.radians(drop))
        for k in range(count):
            z = z0 + pitch * (k + 1)
            a0 = math.radians(60.0 * (k + 1))
            a1 = math.radians(60.0 * (k + 2))
            pts.append((r_mid * math.cos(a0), r_mid * math.sin(a0),
                        z + 0.2 * band))
            pts.append((r_mid * math.cos(a1), r_mid * math.sin(a1),
                        z + 0.8 * band))
    elif style == "tesla":
        # Up one side past a baffle's free end, straight across in the clear
        # band above it, then up the other side: the baffles alternate which
        # side is buried, so that is the whole route. The crossing height has
        # to sit above the free end's top corner and below the next baffle's
        # foot, and those two are only a few mm apart -- a diagonal crossing
        # climbs out of that band and runs into the next plate.
        lean = math.radians(louvre_tilt(p))
        length = (r_bore + 0.2 + louvre_reach(p)) / math.cos(lean)
        top = length * math.sin(lean) + (LOUVRE_THICK / 2.0) * math.cos(lean)
        side = r_bore - 0.8
        for k in range(count):
            z = z0 + pitch * (k + 1) + top + (pitch - top) * 0.5
            buried = 1.0 if k % 2 == 0 else -1.0
            pts.append((-buried * side, 0.0, z))
            pts.append((buried * side, 0.0, z))
    pts.append((0.0, 0.0, d["z_body1"]))
    pts.append((0.0, 0.0, d["z_body1"] + p["cap"] + d["glen_t"] * 0.5))
    pts.append((0.0, 0.0, d["z_body1"] + p["cap"] + d["glen_t"]))
    return pts


def path_length(p):
    pts = vapour_path(p)
    return sum(math.dist(pts[i], pts[i + 1]) for i in range(len(pts) - 1))


def build(**p):
    wall, tol = p["wall"], p["tol"]
    d = dims(p)
    r_body = d["r_body"]

    for key in ("bottom", "top"):
        spec = p[key + "_joint"]
        gender = p[key + "_gender"]
        r_mouth, r_tip, glen = J.joint_dims(spec)
        if gender not in ("male", "female"):
            raise ValueError(f"{key} gender is {gender!r}, not male or female")
        if gender == "male" and r_tip - p["tip_wall"] <= 0.0:
            raise ValueError(f"tip wall {p['tip_wall']:g} mm leaves no bore at "
                             f"a {spec} joint")
    for which in ("bottom", "top"):
        glen, free, inward = end_profile(p, which)
        need = max(abs(d["r_bore"] - inward[0]), abs(d["r_body"] - inward[1]))
        # A transition shallower than the printer's overhang angle becomes a
        # ceiling: the bore would open out over a flat shoulder, and a flat
        # shoulder in a bore is support that cannot be got out again.
        need /= math.tan(math.radians(90.0 - p["overhang"]))
        if p["cap"] < need:
            raise ValueError(
                f"the {which} end needs a {need:.1f} mm transition from its "
                f"{p[which + '_joint']} joint to the {p['bore']:g} mm bore, "
                f"but CAP is {p['cap']:g} mm: a transition shallower than "
                f"{90.0 - p['overhang']:g} deg from horizontal will not print")

    n = t.auto_segments(r_body, tol)
    # The body and the staging are built separately and unioned. The staging is
    # attached to the wall and crosses it by BITE, so the union is a crossing
    # rather than a coincidence; left as separate shells they interpenetrate,
    # and a slicer that fills even-odd prints the overlap as a hole.
    body = t.Mesh("body")
    sweep_body(body, p, n)
    staging = t.Mesh("staging")
    info = stage_column(staging, p, n)
    # The staging is unioned with the body, one call, all pieces at once.
    #
    # One call matters: `csg.union` keeps its intermediates as polygons and
    # repairs once at the end, where chaining meshes welds and repairs between
    # every step and multiplies the work. The three helix arms take 12.5 s and
    # come out as one shell with no gaps; chained, the third union could not
    # finish. Staggering the arms' ends (BAND_STAGGER) is what makes it clean:
    # their caps are flat, and three caps in one plane is the coplanar overlap a
    # BSP resolves worst.
    pieces = [mesh_of(staging, sh) for sh in K.shell_triangles(staging.verts,
                                                              staging.tris)]
    if len(pieces) <= UNION_MAX_PIECES:
        mesh = csg.union(body, *pieces, name="column")
        info["unioned"] = len(pieces)
    else:
        mesh = body
        mesh.merge(staging)
        info["unioned"] = 0
        print(f"  note: {len(pieces)} staging pieces is more than the "
              f"{UNION_MAX_PIECES} this will union in one call; they keep "
              f"their own shells, overlapping the wall at their rims")
    # Unioned only when the staging pieces cross the wall and nothing else.
    # With HELIX_FRAC negative every arm reaches past the axis, so the arms
    # cross each other in a region a millimetre wide as well as crossing the
    # wall -- and a boolean between surfaces that cross in that many places
    # leaves cracks: measured, 190 of them. With the arms starting outside the
    # axis they only cross the wall, and the same union comes out with none.
    # The overlap is the lesser evil, so a self-crossing helix keeps its
    # shells and a slicer that fills by winding unions it.
    return mesh, {"n": n, "dims": d, "staging": info, "tube_r": r_body,
                  "shells": len(mesh.shells()),
                  "cracks": csg.crack_length(mesh.verts, mesh.tris),
                  "bbox": mesh.bbox(), "path_length": path_length(p),
                  "vapour_path": vapour_path(p),
                  "sagitta": r_body * (1.0 - math.cos(math.pi / n))}


def main():
    par, out, args = dict(PARAMS), None, list(sys.argv[1:])
    i = 0
    while i < len(args):
        a = args[i]
        if a == "--angle":
            i += 1
            par["overhang"] = float(args[i])
        elif a.startswith("--angle="):
            par["overhang"] = float(a.split("=", 1)[1])
        elif a == "--elements":
            i += 1
            par["elements"] = float(args[i])
        elif a.startswith("--elements="):
            par["elements"] = float(a.split("=", 1)[1])
        elif a in ("--finger-r", "--finger-tip", "--finger-drop", "--finger-taper"):
            key = "finger_" + a.split("-")[-1]
            i += 1
            par[key] = float(args[i])
        elif a in STYLES:
            par["style"] = a
        else:
            out = a
        i += 1
    if out is None:
        out = t.stamped_filename("column", par)
    mesh, s = build(**par)
    d, info = s["dims"], s["staging"]

    desc = (f"Fractionating column, {par['style']} staging. "
            f"{par['bore']:g} mm bore over {par['body']:g} mm, "
            f"{info['elements']} stagings; bottom {par['bottom_joint']} "
            f"{par['bottom_gender']}, top {par['top_joint']} "
            f"{par['top_gender']}. {par['wall']:g} mm wall.")
    t.write_3mf(out, (mesh.name, mesh), title="Fractionating column",
                description=desc, params=par)

    bb = s["bbox"]
    print(f"wrote {out}")
    print(f"  params            : {t.describe(par)}")
    print(f"  param hash        : {t.params_hash(par)}")
    print(f"  printer           : faces down to {par['overhang']:g} deg from "
          f"horizontal print without support, so that is what the staging was "
          f"shaped to beat")
    print(f"  style             : {info['style']}, {info['elements']} stagings"
          + (f", holes {info['hole'] * 2:g} mm" if "hole" in info else "")
          + (f", {info['arms']} arm(s), {info['turns']:g} turns each leaning "
             f"{info['lean']:.1f} deg" if "turns" in info else "")
          + (f", {info['rows']}x{info['per_row']} fingers "
             f"{info['length']:.1f} mm dropping {info['drop']:g} deg, tips "
             f"{info['tip']:.1f} mm off the axis"
             + (f" (the {info['tip_wanted']:g} mm asked for needs a longer "
                f"pitch than {pitch:.1f} mm)" if info.get("clipped") else "")
             if "rows" in info else "")
          + (f", plates tilted {info['tilt']:g} deg" if "tilt" in info else ""))
    print(f"  body              : {d['r_bore'] * 2:g} mm bore, "
          f"{s['tube_r'] * 2:g} mm OD, "
          f"{d['z_body0']:g} to {d['z_body1']:g} mm")
    print(f"  joints            : bottom {par['bottom_joint']} "
          f"{par['bottom_gender']} at z=0, top {par['top_joint']} "
          f"{par['top_gender']} at z={bb[5]:.1f}")
    print(f"  vapour path       : {s['path_length']:.1f} mm to climb "
          f"{bb[5]:.1f} mm of column "
          f"({s['path_length'] / bb[5]:.2f}x the height)")
    print(f"  round facets      : {s['n']} segments, chord {s['sagitta']:.5f} mm")
    print(f"  bounding box      : {bb[3] - bb[0]:.1f} x {bb[4] - bb[1]:.1f}"
          f" x {bb[5] - bb[2]:.1f} mm")
    print(f"  vertices / tris   : {len(mesh.verts)} / {len(mesh.tris)}")
    print(f"  volume            : {mesh.signed_volume():.2f} mm^3")
    print(f"  shells            : {s['shells']} -- "
          + ("body and staging unioned, so no surface lies inside another"
             if s["shells"] == 1 else
             "the staging keeps its own shells, overlapping the wall at its "
             "rim; a slicer that fills by winding unions it"))
    print(f"  cracks            : {s['cracks']:.2e} mm of edge used once")


if __name__ == "__main__":
    main()

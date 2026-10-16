#!/usr/bin/env python3
"""Generate a scoopula: a long round channel in line with a flat handle.

A scoopula is the powder tool that is neither a spoon nor a spatula. Seen from
above it is a leaf: a narrow channel much longer than it is wide, whose long axis
runs *in line* with a flat bar handle. Seen from the end it is a **U**: a round
trough of nearly constant thin wall, deepest at the middle, whose opening
narrows at both ends -- converging into the shaft at the heel, and running out to
a thin scraping edge, the knife, at the toe.

    SCOOP_L = 62   trough length, heel (x=0) to toe, mm; overall 160 mm
    HND_L   = 98   handle length beyond the heel, mm
    RIM_W   = 14   width of the channel's mouth at the rim, at mid-trough
    DEPTH   = 8    channel depth at the middle, floor to rim, mm
    WALL    = 1.4  wall thickness measured across itself; 3.5 beads
    FLOOR   = 1.6  bed face to channel floor, mm: 4 beads, and on the grid
    LEAN    = 0.36 horizontal per vertical: how far the sides splay (20 deg)
    MW_HEEL = 3.5  half-width of the mouth at the heel, mm
    MW_MID  = 7.0  ... and at the middle, which is where RIM_W is measured
    MW_TOE  = 2.2  ... and at the toe, where the channel runs out to its edge
    BACK_WALL = 3  solid heel behind the channel, mm: over two walls thick
    HND_W   = 7.0  handle width at the joint, tapering to HND_W_TIP = 6.0
    HND_T   = 1.2  handle thickness, mm: three beads, thinner than the floor

The ratio that makes it a scoopula rather than a spoon is 62 mm of channel 14 mm
across, an almost round mouth (14 wide over 8 mm deep). The ratio that
makes it *usable* is the wall: 1.4 mm, three and a half beads, thin enough that
the thing reads as drawn metal and light enough that a 98 mm handle still flexes
like the spring it is.

**The round surface is the inside, and that is what lets it print.** A drawn
steel scoopula is formed over a round mandrel: the surface that carries the
powder is an arc and the outside of the same wall is nearly straight. Give the
outside the arc instead, and the trough standing on its convex back touches the
bed along a tangent line, with a 3 mm shoulder under it falling shallower than
45 deg to horizontal -- the whole scoop wanting support. So the body here is a
wedge: a flat bed band, two straight sides splaying outward at `LEAN`, a flat
rim land. The cutter removes a **U**: an arc across the bottom tangent to two
straight flanks, run up past the rim so the cut crosses the surface rather than
meeting it. That puts every face of the trough's inside looking *up*, and leaves
every face of the outside leaning 20 deg off vertical: support-free, and still a
round trough.

Tangency fixes the arc. Requiring it to touch the floor plane and both flanks at
once gives its radius from `LEAN`, `FLOOR` and the half-width at the rim, so the
channel closes off by itself where the rim comes down -- a U that needs more
height than the wedge has simply is not cut there, and that is what leaves the
heel and the toe solid without any extra geometry. The wall comes out `WALL`
measured *across itself* at every height for the same reason as the trapezoid
did: the horizontal offset is the wall times the secant of the lean, which is the
lesson `make_weighboat.py` records.

Three things about the booleans, all learned the hard way on earlier versions:

* The cutter runs `OVER` = 2 mm *past* the rim, so the cut is a crossing rather
  than two surfaces meeting in one plane. Stopping a cutter tangent to the
  surface it cuts leaves the "trough" a sealed void with a 0.04 mm roof over it:
  one shell, zero gaps, a bounding box and a volume that both looked right, and
  a scoop that could never have held anything. `verify_scoopula.py` casts upward
  from inside the mouth at station after station to prove there is nothing over
  it, and walks the mouth at several heights to prove the channel opens.
* The trough's long axis and the handle's are the *same* axis. The first version
  had the trough extruded across the handle, so the scoop faced sideways down the
  shaft -- which the old verifier agreed with, because it restated the builder's
  own geometry instead of measuring the part.
* The handle is unioned *before* the channel is cut, is thinner than the floor,
  and is narrower than the bed band by a margin where it hides, so no two faces
  of the joint coincide and nothing of the handle can stand up into the channel
  as a fin.

This is the *regular* pattern -- handle in the same plane as the back of the
trough -- not the reverse-curved one, where the bowl is bent down off the shaft.
The bend is what lets a steel one dig into a jar while the grip stays above the
rim, and it would put a downward face on the bed here. Keep it flat and the knife
edge at the toe still starts at the bed plane, which is the part that matters.

On the printer: one piece, no support, and 1200 mm^2 of flat band and handle
on the bed. The thing that will give trouble is the 98 mm of 1.2 mm
bar -- a corner of the grip may lift as it cools, so put a brim along the
handle's length, as `print_notes` already says for the filter. Cheap, and the
shape never needs support.

Usage: python3 make_scoopula.py [output.3mf]
"""

import math
import sys

import memcap
import csglib as csg
import printability as PR
import threelib as t

# --- Parameters ------------------------------------------------------------
SCOOP_L = 62.0      # trough length, heel (x=0) to toe (x=62), mm
HND_L = 98.0        # handle length beyond the heel; overall is 160 mm
RIM_W = 14.0        # width of the channel's mouth at the rim, at mid-trough
DEPTH = 8.0         # channel depth at the middle, floor to rim, mm
WALL = 1.4          # wall thickness measured across itself, at every height
FLOOR = 1.6         # bed face to channel floor, mm: 4 beads and on the grid
LEAN = 0.36         # horizontal per vertical: how far the sides splay outward
MW_HEEL = 3.5       # half-width of the mouth at the heel, mm
MW_MID = 7.0        # ... at the middle: the mouth is RIM_W = 2*MW_MID across
MW_TOE = 2.2        # ... at the toe, where the channel runs out to its edge
HND_W = 7.0         # handle width at the joint with the trough
HND_W_TIP = 6.0     # handle width at its free end
HND_T = 1.2         # handle thickness: 3 beads, thinner than the floor

# Profile stations are fractions of SCOOP_L, so the shape survives a rescale.
W_RISE = 0.22       # the mouth widens from heel to mid over this much length
W_FALL = 0.62       # and narrows to the toe from here on
H_RISE = 0.20       # the rim reaches full height this far past the heel
H_FALL = 0.62       # and starts falling toward the knife edge here
F_RISE = 0.72       # the floor starts lifting into the toe wedge here

OVER = 2.0          # how far the channel cutter stands above the rim
MIN_DEPTH = 1.8     # the channel stops where it would be cut shallower than this
LIFT = 1.6          # how far the floor rises toward the toe
H_LIP = 1.6         # outer height at the toe: the knife edge, 4 beads
H_HEEL = 8.0        # outer height at the heel, just under the middle
GRID = 0.4          # every profile height is quantised to this: 2 layers at 0.2
BURIED = 4.0        # how far inside the heel the handle's end stops
BACK_WALL = 3.0     # solid heel behind the channel: over two walls thick
EDGE = 0.8          # least cover between the handle's side and a channel flank
STATIONS = 41       # sections along the trough, 1.55 mm apart
HND_STATIONS = 17
TAPER_BACK = 24.0   # the handle's taper finishes this far before its buried end
N_ARC = 25          # samples across the channel's round bottom
NOZZLE = PR.NOZZLE_DEFAULT   # 0.4 mm: every wall below is sized in beads
TOL = 0.01          # chord deviation asked of any round surface

PARAMS = dict(scoop_len=SCOOP_L, hnd_l=HND_L, rim_w=RIM_W, depth=DEPTH,
              wall=WALL, floor=FLOOR, lean=LEAN, mw_heel=MW_HEEL,
              mw_mid=MW_MID, mw_toe=MW_TOE, hnd_w=HND_W, hnd_w_tip=HND_W_TIP,
              hnd_t=HND_T, over=OVER, min_depth=MIN_DEPTH, lift=LIFT,
              back_wall=BACK_WALL,
              h_lip=H_LIP, h_heel=H_HEEL, grid=GRID, buried=BURIED,
              stations=STATIONS, n_arc=N_ARC, nozzle=NOZZLE, tol=TOL)


def default_out():
    return t.stamped_filename("scoopula", PARAMS)


# --- profiles --------------------------------------------------------------
def smooth(u):
    """C1 smoothstep, flat at both ends, so a plateau joins without a kink."""
    u = min(1.0, max(0.0, u))
    return u * u * (3.0 - 2.0 * u)


def span(x, lo, mid, hi, rise, fall, length):
    """A profile: `lo` at x=0 rising to `mid` by `rise`, holding, then falling
    to `hi` at x=`length`, from `fall` on."""
    return lo + (mid - lo) * smooth(x / rise) - (mid - hi) * smooth(
        (x - fall) / (length - fall))


def snap(z, grid):
    """Quantise a height onto the layer grid: no feature gets sliced in half.

    `floor`-with-a-half-plus-epsilon, because 1.4 / 0.4 is 3.4999999999999996 in
    binary and would round down a whole 0.2 mm: the floor came out 1.2 mm, which
    is both under the asked thickness and exactly equal to the handle's top, the
    coplanar joint this file elsewhere goes to some trouble to avoid.
    """
    return round(math.floor(z / grid + 0.5 + 1e-9) * grid, 9)


def on_grid(z, grid, name):
    """Refuse a feature height that is not a whole number of layers."""
    if abs(snap(z, grid) - z) > 1e-9:
        raise ValueError(f"{name} {z:g} mm is not a whole number of "
                         f"{grid:g} mm layers -- it is {z / grid:g} -- and "
                         f"snapping it silently would move the geometry by "
                         f"{abs(snap(z, grid) - z) * 1000:.0f} um; ask "
                         f"{math.floor(z / grid) * grid:g} or "
                         f"{math.ceil(z / grid) * grid:g} mm")
    return z


def section(x, z0, z1, h0, h1):
    """A trapezoid at station x: half-width h0 at z0, h1 at z1."""
    return [(x, -h0, z0), (x, h0, z0), (x, h1, z1), (x, -h1, z1)]


def same(a, b, eps=1e-9):
    return all(abs(p - q) <= eps for pa, pb in zip(a, b)
               for p, q in zip(pa, pb))


def dedupe(secs):
    """Drop any section identical to the one before it. A loft cannot turn a
    section into itself; the zero-area run it makes is pure noise."""
    out = []
    for s in secs:
        if not out or not same(s, out[-1]):
            out.append(s)
    return out


def handle_halfwidth(x, p):
    """Half-width of the flat handle at station x, tapering to its free end.

    Zero past the buried end: beyond the joint the trough stands on its own and
    its bed band is set by the channel, not by the handle.
    """
    if x > p["buried"]:
        return 0.0
    tip = 0.5 * p["hnd_w_tip"]
    joint = 0.5 * p["hnd_w"]
    return tip + (joint - tip) * smooth(
        (x + p["hnd_l"]) / (p["hnd_l"] + p["buried"] - TAPER_BACK))


def proportions(p):
    """Everything the section follows from, at the middle of the trough.

    The channel is a U: an arc tangent to the floor plane at `floor` and tangent
    to two straight flanks that lean `lean`, so the flank line is y = bi + lean*z
    and the radius that touches both the floor and a flank is

        rb = (lean*floor + bi) / (sqrt(1+lean^2) - lean)

    `phi` is the half-angle the arc spans from its lowest point: the radius to a
    tangent point stands off the vertical by exactly the flank's lean. `inset` is
    the *horizontal* distance between a flank and the outside face it is offset
    from -- the wall across a leaning face times that secant, so a wall measured
    *across itself* comes out exactly `wall`.
    """
    if abs(p["lean"]) >= 1.0:
        raise ValueError(f"a {p['lean']:g} lean is "
                         f"{math.degrees(math.atan(p['lean'])):.0f} deg off "
                         f"vertical: the outside would fall shallower than 45 "
                         f"deg to horizontal and need support")
    if p["mw_mid"] * 2.0 != p["rim_w"]:
        raise ValueError(f"rim_w {p['rim_w']:g} is not the mouth: mw_mid is "
                         f"{p['mw_mid']:g}, so the mouth is "
                         f"{2 * p['mw_mid']:g} mm across")
    s1 = math.hypot(1.0, p["lean"])
    sec = s1 - p["lean"]
    inset = p["wall"] * s1
    h_mid = p["floor"] + p["depth"]
    aw_mid = p["mw_mid"]
    bi_mid = aw_mid - p["lean"] * h_mid
    if bi_mid <= 0.0:
        raise ValueError(f"a {p['lean']:g} lean under a {h_mid:g} mm rim needs "
                         f"{p['lean'] * h_mid:.2f} mm of half-width and the "
                         f"mouth is only {aw_mid:g}: the flanks cross below "
                         "the rim")
    phi = math.atan2(1.0, p["lean"])
    rb_mid = (p["lean"] * p["floor"] + bi_mid) / sec
    z_t_mid = p["floor"] + rb_mid * (1.0 - math.cos(phi))
    if z_t_mid >= h_mid:
        raise ValueError(f"the round bottom reaches its flanks at "
                         f"{z_t_mid:.2f} mm and the rim is at {h_mid:g}: a U "
                         f"{2 * aw_mid:g} mm wide does not fit under a "
                         f"{p['depth']:g} mm depth; deepen the channel or lean "
                         "the flanks less")
    if rb_mid <= 2.0 * p["nozzle"]:
        raise ValueError(f"the channel's bottom would round over a radius of "
                         f"{rb_mid:.2f} mm, under two {p['nozzle']:g} mm beads; "
                         "the cutter cannot leave a floor that narrow")
    if p["floor"] >= rb_mid:
        raise ValueError(f"a {p['floor']:g} mm floor is thicker than the "
                         f"{rb_mid:.2f} mm radius of the channel's own bottom: "
                         "the arc would have to reach under the bed")
    if h_mid <= p["h_heel"] or h_mid <= p["h_lip"]:
        raise ValueError(f"the middle stands {h_mid:g} mm, the heel "
                         f"{p['h_heel']:g} and the toe {p['h_lip']:g}: a trough "
                         "is deepest in the middle or it is a ramp")
    hw_mid = bi_mid + inset
    rim_out = 2.0 * (hw_mid + p["lean"] * h_mid)
    return dict(inset=inset, sec=sec, s1=s1, phi=phi,
                lean_deg=math.degrees(math.atan(p["lean"])),
                dip=90.0 - math.degrees(math.atan(p["lean"])),
                h_mid=h_mid, aw_mid=aw_mid, bi_mid=bi_mid, rb_mid=rb_mid,
                z_t_mid=z_t_mid, hw_mid=hw_mid, rim_out=rim_out,
                back_w=2.0 * hw_mid, land=0.5 * rim_out - aw_mid)


def ring_at(x, p, pr, h, f, aw):
    """The channel's cutter section at station x, CCW, run past the rim.

    `N_ARC` samples of the round bottom from left to right, then up the right
    flank to the cutter's top, across it, and back down the left flank to close.
    Constant point count, so the loft never has to invent or drop a section.
    """
    L = p["lean"]
    top = h + p["over"]
    bi = aw - L * h
    rb = (L * f + bi) / pr["sec"]
    zc = f + rb
    phi = pr["phi"]
    arc = []
    for i in range(p["n_arc"]):
        a = -phi + 2.0 * phi * i / (p["n_arc"] - 1)
        arc.append((x, rb * math.sin(a), zc - rb * math.cos(a)))
    aw_top = bi + L * top
    return arc + [(x, aw_top, top), (x, -aw_top, top)]


def mouth_area(p, pr, h, f, aw):
    """The channel's cross-section up to the rim, by shoelace: capacity that is
    actually in the part, not the cutter's own volume above it."""
    L = p["lean"]
    poly = [(q[1], q[2]) for q in ring_at(0.0, p, pr, h, f, aw)][:-2]
    poly += [(aw, h), (-aw, h)]                    # flanks up to the rim plane
    return 0.5 * abs(sum(y * z2 - z * y2
                         for (y, z), (y2, z2) in zip(poly, poly[1:] + poly[:1])))


def cuttable(p, pr, x, h, f, aw):
    """Whether a U of that shape actually fits under the rim at this station.

    And not within `back_wall` of the heel: the trough has to be closed at the
    back or the powder simply runs out behind the scoop into the handle.
    """
    if x < p["back_wall"]:
        return False
    if h - f < p["min_depth"] or aw <= p["nozzle"]:
        return False
    bi = aw - p["lean"] * h
    if bi <= 0.0:
        return False
    rb = (p["lean"] * f + bi) / pr["sec"]
    if rb <= 2.0 * p["nozzle"]:
        return False
    z_t = f + rb * (1.0 - math.cos(pr["phi"]))
    return z_t < h - p["nozzle"]


def stations(p, pr):
    """[(x, mouth half-width, rim height, floor height, bed half-width)].

    The mouth and the rim are what the eye reads, so those are the profiles. The
    bed band follows from them -- `inset` outside the flank, except where the
    buried handle needs more room, which is what keeps the joint's faces from
    coinciding. Both heights come back on the layer grid; the floor lifts toward
    the toe, and that is what makes the nose read as an edge and not a brick.
    """
    length = p["scoop_len"]
    out = []
    for i in range(p["stations"]):
        x = i * length / (p["stations"] - 1)
        aw = span(x, p["mw_heel"], pr["aw_mid"], p["mw_toe"],
                  W_RISE * length, W_FALL * length, length)
        h = snap(span(x, p["h_heel"], pr["h_mid"], p["h_lip"],
                      H_RISE * length, H_FALL * length, length), p["grid"])
        f = snap(p["floor"] + p["lift"] * smooth((x - F_RISE * length) /
                                                 ((1.0 - F_RISE) * length)),
                 p["grid"])
        hw = max(aw - p["lean"] * h + pr["inset"],
                 handle_halfwidth(x, p) + max(EDGE, p["nozzle"]) + 0.15)
        out.append((x, aw, h, f, hw))
    return out


def build(**p):
    if p["wall"] < 2.0 * p["nozzle"]:
        raise ValueError(
            f"a {p['wall']:g} mm wall is under two "
            f"{p['nozzle']:g} mm beads; a scoop that carries powder needs "
            f"{2.0 * p['nozzle']:g} mm at least")
    if p["floor"] < 2.0 * p["nozzle"]:
        raise ValueError(f"a {p['floor']:g} mm floor is under two "
                         f"{p['nozzle']:g} mm beads; the band that holds the "
                         f"whole trough down needs {2.0 * p['nozzle']:g} mm")
    if p["hnd_t"] >= p["floor"]:
        raise ValueError(f"a {p['hnd_t']:g} mm handle is not thinner than the "
                         f"{p['floor']:g} mm floor: its top face would either "
                         f"coincide with the channel's floor plane or stand in "
                         f"the channel as a fin")
    if p["h_lip"] < 2.0 * p["nozzle"]:
        raise ValueError(f"a {p['h_lip']:g} mm knife edge is under two "
                         f"{p['nozzle']:g} mm beads of height; the lowest lip "
                         f"that prints as itself is {2.0 * p['nozzle']:g} mm")
    if p["back_wall"] < 2.0 * p["wall"]:
        raise ValueError(f"a {p['back_wall']:g} mm back wall is under two "
                         f"{p['wall']:g} mm walls: the heel the whole trough "
                         f"opens onto wants {2.0 * p['wall']:g} mm of its own")
    if p["mw_heel"] > p["mw_mid"] or p["mw_toe"] > p["mw_mid"]:
        raise ValueError(f"the opening is {p['mw_heel']:g}/{p['mw_mid']:g}/"
                         f"{p['mw_toe']:g} mm of half-width at heel, middle and "
                         "toe: a trough is widest through its middle or it is a "
                         "funnel")
    for name in ("floor", "h_lip", "h_heel"):
        on_grid(p[name], p["grid"], name)
    on_grid(p["floor"] + p["depth"], p["grid"], "the middle's rim height")
    pr = proportions(p)
    phi = pr["phi"]
    need = math.ceil(phi / math.acos(max(-1.0, min(1.0, 1.0 - p["tol"] /
                                                   pr["rb_mid"]))) + 1)
    sag = pr["rb_mid"] * (1.0 - math.cos(phi / (p["n_arc"] - 1)))
    if p["n_arc"] < need:
        raise ValueError(f"{p['n_arc']:g} samples put the channel's round bottom "
                         f"{sag:.3f} mm off its own arc, over the "
                         f"{p['tol']:g} mm this file asks of a curve: "
                         f"{need:g} samples bring it under")
    st = stations(p, pr)

    # --- the solid the channel is cut out of -------------------------------
    # A wedge: flat bed band, sides splaying outward at LEAN, flat top at the rim
    # height. Nothing on it faces down more steeply than the lean allows.
    body = t.Mesh("body")
    t.loft(body, dedupe([section(x, 0.0, h, hw, hw + p["lean"] * h)
                         for x, aw, h, f, hw in st]))
    if body.signed_volume() <= 0.0:
        raise ValueError("the body loft came out inside-out")

    # --- the channel -------------------------------------------------------
    rings = dedupe([ring_at(x, p, pr, h, f, aw)
                    for x, aw, h, f, hw in st
                    if cuttable(p, pr, x, h, f, aw)])
    if len(rings) < 2:
        raise ValueError(f"nothing to cut: only {len(rings)} of "
                         f"{len(st)} stations have room for a U "
                         f"{2 * pr['rb_mid']:.2f} mm across the bottom under a "
                         f"{pr['h_mid']:g} mm rim over a {p['floor']:g} mm "
                         f"floor; min_depth {p['min_depth']:g} has to be under "
                         f"depth {p['depth']:g} minus floor {p['floor']:g}")
    cavity = t.Mesh("cavity")
    t.loft(cavity, rings)
    if cavity.signed_volume() <= 0.0:
        raise ValueError("the cavity loft came out inside-out")

    # --- the handle --------------------------------------------------------
    # On the bed with the trough, in line with its axis, buried under the solid
    # heel, and narrower than the bed band by a margin all the way along.
    hx0, hx1 = -p["hnd_l"], p["buried"]
    handle = t.Mesh("handle")
    t.loft(handle, dedupe([section(
        hx0 + i * (hx1 - hx0) / (HND_STATIONS - 1), 0.0, p["hnd_t"],
        handle_halfwidth(hx0 + i * (hx1 - hx0) / (HND_STATIONS - 1), p),
        handle_halfwidth(hx0 + i * (hx1 - hx0) / (HND_STATIONS - 1), p))
        for i in range(HND_STATIONS)]))
    if handle.signed_volume() <= 0.0:
        raise ValueError("the handle loft came out inside-out")

    cover = min(hw - handle_halfwidth(x, p) for x, aw, h, f, hw in st
                if x <= p["buried"] + 1e-9)
    if cover < max(EDGE, p["nozzle"]):
        raise ValueError(f"the handle has only {cover:.2f} mm of cover under "
                         f"the trough where it hides, needs "
                         f"{max(EDGE, p['nozzle']):g} mm so that no face of the "
                         "joint coincides with a face of the channel")

    # Union first, cut last: anything of the handle above the floor is removed
    # by the channel rather than left standing in it as a fin.
    mesh = csg.combine([("start", body), ("union", handle),
                        ("subtract", cavity)], name="scoopula")

    # Capacity in closed form, taken to the rim and not to the cutter's top: the
    # U's area at each open station, trapezoid-integrated along x. Derived here
    # so the reported figure is not the mesh describing itself.
    xs, area = [], []
    for x, aw, h, f, hw in st:
        if not cuttable(p, pr, x, h, f, aw):
            continue
        xs.append(x)
        area.append(mouth_area(p, pr, h, f, aw))
    capacity = sum(0.5 * (area[i] + area[i + 1]) * (xs[i + 1] - xs[i])
                   for i in range(len(xs) - 1))

    return mesh, {"inset": pr["inset"], "lean_deg": pr["lean_deg"],
                  "dip": pr["dip"], "hw_mid": pr["hw_mid"],
                  "h_mid": pr["h_mid"], "rb_mid": pr["rb_mid"],
                  "z_t_mid": pr["z_t_mid"], "sag": sag,
                  "rim_out": pr["rim_out"],
                  "back_w": pr["back_w"], "land": pr["land"],
                  "aw_mid": pr["aw_mid"], "cover": cover,
                  "cavity_from": xs[0], "cavity_to": xs[-1],
                  "capacity": capacity, "n_open": len(xs),
                  "shells": len(mesh.shells()),
                  "cracks": csg.crack_length(mesh.verts, mesh.tris),
                  "gaps": len(csg.real_gaps(mesh.verts, mesh.tris)),
                  "bbox": mesh.bbox(),
                  "bed": PR.bed_contact(mesh.verts, mesh.tris)}


def main():
    memcap.cap()
    out = sys.argv[1] if len(sys.argv) > 1 else default_out()
    mesh, s = build(**PARAMS)
    p = PARAMS
    desc = (f"Scoopula: a {p['scoop_len']:g} mm round channel {p['rim_w']:g} mm "
            f"across the mouth and {p['depth']:g} mm deep with a {p['wall']:g} "
            f"mm wall, in line with a {p['hnd_l']:g} mm flat handle.")
    t.write_3mf(out, (mesh.name, mesh), title="Scoopula",
                description=desc, params=p)

    bb = s["bbox"]
    print(f"wrote {out}")
    print(f"  params            : {t.describe(p)}")
    print(f"  param hash        : {t.params_hash(p)}")
    print(f"  channel           : {2 * s['aw_mid']:.2f} mm mouth, "
          f"{p['depth']:g} mm rim to floor, round bottom R="
          f"{s['rb_mid']:.2f} mm tangent to flanks leaning "
          f"{s['lean_deg']:.1f} deg off vertical, {s['n_open']} stations open")
    print(f"  wall and floor    : {p['wall']:g} mm across the wall at every "
          f"height (cut {s['inset']:.3f} mm inside the face), {p['floor']:g} mm "
          f"of floor under the channel, flat rim land {s['land']:.2f} mm a side")
    print(f"  outside           : {s['rim_out']:.2f} mm over the rim, "
          f"{s['back_w']:.2f} mm over the flat bed band at the middle, "
          f"{2 * (p['mw_toe'] + s['inset']) - 2 * p['lean'] * p['h_lip']:.2f} mm "
          f"across the nose")
    print(f"  open from         : x={s['cavity_from']:.1f} to "
          f"x={s['cavity_to']:.1f} mm; {s['cavity_from']:.1f} mm of heel and "
          f"{p['scoop_len'] - s['cavity_to']:.1f} mm of solid toe wedge")
    print(f"  knife edge        : {p['h_lip']:g} mm high at the toe, whose "
          f"lowest edge is the bed plane, so it enters a pile at floor level")
    print(f"  handle            : {p['hnd_l']:g} mm beyond the heel, "
          f"{p['hnd_w']:g} mm at the joint tapering to {p['hnd_w_tip']:g} mm, "
          f"{p['hnd_t']:g} mm thick, on the trough's own axis, buried "
          f"{p['buried']:.1f} mm with {s['cover']:.2f} mm of cover")
    print(f"  size              : {bb[3] - bb[0]:.1f} x "
          f"{bb[4] - bb[1]:.1f} x {bb[5] - bb[2]:.1f} mm")
    print(f"  volume            : {mesh.signed_volume():.1f} mm^3 of material")
    print(f"  capacity          : {s['capacity'] / 1000.0:.2f} ml to the rim")
    print(f"  bed contact       : {s['bed']:.1f} mm^2")
    print(f"  print orientation : flat on its back: it stands on a strip "
          f"{s['back_w']:.1f} mm wide, every outside face leans "
          f"{s['lean_deg']:.0f} deg off vertical ({s['dip']:.1f} deg from "
          f"horizontal) and the whole round inside looks up. No support")
    print(f"  vertices / tris   : {len(mesh.verts)} / {len(mesh.tris)}")
    print(f"  shells            : {s['shells']} (body unioned with the handle, "
          f"channel cut out)")
    print(f"  gaps              : {s['gaps']} real, {s['cracks']:.2e} mm of "
          f"edge used once")


if __name__ == "__main__":
    main()

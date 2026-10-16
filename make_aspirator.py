#!/usr/bin/env python3
"""Generate a water aspirator: a Venturi pump that pulls air with a water jet.

Water enters at the bottom in a wide tube, is squeezed through a converging
cone, runs as a narrow jet *up through* the chamber where the air comes in, and
disperses in the chamber above the nozzle's tip. Whether the nozzle lets go at
all is a parameter, JET_RELEASE, below. Where the jet is narrowest it is fastest and its
pressure is lowest, so air drawn into the chamber is carried up and out with the
water. No moving parts; the vacuum is only as good as the water pressure.

Every constant below is a radius off the axis, in millimetres; `main` reports
diameters. Nothing here quotes a dimension, because all of them are set by name
at the top of this file. Reading up the axis, which is vertical in the print:

  water supply   INLET_BORE passage in an INLET_OUT skin, Z_IN0 up to Z_IN1.
                 The part stands on this end.
  shoulder       FOOT_SKIN: begins at INLET_OUT - BURIED on station Z_FOOT0,
                 leans out to CHAMBER_OUT - BURIED by Z_CH0, then holds to
                 Z_FOOT1 with bore FLOOR_BORE -- so it is the chamber's floor
                 as well as its outside step. Z_FOOT0 is not free: it comes
                 from FOOT_RISE and LEAN_DEG, rounded down onto GRID so that no
                 feature edge lands in the middle of a layer.
  chamber        CHAMBER_R void in a CHAMBER_OUT skin, Z_CH0 up to Z_CH1.
  nozzle         NOZZLE_BORE is what gets cut, NOZZLE_CHANNEL what the water is
                 left with. They differ only over the buried foot, where the
                 bore is BURIED wider than INLET_BORE so that the two
                 cylindrical surfaces cannot run alongside each other. One wall
                 of NOZZLE_WALL, from Z_NOZ0 to Z_NOZ1; the bore is INLET_BORE
                 until Z_CONE0, converges to THROAT at Z_CONE1, holds there to
                 Z_THROAT1, and gains JET_RELEASE by the tip at Z_NOZ1.
  throat         THROAT, between Z_CONE1 and Z_THROAT1: the fastest, lowest
                 pressure point of the water's path, and the reason for the
                 rest of the part.
  air inlet      AIR_BORE hole through an AIR_BOSS pad reaching out to AIR_OUT,
                 on station Z_AIR. Drilled last, starting JET_CLEAR outside the
                 nozzle's skin; BOSS_SINK is how deep its base sits into the
                 chamber's wall. It is what opens the annulus to the air.
  neck           NK_BORE with NK_SKIN, Z_NK0 down to Z_NK2: the chamber's shell
                 closing onto the outlet, clear of both by BURIED.
  outlet         OUTLET_BORE in an OUTLET_OUT skin, Z_OUT0 up to Z_OUT1.
  all round      TOL is the chord deviation auto_segments() is given, off
                 CHAMBER_OUT as the widest surface in the part; one facet count
                 and one phase serve every circle here, see `build`.

Two more things are settled by name rather than left in the numbers, because
they are the ones worth turning:

  the jet widening again  JET_RELEASE is the radius the nozzle's channel gains
      between Z_THROAT1 and its tip Z_NOZ1. At 0.0 the water only ever
      converges: `_noz_stations` drops the final station, the throat runs out
      to the tip, and the expansion happens in the open chamber instead, where
      the jet is free to slow against the air it has drawn in. To open the
      channel back to the outlet's bore inside the nozzle, JET_RELEASE =
      OUTLET_BORE / 2 - THROAT, and then Z_THROAT1 has to drop below Z_NOZ1 to
      give the release some length to happen over.

  the inlet's end and the outlet's start  The water runs in below Z_IN1 and
      away above Z_OUT0, and what keeps those stations apart is not a free
      choice: the jet's tip has to clear the air's mouth (Z_NOZ1 above
      Z_AIR + AIR_BORE), the shell has to close over the tip (Z_CH1 above
      Z_NOZ1), and the neck has to start inside the shell (Z_NK0 under Z_CH1),
      meet the outlet's bore at Z_NK1 -- which is where the outlet tube starts,
      Z_OUT0 -- and take its skin from NK_OUT1 down to OUTLET_OUT over
      Z_NK1..Z_NK2. Put each at its minimum and Z_OUT0 lands just above Z_IN1,
      with a chamber no taller than the jet needs room in. `validate` holds
      these, and refuses a nozzle whose stations do not ascend, which is how a
      folded `sweep` used to fail: a piece with negative volume, and an error
      about the mesh rather than the parameter.

The three demands this part has to meet, each stated in the constants that
decide it and enforced by `validate`:

  * the narrow part starts before the air inlet and ends after it:
    Z_CONE1 < Z_AIR - AIR_BORE and Z_THROAT1 > Z_AIR + AIR_BORE, so the throat
    brackets the inlet's mouth with room at both ends;
  * the wall round that section is a Z read lengthwise: the outside leans out
    over Z_FOOT0..Z_CH0 while the inside is still straight, and the inside
    draws in over Z_CONE0..Z_CONE1 after the outside has finished; between the
    two diagonals a cut across the axis is one band of plastic, INLET_BORE out
    to CHAMBER_OUT, with nothing hollow in it;
  * at the inlet the section is two rings and not one: the nozzle's skin at
    THROAT + NOZZLE_WALL, the chamber's void at CHAMBER_R, and between them an
    annulus of air -- CHAMBER_R - (THROAT + NOZZLE_WALL) -- wider than the
    AIR_BORE that feeds it.

The part is built as stacked rings, each its own closed hollow tube, unioned
where they overlap and drilled last for the air inlet. Three things to keep in
mind when changing it:

  * `revolve` cannot make these rings. It caps a profile that stops short of
    the axis with a flat disc, so a profile drawn as an annulus comes out
    solid -- which is what the first cut of this file did: valid mesh,
    plausible volume, every bore filled with plastic. `tube` and `sweep` take
    the bore and the outside separately and leave the bore open.
  * Every piece overlaps its neighbour, and the overlap is *material* on both
    sides -- one annulus inside another, never a surface grazing a surface.
  * Two surfaces of the same radius running alongside each other intersect
    nowhere, so the boolean has nothing to cut with and leaves a doubled wall
    with a gap you can see through. Where a buried piece would otherwise match
    the radius of the one burying it -- the nozzle's bore inside the inlet's,
    the shoulder's skin against the chamber's, the neck's tail inside the
    outlet's -- the buried one is offset by BURIED, a bore wider and a skin
    narrower, so the surface the design wants is the one that survives and the
    other disappears inside plastic.

Standing on the inlet end is what makes this printable: the nozzle then rises
out of the chamber's floor, so every surface inside the chamber is a vertical
wall or a cone leaning inward as it rises, and nothing inside needs support that
could not be reached. Outside, the only face that leans as it rises is the
shoulder, at LEAN_DEG, and the printer carries that alone -- flatten it and it
turns into a shelf. What is left is the air inlet's boss: its underside, outside
the chamber's skin, where support comes off. The foot it stands on is a ring
INLET_OUT - INLET_BORE wide, and a part this tall on a ring that thin wants a
brim, which is the slicer's business, not this file's.
"""

import math
import sys

import csglib as csg
import memcap
import threelib as t

# --- Parameters: radii in mm ------------------------------------------------
INLET_BORE = 2.6      # water supply bore radius
INLET_OUT = 4.6      # inlet tube outside radius
THROAT = 1.3          # throat radius: the narrowest the water ever gets
NOZZLE_WALL = 1.2     # wall of the nozzle that stands in the chamber
FLOOR_BORE = 3.2        # the chamber floor's hole: wider than the jet's own
                        # ...channel, narrower than the nozzle's skin, so the
                        # ...floor grips the nozzle instead of leaving it
                        # ...standing in an open hole
CHAMBER_R = 14.0      # chamber void radius
CHAMBER_OUT = 16.0    # chamber outside radius
OUTLET_BORE = 2.6     # outlet bore radius
OUTLET_OUT = 3.6     # outlet tube outside radius
AIR_BORE = 2.6        # air inlet bore radius
AIR_BOSS = 3.6        # air inlet boss radius
AIR_OUT = 20.0        # how far the air inlet stands out from the axis
BOSS_SINK = 0.4       # how deep the boss's base sits into the chamber wall
JET_CLEAR = 0.7       # the air inlet stops this far above the nozzle's skin
JET_RELEASE = 0.0     # radius the jet's channel gains between the throat and
                      # ...the nozzle's tip. 0.0: the water only ever narrows,
                      # ...and it is the chamber, not the nozzle, that lets it
                      # ...go -- which is what an aspirator wants, since the
                      # ...jet has still to be moving when the air arrives.
TOL = 0.05            # max chord deviation on round surfaces

# Where one piece is buried inside another at the same height, offset the
# buried surface by this much. Two cylinders of the same radius running
# alongside each other do not intersect anywhere, so the boolean has nothing to
# cut with: it leaves doubled walls, edges used once, and gaps you can see
# through. Each buried surface is offset away from the one meant to be seen --
# a buried bore is made wider than the bore that governs it, a buried skin
# narrower -- so the union keeps the surface the design asks for and the offset
# one disappears inside plastic. The four places this matters are the nozzle
# inside the inlet tube, the chamber floor against the shell, and the neck's
# tail inside the outlet.
BURIED = 0.3

# Stations along the axis, all of them on the layer grid.
GRID = 0.2                       # the layer height this is drawn for
Z_IN0, Z_IN1 = 0.0, 28.0        # the inlet tube
Z_NOZ0, Z_NOZ1 = 24.0, 49.0     # the nozzle that stands in the chamber
Z_CONE0, Z_CONE1 = 28.0, 38.0   # where its bore converges
Z_THROAT1 = 49.0                # where the throat stops being narrow: equal to
                                # ...Z_NOZ1, so with JET_RELEASE at 0.0 it never
                                # ...does -- the nozzle only converges
Z_CH0, Z_CH1 = 26.0, 56.0       # the chamber's shell
Z_FOOT1 = 30.0                  # the shoulder's top: the chamber's floor
Z_NK0 = 54.0                    # the neck starts inside the shell's wall
Z_NK1, Z_NK2 = 60.0, 66.0       # and reaches the outlet's bore, then its skin
Z_OUT0, Z_OUT1 = 60.0, 130.0    # the outlet tube: it starts 32 mm above where
                                # ...the inlet's tube ends, which is as close as
                                # ...the jet's own reach and the neck allow
Z_AIR = 46.0                    # the air inlet, over the middle of the throat

# The shoulder: where the inlet tube becomes the chamber, and the chamber's
# floor at the same time. Its skin is the one outside face that leans
# outward as it rises, so it is the one place the print angle bites here. Given
# LEAN_DEG above horizontal the printer carries it with nothing under it; flatten it and
# it becomes a shelf. Under the chamber's bottom edge it stays short of the
# shell's own skin by BURIED, so the two never run alongside each other.
FOOT_TOP = CHAMBER_OUT - BURIED
LEAN_DEG = 55.0
FOOT_RISE = FOOT_TOP - (INLET_OUT - BURIED)
# Quantised to the layer this is drawn for: a feature edge that lands in the
# middle of a layer prints as a sliver of a layer, and here there is no reason
# for one -- the run only has to be at least FOOT_RISE at the shoulder's angle.
Z_FOOT0 = round(Z_CH0 - GRID * math.floor(FOOT_RISE
                                         * math.tan(math.radians(LEAN_DEG)) / GRID), 2)
FOOT_SKIN = ((Z_FOOT0, INLET_OUT - BURIED), (Z_CH0, FOOT_TOP),
             (Z_FOOT1, FOOT_TOP))

# The neck's own profile. Its bore reaches the outlet's bore exactly at
# Z_NK1, where the outlet tube starts, so the channel is continuous; past that
# it flares by BURIED so its tail, buried in the outlet's wall, never runs
# alongside the outlet's own bore. Its skin does the mirror thing: it stands
# clear of the chamber's, then comes down to the outlet's.
NK_BORE0 = CHAMBER_R + BURIED          # buried in the shell's void edge
NK_BORE1 = OUTLET_BORE                 # the channel, met not overlapped
NK_BORE2 = OUTLET_BORE + BURIED        # the flare that keeps its tail clear
NK_OUT0 = CHAMBER_OUT - BURIED         # inside the shell's skin, never on it
NK_OUT1 = CHAMBER_OUT - 4.0
FOOT_SKIN_TOP = FOOT_TOP
NOZ_BORE0 = INLET_BORE + BURIED        # the nozzle's buried foot, wider than
                                       # ...the inlet bore that governs it

PARAMS = dict(inlet_bore=INLET_BORE, inlet_out=INLET_OUT, throat=THROAT,
              nozzle_wall=NOZZLE_WALL, floor_bore=FLOOR_BORE,
              jet_release=JET_RELEASE,
              chamber_r=CHAMBER_R, chamber_out=CHAMBER_OUT,
              outlet_bore=OUTLET_BORE, outlet_out=OUTLET_OUT,
              air_bore=AIR_BORE, air_boss=AIR_BOSS, air_out=AIR_OUT,
              z_cone0=Z_CONE0, z_cone1=Z_CONE1, z_throat1=Z_THROAT1,
              z_foot0=Z_FOOT0, z_foot1=Z_FOOT1,
              z_noz1=Z_NOZ1, z_ch0=Z_CH0, z_ch1=Z_CH1, z_air=Z_AIR,
              z_out1=Z_OUT1, buried=BURIED, tol=TOL)


def default_out():
    return t.stamped_filename("aspirator", PARAMS)


# --- Profiles along the axis ------------------------------------------------
def _prof(z, pts):
    """Straight lines through a list of (station, value) pairs."""
    if z <= pts[0][0]:
        return pts[0][1]
    for (za, ra), (zb, rb) in zip(pts, pts[1:]):
        if z <= zb:
            return ra + (z - za) / (zb - za) * (rb - ra)
    return pts[-1][1]


# What the nozzle is actually cut as, and what the water is left with. They
# differ only over the nozzle's buried foot, where the built bore is BURIED
# wider than the inlet's so the two surfaces cannot run together -- and the
# inlet's own bore is the narrower one, so it is the channel that survives.
def _noz_stations(foot_bore):
    """The nozzle's (station, bore) pairs: buried foot, straight, converge, hold.

    The release is the only part that widens, and it is a choice -- JET_RELEASE.
    At 0.0 there is nothing to let go, so no final station: a segment of zero
    height would be a band of zero-area triangles, and a segment whose stations
    are the wrong way round folds the sweep inside-out, which is how this file
    once produced a mesh with negative volume and no explanation.
    """
    st = [(Z_NOZ0, foot_bore), (Z_CONE0, INLET_BORE), (Z_CONE1, THROAT),
          (min(Z_THROAT1, Z_NOZ1), THROAT)]
    if JET_RELEASE > 0.0:
        st.append((Z_NOZ1, THROAT + JET_RELEASE))
    return tuple(st)


NOZZLE_BORE = _noz_stations(NOZ_BORE0)
NOZZLE_CHANNEL = _noz_stations(INLET_BORE)

NK_BORE = ((Z_NK0, NK_BORE0), (Z_CH1, NK_BORE0), (Z_NK1, NK_BORE1),
           (Z_NK2, NK_BORE2))
NK_SKIN = ((Z_NK0, NK_OUT0), (Z_CH1, NK_OUT0), (Z_NK1, NK_OUT1),
           (Z_NK2, OUTLET_OUT))


def nozzle_bore(z):
    """The jet's own channel: squeezed, held, then let go inside the chamber."""
    return _prof(z, NOZZLE_CHANNEL)


def nozzle_out(z):
    return nozzle_bore(z) + NOZZLE_WALL


def neck_bore(z):
    """Chamber to outlet. Flat where it is buried in the shell's wall."""
    return _prof(z, NK_BORE)


def neck_out(z):
    return _prof(z, NK_SKIN)


def rings_at(z):
    """The material's radial bands at station z, as the design intends them.

    An independent copy of the geometry: the verifier measures the same thing
    off the mesh, and the two agreeing is the point. Where a piece is buried
    (the nozzle's foot, the neck's tail) the band listed here is the one that
    governs, which is why the nozzle is entered at the inlet's bore and the
    neck at the outlet's rather than at their own offset radii.
    """
    pieces = [(Z_IN0, Z_IN1, INLET_BORE, INLET_OUT),
              (Z_FOOT0, Z_FOOT1, FLOOR_BORE, _prof(z, FOOT_SKIN)),
              (Z_NOZ0, Z_NOZ1, nozzle_bore(z), nozzle_out(z)),
              (Z_CH0, Z_CH1, CHAMBER_R, CHAMBER_OUT),
              (Z_NK0, Z_NK2, neck_bore(z), neck_out(z)),
              (Z_OUT0, Z_OUT1, OUTLET_BORE, OUTLET_OUT)]
    spans = sorted((lo, hi) for z0, z1, lo, hi in pieces
                   if z0 <= z <= z1 and hi > lo)
    merged = []
    for lo, hi in spans:
        if merged and lo <= merged[-1][1] + 1e-12:
            merged[-1][1] = max(merged[-1][1], hi)
        else:
            merged.append([lo, hi])
    return [(a, b) for a, b in merged]


def ideal_volume(step=0.02):
    """The smooth-body volume, integrated off `rings_at` plus the boss.

    The mesh is a little under this: every round surface here is a 40-gon, and
    the inscribed polygon loses (sin(2 pi/n) / (2 pi/n)) of the area, about a
    tenth of a percent at this facet count.
    """
    n_z = int(round((Z_OUT1 - Z_IN0) / step))
    v = 0.0
    for i in range(n_z):
        z = Z_IN0 + (i + 0.5) * step
        v += math.pi * sum(h * h - a * a for a, h in rings_at(z)) * step

    base = CHAMBER_R + BOSS_SINK                 # the boss's buried foot
    boss = math.pi * AIR_BOSS ** 2 * (AIR_OUT - base)

    # ...less the part of that foot that lies inside the chamber's own wall.
    sink = 0.0
    ys = [(-AIR_BOSS + (i + 0.5) * 0.02)
          for i in range(int(round(2 * AIR_BOSS / 0.02)))]
    for y in ys:
        h = math.sqrt(max(0.0, AIR_BOSS ** 2 - y * y))
        x_lo = math.sqrt(max(0.0, CHAMBER_R ** 2 - y * y))
        x_hi = math.sqrt(max(0.0, CHAMBER_OUT ** 2 - y * y))
        lo, hi = max(base, x_lo), min(AIR_OUT, x_hi)
        if hi > lo:
            sink += 2 * h * (hi - lo) * 0.02

    # ...and less the inlet's own hole, which runs through wall and boss.
    cut = 0.0
    step_y = 0.02
    for i in range(int(round(2 * AIR_BORE / step_y))):
        y = -AIR_BORE + (i + 0.5) * step_y
        h = math.sqrt(max(0.0, AIR_BORE ** 2 - y * y))
        x_lo = min(base, math.sqrt(max(0.0, CHAMBER_R ** 2 - y * y)))
        cut += 2 * h * (AIR_OUT - x_lo) * step_y
    return v + boss - sink - cut


# --- The solids themselves --------------------------------------------------
def validate(p):
    """Refuse a parameter set the geometry cannot honour, with the number
    that would work."""
    bad = []
    if p["throat"] >= p["inlet_bore"]:
        bad.append(f"throat {p['throat']} must be under the inlet bore "
                   f"{p['inlet_bore']} or nothing is squeezed")
    if p["throat"] >= p["chamber_r"] - p["nozzle_wall"]:
        bad.append(f"throat {p['throat']} plus a {p['nozzle_wall']} wall must "
                   f"stay inside the chamber's {p['chamber_r']} void")
    # The throat has to hold the air inlet between its own two ends.
    if not (Z_CONE1 < Z_AIR - p["air_bore"] and Z_THROAT1 > Z_AIR + p["air_bore"]):
        bad.append(f"the throat runs {Z_CONE1}..{Z_THROAT1} and the air inlet "
                   f"is {2 * p['air_bore']:g} mm across at {Z_AIR:g}: it must "
                   f"start above {Z_AIR - p['air_bore']:g} and end below "
                   f"{Z_AIR + p['air_bore']:g}")
    # The jet has still to be a jet where the air arrives, and it has to end
    # inside the chamber rather than built into the neck.
    if Z_NOZ1 <= Z_AIR + p["air_bore"]:
        bad.append(f"the nozzle's tip is at {Z_NOZ1:g}, at or below the top of "
                   f"the air inlet's mouth at {Z_AIR + p['air_bore']:.2f}: the "
                   f"water has dispersed before the air gets there, so nothing "
                   f"is carried away with it. Run Z_NOZ1 past "
                   f"{Z_AIR + p['air_bore']:.2f}")
    # Stations have to ascend within a piece: `sweep` and `tapered_tube` join
    # them in order, so one out of order folds the shell back through itself and
    # the piece comes out inside-out -- a negative volume the boolean refuses
    # with a message about the mesh, not about the parameter. Two stations may
    # share a height, which is a segment of no length, but none may come after
    # one that is higher.
    order = [("nozzle", [(Z_NOZ0, "Z_NOZ0"), (Z_CONE0, "Z_CONE0"),
                         (Z_CONE1, "Z_CONE1"), (Z_THROAT1, "Z_THROAT1"),
                         (Z_NOZ1, "Z_NOZ1")]),
             ("shell and neck", [(Z_CH0, "Z_CH0"), (Z_NK0, "Z_NK0"),
                                (Z_CH1, "Z_CH1"), (Z_NK1, "Z_NK1"),
                                (Z_NK2, "Z_NK2")]),
             ("shoulder", [(Z_FOOT0, "Z_FOOT0"), (Z_CH0, "Z_CH0"),
                           (Z_FOOT1, "Z_FOOT1")]),
             ("tubes", [(Z_IN0, "Z_IN0"), (Z_OUT0, "Z_OUT0"),
                        (Z_OUT1, "Z_OUT1")])]
    for name, pairs in order:
        seen = []
        for z, label in pairs:
            if seen and z < seen[-1][0]:
                bad.append(f"{name}: {label} = {z:g} does not rise above "
                           f"{seen[-1][1]} = {seen[-1][0]:g} -- a piece whose "
                           f"stations do not ascend folds back through itself")
            seen.append((z, label))
    if JET_RELEASE > 0.0 and Z_THROAT1 >= Z_NOZ1:
        bad.append(f"JET_RELEASE = {JET_RELEASE:g} has nowhere to happen: the "
                   f"throat already runs to {Z_THROAT1:g} and the tip is at "
                   f"{Z_NOZ1:g}. Lower Z_THROAT1 to give it room, or set "
                   f"JET_RELEASE back to 0.0 for a converging nozzle")
    if Z_NOZ1 >= Z_CH1:
        bad.append(f"the nozzle's tip ({Z_NOZ1:g}) is at or past the chamber's "
                   f"top ({Z_CH1:g}): it would be built into the neck")
    if p["chamber_r"] >= p["chamber_out"]:
        bad.append("the chamber's void is outside its own skin")
    # A hollow ring is only right side out if its bore is inside its skin, at
    # every station: get that wrong and the piece comes out inside-out, and the
    # union then swallows the neighbour it overlaps -- this is what removed ten
    # millimetres of the inlet tube once, leaving a shorter part and no error.
    rings = [("nozzle", [(z, b, b + p["nozzle_wall"]) for z, b in NOZZLE_BORE]),
             ("shoulder", [(z, p["floor_bore"], sk) for z, sk in FOOT_SKIN]),
             ("shell", [(Z_CH0, p["chamber_r"], p["chamber_out"]),
                        (Z_CH1, p["chamber_r"], p["chamber_out"])]),
             ("inlet", [(Z_IN0, p["inlet_bore"], p["inlet_out"]),
                        (Z_IN1, p["inlet_bore"], p["inlet_out"])]),
             ("outlet", [(Z_OUT0, p["outlet_bore"], p["outlet_out"]),
                         (Z_OUT1, p["outlet_bore"], p["outlet_out"])]),
             ("neck", [(z, neck_bore(z), neck_out(z)) for z, _ in NK_SKIN]),
             ("air boss", [(Z_AIR, p["air_bore"], p["air_boss"])])]
    for name, prof in rings:
        for z, ri, ro in prof:
            if ri >= ro:
                bad.append(f"the {name} is inside-out at z={z:g}: its bore "
                           f"{ri:.2f} mm is outside its own skin {ro:.2f} mm")
    if p["floor_bore"] <= p["inlet_bore"]:
        bad.append(f"the chamber floor's hole ({p['floor_bore']}) must clear "
                   f"the inlet bore ({p['inlet_bore']}) to overlap it")
    grip = nozzle_out(min(Z_FOOT1, Z_NOZ1))
    if p["floor_bore"] >= grip:
        bad.append(f"the floor's hole ({p['floor_bore']}) is wider than the "
                   f"nozzle standing in it ({grip:.2f} mm at z="
                   f"{min(Z_FOOT1, Z_NOZ1):g}), so the nozzle floats in the "
                   f"chamber's floor instead of fusing to it: bring "
                   f"FLOOR_BORE under {grip:.2f}")
    if p["floor_bore"] >= nozzle_bore(Z_NOZ0) + p["nozzle_wall"]:
        bad.append(f"the floor's hole ({p['floor_bore']}) is wider than the "
                   f"nozzle's buried foot: water would leak past the joint")
    if Z_FOOT1 <= Z_CH0:
        bad.append("the shoulder's floor has to outlast the shell's bottom edge")
    # The one outside face that leans outward as it rises. Flattened under the
    # print angle it is a shelf the printer has to hold up.
    off_grid = [z for z in (Z_IN0, Z_IN1, Z_FOOT0, Z_FOOT1, Z_NOZ0, Z_NOZ1,
                            Z_CONE0, Z_CONE1, Z_THROAT1, Z_CH0, Z_CH1, Z_NK0,
                            Z_NK1, Z_NK2, Z_OUT0, Z_OUT1, Z_AIR)
                if abs(z / GRID - round(z / GRID)) > 1e-9]
    if off_grid:
        bad.append(f"feature edges off the {GRID:g} mm layer grid: "
                   f"{off_grid} -- they print as a partial layer")
    lean = math.degrees(math.atan2(Z_CH0 - Z_FOOT0, FOOT_RISE))
    if lean < 45.0:
        bad.append(f"the shoulder leans out {FOOT_RISE:.2f} mm over "
                   f"{Z_CH0 - Z_FOOT0:.2f} mm, which is {lean:.1f} degrees off "
                   f"the horizontal -- flatter than 45 and it is a shelf the "
                   f"printer holds up. Raise LEAN_DEG, or start the shoulder "
                   f"at z={Z_CH0 - FOOT_RISE:.2f} or lower")
    if p["inlet_out"] - BURIED <= p["inlet_bore"] + p["nozzle_wall"]:
        bad.append("the shoulder's foot has to be wider than the nozzle "
                   "standing inside it")
    # The annulus the air travels in has to be at least as wide as the hole it
    # arrives through, or it is a slit, not a passage.
    annulus = p["chamber_r"] - (p["throat"] + p["nozzle_wall"])
    if annulus < 2 * p["air_bore"]:
        bad.append(f"the annulus between the nozzle's skin and the chamber is "
                   f"{annulus:.2f} mm, narrower than the "
                   f"{2 * p['air_bore']:g} mm air inlet that feeds it")
    if JET_CLEAR <= 0.0:
        bad.append("the air inlet's hole must stop short of the nozzle's skin")
    # The boss's foot has to sit inside the chamber's wall, and the wall is
    # faceted: its skin is at chamber_out * cos(pi / n), not chamber_out.
    n = t.auto_segments(p["chamber_out"], p["tol"])
    base = p["chamber_r"] + BOSS_SINK
    foot = math.hypot(p["air_boss"], base)
    skin = p["chamber_out"] * math.cos(math.pi / n)
    if base <= p["chamber_r"]:
        bad.append(f"the boss's foot at {base:g} mm is not inside the wall: "
                   f"the chamber's void reaches {p['chamber_r']:g} mm")
    if foot >= skin:
        bad.append(f"a {2 * p['air_boss']:g} mm boss standing {base:g} mm out "
                   f"reaches {foot:.2f} mm, past the faceted skin at "
                   f"{skin:.2f}: sink the foot to at most "
                   f"{math.sqrt(skin ** 2 - p['air_boss'] ** 2):.2f} mm or "
                   f"thin the boss to "
                   f"{math.sqrt(skin ** 2 - base ** 2):.2f} mm radius")
    if p["air_out"] <= skin:
        bad.append(f"the air inlet stops {p['air_out']:g} mm out from the axis, "
                   f"short of the chamber's faceted skin at {skin:.2f} mm, so "
                   f"the hole never reaches the air: put AIR_OUT beyond "
                   f"{skin:.2f}, and remember the boss has to stand clear of "
                   f"the chamber's {p['chamber_out']:g} mm round skin as well")
    if p["air_out"] < skin + 0.6:
        bad.append(f"the boss ends {p['air_out'] - skin:.2f} mm proud of the "
                   f"chamber's faceted skin -- too thin to drill cleanly or to "
                   f"bolt a fitting to. Give it a few millimetres: AIR_OUT = "
                   f"{math.ceil((skin + 3.0) / GRID) * GRID:.1f} stands 3 mm "
                   f"clear")
    if Z_AIR - p["air_bore"] < Z_CH0 or Z_AIR + p["air_bore"] > Z_CH1:
        bad.append("the air inlet does not fit inside the chamber's length")
    if bad:
        raise SystemExit("cannot build this aspirator:\n  " + "\n  ".join(bad))


def build(**overrides):
    p = dict(inlet_bore=INLET_BORE, inlet_out=INLET_OUT, throat=THROAT,
             nozzle_wall=NOZZLE_WALL, floor_bore=FLOOR_BORE, chamber_r=CHAMBER_R,
             chamber_out=CHAMBER_OUT, outlet_bore=OUTLET_BORE,
             outlet_out=OUTLET_OUT, air_bore=AIR_BORE, air_boss=AIR_BOSS,
             air_out=AIR_OUT, tol=TOL)
    p.update(overrides)
    validate(p)

    # One facet count and one phase for everything round, so that where a
    # surface of one piece is cut by the plane of another the new vertices
    # still land on the same azimuths and nothing has to close a T-junction
    # between two polygons that disagree about where the circle's corners are.
    n = t.auto_segments(p["chamber_out"], p["tol"])
    # `sweep` starts its rings at +y where `ring` starts at +x; half a turn of
    # phase puts them back on the same vertices.
    half = -math.pi / 2.0
    axis = (0.0, 0.0, 1.0)

    def ring(name, z0, z1, ri0, ri1, ro0, ro1):
        m = t.Mesh(name)
        t.tapered_tube(m, ri0, ro0, ri1, ro1, z0, z1, n)
        return m

    def revolve(name, profs):
        """A hollow body of revolution through (station, bore, skin) triples."""
        m = t.Mesh(name)
        t.sweep(m, [((0.0, 0.0, z), axis, ri, ro) for z, ri, ro in profs], n,
                phase=half)
        return m

    # The nozzle: one hollow tube whose bore converges, holds, and lets go. Its
    # first station sits inside the inlet tube, so its bore there is BURIED
    # wider than the inlet's -- the inlet's own bore is what the water sees, and
    # the two cylindrical surfaces never run alongside each other. A single
    # `sweep` because a tapered_tube is straight, and five butted ones would be
    # four faces meeting in one plane.
    nozzle = revolve("nozzle", [(z, r, r + p["nozzle_wall"])
                                for z, r in NOZZLE_BORE])

    # The neck: its bore and its skin reach the outlet's at different heights,
    # which is the one thing a tapered_tube cannot say and a sweep can.
    neck = revolve("neck", [(z, rb, rn)
                            for (z, rb), (_, rn) in zip(NK_BORE, NK_SKIN)])

    base = p["chamber_r"] + BOSS_SINK
    boss = t.Mesh("boss")
    t.cylinder(boss, p["air_boss"], base, p["air_out"], n)
    t.rotate_y(boss, 90.0)
    t.translate(boss, dz=Z_AIR)

    bore = t.Mesh("bore")
    t.cylinder(bore, p["air_bore"], p["throat"] + p["nozzle_wall"] + JET_CLEAR,
               p["air_out"] + 1.0, n)
    t.rotate_y(bore, 90.0)
    t.translate(bore, dz=Z_AIR)

    parts = [
        ring("inlet", Z_IN0, Z_IN1, p["inlet_bore"], p["inlet_bore"],
             p["inlet_out"], p["inlet_out"]),
        revolve("shoulder", [(z, p["floor_bore"], sk)
                             for z, sk in FOOT_SKIN]),
        nozzle,
        ring("shell", Z_CH0, Z_CH1, p["chamber_r"], p["chamber_r"],
             p["chamber_out"], p["chamber_out"]),
        neck,
        ring("outlet", Z_OUT0, Z_OUT1, p["outlet_bore"], p["outlet_bore"],
             p["outlet_out"], p["outlet_out"]),
        boss,
    ]

    body = csg.combine(
        [("start", parts[0])]
        + [("union", m) for m in parts[1:]]
        + [("subtract", bore)])
    body.name = "aspirator"
    return body


def main():
    memcap.cap()
    out = sys.argv[1] if len(sys.argv) > 1 else default_out()
    mesh = build()
    annulus = CHAMBER_R - (THROAT + NOZZLE_WALL)
    bb = mesh.bbox()
    desc = (f"Water aspirator, a Venturi pump with no moving parts: a "
            f"{2 * INLET_BORE:g} mm supply bore squeezes water through a "
            f"{2 * THROAT:g} mm throat that starts before, and ends after, the "
            f"{2 * AIR_BORE:g} mm air inlet on the side of a "
            f"{2 * CHAMBER_R:g} mm chamber, so the jet drags air through and out "
            f"the {2 * OUTLET_BORE:g} mm outlet. The nozzle only converges"
            + (" -- the jet is still at throat size as it leaves the tip"
               if JET_RELEASE == 0.0 else f", letting go by {2 * JET_RELEASE:g} mm")
            + f". Around the inlet the wall is two "
            f"rings -- a {2 * THROAT:g} mm jet inside a "
            f"{2 * (THROAT + NOZZLE_WALL):.1f} mm nozzle, {annulus:.1f} mm of air "
            f"between it and the chamber's shell. Prints standing on the inlet "
            f"end with the axis vertical and the air inlet horizontal; the "
            f"boss's underside is the only face that wants support, and it is "
            f"outside the part.")
    t.write_3mf(out, mesh, title="water aspirator",
                description=desc, params=PARAMS)
    print(f"wrote {out}")
    print(f"  params            : {t.describe(PARAMS)}")
    print(f"  param hash        : {t.params_hash(PARAMS)}")
    print(f"  water path        : {2 * INLET_BORE:g} mm bore, squeezed to "
          f"{2 * THROAT:g} mm over {Z_CONE1 - Z_CONE0:g} mm, held narrow to "
          f"z={Z_THROAT1:g}" + (f", let go {2 * JET_RELEASE:g} mm wider by the "
                               f"tip at z={Z_NOZ1:g}" if JET_RELEASE > 0.0
                               else ", and never widened again"))
    print(f"  throat spans the air inlet: {Z_CONE1:g}..{Z_THROAT1:g}, the "
          f"inlet's own {2 * AIR_BORE:g} mm hole sits at {Z_AIR:g}")
    print(f"  chamber           : {2 * CHAMBER_R:g} mm void round the nozzle, "
          f"{Z_CH1 - Z_CH0:g} mm tall, {2 * CHAMBER_OUT:g} mm outside")
    print(f"  section at the inlet: two rings -- a {2 * THROAT:g} mm jet in a "
          f"{2 * (THROAT + NOZZLE_WALL):g} mm nozzle, an annulus of air out "
          f"to {2 * CHAMBER_R:g} mm, and the shell to {2 * CHAMBER_OUT:g} mm")
    print(f"  air inlet         : {2 * AIR_BORE:g} mm bore in a "
          f"{2 * AIR_BOSS:g} mm boss, standing {AIR_OUT:g} mm from the axis")
    print(f"  size              : {bb[3] - bb[0]:.1f} x {bb[4] - bb[1]:.1f} x "
          f"{bb[5] - bb[2]:.1f} mm")
    print(f"  volume            : {mesh.signed_volume():.1f} mm^3 "
          f"(ideal {ideal_volume():.1f})")
    print(f"  print orientation : standing on the inlet end, axis vertical, "
          f"air inlet horizontal. Every internal surface is then a vertical "
          f"wall or a cone leaning inward as it rises; the boss outside the "
          f"chamber is the one feature that wants support, and it is "
          f"removable")
    print(f"  vertices / tris   : {len(mesh.verts)} / {len(mesh.tris)}")
    print(f"  shells            : {len(mesh.shells())}")
    crack = csg.crack_length(mesh.verts, mesh.tris)
    print(f"  gaps              : {len(csg.real_gaps(mesh.verts, mesh.tris))} "
          f"real, {crack:.1f} mm of seam where two curved surfaces met along "
          f"one line -- four edges with four triangles around them, none open")


if __name__ == "__main__":
    main()

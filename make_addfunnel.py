#!/usr/bin/env python3
"""Pressure-equalising liquid addition funnel.

The storage cone is printed upright.  Liquid leaves through a coaxial
delivery tube; the larger pressure-equalising tube surrounds it and opens
into the headspace above the liquid.  The stopcock key and storage stopper
are separate 3MF objects so they can be printed in their own stable poses.

Usage: python3 make_addfunnel.py [output.3mf]
"""

import math
import os
import sys

import csglib as csg
import memcap
import threelib as t

P = {
    "funnel_z0": 19.5, "funnel_z1": 88.0,
    "funnel_rim_r": 24.0, "funnel_rim_out": 26.0,
    "funnel_bore0": 3.8, "funnel_out0": 9.2,
    "barrel_r": 8.0, "barrel_half": 8.0, "barrel_zc": 13.0,
    "bore_narrow": 4.4, "bore_wide": 5.2,
    "key_clear": 0.08, "key_tip_x": -7.6, "key_end_x": 7.6,
    "port_half": 2.0, "sleeve_out": 4.0, "top_sleeve_out": 3.7,
    "delivery_top": 7.4, "delivery_tip": -80.0,
    "eq_r_in": 5.4, "eq_r_out": 6.7,
    "port_z": 74.0, "bend_r": 12.0, "transition_r": 12.0,
    "transition_angle": 35.0, "transition_top": -18.0,
    "eq_tip": -80.0,
    "stopper_z0": 48.0, "stopper_z1": 88.0,
    "stopper_clear": 0.15, "flange_r": 27.0, "flange_t": 3.0,
    "knob_af": 5.0, "knob_len": 9.4,
    "plate_gap": 60.0,
    "facet_tol": 0.01,
}


def default_out():
    return t.stamped_filename("addfunnel", P)


def funnel_r_in(p, z):
    return p["funnel_bore0"] + (p["funnel_rim_r"] - p["funnel_bore0"]) * (
        (z - p["funnel_z0"]) / (p["funnel_z1"] - p["funnel_z0"]))


def funnel_r_out(p, z):
    return p["funnel_out0"] + (p["funnel_rim_out"] - p["funnel_out0"]) * (
        (z - p["funnel_z0"]) / (p["funnel_z1"] - p["funnel_z0"]))


def bore_r(p, x):
    f = p["barrel_half"]
    return p["bore_narrow"] + (p["bore_wide"] - p["bore_narrow"]) * (
        (x + f) / (2 * f))


def key_r(p, x):
    return bore_r(p, x) - p["key_clear"]


def eq_port_x(p):
    return funnel_r_in(p, p["port_z"]) - 2.0


def eq_centerline(p):
    """Return independent station data for the equaliser centreline."""
    x0 = eq_port_x(p)
    rb = p["bend_r"]
    rs = p["transition_r"]
    ang = math.radians(p["transition_angle"])
    zt = p["transition_top"]
    stations = []

    def add(pos, tangent):
        stations.append((tuple(pos), tuple(tangent)))

    # Horizontal opening turns through a 90 degree radius into the descent.
    for i in range(23):
        a = math.pi * i / 44
        add((x0 + rb * math.sin(a), 0.0, p["port_z"] - rb + rb * math.cos(a)),
            (math.cos(a), 0.0, -math.sin(a)))
    r_desc = x0 + rb
    add((r_desc, 0.0, zt), (0.0, 0.0, -1.0))

    # S transition uses two circular arcs separated by a straight section.
    c1 = (r_desc - rs, 0.0, zt)
    for i in range(1, 8):
        a = ang * i / 8
        add((c1[0] + rs * math.cos(a), 0.0, zt - rs * math.sin(a)),
            (-math.sin(a), 0.0, -math.cos(a)))
    p1 = (r_desc - rs * (1 - math.cos(ang)), 0.0, zt - rs * math.sin(ang))
    direction = (-math.sin(ang), 0.0, -math.cos(ang))
    length = (r_desc - 2 * rs * (1 - math.cos(ang))) / math.sin(ang)
    p2 = tuple(p1[i] + direction[i] * length for i in range(3))
    add(p1, direction)
    add(p2, direction)
    radial = (-math.cos(ang), 0.0, math.sin(ang))
    c2 = (p2[0] + rs * math.cos(ang), 0.0, p2[2] - rs * math.sin(ang))
    for i in range(1, 8):
        a = ang * i / 8
        pos = tuple(c2[j] + rs * (radial[j] * math.cos(a) + direction[j] * math.sin(a))
                    for j in range(3))
        tangent = tuple(-math.sin(a) * radial[j] + math.cos(a) * direction[j]
                        for j in range(3))
        add(pos, tangent)
    add((0.0, 0.0, p2[2] - rs * math.sin(ang)), (0.0, 0.0, -1.0))
    add((0.0, 0.0, p["eq_tip"]), (0.0, 0.0, -1.0))
    return stations, p2[2] - rs * math.sin(ang)


def body_mesh(p):
    return t.tapered_tube(
        t.Mesh(), p["funnel_bore0"], p["funnel_out0"],
        p["funnel_rim_r"], p["funnel_rim_out"],
        p["funnel_z0"], p["funnel_z1"], t.auto_segments(26, p["facet_tol"]))


def stopcock(p):
    mesh = t.Mesh()
    wall = t.tapered_tube(
        mesh, p["bore_narrow"], p["barrel_r"], p["bore_wide"], p["barrel_r"],
        -p["barrel_half"], p["barrel_half"], t.auto_segments(8, p["facet_tol"]))
    wall = t.translate(t.rotate_y(wall, 90), 0, 0, p["barrel_zc"])
    swath = t.box(t.Mesh(), -p["port_half"], -p["barrel_r"] - 0.5,
                  p["barrel_zc"] - p["barrel_r"] - 0.5,
                  p["port_half"], p["barrel_r"] + 0.5,
                  p["barrel_zc"] + p["barrel_r"] + 0.5)
    top = t.tube(t.Mesh(), p["port_half"], p["top_sleeve_out"], 19.4, 26.0,
                 t.auto_segments(4, p["facet_tol"]))
    delivery = t.tube(t.Mesh(), p["port_half"], p["sleeve_out"],
                      p["delivery_tip"], p["delivery_top"],
                      t.auto_segments(4, p["facet_tol"]))
    return csg.combine([
        ("start", wall), ("subtract", swath), ("union", top), ("union", delivery)])


def equalizer(p, mesh):
    stations, eq_start = eq_centerline(p)
    rings = [((x, y, z), (tx, ty, tz), p["eq_r_in"], p["eq_r_out"])
             for (x, y, z), (tx, ty, tz) in stations]
    return t.sweep(mesh, rings, t.auto_segments(6.7, p["facet_tol"]),
                   ref=(0, 1, 0)), eq_start


def assembly(p):
    mesh = body_mesh(p)
    eq_mesh, _ = equalizer(p, t.Mesh())
    return csg.combine([
        ("start", mesh), ("union", stopcock(p)), ("union", eq_mesh)])


def key(p):
    mesh = t.Mesh()
    plug = t.revolve(mesh, [
        (0.0, p["key_tip_x"]), (key_r(p, p["key_tip_x"]), p["key_tip_x"]),
        (key_r(p, p["key_end_x"]), p["key_end_x"]), (0.0, p["key_end_x"]),
    ], t.auto_segments(8, p["facet_tol"]))
    half = p["port_half"]
    slot = t.box(t.Mesh(), -(key_r(p, p["key_end_x"]) + 1), -half, -half,
                 key_r(p, p["key_end_x"]) + 1, half, half)
    knob_r = p["knob_af"] / (2 * math.cos(math.pi / 6))
    knob_poly = [(knob_r * math.cos(i * math.pi / 3 + math.pi / 6),
                  knob_r * math.sin(i * math.pi / 3 + math.pi / 6))
                 for i in range(6)]
    knob = t.prism(t.Mesh(), knob_poly, p["key_end_x"],
                   p["key_end_x"] + p["knob_len"])
    key = csg.combine([("start", plug), ("subtract", slot), ("union", knob)])
    # The build frame already stands the key on its tip: axis vertical, slot
    # horizontal (it turns vertical about y when the key is inserted).
    return t.translate(key, p["plate_gap"], 0, -p["key_tip_x"])


def stopper(p):
    mesh = t.Mesh()
    z0, z1 = p["stopper_z0"], p["stopper_z1"]
    profile = [(0.0, z0), (funnel_r_in(p, z0) - p["stopper_clear"], z0),
               (funnel_r_in(p, z1) - p["stopper_clear"], z1), (0.0, z1)]
    plug = t.revolve(mesh, profile, t.auto_segments(26, p["facet_tol"]))
    flange = t.cylinder(t.Mesh(), p["flange_r"], z1, z1 + p["flange_t"],
                        t.auto_segments(26, p["facet_tol"]))
    stop = csg.combine([("start", plug), ("union", flange)])
    stop = t.translate(t.rotate_x(stop, 180), 0, 0, z1 + p["flange_t"])
    return t.translate(stop, -p["plate_gap"], 0, 0)


def build(**p):
    params = dict(P)
    params.update(p)
    return [
        ("assembly", assembly(params)),
        ("key", key(params)),
        ("stopper", stopper(params)),
    ]


def main():
    memcap.cap()
    out = sys.argv[1] if len(sys.argv) > 1 else default_out()
    t.write_3mf(out, build(), params=P)
    print(f"wrote {out} ({os.path.getsize(out)} bytes)")


if __name__ == "__main__":
    main()

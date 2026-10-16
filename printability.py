"""printability -- sanity-check a mesh against printer process limits.

These are heuristics for catching gross problems *before* slicing, not a
substitute for slicing. Every function returns plain data so it can be used
in a test suite as well as printed.

The governing idea: a mesh is not a solid until a nozzle with a given
diameter has had to build it, one layer at a time. Features smaller than the
nozzle cannot be represented, and features not aligned to the layer grid get
staircased.

Limitations worth knowing:
* `lateral_thickness` measures each shell's XY bounding box, which is only
  meaningful for box-like solids. A tube's bbox spans its whole diameter and
  hides the wall thickness entirely.
* `overhangs` looks at triangles independently, so a curved underside reports
  as many small steep faces even where a shell-based slicer would bridge them.
"""

from __future__ import annotations

import math
import os

__all__ = [
    "NOZZLE_DEFAULT", "LAYER_DEFAULT", "OVERHANG_DEG_DEFAULT", "OVERHANG_DEG",
    "triangle_normals", "overhangs", "lateral_thickness", "layer_alignment",
    "bed_contact", "assess",
]

NOZZLE_DEFAULT = 0.4      # mm, common FDM extrusion width
LAYER_DEFAULT = 0.2       # mm
OVERHANG_DEG_DEFAULT = 45.0

# The shallowest face that still prints unsupported, in degrees from the
# horizontal. This is a property of the printer and the plastic, not of the
# part: a well-tuned machine copes with far shallower overhangs than the usual
# 45 degree rule, and one that is running cold copes with less. Set it in the
# environment (`OVERHANG_ANGLE=20`) or pass `max_angle` per call; every
# verifier here takes `--angle` and reports which value it used, because a
# check that quietly assumed 45 would pass a part that will not print.
OVERHANG_DEG = float(os.environ.get("OVERHANG_ANGLE", OVERHANG_DEG_DEFAULT))
# n.z at or below this counts as horizontal. -0.98 is ~11 degrees off flat.
FLAT_CUTOFF_DEFAULT = -0.98


def triangle_normals(verts, tris):
    """(unit normal, area, lowest z, highest z) for each triangle."""
    out = []
    for a, b, c in tris:
        ax, ay, az = verts[a]
        bx, by, bz = verts[b]
        cx, cy, cz = verts[c]
        ux, uy, uz = bx - ax, by - ay, bz - az
        vx, vy, vz = cx - ax, cy - ay, cz - az
        nx = uy * vz - uz * vy
        ny = uz * vx - ux * vz
        nz = ux * vy - uy * vx
        ln = math.sqrt(nx * nx + ny * ny + nz * nz)
        if ln == 0.0:
            out.append(((0.0, 0.0, 0.0), 0.0, min(az, bz, cz),
                        max(az, bz, cz)))
            continue
        out.append(((nx / ln, ny / ln, nz / ln), 0.5 * ln,
                    min(az, bz, cz), max(az, bz, cz)))
    return out


def overhangs(verts, tris, layer_height=LAYER_DEFAULT,
              max_angle=None, ignore_below=None,
              flat_cutoff=FLAT_CUTOFF_DEFAULT, exposed_only=False,
              shells=None, probe=0.15):
    """Downward faces shallower than `max_angle` from horizontal (need support).

    A vertical wall has a horizontal normal and needs nothing. A flat ceiling
    has n.z == -1 and needs full support, so it is always reported. The limit
    is `max_angle`, in degrees from the horizontal, and it defaults to
    `OVERHANG_DEG` -- 45 unless the environment says otherwise, because what a
    printer can bridge is the printer's business.

    Faces lying entirely at or below the first layer are skipped -- those sit
    on the bed and are supported by definition. The test is on the face's
    HIGHEST point: a tall wall that merely touches the bed at its base is not
    supported by it, and testing the lowest point hides every such overhang.

    Each face is tagged "sloped" or "flat". The angle rule only means anything
    for a sloped face; a horizontal down-face carries no angle information at
    all. A lattice bar bridging a 0.5 mm gap and a 70 mm ceiling share the
    same -Z normal, so callers must treat the two kinds differently.

    `exposed_only` drops the faces that are buried in the material. Parts here
    are built by overlapping shells, so plenty of their faces are internal --
    a cap where one tube stops inside another -- and the slicer never sees
    them: the union of the shells is what prints. A face counts as exposed
    when stepping off it along its own normal lands outside every shell, which
    is true of both the outside skin and the walls of a passage, and those are
    the ones a printer has to carry.
    """
    max_angle = OVERHANG_DEG if max_angle is None else max_angle
    cutoff = -math.cos(math.radians(max_angle))
    ignore_below = layer_height if ignore_below is None else ignore_below
    boxes = None
    if exposed_only:
        import checks
        if shells is None:
            shells = checks.shell_triangles(verts, tris)
        boxes = checks.shell_boxes(verts, tris, shells)

    out = []
    for (n, area, zmin, zmax), tri in zip(triangle_normals(verts, tris), tris):
        if area <= 0 or n[2] > cutoff or zmax <= ignore_below + 1e-9:
            continue
        cx = cy = cz = 0.0
        for v in tri:
            cx += verts[v][0] / 3.0
            cy += verts[v][1] / 3.0
            cz += verts[v][2] / 3.0
        if exposed_only:
            step = (cx + probe * n[0], cy + probe * n[1], cz + probe * n[2])
            if checks.inside(verts, tris, step, shells=shells, boxes=boxes):
                continue
        out.append({"c": (cx, cy, cz), "n": n, "area": area, "z": zmin,
                    "kind": "flat" if n[2] <= flat_cutoff else "sloped",
                    "angle": math.degrees(math.acos(max(-1.0, min(1.0, n[2]))))})
    out.sort(key=lambda d: d["z"])
    return out


def lateral_thickness(verts, shells, nozzle=NOZZLE_DEFAULT, marginal=1.5):
    """Classify shells by lateral thickness against the nozzle.

    Two very different regimes, easy to conflate:

    * **narrower than the nozzle** -- cannot be extruded at all; the nozzle
      is physically wider than the feature it is asked to fill.
    * **up to ~1.5 nozzle diameters** -- prints as a single bead. This is
      normal and usually exactly what you want for a screen, grid, lattice
      or light fit. The real caveats are that it will not carry load, and its
      actual width follows your extrusion calibration, not the nominal value.

    Only meaningful for box-like shells -- see the module docstring.
    """
    out = []
    for sh in shells:
        x0, y0, _, x1, y1, _ = sh["bbox"]
        w = min(x1 - x0, y1 - y0)
        if w >= nozzle * marginal:
            continue
        ratio = w / nozzle
        out.append({"width": w, "ratio": ratio, "bbox": sh["bbox"],
                    "n_verts": sh["n"],
                    "tier": "impossible" if ratio < 1.0 else "single-bead"})
    out.sort(key=lambda d: d["width"])
    return out


def layer_alignment(verts, layer_height=LAYER_DEFAULT, tol=0.02):
    """Distinct z values that do not sit on the layer grid.

    Anything not a multiple of the layer height gets sliced across, which
    rounds it and makes the real dimension depend on the layer height.
    """
    bad = []
    for z in sorted({round(v[2], 9) for v in verts}):
        k = z / layer_height
        off = abs(k - round(k)) * layer_height
        if off > tol and z > tol:
            bad.append({"z": z, "off_by": off,
                        "nearest_layers": (math.floor(k), math.ceil(k))})
    return bad


def bed_contact(verts, tris, layer_height=LAYER_DEFAULT, tol=1e-6):
    """Total downward-facing area sitting on the first layer."""
    total = 0.0
    for (n, area, zmin, _zmax) in triangle_normals(verts, tris):
        if area > 0 and n[2] < -0.999 and zmin <= layer_height + tol:
            total += area
    return total


def assess(verts, tris, shells=None, nozzle=NOZZLE_DEFAULT,
           layer_height=LAYER_DEFAULT, overhang_deg=None,
           marginal=1.5):
    """Run every check; returns (findings, summary). Findings are dicts.

    `warn` means the model will not print as drawn. `info` means it prints
    but carries a caveat worth knowing -- do not treat every info as a defect.
    """
    if shells is None:
        import checks
        shells = checks.shells(verts, tris)

    overhang_deg = OVERHANG_DEG if overhang_deg is None else overhang_deg
    findings = []
    ov = overhangs(verts, tris, layer_height, overhang_deg)
    sloped = [o for o in ov if o["kind"] == "sloped"]
    flat = [o for o in ov if o["kind"] == "flat"]
    if sloped:
        findings.append({
            "level": "warn", "check": "overhang",
            "msg": f"{len(sloped)} faces steeper than {overhang_deg:g} deg need "
                   f"support ({sum(o['area'] for o in sloped):.1f} mm^2, "
                   f"lowest at z={sloped[0]['z']:.2f})",
        })
    if flat:
        findings.append({
            "level": "info", "check": "ceiling",
            "msg": f"{len(flat)} horizontal down-facing faces "
                   f"({sum(o['area'] for o in flat):.1f} mm^2) -- on a lattice "
                   f"these bridge short gaps, but a wide flat ceiling needs "
                   f"support; the angle rule cannot tell the two apart",
        })

    thin = lateral_thickness(verts, shells, nozzle, marginal)
    impossible = [x for x in thin if x["tier"] == "impossible"]
    single = [x for x in thin if x["tier"] == "single-bead"]
    if impossible:
        findings.append({
            "level": "warn", "check": "thin",
            "msg": f"{len(impossible)} shells are narrower than the "
                   f"{nozzle:g} mm nozzle (thinnest {impossible[0]['width']:.3f} mm) "
                   f"-- cannot be extruded at all",
        })
    if single:
        findings.append({
            "level": "info", "check": "thin",
            "msg": f"{len(single)} shells print as a single {nozzle:g} mm bead "
                   f"(thinnest {single[0]['width']:.3f} mm = "
                   f"{single[0]['ratio']:.2f}x nozzle); fine for a screen or grid, "
                   f"but not load-bearing and not dimensionally exact",
        })

    misaligned = layer_alignment(verts, layer_height)
    if misaligned:
        findings.append({
            "level": "warn", "check": "layers",
            "msg": f"{len(misaligned)} z values off the {layer_height:g} mm layer "
                   f"grid, e.g. z={misaligned[0]['z']:.3f} "
                   f"(off by {misaligned[0]['off_by']:.3f} mm)",
        })

    contact = bed_contact(verts, tris, layer_height)
    if contact < 20.0:
        findings.append({
            "level": "info", "check": "bed",
            "msg": f"only {contact:.1f} mm^2 on the build plate -- "
                   f"weak first-layer adhesion",
        })
    return findings, {"overhang_area": sum(o["area"] for o in ov),
                      "n_overhang_faces": len(ov),
                      "n_thin_impossible": len(impossible),
                      "n_thin_single_bead": len(single),
                      "thinnest": thin[0]["width"] if thin else None,
                      "thinnest_ratio": thin[0]["ratio"] if thin else None,
                      "n_misaligned": len(misaligned),
                      "bed_contact": contact}

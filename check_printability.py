#!/usr/bin/env python3
"""Run printability checks against a generated 3MF."""

import sys

import checks as K
import printability as P

NOZZLE = float(sys.argv[2]) if len(sys.argv) > 2 else P.NOZZLE_DEFAULT
LAYER = float(sys.argv[3]) if len(sys.argv) > 3 else P.LAYER_DEFAULT

path = sys.argv[1]
objs, _, _ = K.load_3mf(path)
verts, tris = objs[0]["verts"], objs[0]["tris"]
shells = K.shells(verts, tris)

findings, summary = P.assess(verts, tris, shells, nozzle=NOZZLE,
                             layer_height=LAYER)

print(f"== {path} == nozzle {NOZZLE:g} mm, layer {LAYER:g} mm ==")
print(f"  shells {len(shells)}  verts {len(verts)}  tris {len(tris)}")
if findings:
    warns = [f for f in findings if f["level"] == "warn"]
    infos = [f for f in findings if f["level"] == "info"]
    for f in warns + infos:
        print(f"  [{f['level'].upper():5}] {f['check']}: {f['msg']}")
else:
    warns = []
    print("  no findings")
th = summary["thinnest"]
print(f"  overhang {summary['n_overhang_faces']} faces / "
      f"{summary['overhang_area']:.1f} mm^2"
      + (f" | thinnest {th:.3f} mm ({summary['thinnest_ratio']:.2f}x nozzle)"
         if th is not None else "")
      + f" | bed contact {summary['bed_contact']:.0f} mm^2")
sys.exit(1 if warns else 0)

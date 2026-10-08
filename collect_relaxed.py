#!/usr/bin/env python3
"""Copy the DZVP-relaxed ribbons (runs/<name>/relaxed.xyz) into structures_relaxed/ for build_wet.py.

    python3 collect_relaxed.py    # then:  python3 build_wet.py && python3 make_cp2k.py
"""
import shutil
from pathlib import Path

HERE = Path(__file__).resolve().parent
(HERE / "structures_relaxed").mkdir(exist_ok=True)
for name in ("rib_y_o00_si", "rib_y_o00_al", "rib_x_o00_si"):
    src = HERE / "runs" / name / "relaxed.xyz"
    if src.exists():
        shutil.copy(src, HERE / "structures_relaxed" / f"{name}.xyz"); print("collected", name)
    else:
        print("MISSING", name, "(run its optimisation first)")

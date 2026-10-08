#!/usr/bin/env python3
"""Write CP2K inputs for every structure in ../structures -> cp2k/runs/<name>/{start.xyz,opt.inp,sp.inp}.

    python3 make_cp2k.py

opt: GEO_OPT (LBFGS, PBE, DZVP-MOLOPT-SR-GTH, 400 Ry), ions only, cell fixed.
sp : ENERGY on opt's last geometry (relaxed.xyz), DZVP-MOLOPT-SR-GTH, 400 Ry,
     Hirshfeld + Mulliken charges.  Cells are 3D periodic with vacuum.
"""
import json
import shutil
from pathlib import Path

HERE = Path(__file__).resolve().parent
KIND = {"H": "q1", "O": "q6", "Al": "q3", "Si": "q4", "Mg": "q10", "Na": "q9"}


def kinds(basis, elements):
    return "\n".join(f"    &KIND {e}\n      BASIS_SET {basis}\n      POTENTIAL GTH-PBE-{q}\n    &END KIND"
                     for e, q in KIND.items() if e in elements)


def inp(name, cell, mode, robust=False, elements=tuple(KIND)):
    opt = mode == "opt"
    basis = "DZVP-MOLOPT-SR-GTH"          # same level for geometry and charges (SZV gave Si-O 1.73 A, O-H 1.03 A: far too long)
    cutoff, rel = 400, 50
    coord = "start.xyz" if opt else "relaxed.xyz"
    printblk = "" if opt else """
    &PRINT
      &HIRSHFELD
        SELF_CONSISTENT .FALSE.
        SHAPE_FUNCTION DENSITY
        REFERENCE_CHARGE ATOMIC
      &END HIRSHFELD
      &MULLIKEN ON
      &END MULLIKEN
    &END PRINT"""
    motion = """
&MOTION
  &GEO_OPT
    TYPE MINIMIZATION
    OPTIMIZER LBFGS
    MAX_ITER 100
    MAX_FORCE 1.0E-3
  &END GEO_OPT
&END MOTION
""" if opt else ""
    eps = "1.0E-5" if opt else "5.0E-6"
    if robust:   # fallback when OT/DIIS stalls (edge states, poor starting geometry): diagonalisation + Broyden + smearing
        scf_block = f"""    &SCF
      SCF_GUESS ATOMIC
      EPS_SCF {eps}
      MAX_SCF 400
      ADDED_MOS 40
      &DIAGONALIZATION
        ALGORITHM STANDARD
      &END DIAGONALIZATION
      &MIXING
        METHOD BROYDEN_MIXING
        ALPHA 0.15
        BETA 1.5
        NBROYDEN 8
      &END MIXING
      &SMEAR
        METHOD FERMI_DIRAC
        ELECTRONIC_TEMPERATURE [K] 300.0
      &END SMEAR
    &END SCF"""
    else:
        scf_block = f"""    &SCF
      SCF_GUESS ATOMIC
      EPS_SCF {eps}
      MAX_SCF 300
      &OT
        PRECONDITIONER FULL_SINGLE_INVERSE
        MINIMIZER DIIS
      &END OT
    &END SCF"""
    return f"""&GLOBAL
  PROJECT {mode}
  RUN_TYPE {"GEO_OPT" if opt else "ENERGY"}
  PRINT_LEVEL LOW
&END GLOBAL
{motion}
&FORCE_EVAL
  METHOD Quickstep
  &DFT
    BASIS_SET_FILE_NAME BASIS_MOLOPT
    POTENTIAL_FILE_NAME GTH_POTENTIALS
    &MGRID
      CUTOFF {cutoff}
      REL_CUTOFF {rel}
    &END MGRID
    &QS
      EPS_DEFAULT 1.0E-10
    &END QS
{scf_block}
    &XC
      &XC_FUNCTIONAL PBE
      &END XC_FUNCTIONAL
    &END XC{printblk}
  &END DFT
  &SUBSYS
    &CELL
      ABC {cell[0]:.5f} {cell[1]:.5f} {cell[2]:.5f}
      PERIODIC XYZ
    &END CELL
    &TOPOLOGY
      COORD_FILE_NAME {coord}
      COORD_FILE_FORMAT XYZ
    &END TOPOLOGY
{kinds(basis, elements)}
  &END SUBSYS
&END FORCE_EVAL
"""


def write_set(src_dir, out_dir, opt):
    for xyz in sorted(src_dir.glob("*.xyz")):
        name = xyz.stem
        cell = json.loads(xyz.with_suffix(".json").read_text())["cell"]
        elements = {l.split()[0] for l in xyz.read_text().splitlines()[2:] if l.strip()}
        d = HERE / out_dir / name
        d.mkdir(parents=True, exist_ok=True)
        warm = HERE / "warmstart" / f"{name}.xyz"
        shutil.copy(warm if (opt and warm.exists()) else xyz, d / "start.xyz")   # warm start from an earlier relaxed geometry
        if opt:
            (d / "opt.inp").write_text(inp(name, cell, "opt", elements=elements))
            (d / "opt_robust.inp").write_text(inp(name, cell, "opt", robust=True, elements=elements))
        else:
            (d / "NO_OPT").write_text("single point on the given geometry\n")
        (d / "sp.inp").write_text(inp(name, cell, "sp", elements=elements))
        (d / "sp_robust.inp").write_text(inp(name, cell, "sp", robust=True, elements=elements))
        print("wrote", d)


def main():
    write_set(HERE / "structures", "runs", True)
    wet = HERE / "structures_wet"
    if wet.exists():
        write_set(wet, "runs_wet", False)


if __name__ == "__main__":
    main()

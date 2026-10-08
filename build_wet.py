#!/usr/bin/env python3
"""Put water (and optionally Mg-substitution + Na+) around the DFT-relaxed edge ribbons.

    python3 build_wet.py            # needs numpy + the `lammps` python module (any LAMMPS with KSPACE, MOLECULE)

For each relaxed ribbon (structures_relaxed/<rib>.xyz, cell in structures/<rib>.json):
  1. fill the vacuum gap between the two edges (the periodic image of the ribbon) with ~1 g/cm3 SPC water,
     confined in z to the layer thickness +- 1 A (so the water sees the edges, not a free surface),
  2. variant 'Na': one interior octahedral Al -> Mg (layer charge -1, same 1/8 density as the
     Si8Al3.5Mg0.5 montmorillonite model) and one Na+ in the water,
  3. equilibrate the water/ion with a frozen ribbon (classical ClayFF/SPC charges, LAMMPS NVT 300 K, 12 ps),
  4. write 2 snapshots (6 ps, 12 ps) as structures_wet/wet_<rib>_<w|wNa><k>.{xyz,json}.
The ribbon geometry is the DFT one, unchanged (the Mg site is NOT re-relaxed; Hirshfeld on Mg is only
approximate). The classical model is used ONLY to place water; the charges come from the CP2K single points.
Atom order in every xyz: ribbon atoms (as in the dry file), then Na, then water O H H.
"""
import json
import sys
from pathlib import Path

import numpy as np
from lammps import lammps

HERE = Path(__file__).resolve().parent
RIBBONS = ["rib_y_o00_si", "rib_y_o00_al", "rib_x_o00_si"]
EPS_SIG = {"st": (1.84e-6, 3.302), "ao": (1.33e-6, 4.2712), "mgo": (9.0e-7, 5.2643), "ob": (0.15540153, 3.1655),
           "oh": (0.15540153, 3.1655), "ho": (0.0, 0.0), "ow": (0.15540153, 3.1655), "hw": (0.0, 0.0),
           "na": (0.02909167, 2.610333)}
TYPES = ["st", "ao", "mgo", "ob", "oh", "ho", "ow", "hw", "na"]       # ids 1..9
MASS = {"st": 28.09, "ao": 26.98, "mgo": 24.31, "ob": 15.9994, "oh": 15.9994, "ho": 1.008, "ow": 15.9994, "hw": 1.008, "na": 22.99}
CHG = {"st": 2.1, "ao": 1.575, "mgo": 1.36, "ob": -1.05, "oh": -0.95, "ho": 0.425, "ow": -0.82, "hw": 0.41, "na": 1.0}
N_PER_A3 = 0.0284                                                      # molecules / A^3 at 1 g/cm3


def read_xyz(p):
    L = Path(p).read_text().splitlines()[2:]
    return [l.split()[0] for l in L if l.strip()], np.array([[float(x) for x in l.split()[1:4]] for l in L if l.strip()])


def mic(d, box):
    return d - box * np.round(d / box)


def classify(el, pos, box):
    """ClayFF-like type per ribbon atom (O by number of H neighbours; cap H of aquo groups = water-like)."""
    n = len(el)
    d = mic(pos[:, None] - pos[None], box)
    r = np.linalg.norm(d, axis=2)
    nH = np.array([sum(1 for j in range(n) if el[j] == "H" and r[i, j] < 1.25) if el[i] == "O" else 0 for i in range(n)])
    t = []
    for i, e in enumerate(el):
        if e == "Si": t.append("st")
        elif e == "Al": t.append("ao")
        elif e == "O": t.append("ob" if nH[i] == 0 else "oh" if nH[i] == 1 else "ow")
        else:
            o = next(j for j in range(n) if el[j] == "O" and r[i, j] < 1.25)
            t.append("ho" if nH[o] == 1 else "hw")
    return t


def random_rotation(rng):
    q = rng.normal(size=4); q /= np.linalg.norm(q)
    a, b, c, d = q
    return np.array([[a*a+b*b-c*c-d*d, 2*(b*c-a*d), 2*(b*d+a*c)], [2*(b*c+a*d), a*a-b*b+c*c-d*d, 2*(c*d-a*b)],
                     [2*(b*d-a*c), 2*(c*d+a*b), a*a-b*b-c*c+d*d]])


WATER = np.array([[0.0, 0.0, 0.0], [0.8165, 0.5773, 0.0], [-0.8165, 0.5773, 0.0]])   # SPC geometry (1.0 A, 109.47)


def build_case(rib, with_na, rng):
    el, pos = read_xyz(HERE / "structures_relaxed" / f"{rib}.xyz")
    box = np.array(json.loads((HERE / "structures" / f"{rib}.json").read_text())["cell"])
    ax = 1 if rib.startswith("rib_y") else 0
    el = list(el)
    types = classify(el, pos, box)
    mg_index = None
    if with_na:                                         # interior Al -> Mg (closest to the ribbon centre along the edge normal)
        centre = 0.5 * (pos[:, ax].min() + pos[:, ax].max())
        al = [i for i, e in enumerate(el) if e == "Al"]
        mg_index = min(al, key=lambda i: abs(pos[i, ax] - centre))
        el[mg_index] = "Mg"; types[mg_index] = "mgo"
    lo, hi = pos[:, ax].min(), pos[:, ax].max()
    zlo, zhi = pos[:, 2].min() - 1.0, pos[:, 2].max() + 1.0
    glo, ghi = hi + 2.0, lo + box[ax] - 2.0              # water gap along the edge normal
    other = 1 - ax
    vol = box[other] * (ghi - glo) * (zhi - zlo)
    nw = int(round(N_PER_A3 * vol))
    placed = []                                          # (O, H, H) positions
    allpos = pos.copy()
    extra = []
    tries = 0
    n_objects = nw + (1 if with_na else 0)
    kinds = ["na"] * (1 if with_na else 0) + ["w"] * nw
    for kind in kinds:
        tries = 0
        while True:
            tries += 1
            assert tries < 50000, "could not place water"
            c = np.zeros(3); c[ax] = rng.uniform(glo, ghi); c[other] = rng.uniform(0, box[other]); c[2] = rng.uniform(zlo + 1.3, zhi - 1.3)
            geo = (c[None] if kind == "na" else c + WATER @ random_rotation(rng).T)
            cur = np.vstack([allpos] + [np.array(e) for e in extra]) if extra else allpos
            dmin = min(np.linalg.norm(mic(cur - g, box), axis=1).min() for g in geo)
            if dmin > (2.0 if kind == "na" else 1.9):
                extra.append(geo.tolist()); break
    return dict(rib=rib, el=el, pos=pos, types=types, box=box, ax=ax, with_na=with_na, mg_index=mg_index,
                extra=extra, nw=nw, zlim=(zlo, zhi))


def write_data(case, path):
    n_r = len(case["el"])
    atoms, bonds, angles = [], [], []
    for t, p in zip(case["types"], case["pos"]):
        atoms.append((0, t, p))
    mol = 1
    nxt = n_r
    for g in case["extra"]:
        g = np.array(g)
        if len(g) == 1:
            mol += 1; atoms.append((mol, "na", g[0]))
        else:
            mol += 1
            atoms += [(mol, "ow", g[0]), (mol, "hw", g[1]), (mol, "hw", g[2])]
            i = len(atoms)
            bonds += [(1, i - 2, i - 1), (1, i - 2, i)]
            angles.append((1, i - 1, i - 2, i))
    q = np.array([CHG[t] for _, t, _ in atoms])
    # the classical ClayFF charges of a cut ribbon are not exactly neutral: spread the residual over the ribbon
    resid = q.sum()
    q[:n_r] -= resid / n_r
    L = case["box"]
    out = ["LAMMPS data: ribbon + water", "", f"{len(atoms)} atoms", f"{len(bonds)} bonds", f"{len(angles)} angles", "",
           f"{len(TYPES)} atom types", "1 bond types", "1 angle types", "",
           f"0 {L[0]:.6f} xlo xhi", f"0 {L[1]:.6f} ylo yhi", f"-10 {L[2] + 10:.6f} zlo zhi", "", "Masses", ""]
    out += [f"{i + 1} {MASS[t]}" for i, t in enumerate(TYPES)]
    out += ["", "Atoms # full", ""]
    for i, ((m, t, p), qq) in enumerate(zip(atoms, q), 1):
        out.append(f"{i} {m} {TYPES.index(t) + 1} {qq:.6f} {p[0]:.6f} {p[1]:.6f} {p[2]:.6f}")
    out += ["", "Bonds", ""] + [f"{i} {b[0]} {b[1]} {b[2]}" for i, b in enumerate(bonds, 1)]
    out += ["", "Angles", ""] + [f"{i} {a[0]} {a[1]} {a[2]} {a[3]}" for i, a in enumerate(angles, 1)]
    Path(path).write_text("\n".join(out) + "\n")
    return n_r, resid


def run_md(case, tag, seeds, snap_ps=(6, 12)):
    data = HERE / f"_tmp_{tag}.data"
    n_r, resid = write_data(case, data)
    lmp = lammps(cmdargs=["-log", "none", "-screen", "none"])
    c = lmp.command
    for l in ["units real", "atom_style full", "boundary p p p", "bond_style harmonic", "angle_style harmonic",
              "pair_style lj/cut/coul/long 10.0", "special_bonds lj/coul 0 0 0", f"read_data {data}"]:
        c(l)
    for i, ti in enumerate(TYPES, 1):
        for j, tj in enumerate(TYPES, 1):
            if j >= i:
                e = (EPS_SIG[ti][0] * EPS_SIG[tj][0]) ** 0.5
                s = 0.5 * (EPS_SIG[ti][1] + EPS_SIG[tj][1])
                c(f"pair_coeff {i} {j} {e:.8g} {s:.8g}")
    c("bond_coeff 1 1000 1.0"); c("angle_coeff 1 100 109.47")
    c("kspace_style ewald 1.0e-5"); c("neighbor 2.0 bin"); c("timestep 1.0")
    n_r_ = n_r
    c(f"group ribbon id 1:{n_r_}"); c("group mobile subtract all ribbon")
    c("group water type 7 8")
    c("thermo 2000"); c("thermo_style custom step temp pe")
    zlo, zhi = case["zlim"]
    c(f"region slab block INF INF INF INF {zlo - 0.5} {zhi + 0.5} units box")
    c("fix wall mobile wall/region slab harmonic 10.0 1.0 1.0")
    # relax the placed water/ion first (ribbon frozen: it is never integrated), flexible water
    c("fix frz ribbon setforce 0 0 0"); c("min_style cg"); c("minimize 1e-4 1e-6 500 5000"); c("unfix frz")
    c("velocity mobile create 300 %d dist gaussian" % seeds)
    c("fix shk water shake 1e-6 20 0 b 1 a 1")
    c("fix nvt mobile nvt temp 300 300 100")
    snaps = []
    steps_done = 0
    for ps in snap_ps:
        n = int(ps * 1000) - steps_done
        c(f"run {n}"); steps_done += n
        snaps.append((ps, np.array(lmp.gather_atoms("x", 1, 3)).reshape(-1, 3).copy(), None))
    lmp.close()
    return snaps


def unwrap_molecules(x, case):
    """Make every water molecule whole (O is the reference) and keep Na/ribbon as is."""
    box = case["box"]
    n_r = len(case["el"])
    x = x.copy()
    # ribbon atoms are frozen at their original coordinates
    x[:n_r] = case["pos"]
    k = n_r + (1 if case["with_na"] else 0)
    while k < len(x):
        o = x[k]
        for h in (k + 1, k + 2):
            x[h] = o + mic(x[h] - o, box)
        k += 3
    return x


def write_snapshots(case, snaps, name_base):
    n_r = len(case["el"])
    els = list(case["el"]) + (["Na"] if case["with_na"] else []) + ["O", "H", "H"] * case["nw"]
    out = HERE / "structures_wet"
    out.mkdir(exist_ok=True)
    for k, (ps, x, _) in enumerate(snaps, 1):
        xs = unwrap_molecules(x, case)
        name = f"{name_base}{k}"
        (out / f"{name}.xyz").write_text(f"{len(els)}\n{name} t={ps}ps nw={case['nw']} na={int(case['with_na'])}\n" + "".join(
            f"{e} {p[0]:.6f} {p[1]:.6f} {p[2]:.6f}\n" for e, p in zip(els, xs)))
        (out / f"{name}.json").write_text(json.dumps(dict(cell=case["box"].tolist(), pbc=[True, True, True], ribbon=case["rib"],
                                                          n_ribbon=n_r, n_water=case["nw"], na=case["with_na"],
                                                          mg_index=case["mg_index"], snapshot_ps=ps)))
        print("wrote", name, len(els), "atoms")


def main():
    rng = np.random.default_rng(7)
    for rib in RIBBONS:
        for with_na, tag in ((False, "w"), (True, "wNa")):
            case = build_case(rib, with_na, rng)
            print(rib, tag, "waters", case["nw"], flush=True)
            snaps = run_md(case, f"{rib}_{tag}", seeds=int(rng.integers(1, 10**6)))
            write_snapshots(case, snaps, f"wet_{rib}_{tag}")
    for p in HERE.glob("_tmp_*.data"):
        p.unlink()


if __name__ == "__main__":
    main()

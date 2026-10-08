#!/usr/bin/env python3
"""Build pyrophyllite edge ribbons (periodic along one axis) for DFT charge labels.

Bulk layer: ClayCode C2000 unit cell (Si8Al4O20(OH)4, 40 atoms, a=5.16 b=8.966 c=9.347 A).
A ribbon is a rectangular piece of one layer, cut along x or y. Bonds cut by the
cut are healed with the *bulk positions* of the missing ligands:
  - a metal that lost O ligands gets cap O atoms at the missing sites,
  - an O that lost a metal neighbour gets H on the missing-metal side,
  - metals left with < 3 (Si) / < 4 (Al) ligands are removed (iterate).
Hydrogens are then distributed over the edge O atoms so the cell is neutral with
formal charges (Si+4 Al+3 O-2 H+1). Different distributions are the protonation
variants: 'si' (SiOH first, Al caps end up bare OH/O), 'al' (Al-OH2 first),
'mix' (alternate). Geometry is crude; the CP2K pre-relaxation fixes it.

    python build_edges.py     # writes structures/<name>.xyz + .json (cell, pbc)
"""
import json
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
GRO = HERE / "C2000.gro"  # ClayCode CD21 pyrophyllite cell (MIT, see ClayCode.LICENSE.txt)
CUT = {("Si", "O"): 2.0, ("Al", "O"): 2.4, ("O", "H"): 1.2}
NEED = {"Si": 3, "Al": 4}  # minimum ligands to keep a metal
BOND = {"Si": 1.62, "Al": 1.92}
OH = 0.97


def read_gro(path):
    lines = path.read_text().splitlines()
    n = int(lines[1])
    names, pos = [], []
    for l in lines[2:2 + n]:
        nm = l[10:15].strip()
        names.append("Si" if nm.startswith("ST") else "Al" if nm.startswith("AO") else
                     "H" if nm.startswith("HO") else "O")
        pos.append([float(l[20:28]), float(l[28:36]), float(l[36:44])])
    box = np.array([float(x) for x in lines[2 + n].split()]) * 10
    pos = np.array(pos) * 10 % box
    return names, pos, box


def supercell(names, pos, box, nx, ny):
    N, P = [], []
    for i in range(nx):
        for j in range(ny):
            N += names
            P.append(pos + np.array([i * box[0], j * box[1], 0.0]))
    return N, np.vstack(P), np.array([nx * box[0], ny * box[1], box[2]])


def neighbours(names, pos, cell):
    """Bulk (periodic in x,y) neighbour table with minimum-image shift vectors."""
    n = len(names)
    nb = {i: [] for i in range(n)}
    for i in range(n):
        d = pos - pos[i]
        d[:, :2] -= np.round(d[:, :2] / cell[:2]) * cell[:2]
        r = np.linalg.norm(d, axis=1)
        for j in range(n):
            if j == i:
                continue
            key = (names[i], names[j])
            if key in CUT and r[j] < CUT[key]:
                nb[i].append((j, d[j]))
            elif key[::-1] in CUT and r[j] < CUT[key[::-1]]:
                nb[i].append((j, d[j]))
    return nb


def unit(v):
    return v / np.linalg.norm(v)


def build(nx, ny, axis, offset, variant, width_cells):
    """axis: 'y' -> ribbon periodic in x, cut normal to y; 'x' -> periodic in y."""
    n0, p0, b0 = read_gro(GRO)
    names, pos, cell = supercell(n0, p0, b0, nx, ny)
    nb = neighbours(names, pos, cell)
    ax = 1 if axis == "y" else 0
    L = cell[ax]
    w = width_cells * b0[ax]
    lo = offset * b0[ax]
    # keep window along the cut axis
    keep = np.array([(pos[i, ax] - lo) % L < w for i in range(len(names))])
    # iterative metal pruning
    changed = True
    while changed:
        changed = False
        for i in range(len(names)):
            if keep[i] and names[i] in NEED:
                if sum(keep[j] for j, _ in nb[i]) < NEED[names[i]]:
                    keep[i] = False
                    changed = True
        for i in range(len(names)):  # O (and its H) with no kept metal go too
            if keep[i] and names[i] == "O" and not any(keep[j] for j, _ in nb[i] if names[j] in NEED):
                keep[i] = False
                changed = True
        for i in range(len(names)):
            if keep[i] and names[i] == "H" and not any(keep[j] for j, _ in nb[i]):
                keep[i] = False
                changed = True

    atoms = [(names[i], pos[i].copy()) for i in range(len(names)) if keep[i]]
    idx = {i: k for k, i in enumerate(i for i in range(len(names)) if keep[i])}
    edge_O = []  # (atom index in list, 'Si'|'Al', direction toward lost neighbour)
    # (1) cap O on metals with missing ligands
    for i in range(len(names)):
        if not keep[i] or names[i] not in NEED:
            continue
        for j, d in nb[i]:
            if names[j] == "O" and not keep[j]:
                atoms.append(("O", pos[i] + d))
                edge_O.append((len(atoms) - 1, names[i], unit(d)))
    # (2) kept O that lost a metal -> edge O (pointing toward lost metal)
    for i in range(len(names)):
        if not keep[i] or names[i] != "O":
            continue
        lost = [d for j, d in nb[i] if names[j] in NEED and not keep[j]]
        have = [names[j] for j, _ in nb[i] if names[j] in NEED and keep[j]]
        if lost:
            kind = "Si" if "Si" in have else "Al" if "Al" in have else "Si"
            edge_O.append((idx[i], kind, unit(sum(unit(d) for d in lost))))
    # formal-charge neutrality -> hydrogens to add
    nSi = sum(a[0] == "Si" for a in atoms)
    nAl = sum(a[0] == "Al" for a in atoms)
    nO = sum(a[0] == "O" for a in atoms)
    nH_have = sum(a[0] == "H" for a in atoms)
    # structural OH whose O was cut away keeps its H only if O is kept (handled above)
    nH_need = 2 * nO - 4 * nSi - 3 * nAl
    extra = nH_need - nH_have
    # do existing H on edge O atoms count? hydroxyl H are bound to kept O -> already in nH_have
    def nh(k):  # hydrogens already bound to atom k
        return sum(1 for a in atoms if a[0] == "H" and np.linalg.norm(a[1] - atoms[k][1]) < 1.2)
    cap = {e[0]: (1 if e[1] == "Si" else 2) - nh(e[0]) for e in edge_O}
    sites = {"Si": [e for e in edge_O if e[1] == "Si" and cap[e[0]] > 0],
             "Al": [e for e in edge_O if e[1] == "Al" and cap[e[0]] > 0]}
    room = sum(cap[e[0]] for e in sites["Si"] + sites["Al"])
    if extra < 0 or extra > room:
        raise ValueError(f"cannot neutralise: need {extra} extra H, room {room}")
    s, a = sites["Si"], sites["Al"]
    if variant == "si":
        order = s + a + a                    # SiOH first, then Al caps (second pass -> aquo)
    elif variant == "al":
        order = a + a + s                    # Al aquo first
    else:
        order = [x for pair in zip(s + [None] * len(a), a + [None] * len(s)) for x in pair if x] + a
    plan = []  # edge-O entries receiving one H each
    left = extra
    for e in order:
        if left <= 0:
            break
        if sum(1 for p in plan if p[0] == e[0]) < cap[e[0]]:
            plan.append(e)
            left -= 1
    if left > 0:
        raise ValueError(f"{left} protons left over")
    counts = {}
    for e in plan:
        k = counts.get(e[0], 0) + nh(e[0])  # existing H occupy the first slot
        counts[e[0]] = counts.get(e[0], 0) + 1
        o = atoms[e[0]][1]
        d = e[2]
        # first H along the direction of the (lost) neighbour; second H tilted ~105 deg
        perp = unit(np.cross(d, [0.0, 0.0, 1.0] if abs(d[2]) < 0.9 else [1.0, 0.0, 0.0]))
        v = d if k == 0 else unit(d * np.cos(np.radians(105)) + perp * np.sin(np.radians(105)))
        atoms.append(("H", o + OH * v))
    # vacuum: shift ribbon, cell = nx*a (periodic) x width + vacuum, z gets vacuum
    P = np.array([a[1] for a in atoms])
    out_cell = cell.copy()
    out_cell[ax] = w + 14.0
    out_cell[2] = b0[2] + 8.0
    P[:, ax] = ((P[:, ax] - lo + 4.0) % L) - 4.0 + 7.0  # caps may sit just outside the window
    P[:, 2] += 4.0
    pbc = [True, False, False] if axis == "y" else [False, True, False]
    el = [a[0] for a in atoms]
    return el, P, out_cell, pbc, dict(nSi=nSi, nAl=nAl, nO=nO, nH=sum(e == "H" for e in el),
                                      edge_O=len(edge_O), extra_H=extra, variant=variant)


def write(name, el, P, cell, pbc, info):
    with open(HERE / "structures" / f"{name}.xyz", "w") as f:
        f.write(f"{len(el)}\n{name} {json.dumps(info)}\n")
        for e, p in zip(el, P):
            f.write(f"{e} {p[0]:.6f} {p[1]:.6f} {p[2]:.6f}\n")
    (HERE / "structures" / f"{name}.json").write_text(json.dumps(dict(cell=cell.tolist(), pbc=pbc, **info)))
    print(name, len(el), "atoms", info)


def main():
    # edge normal to y, ribbon periodic along x (1 cell wide): two cut offsets x three protonation variants
    for off in (0.0, 0.37):
        for var in ("si", "al"):
            write(f"rib_y_o{int(off * 100):02d}_{var}", *build(1, 3, "y", off, var, 2.0))
    # held-out geometry family: edge normal to x, ribbon periodic along y
    for var in ("si", "al"):
        write(f"rib_x_o00_{var}", *build(3, 1, "x", 0.0, var, 2.0))


if __name__ == "__main__":
    main()

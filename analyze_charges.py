#!/usr/bin/env python3
"""test3 step 3: learn SOAP -> DFT (Hirshfeld) charge, validate on held-out edge geometries.

Data: bulk pyrophyllite cell + four y-normal ribbons (train) ; two x-normal ribbons
(held out: a different edge orientation never seen in training). Also reports a
leave-one-ribbon-out score within the training family.

Baselines compared on the same atoms:
  fixed-type : mean *bulk* Hirshfeld charge per ClayFF-like class (Si, Al, O with H
               -> 'oh', O without H -> 'ob', H), i.e. the "human-assigned per-type" idea;
  SOAP-KRR   : same model family as test2 (polynomial-2 kernel on normalised SOAP).
Edge zone = atoms within 4.5 A of the ribbon edge (the only place the types are doubtful).

    python analyze_charges.py   # after runs/<name>/sp.out exist for all structures
"""
import json
import re
from pathlib import Path

import numpy as np
from ase import Atoms
from dscribe.descriptors import SOAP
from sklearn.kernel_ridge import KernelRidge

HERE = Path(__file__).resolve().parent
CLAYFF = {"Si": 2.1, "Al": 1.575, "H": 0.425, "ob": -1.05, "oh": -0.95}
TRAIN = ["bulk", "rib_y_o00_si", "rib_y_o00_al", "rib_y_o37_si", "rib_y_o37_al"]
TEST = ["rib_x_o00_si", "rib_x_o00_al"]
EDGE_ZONE = 4.5


def read_xyz(path):
    L = Path(path).read_text().splitlines()[2:]
    return [l.split()[0] for l in L], np.array([[float(x) for x in l.split()[1:4]] for l in L])


def hirshfeld(out, n):
    """Parse 'Hirshfeld Charges' block; returns net charges (len n)."""
    T = Path(out).read_text().splitlines()
    i = next(k for k, l in enumerate(T) if "Hirshfeld Charges" in l)
    q = []
    for l in T[i + 3:]:
        f = l.split()
        if len(f) >= 6 and f[0].isdigit():
            q.append(float(f[-1]))
        elif q and not f:
            break
    return np.array(q[:n])


def load(name):
    el, pos = read_xyz(HERE / "runs" / name / "relaxed.xyz")
    meta = json.loads((HERE / "structures" / f"{name}.json").read_text())
    q = hirshfeld(HERE / "runs" / name / "sp.out", len(el))
    assert len(q) == len(el), (name, len(q), len(el))
    return dict(name=name, el=el, pos=pos, cell=np.array(meta["cell"]), pbc=meta["pbc"], q=q)


def classes(el, pos, cell):
    """ClayFF-like class: O with an H within 1.2 A -> 'oh', else 'ob'."""
    d = pos[:, None] - pos[None]
    d -= np.round(d / cell) * cell
    r = np.linalg.norm(d, axis=2)
    out = []
    for i, e in enumerate(el):
        if e == "O":
            out.append("oh" if any(el[j] == "H" and r[i, j] < 1.2 for j in range(len(el))) else "ob")
        else:
            out.append(e)
    return np.array(out)


def edge_mask(s):
    if s["name"] == "bulk":
        return np.zeros(len(s["el"]), bool)
    ax = s["pbc"].index(False) if False in s["pbc"] else 1
    ax = next(k for k in (0, 1) if not s["pbc"][k])
    c = s["pos"][:, ax]
    metal = np.array([e in ("Si", "Al") for e in s["el"]])
    lo, hi = c[metal].min(), c[metal].max()
    return (c - lo < EDGE_ZONE) | (hi - c < EDGE_ZONE)


def ready(n):
    p = HERE / "runs" / n / "sp.out"
    return p.exists() and "PROGRAM ENDED" in p.read_text()


def main():
    global TRAIN, TEST
    TRAIN, TEST = [n for n in TRAIN if ready(n)], [n for n in TEST if ready(n)]
    print("structures: train", TRAIN, "test", TEST)
    soap = SOAP(species=["H", "O", "Al", "Si"], r_cut=5.0, n_max=4, l_max=3, sigma=0.4, periodic=True)
    S = {n: load(n) for n in TRAIN + TEST}
    for s in S.values():
        X = soap.create(Atoms(s["el"], positions=s["pos"] % s["cell"], cell=s["cell"], pbc=True), n_jobs=1)
        s["X"] = X / np.linalg.norm(X, axis=1, keepdims=True)
        s["cls"] = classes(s["el"], s["pos"], s["cell"])
        s["edge"] = edge_mask(s)

    # baseline: bulk mean per class
    b = S["bulk"]
    table = {c: float(b["q"][b["cls"] == c].mean()) for c in set(b["cls"])}
    scale = {c: CLAYFF[c] / table[c] for c in table}
    print("bulk Hirshfeld per class vs ClayFF:")
    for c in sorted(table):
        print(f"  {c:3s} DFT {table[c]:+.3f}   ClayFF {CLAYFF[c]:+.3f}   ratio {scale[c]:.2f}")

    def fit(names):
        X = np.vstack([S[n]["X"] for n in names])
        y = np.concatenate([S[n]["q"] for n in names])
        return KernelRidge(kernel="polynomial", degree=2, gamma=1.0, coef0=0.0, alpha=1e-4).fit(X, y)

    def score(model, n):
        s = S[n]
        base = np.array([table[c] for c in s["cls"]])
        pred = model.predict(s["X"])
        r = {}
        for tag, m in (("all", np.ones(len(base), bool)), ("edge", s["edge"])):
            if m.any():
                r[tag] = dict(n=int(m.sum()), mae_fixed=float(np.abs(base - s["q"])[m].mean()),
                              mae_soap=float(np.abs(pred - s["q"])[m].mean()),
                              max_fixed=float(np.abs(base - s["q"])[m].max()),
                              max_soap=float(np.abs(pred - s["q"])[m].max()))
        return r, pred

    results = {}
    print("\n== held-out edge orientation (train: bulk + y-ribbons, test: x-ribbons) ==")
    model = fit(TRAIN)
    for n in TEST:
        results[n], _ = score(model, n)
        for tag, v in results[n].items():
            print(f"{n:14s} {tag:4s} n={v['n']:3d}  MAE fixed {v['mae_fixed']:.4f}  SOAP {v['mae_soap']:.4f}"
                  f"   max err fixed {v['max_fixed']:.3f}  SOAP {v['max_soap']:.3f}")
    print("\n== leave-one-ribbon-out within the y-family ==")
    for n in TRAIN[1:]:
        m = fit([t for t in TRAIN if t != n])
        results[n], _ = score(m, n)
        for tag, v in results[n].items():
            print(f"{n:14s} {tag:4s} n={v['n']:3d}  MAE fixed {v['mae_fixed']:.4f}  SOAP {v['mae_soap']:.4f}"
                  f"   max err fixed {v['max_fixed']:.3f}  SOAP {v['max_soap']:.3f}")

    # final model on everything, expressed on the ClayFF scale for MD use
    if not TEST:
        return
    full = fit(TRAIN + TEST)
    elem_scale = {e: np.mean([scale[c] for c in scale if c.startswith(e.lower()) or c == e])
                  for e in ("Si", "Al", "H", "O")}
    (HERE / "results").mkdir(exist_ok=True)
    np.savez(HERE / "results" / "dft_charges.npz", **{n: S[n]["q"] for n in S},
             table=json.dumps(table), scale=json.dumps(scale))
    (HERE / "results" / "dft_validation.json").write_text(json.dumps(
        dict(bulk_class_charge=table, clayff=CLAYFF, ratio_clayff_over_dft=scale, results=results,
             edge_zone_A=EDGE_ZONE, label="Hirshfeld (CP2K, PBE, DZVP-MOLOPT-SR-GTH)"), indent=2))
    import pickle
    pickle.dump(dict(model=full, scale=scale, elem_scale=elem_scale, table=table), open(HERE / "results" / "dft_model.pkl", "wb"))


if __name__ == "__main__":
    main()

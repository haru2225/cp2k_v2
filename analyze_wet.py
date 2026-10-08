#!/usr/bin/env python3
"""How much does water (and Mg/Na) change the Hirshfeld charges of the ribbon atoms?

    python3 analyze_wet.py        # needs runs/<rib>/sp.out (dry) and runs_wet/<case>/sp.out ; numpy only

For every wet case the first n_ribbon atoms are the same ribbon atoms as in the dry (relaxed) calculation, so
dq_i = q_wet,i - q_dry,i is a per-atom polarisation / charge-transfer response to water (+ Na+).
For the 'Na' variant one Al is replaced by Mg and the layer carries -1 e, so its dq also contains the effect of
the substitution (not separable here); the water-only cases are the clean polarisation test.
"""
import json
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent


def hirshfeld(out):
    T = Path(out).read_text().splitlines()
    i = next(k for k, l in enumerate(T) if "Hirshfeld Charges" in l)
    el, q = [], []
    for l in T[i + 3:]:
        f = l.split()
        if len(f) >= 6 and f[0].isdigit():
            el.append(f[1]); q.append(float(f[-1]))
        elif q and not f:
            break
    return el, np.array(q)


def main():
    rows = []
    for d in sorted((HERE / "runs_wet").glob("*/")):
        if not (d / "sp.out").exists():
            continue
        meta = json.loads((HERE / "structures_wet" / f"{d.name}.json").read_text())
        rib, nr = meta["ribbon"], meta["n_ribbon"]
        dry_out = HERE / "runs" / rib / "sp.out"
        if not dry_out.exists():
            print(d.name, ": dry reference missing"); continue
        el_w, q_w = hirshfeld(d / "sp.out")
        el_d, q_d = hirshfeld(dry_out)
        dq = q_w[:nr] - q_d
        sel = np.ones(nr, bool)
        if meta["na"]:
            sel[meta["mg_index"]] = False          # the substituted site changes element
        k = nr + (1 if meta["na"] else 0)
        qO = q_w[k::3]; qH = np.concatenate([q_w[k + 1::3], q_w[k + 2::3]])
        print(f"{d.name}: ribbon total charge dry {q_d.sum():+.3f} -> wet {q_w[:nr].sum():+.3f} e"
              f" | water O {qO.mean():+.3f} H {qH.mean():+.3f} (net/molecule {(qO.mean()+2*qH.mean()):+.3f})"
              + (f" | Na {q_w[nr]:+.3f}" if meta["na"] else ""))
        for e in ("Si", "Al", "O", "H"):
            m = np.array([x == e for x in el_d]) & sel
            if m.any():
                print(f"   {e:2s} dq mean {dq[m].mean():+.4f}  rms {np.sqrt((dq[m]**2).mean()):.4f}  max|dq| {np.abs(dq[m]).max():.4f}")
        rows.append((d.name, dq[sel]))
    if not rows:
        print("no finished runs_wet/*/sp.out found")


if __name__ == "__main__":
    main()

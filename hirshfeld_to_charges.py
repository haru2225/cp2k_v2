#!/usr/bin/env python3
"""runs/<name>/sp.out (CP2K Hirshfeld block) -> runs/<name>/charges.dat  (atom order = structures/<name>.xyz).

    python3 hirshfeld_to_charges.py <name>
"""
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent


def main(name):
    d = HERE / "runs" / name
    lines = (d / "sp.out").read_text().splitlines()
    i = next(k for k, l in enumerate(lines) if "Hirshfeld Charges" in l)
    rows = []
    for l in lines[i + 3:]:
        f = l.split()
        if len(f) >= 6 and f[0].isdigit():
            rows.append((int(f[0]) - 1, f[1], float(f[-1])))
        elif rows and not f:
            break
    n = int((d / "start.xyz").read_text().splitlines()[0])
    assert len(rows) == n, (len(rows), n)
    (d / "charges.dat").write_text("# index element net_charge(Hirshfeld)\n" + "".join(
        f"{i} {e} {q:.5f}\n" for i, e, q in rows))
    print(f"{name}: total charge {sum(q for _, _, q in rows):+.4f} e over {n} atoms")


if __name__ == "__main__":
    main(sys.argv[1])

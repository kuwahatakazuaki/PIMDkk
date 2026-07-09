#!/usr/bin/env bash
set -euo pipefail

script_dir=$(cd "$(dirname "$0")" && pwd)
pimd8=${PIMD8:-"$script_dir/../.."}

cd "$script_dir"
rm -rf Result Scr std.out std.err
mkdir -p Result Scr

if [ ! -x "$pimd8/pimd.exe" ]; then
  make -C "$pimd8" pimd.exe >/dev/null
fi

"$pimd8/pimd.exe" >/dev/null

python3 - <<'PY'
from math import isfinite
from pathlib import Path

path = Path("Result/pressure.dat")
if not path.exists():
    raise SystemExit("FAIL: Result/pressure.dat was not created")

lines = path.read_text().splitlines()
if not lines or not lines[0].startswith("# 1step 2time_fs"):
    raise SystemExit("FAIL: pressure.dat header is missing or unexpected")

rows = [line.split() for line in lines[1:] if line.strip()]
if len(rows) != 3:
    raise SystemExit(f"FAIL: expected 3 pressure rows for steps 0,1,2; got {len(rows)}")

expected_steps = [0, 1, 2]
for row, expected_step in zip(rows, expected_steps):
    if len(row) != 10:
        raise SystemExit(f"FAIL: expected 10 columns, got {len(row)} in: {' '.join(row)}")
    step = int(row[0])
    values = [float(x) for x in row[1:]]
    if step != expected_step:
        raise SystemExit(f"FAIL: expected step {expected_step}, got {step}")
    if not all(isfinite(x) for x in values):
        raise SystemExit(f"FAIL: non-finite pressure row: {' '.join(row)}")

    p_cv = float(row[2])
    p_prim = float(row[3])
    volume = float(row[5])
    a_len, b_len, c_len = map(float, row[6:9])

    if abs(p_cv - p_prim) > max(1.0e-7, 1.0e-12 * abs(p_cv)):
        raise SystemExit(f"FAIL: Nbead=1 should give P_cv == P_prim; got {p_cv} vs {p_prim}")
    if abs(volume - 1000.0) > 1.0e-8:
        raise SystemExit(f"FAIL: expected volume 1000 A^3, got {volume}")
    if max(abs(a_len - 10.0), abs(b_len - 10.0), abs(c_len - 10.0)) > 1.0e-10:
        raise SystemExit(f"FAIL: expected 10 A cell lengths, got {a_len}, {b_len}, {c_len}")

print("PASS")
PY

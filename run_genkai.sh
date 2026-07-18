#!/bin/bash
#PJM -L rscgrp='a-batch'
#PJM -L vnode-core=16
#PJM -L elapse=1:00:00
#PJM -o log.out
#PJM -e err.out

input="input.inp"
dir_pimd="/home/pj24003139/share/PIMDprogram"  # PIMD プログラム一式の配置先
dir_prog="$dir_pimd/PIMDkk"                    # PIMD プログラム本体の配置先
dir_exe="$dir_pimd/bin2"                       # 実行ファイルの配置先
dir_result=$(grep -A1 '$address_result' "$input" | tail -1)
dir_scratch=$(grep -A1 '$address_scr' "$input" | tail -1)
Lrestart=$(grep -A1 '$Lrestart' "$input" | tail -1)
Isimulation=$(grep -A1 '$Isimulation' "$input" | tail -1)

# Set the MACE model when using MACE.
# export MACE_PYTHON_DIR="$pimddir/PIMDkk/MACE"
export MACE_PYTHON_DIR="$dir_prog/MACE"  # MACE Python モジュール、変更不要。
export MACE_MODEL="medium.MACE.model"    # 使用する MACE モデルを設定する。

# Make user-installed Python packages and Matplotlib cache available in the job.
python_user_site=$(python3 -m site --user-site 2>/dev/null)
if [ -n "$python_user_site" ]; then
  export PYTHONPATH="${python_user_site}${PYTHONPATH:+:${PYTHONPATH}}"
fi
export MPLCONFIGDIR="${dir_scratch%/}/.matplotlib"
mkdir -p "$MPLCONFIGDIR"

vnode_core=${PJM_VNODE_CORE:-1}

# 2コア以上では MPI 並列版を、1コアでは逐次版を実行する。
if [ "$vnode_core" -gt 1 ]; then
  if ! cp "$dir_exe/pimd_mpi.exe" ./pimd_mpi.exe; then
    echo "Error: failed to copy $dir_exe/pimd_mpi.exe" >&2
    exit 1
  fi

  mpirun -np "$vnode_core" ./pimd_mpi.exe
else
  if ! cp "$dir_exe/pimd.exe" ./pimd.exe; then
    echo "Error: failed to copy $dir_exe/pimd.exe" >&2
    exit 1
  fi

  ./pimd.exe
fi

exit 0

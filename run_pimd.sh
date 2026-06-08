#!/bin/bash
input="input.inp"
output="err.out"
dir_pimd="/Users/kuwahatakazuaki/Desktop/PIMD/NewPIMD"
# dir_pimd="/home/pj24003139/share/PIMDprogram"
dir_prog="$dir_pimd/PIMD8"
dir_exe="$dir_pimd/bin"
dir_result=$(grep -A1 '$address_result' "$input" | tail -1)
dir_scratch=$(grep -A1 '$address_scr' "$input" | tail -1)
Nproc=$(grep -A1 '$The Number of Proccesors' "$input" | tail -1)
Lrestart=$(grep -A1 '$Lrestart' "$input" | tail -1)
Isimulation=$(grep -A1 '$Isimulation' "$input" | tail -1)

# Set the MACE model when using MACE.
export MACE_PYTHON_DIR="$dir_prog/MACE"
export MACE_MODEL="$MACE_PYTHON_DIR/small_mace_model.model"

export TIME="User: %U, System: %S, Elapsed: %E
CPU: %P, Resident size: %M KB
Share Text: %X, Unshare Data: %D Kb
Input: %I, Output: %O
Pagefaults, Major: %F, Minor: %R, Swaps: %W\n"


# Checking the setting
if [ "$dir_result" = "" ];then
  echo " ERROR!!! Result directory is missing"
  exit 1
elif [ ${Lrestart//./} == "F" ]; then
  echo -e " +++ Cleaning $dir_result\n"
  rm -f $dir_result/*
fi

if [ "$dir_scratch" = "" ];then
  echo " ERROR!!! Scratch directory is missing"
  exit 1
fi
echo -e " +++ Result directory is \n  $dir_result\n"
echo -e " +++ Scratch directory is \n  $dir_scratch\n"

# if [ "$1" = "" ];then
#   Nproc=`grep -A1 "The Number of Proccesors" $input |tail -1`
# else
#   Nproc="$1"
# fi
# echo " +++ Nproc is $Nproc +++"

if [ "$vnode_core" -gt 1 ]; then
  mpirun -np $Nproc "$pimddir/pimd_mpi.exe" > $output
else
  "$dir_exe/pimd.exe" > $output
fi

exit 0


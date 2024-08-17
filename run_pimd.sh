#!/bin/bash
input="input.inp"
output="output"
machinefile="mpihosts"
pimddir="/Users/kuwahatakazuaki/Desktop/PIMD/NewPIMD/PIMD6"
dir_result=`grep -A1 "address_result" $input |tail -1`
dir_scratch=`grep -A1 "address_scr" $input |tail -1`
Lrestart=`grep -A1 "Lrestart" $input |tail -1`
simulation=`grep -A1 "Simulation" $input |tail -1`

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

if [ "$1" = "" ];then
  Nproc=`grep -A1 "The Number of Proccesors" $input |tail -1`
else
  Nproc="$1"
fi
echo " +++ Nproc is $Nproc +++"

# nohup mpirun -np $Nproc $pimddir/pimd.exe > $output  &
nohup $pimddir/pimd.exe > $output &

exit 0




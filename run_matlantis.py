#!/usr/bin/env python3
from ase import Atoms
from ase.io import read
# from ase.calculators.emt import EMT
import os
import numpy as np
import sys
# gRPC のログレベルを設定
os.environ["GRPC_VERBOSITY"] = "ERROR"
os.environ["GRPC_TRACE"] = ""  # トレースログも無効化

# === For Matlantis ===
import pfp_api_client
from pfp_api_client.pfp.calculators.ase_calculator import ASECalculator
from pfp_api_client.pfp.estimator import Estimator, EstimatorCalcMode

# === Choose the EstimatorCalcMode from CRYSTAL, CRYSTAL_U0, CRYSTAL_PLUS_D3, MOLECULE ===
estimator = Estimator(calc_mode=EstimatorCalcMode.MOLECULE)
calculator = ASECalculator(estimator)
# === For Matlantis ===

# # === For Effective Medium Theory ===
# calculator = EMT()
# # === For Effective Medium Theory ===

Fforce = 'forces.out'
Fenergy = 'energy.out'

# (path_dir, Ista, Iend) = sys.argv[1:]
(path_dir, Ista, Iend, Lperi) = sys.argv[1:]
Ista = int(Ista)
Iend = int(Iend)

if Lperi == "T":
    extension = "vasp"
elif Lperi == "F":
    extension = "xyz"
else:
    print("ERROR!!! Bad statement for periodic condition in run_matlantis.py")
    sys.exit()

# if os.path.isfile(Fforce):
#     os.remove(Fforce)

all_forces = np.array([])
all_energy = np.array([])

for i in range(Ista,Iend+1):
    path_inp=''
    Fname = 'str{0:05}.'.format(i) + extension
    path_inp = path_dir + Fname

    atoms = read(path_inp)
    atoms.calc = calculator

    energy = atoms.get_total_energy()
    forces = atoms.get_forces()
    # charges = atoms.get_charges()

    if i == Ista:
        all_forces = forces
    else:
        all_forces = np.concatenate([all_forces,forces])
    all_energy = np.append(all_energy,energy)


np.savetxt(Fforce,all_forces)
np.savetxt(Fenergy,all_energy)



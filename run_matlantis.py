#!/usr/bin/env python3
from ase import Atoms
from ase.io import read
from ase.calculators.emt import EMT
import os
import numpy as np
import sys

# === For Matlantis === 
# import pfp_api_client
# from pfp_api_client.pfp.calculators.ase_calculator import ASECalculator
# from pfp_api_client.pfp.estimator import Estimator, EstimatorCalcMode

# estimator = Estimator()
# calculator = ASECalculator(estimator)
# === For Matlantis === 

# === For Effective Medium Theory === 
calculator = EMT()
# === For Effective Medium Theory === 

Fforce = 'forces.out'
Fenergy = 'energy.out'

(path_dir, Ista, Iend) = sys.argv[1:]
Ista = int(Ista)
Iend = int(Iend)

if os.path.isfile(Fforce):
    os.remove(Fforce)

all_forces = np.array([])
all_energy = np.array([])

for i in range(Ista,Iend+1):
    path_inp=''
    Fname = 'str{0:05}.xyz'.format(i)
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
# print(all_energy)



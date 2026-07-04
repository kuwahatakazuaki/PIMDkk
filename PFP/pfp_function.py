import numpy as np
from ase import Atoms
from pfp_api_client.pfp.calculators.ase_calculator import ASECalculator
from pfp_api_client.pfp.estimator import Estimator, EstimatorCalcMode

_CALCULATOR_CACHE = {}

_SUPPORTED_CALC_MODES = ("MOLECULE", "CRYSTAL", "CRYSTAL_U0", "CRYSTAL_PLUS_D3")


def calculate_energy_and_forces(atomic_numbers, positions, cell, pbc, calc_mode):
    """
    Calculate potential energy and forces of a structure using Matlantis PFP.

    Parameters
    ----------
    atomic_numbers : list or numpy array
        Atomic numbers of the atoms.
    positions : numpy array
        Cartesian coordinates of the atoms, shape (n_atoms, 3), in Angstrom.
    cell : numpy array
        Lattice vectors, shape (3, 3), in Angstrom. Ignored when pbc is False.
    pbc : bool
        Periodic boundary condition flag.
    calc_mode : str
        Name of an EstimatorCalcMode member (e.g. "MOLECULE", "CRYSTAL").

    Returns
    -------
    energy : float
        Potential energy in eV.
    forces : numpy array
        Forces acting on each atom in eV/Angstrom, shape (n_atoms, 3).
    """
    atomic_numbers = np.array(atomic_numbers)
    positions = np.array(positions)

    atoms = Atoms(numbers=atomic_numbers, positions=positions)
    atoms.pbc = bool(pbc)
    if pbc:
        atoms.cell = np.array(cell, dtype=float)

    atoms.calc = _get_calculator(calc_mode)

    energy = atoms.get_potential_energy()
    forces = atoms.get_forces()

    return energy, forces


def _get_calculator(calc_mode):
    if calc_mode in _CALCULATOR_CACHE:
        return _CALCULATOR_CACHE[calc_mode]

    if calc_mode not in _SUPPORTED_CALC_MODES:
        raise ValueError(
            f"Unknown PFP calc_mode: {calc_mode!r}. "
            f"Supported modes are: {', '.join(_SUPPORTED_CALC_MODES)}"
        )

    estimator = Estimator(calc_mode=EstimatorCalcMode[calc_mode])
    calculator = ASECalculator(estimator)

    _CALCULATOR_CACHE[calc_mode] = calculator
    return calculator

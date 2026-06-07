import torch
from mace.calculators import MACECalculator
from ase import Atoms
import numpy as np
import warnings
import os

_MACE_PATCHED = False
_CALCULATOR_CACHE = {}


def _ensure_mace_compatibility():
    global _MACE_PATCHED
    if _MACE_PATCHED:
        return

    try:
        from e3nn.util.codegen import _mixin as _codegen_mixin
    except Exception:  # pragma: no cover - e3nn not available
        return

    original_setstate = getattr(_codegen_mixin.CodeGenMixin, "__setstate__", None)
    if not callable(original_setstate):
        _MACE_PATCHED = True
        return

    def _compat_setstate(self, state):
        if "__codegen__" in state and isinstance(state["__codegen__"], dict):
            converted = {}
            changed = False
            for key, val in state["__codegen__"].items():
                if isinstance(val, bytes):
                    converted[key] = ("torchscript", val)
                    changed = True
                else:
                    converted[key] = val
            if changed:
                state = dict(state)
                state["__codegen__"] = converted
        return original_setstate(self, state)

    _codegen_mixin.CodeGenMixin.__setstate__ = _compat_setstate

    try:
        from e3nn.o3 import _spherical_harmonics as _sh
        from e3nn.nn._activation import Activation as _Activation
    except Exception:  # pragma: no cover - optional dependencies missing
        _MACE_PATCHED = True
        return

    def _post_patch(model):
        for module in model.modules():
            if module.__class__.__name__ == "SphericalHarmonics" and not hasattr(module, "sph_func"):
                module.sph_func = _sh._spherical_harmonics
            if isinstance(module, _Activation) and not hasattr(module, "paths"):
                module.paths = [(mul, ir, act) for (mul, ir), act in zip(module.irreps_in, module.acts)]

    _codegen_mixin._mace_post_patch = _post_patch  # type: ignore[attr-defined]
    _MACE_PATCHED = True

warnings.filterwarnings("ignore")
os.environ['PYTHONWARNINGS'] = 'ignore'


def calculate_energy_and_forces(atomic_numbers, positions, model_path, cell=None, pbc=True, device="cpu"):
    """
    Calculate potential energy and forces of a molecule using MACE.
    
    Parameters:
    -----------
    atomic_numbers : list or numpy array
        Atomic numbers of the atoms in the molecule
    positions : numpy array
        Cartesian coordinates of the atoms in the molecule, shape (n_atoms, 3)
    
    Returns:
    --------
    energy : float
        Potential energy of the molecule in eV
    forces : numpy array
        Forces acting on each atom in eV/Å, shape (n_atoms, 3)
    """
    # Convert inputs to appropriate formats if needed
    atomic_numbers = np.array(atomic_numbers)
    positions = np.array(positions)
    
    # Create an ASE Atoms object
    atoms = Atoms(numbers=atomic_numbers, positions=positions)
    if cell is not None:
        atoms.cell = np.array(cell, dtype=float)
        atoms.pbc = pbc

    calculator = _get_calculator(model_path, positions, device)

    # Attach the calculator to the atoms
    atoms.calc = calculator
    
    # Calculate energy and forces
    energy = atoms.get_potential_energy()
    forces = atoms.get_forces()
    
    return energy, forces


def _get_calculator(model_path, positions, device):
    cache_key = (os.path.abspath(model_path), device)
    if cache_key in _CALCULATOR_CACHE:
        return _CALCULATOR_CACHE[cache_key]

    try:
        try:
            calculator = MACECalculator(model_paths=[model_path], device=device)
        except TypeError:
            calculator = MACECalculator(model_path, device=device)
        post_patch = getattr(__import__("e3nn.util.codegen", fromlist=["_mixin"])._mixin, "_mace_post_patch", None)  # type: ignore[attr-defined]
        if callable(post_patch):
            for model in getattr(calculator, "models", []):
                post_patch(model)
    except Exception as exc:
        raise RuntimeError(f"MACE model could not be loaded: {model_path}") from exc

    _CALCULATOR_CACHE[cache_key] = calculator
    return calculator


class _ZeroCalculator:
    def __init__(self, positions):
        self.positions = positions

    def get_potential_energy(self, atoms=None):
        return 0.0

    def get_forces(self, atoms=None):
        return np.zeros_like(np.array(self.positions), dtype=float)

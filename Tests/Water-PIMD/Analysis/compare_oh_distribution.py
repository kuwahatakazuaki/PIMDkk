#!/usr/bin/env python3
"""Compare SPC/F2 reference O-H densities with CL and QM sampled O-H distances."""

from __future__ import annotations

import os
from pathlib import Path
import shutil
import subprocess

_CACHE_DIR = Path(".compare_oh_distribution_cache")
_CACHE_DIR.mkdir(exist_ok=True)
(_CACHE_DIR / "matplotlib").mkdir(exist_ok=True)
(_CACHE_DIR / "xdg").mkdir(exist_ok=True)
os.environ.setdefault("MPLCONFIGDIR", str(_CACHE_DIR / "matplotlib"))
os.environ.setdefault("XDG_CACHE_HOME", str(_CACHE_DIR / "xdg"))

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

TRAJECTORIES = (
    ("CL", Path("../CL/Result/coor.xyz"), "#4C78A8"),
    ("QM", Path("../QM/Result/coor.xyz"), "#F58518"),
)                                                         # Input trajectories in XYZ format.
EQUIL_STEPS = 501                                           # Frames before this step are discarded.
BINS = 50                                                   # Number of histogram bins for sampled distances.
R_MIN_ANG = 0.60                                            # Plot range minimum in Angstrom.
R_MAX_ANG = 1.40                                            # Plot range maximum in Angstrom.
NUM_POINTS = 400                                            # Number of grid points for the theoretical curve.
OUTPUT_PATH = Path("oh_distribution.png")                   # Output comparison figure.
THEORY_DATA_PATH = Path("distribution_exact.dat")           # Output table for theoretical densities.
CLASSICAL_TEMPERATURE_K = 300.0                             # Classical Boltzmann curve temperature in Kelvin.

AU_LENGTH = 0.529177249e-10
AU_ENERGY = 4.3597482e-18
AU2ANG = AU_LENGTH * 1.0e10
AU_MASS = 9.1093897e-31
AMU = 1.66053906660e-27
KB_AU_PER_K = 3.166811563e-6

# SPC/F2 O-H constants copied from force_spcf.
RHO_W = 2.361e10 * AU_LENGTH
D_W = 0.708e-18 / AU_ENERGY
B_OH = 1.000e-10 / AU_LENGTH
MASS_O_AU = 15.999 * AMU / AU_MASS
MASS_H_AU = 1.00784 * AMU / AU_MASS


def ang_to_au(value_ang: np.ndarray | float) -> np.ndarray | float:
    return np.asarray(value_ang) / AU2ANG


def reduced_mass_oh() -> float:
    return MASS_O_AU * MASS_H_AU / (MASS_O_AU + MASS_H_AU)


def omega_oh() -> float:
    force_constant = 2.0 * RHO_W * RHO_W * D_W
    return np.sqrt(force_constant / reduced_mass_oh())


def psi0_oh(r_oh_au: np.ndarray) -> np.ndarray:
    mu = reduced_mass_oh()
    omega = omega_oh()
    prefactor = (mu * omega / np.pi) ** 0.25
    return prefactor * np.exp(-0.5 * mu * omega * (r_oh_au - B_OH) ** 2)


def psi0_sq_angstrom(r_oh_ang: np.ndarray) -> np.ndarray:
    # psi0_oh is normalized in Bohr, so divide by sqrt(AU2ANG) before squaring.
    psi0_ang = psi0_oh(ang_to_au(r_oh_ang)) / np.sqrt(AU2ANG)
    return psi0_ang * psi0_ang


def classical_density_angstrom(r_oh_ang: np.ndarray, temperature_k: float) -> np.ndarray:
    beta = 1.0 / (KB_AU_PER_K * temperature_k)
    r_oh_au = ang_to_au(r_oh_ang)
    curvature = beta * RHO_W * RHO_W * D_W
    density_au = np.sqrt(curvature / np.pi) * np.exp(-curvature * (r_oh_au - B_OH) ** 2)
    return density_au / AU2ANG


def read_frames(path: Path) -> list[tuple[int, list[list[str]]]]:
    frames = []
    with path.open() as f:
        while True:
            line = f.readline()
            if not line:
                break

            natom_total = int(line)
            step = int(f.readline().strip())
            atoms = [f.readline().split() for _ in range(natom_total)]
            if step < EQUIL_STEPS:
                continue
            frames.append((step, atoms))

    return frames


def oh_distances(frames: list[tuple[int, list[list[str]]]]) -> np.ndarray:
    distances = []

    for _, atoms in frames:
        if len(atoms) % 3 != 0:
            raise ValueError(f"Frame atom count {len(atoms)} is not divisible by 3.")

        for i in range(0, len(atoms), 3):
            triplet = atoms[i:i + 3]
            elements = [atom[0] for atom in triplet]
            if elements != ["O", "H", "H"]:
                raise ValueError(
                    f"Unexpected atom order {elements} at atoms {i + 1}-{i + 3}; expected O, H, H."
                )

            coords = np.array([[float(x), float(y), float(z)] for _, x, y, z in triplet])
            o = coords[0]
            h1 = coords[1]
            h2 = coords[2]
            distances.append(np.linalg.norm(h1 - o))
            distances.append(np.linalg.norm(h2 - o))

    return np.array(distances)


def build_histogram(distances: np.ndarray, bins: int) -> tuple[np.ndarray, np.ndarray]:
    density, edges = np.histogram(distances, bins=bins, density=True)
    centers = 0.5 * (edges[:-1] + edges[1:])
    return centers, density


def plot_comparison(
    r_ang: np.ndarray,
    psi0_sq_ang: np.ndarray,
    classical_density_ang: np.ndarray,
    sampled_curves: list[tuple[str, str, np.ndarray, np.ndarray]],
    path: Path,
) -> None:
    fig, ax = plt.subplots(figsize=(7.0, 4.8))

    for label, color, hist_centers, hist_density in sampled_curves:
        ax.plot(
            hist_centers,
            hist_density,
            color=color,
            linewidth=2.0,
            label=f"{label} sampled O-H histogram",
        )
    ax.plot(
        r_ang,
        psi0_sq_ang,
        color="tab:red",
        linewidth=2.2,
        label=r"quantum $\psi_0(r)^2$",
    )
    ax.plot(
        r_ang,
        classical_density_ang,
        color="black",
        linewidth=2.0,
        linestyle="--",
        label=f"classical Boltzmann ({CLASSICAL_TEMPERATURE_K:.1f} K)",
    )
    ax.set_xlim(R_MIN_ANG, R_MAX_ANG)
    ax.set_xlabel("O-H distance / Angstrom")
    ax.set_ylabel(r"Probability density / Angstrom$^{-1}$")
    ax.set_title("O-H distribution: CL, QM, quantum, and classical SPC/F2")
    ax.grid(True, alpha=0.3)
    ax.legend(loc="upper right")

    fig.tight_layout()
    fig.savefig(path, dpi=200)
    plt.close(fig)


def open_output(path: Path) -> None:
    if shutil.which("open") is None:
        return
    subprocess.run(["open", str(path)], check=False, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)


def main() -> None:
    if R_MIN_ANG >= R_MAX_ANG:
        raise SystemExit("R_MIN_ANG must be smaller than R_MAX_ANG")
    if NUM_POINTS < 2:
        raise SystemExit("NUM_POINTS must be at least 2")
    if BINS < 1:
        raise SystemExit("BINS must be at least 1")

    r_ang = np.linspace(R_MIN_ANG, R_MAX_ANG, NUM_POINTS)
    psi0_sq_ang = psi0_sq_angstrom(r_ang)
    classical_density_ang = classical_density_angstrom(r_ang, CLASSICAL_TEMPERATURE_K)

    theory_data = np.column_stack((r_ang, psi0_sq_ang, classical_density_ang))
    np.savetxt(
        THEORY_DATA_PATH,
        theory_data,
        header=(
            "r_oh_angstrom "
            "psi0_sq_angstrom_min1 "
            f"classical_density_angstrom_min1_T{CLASSICAL_TEMPERATURE_K:.1f}K"
        ),
    )

    summaries = []
    sampled_curves = []
    histogram_paths = []
    for label, xyz_path, color in TRAJECTORIES:
        frames = read_frames(xyz_path)
        if not frames:
            raise ValueError(f"No frames found at or after step {EQUIL_STEPS} in {xyz_path}.")

        distances = oh_distances(frames)
        hist_centers, hist_density = build_histogram(distances, BINS)
        hist_data_path = Path(f"distribution_simulation_{label.lower()}.dat")
        histogram_data = np.column_stack((hist_centers, hist_density))
        np.savetxt(
            hist_data_path,
            histogram_data,
            header="bin_center_angstrom probability_density",
        )

        sampled_curves.append((label, color, hist_centers, hist_density))
        histogram_paths.append(hist_data_path)
        summaries.append((label, xyz_path, frames, distances))

    plot_comparison(
        r_ang,
        psi0_sq_ang,
        classical_density_ang,
        sampled_curves,
        OUTPUT_PATH,
    )
    open_output(OUTPUT_PATH)

    print(f"equilibration cutoff: {EQUIL_STEPS}")
    print(f"classical temperature: {CLASSICAL_TEMPERATURE_K:.1f} K")
    for label, xyz_path, frames, distances in summaries:
        print(f"{label} source: {xyz_path}")
        print(f"{label} frames: {len(frames)}")
        print(f"{label} steps: {frames[0][0]} -> {frames[-1][0]}")
        print(f"{label} samples: {len(distances)}")
        print(f"{label} mean: {distances.mean():.6f} Angstrom")
        print(f"{label} std: {distances.std():.6f} Angstrom")
    print(f"saved: {OUTPUT_PATH}")
    print(f"saved: {THEORY_DATA_PATH}")
    for hist_path in histogram_paths:
        print(f"saved: {hist_path}")


if __name__ == "__main__":
    main()

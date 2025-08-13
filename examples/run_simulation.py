"""
Example: 2D TEz FDTD with CPML and a Gaussian Hz source.

Frames-only headless mode (default in Docker): saves PNGs to outputs/.
"""

import os
import sys

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
PROJECT_ROOT = os.path.dirname(SCRIPT_DIR)
if PROJECT_ROOT not in sys.path:
    sys.path.insert(0, PROJECT_ROOT)

import numpy as np
from fdtd_solver import FDTDSolverTE
from boundary_conditions import CPML, CPMLParams
from source import GaussianPulseSource
from visualization import animate_hz
from utils import compute_cfl_dt


def main():
    nx, ny = 200, 200
    dx = dy = 1.0e-3
    cfl = 0.7
    dt = compute_cfl_dt(dx, dy, cfl=cfl)

    total_steps = 600
    pml_thickness = 20

    solver = FDTDSolverTE(nx=nx, ny=ny, dx=dx, dy=dy, dt=dt, cfl=cfl)
    cpml_params = CPMLParams(
        thickness=pml_thickness,
        order_m=3,
        R_target=1e-6,
        kappa_max=4.0,
        alpha_max=0.02,
    )
    cpml = CPML(nx=nx, ny=ny, dx=dx, dy=dy, dt=solver.dt, params=cpml_params)

    src = GaussianPulseSource(i0=nx // 2, j0=ny // 2, t0=40.0, spread=12.0, amplitude=0.5)

    def step_fn(frame_idx: int):
        solver.step(cpml=cpml, source=src, time_step_index=frame_idx)

    def get_hz():
        return solver.Hz

    headless = os.environ.get("HEADLESS", "1") == "1"
    if headless:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt

        out_dir = os.path.join(PROJECT_ROOT, "outputs")
        os.makedirs(out_dir, exist_ok=True)
        save_every = 25
        for n in range(total_steps):
            step_fn(n)
            if n % save_every == 0 or n == total_steps - 1:
                hz = get_hz()
                extent = [0.0, nx * dx, 0.0, ny * dy]
                plt.figure(figsize=(6, 5))
                im = plt.imshow(hz.T, origin="lower", extent=extent, cmap="RdBu_r", interpolation="nearest", aspect="auto")
                plt.colorbar(im, label="Hz (A/m)")
                plt.xlabel("x [m]")
                plt.ylabel("y [m]")
                plt.title(f"Hz at step {n}")
                plt.tight_layout()
                out_path = os.path.join(out_dir, f"hz_step_{n:05d}.png")
                plt.savefig(out_path, dpi=120)
                plt.close()
        print(f"Headless run complete. Snapshots saved to: {out_dir}")
    else:
        out_dir = os.path.join(PROJECT_ROOT, "outputs")
        os.makedirs(out_dir, exist_ok=True)
        animate_hz(
            step_fn=step_fn,
            get_hz_fn=get_hz,
            nx=nx,
            ny=ny,
            dx=dx,
            dy=dy,
            n_frames=total_steps,
            interval_ms=20,
            title="2D FDTD TEz with CPML (Hz)",
        )
        import matplotlib.pyplot as plt
        plt.show()


if __name__ == "__main__":
    main()
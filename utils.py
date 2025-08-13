"""
Utility constants and helper functions for the FDTD project.
"""

from __future__ import annotations
import math


EPS0 = 8.854187817e-12
MU0 = 4.0e-7 * math.pi
C0 = 1.0 / math.sqrt(EPS0 * MU0)
ETA0 = math.sqrt(MU0 / EPS0)


def compute_cfl_dt(dx: float, dy: float, cfl: float = 0.7) -> float:
    if dx <= 0.0 or dy <= 0.0:
        raise ValueError("dx and dy must be positive.")
    return cfl / (C0 * math.sqrt((1.0 / dx) ** 2 + (1.0 / dy) ** 2))


def validate_parameters(nx: int, ny: int, dx: float, dy: float, dt: float, pml_thickness: int) -> None:
    if nx < 10 or ny < 10:
        raise ValueError("Grid size too small. Use at least 10x10.")
    if dx <= 0.0 or dy <= 0.0:
        raise ValueError("dx and dy must be positive.")
    if dt <= 0.0:
        raise ValueError("dt must be positive.")
    if pml_thickness < 1:
        raise ValueError("PML thickness must be >= 1 cell.")



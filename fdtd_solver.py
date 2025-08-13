"""
2D FDTD solver for electromagnetic wave propagation in TE (Hz, Ex, Ey) mode.

This module implements:
- Yee grid staggering for Ex, Ey, Hz
- Leapfrog time stepping (E updated at n+1/2, H at n+1)
- Vacuum update equations using ε0 and μ0
- Optional CPML absorbing boundary updates (delegated to boundary_conditions.CPML)

Field staggering (indices are array indices):
- Ex[i, j+1/2] -> stored as Ex[nx, ny+1]
- Ey[i+1/2, j] -> stored as Ey[nx+1, ny]
- Hz[i+1/2, j+1/2] -> stored as Hz[nx, ny]

Update equations in vacuum (discrete, centered differences):
- Ex^{n+1/2} = Ex^{n-1/2} + (dt/ε0) * ( (Hz^n[i, j] - Hz^n[i, j-1]) / dy )
- Ey^{n+1/2} = Ey^{n-1/2} - (dt/ε0) * ( (Hz^n[i, j] - Hz^n[i-1, j]) / dx )
- Hz^{n+1}   = Hz^{n}   + (dt/μ0) * ( (Ey^{n+1/2}[i+1, j] - Ey^{n+1/2}[i, j]) / dx
                                      - (Ex^{n+1/2}[i, j+1] - Ex^{n+1/2}[i, j]) / dy )

In the CPML case, the update is handled by boundary_conditions.CPML which augments the
derivatives with auxiliary ψ terms and κ stretching to produce near-zero reflections.
"""

from __future__ import annotations
from typing import Optional

import numpy as np

from utils import EPS0, MU0, compute_cfl_dt, validate_parameters


class FDTDSolverTE:
    """
    Finite-Difference Time-Domain solver for 2D TEz fields: Ex, Ey, Hz.

    Attributes:
    - nx, ny: number of Yee cells along x and y
    - dx, dy: spatial steps [m]
    - dt: time step [s]
    - Ex: shape (nx, ny+1)
    - Ey: shape (nx+1, ny)
    - Hz: shape (nx, ny)
    """

    def __init__(
        self,
        nx: int,
        ny: int,
        dx: float,
        dy: float,
        dt: Optional[float] = None,
        cfl: float = 0.7,
        seed: Optional[int] = None,
    ) -> None:
        if seed is not None:
            np.random.seed(seed)

        if dt is None:
            dt = compute_cfl_dt(dx, dy, cfl=cfl)

        validate_parameters(nx=nx, ny=ny, dx=dx, dy=dy, dt=dt, pml_thickness=10)

        self.nx = nx
        self.ny = ny
        self.dx = dx
        self.dy = dy
        self.dt = dt

        # Yee-staggered field arrays
        self.Ex = np.zeros((nx, ny + 1), dtype=np.float64)
        self.Ey = np.zeros((nx + 1, ny), dtype=np.float64)
        self.Hz = np.zeros((nx, ny), dtype=np.float64)

    def reset_fields(self) -> None:
        self.Ex.fill(0.0)
        self.Ey.fill(0.0)
        self.Hz.fill(0.0)

    # -----------------------------
    # Vacuum updates (no boundaries)
    # -----------------------------
    def update_e_fields_vacuum(self) -> None:
        """
        Update Ex, Ey using standard Yee FDTD update in vacuum (no PML).
        """
        # dHz/dy for Ex (shape nx x (ny+1)). Zero at j=0 and j=ny
        dHz_dy = np.zeros_like(self.Ex)
        dHz_dy[:, 1:-1] = self.Hz[:, 1:] - self.Hz[:, :-1]

        # dHz/dx for Ey (shape (nx+1) x ny). Zero at i=0 and i=nx
        dHz_dx = np.zeros_like(self.Ey)
        dHz_dx[1:-1, :] = self.Hz[1:, :] - self.Hz[:-1, :]

        self.Ex += (self.dt / EPS0) * (dHz_dy / self.dy)
        self.Ey -= (self.dt / EPS0) * (dHz_dx / self.dx)

    def update_h_field_vacuum(self) -> None:
        """
        Update Hz using standard Yee FDTD update in vacuum (no PML).
        """
        # dEy/dx (nx x ny)
        dEy_dx = self.Ey[1:, :] - self.Ey[:-1, :]
        # dEx/dy (nx x ny)
        dEx_dy = self.Ex[:, 1:] - self.Ex[:, :-1]

        self.Hz += (self.dt / MU0) * (dEy_dx / self.dx - dEx_dy / self.dy)

    def step(
        self,
        cpml=None,
        source=None,
        time_step_index: Optional[int] = None,
    ) -> None:
        """
        Advance one leapfrog step.

        - If cpml is provided (boundary_conditions.CPML), PML-aware updates are used.
        - Otherwise, vacuum updates are used.
        - If a `source` is provided (source.GaussianPulseSource), inject into Hz each step.

        Parameters:
        - cpml: instance of CPML or None
        - source: instance providing .inject(Hz, n)
        - time_step_index: current integer time step (for time-dependent sources)
        """
        if cpml is None:
            self.update_e_fields_vacuum()
            self.update_h_field_vacuum()
        else:
            cpml.update_e(self.Ex, self.Ey, self.Hz)
            cpml.update_h(self.Hz, self.Ex, self.Ey)

        if source is not None and time_step_index is not None:
            source.inject(self.Hz, time_step_index)



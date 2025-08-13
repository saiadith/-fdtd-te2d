"""
Perfectly Matched Layer (PML) absorbing boundary conditions using CPML for TEz.
"""

from __future__ import annotations
from dataclasses import dataclass
import numpy as np

from utils import EPS0, MU0, ETA0


@dataclass
class CPMLParams:
    thickness: int = 20
    order_m: int = 3
    R_target: float = 1e-6
    kappa_max: float = 4.0
    alpha_max: float = 0.02


class CPML:
    def __init__(self, nx: int, ny: int, dx: float, dy: float, dt: float, params: CPMLParams = CPMLParams()) -> None:
        self.nx, self.ny = nx, ny
        self.dx, self.dy = dx, dy
        self.dt = dt
        self.p = params

        self.kappa_x_e, self.sigma_x_e, self.alpha_x_e = self._make_profiles(nx + 1, params.thickness, dx, True)
        self.kappa_y_e, self.sigma_y_e, self.alpha_y_e = self._make_profiles(ny + 1, params.thickness, dy, True)
        self.kappa_x_h, self.sigma_x_h, self.alpha_x_h = self._make_profiles(nx, params.thickness, dx, False)
        self.kappa_y_h, self.sigma_y_h, self.alpha_y_h = self._make_profiles(ny, params.thickness, dy, False)

        self.b_x_e, self.a_x_e = self._coeffs(self.sigma_x_e, self.kappa_x_e, self.alpha_x_e, EPS0)
        self.b_y_e, self.a_y_e = self._coeffs(self.sigma_y_e, self.kappa_y_e, self.alpha_y_e, EPS0)
        self.b_x_h, self.a_x_h = self._coeffs(self.sigma_x_h, self.kappa_x_h, self.alpha_x_h, MU0)
        self.b_y_h, self.a_y_h = self._coeffs(self.sigma_y_h, self.kappa_y_h, self.alpha_y_h, MU0)

        self.psi_ex_y = np.zeros((nx, ny + 1), dtype=np.float64)
        self.psi_ey_x = np.zeros((nx + 1, ny), dtype=np.float64)
        self.psi_hz_x = np.zeros((nx, ny), dtype=np.float64)
        self.psi_hz_y = np.zeros((nx, ny), dtype=np.float64)

    def _sigma_max(self, cell: float) -> float:
        delta = self.p.thickness * cell
        return - (self.p.order_m + 1) * np.log(self.p.R_target) / (2.0 * ETA0 * delta)

    def _make_profiles(self, num: int, thickness: int, cell: float, is_e: bool):
        sigma = np.zeros(num)
        kappa = np.ones(num)
        alpha = np.zeros(num)
        half = 0.5 if is_e else 0.0
        smax = self._sigma_max(cell)
        kmax = self.p.kappa_max
        amax = self.p.alpha_max
        m = self.p.order_m

        for i in range(min(thickness, num)):
            x = (thickness - (i + half)) / thickness
            if x > 0:
                sigma[i] += smax * (x ** m)
                kappa[i] += (kmax - 1.0) * (x ** m)
                alpha[i] += amax * (1.0 - x)
        for i in range(num - 1, max(num - thickness - 1, -1), -1):
            xr = (num - 1 - i) + half
            x = (thickness - xr) / thickness
            if x > 0:
                sigma[i] += smax * (x ** m)
                kappa[i] += (kmax - 1.0) * (x ** m)
                alpha[i] += amax * (1.0 - x)
        return kappa, sigma, alpha

    def _coeffs(self, sigma: np.ndarray, kappa: np.ndarray, alpha: np.ndarray, xi: float):
        with np.errstate(divide="ignore", invalid="ignore"):
            b = np.exp(-((sigma / kappa) + alpha) * self.dt / xi)
            denom = sigma * kappa + (kappa ** 2) * alpha
            a = np.zeros_like(b)
            mask = denom > 0
            a[mask] = (sigma[mask] * (b[mask] - 1.0)) / denom[mask]
        interior = (sigma == 0.0) & (alpha == 0.0)
        b[interior] = 1.0
        a[interior] = 0.0
        return b, a

    def update_e(self, Ex: np.ndarray, Ey: np.ndarray, Hz: np.ndarray) -> None:
        dHz_dy = np.zeros_like(Ex)
        dHz_dy[:, 1:-1] = Hz[:, 1:] - Hz[:, :-1]
        dHz_dx = np.zeros_like(Ey)
        dHz_dx[1:-1, :] = Hz[1:, :] - Hz[:-1, :]

        self.psi_ex_y = self.b_y_e[np.newaxis, :] * self.psi_ex_y + self.a_y_e[np.newaxis, :] * dHz_dy
        Ex += (self.dt / EPS0) * ((dHz_dy / self.dy) / self.kappa_y_e[np.newaxis, :] + self.psi_ex_y / self.dy)

        self.psi_ey_x = self.b_x_e[:, np.newaxis] * self.psi_ey_x + self.a_x_e[:, np.newaxis] * dHz_dx
        Ey -= (self.dt / EPS0) * ((dHz_dx / self.dx) / self.kappa_x_e[:, np.newaxis] + self.psi_ey_x / self.dx)

    def update_h(self, Hz: np.ndarray, Ex: np.ndarray, Ey: np.ndarray) -> None:
        dEy_dx = Ey[1:, :] - Ey[:-1, :]
        dEx_dy = Ex[:, 1:] - Ex[:, :-1]

        self.psi_hz_x = self.b_x_h[:, np.newaxis] * self.psi_hz_x + self.a_x_h[:, np.newaxis] * dEy_dx
        self.psi_hz_y = self.b_y_h[np.newaxis, :] * self.psi_hz_y + self.a_y_h[np.newaxis, :] * dEx_dy

        Hz += (self.dt / MU0) * (
            (dEy_dx / self.dx) / self.kappa_x_h[:, np.newaxis] + self.psi_hz_x / self.dx
            - (dEx_dy / self.dy) / self.kappa_y_h[np.newaxis, :] - self.psi_hz_y / self.dy
        )



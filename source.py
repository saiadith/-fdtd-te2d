"""
Source utilities for injecting time-dependent fields (Gaussian pulse).
"""

from __future__ import annotations
from dataclasses import dataclass
import numpy as np


@dataclass
class GaussianPulseSource:
    i0: int
    j0: int
    t0: float
    spread: float
    amplitude: float = 0.5

    def value(self, n: int) -> float:
        return float(self.amplitude * np.exp(-((n - self.t0) ** 2) / (2.0 * (self.spread ** 2))))

    def inject(self, Hz: np.ndarray, n: int) -> None:
        if 0 <= self.i0 < Hz.shape[0] and 0 <= self.j0 < Hz.shape[1]:
            Hz[self.i0, self.j0] += self.value(n)



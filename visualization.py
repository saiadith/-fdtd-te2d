"""
Visualization helpers for 2D TEz FDTD simulations.
"""

from __future__ import annotations
from typing import Optional, Callable

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation


def show_hz_field(Hz: np.ndarray, dx: float, dy: float, title: str = "Hz field", vmin: Optional[float] = None, vmax: Optional[float] = None) -> None:
    nx, ny = Hz.shape
    extent = [0.0, nx * dx, 0.0, ny * dy]
    plt.figure(figsize=(6, 5))
    im = plt.imshow(Hz.T, origin="lower", extent=extent, cmap="RdBu_r", vmin=vmin, vmax=vmax, interpolation="nearest", aspect="auto")
    plt.colorbar(im, label="Hz (A/m)")
    plt.xlabel("x [m]")
    plt.ylabel("y [m]")
    plt.title(title)
    plt.tight_layout()
    plt.show()


def animate_hz(
    step_fn: Callable[[int], None],
    get_hz_fn: Callable[[], np.ndarray],
    nx: int,
    ny: int,
    dx: float,
    dy: float,
    n_frames: int,
    interval_ms: int = 30,
    vmin: Optional[float] = None,
    vmax: Optional[float] = None,
    title: str = "FDTD TEz: Hz propagation",
    save_path: Optional[str] = None,
) -> FuncAnimation:
    fig, ax = plt.subplots(figsize=(7, 5))
    extent = [0.0, nx * dx, 0.0, ny * dy]
    hz = get_hz_fn()
    if vmin is None or vmax is None:
        abs_max = np.max(np.abs(hz))
        if abs_max == 0:
            abs_max = 1.0
        vmin = -abs_max
        vmax = abs_max
    im = ax.imshow(hz.T, origin="lower", extent=extent, cmap="RdBu_r", vmin=vmin, vmax=vmax, interpolation="nearest", aspect="auto")
    fig.colorbar(im, ax=ax, label="Hz (A/m)")
    ax.set_xlabel("x [m]")
    ax.set_ylabel("y [m]")
    ax.set_title(title)

    def _update(frame: int):
        step_fn(frame)
        im.set_data(get_hz_fn().T)
        return [im]

    anim = FuncAnimation(fig, _update, frames=n_frames, interval=interval_ms, blit=True)
    plt.tight_layout()
    if save_path is not None:
        anim.save(save_path, dpi=120, writer="ffmpeg")
    return anim



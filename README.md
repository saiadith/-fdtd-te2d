## FDTD-TE2D: Modular 2D FDTD (TE mode) with CPML in Python

A clean, modular implementation of a 2D Finite-Difference Time-Domain (FDTD) solver for electromagnetic wave propagation in TE mode (fields: `Hz, Ex, Ey`). Includes Yee grid/leapfrog time stepping, CPML absorbing boundaries, a Gaussian pulse source, and visualization utilities.

### Features
- Yee grid TEz (`Hz, Ex, Ey`) with leapfrog time stepping
- Vacuum updates using ε0 and μ0
- CPML absorbing boundaries (configurable thickness and grading)
- Gaussian pulse source injection
- Headless frame rendering (PNG snapshots); GUI animation optional when run natively

### Repository structure
- `fdtd_solver.py` — core solver (Yee grid, updates, step loop)
- `boundary_conditions.py` — CPML implementation and parameters
- `source.py` — Gaussian pulse source
- `visualization.py` — colormap and animation helpers
- `utils.py` — constants (ε0, μ0, c), CFL, validation
- `examples/run_simulation.py` — runnable example (parameters, loop, saving)
- `Dockerfile`, `.dockerignore`, `requirements.txt`

### Quick start (Docker, recommended)
Frames-only by default; images are written to `outputs/` on your host.
```powershell
cd C:\Users\saiad\OneDrive\Desktop\FTDT
docker build -t fdtd-te .
mkdir outputs
docker run --rm -e HEADLESS=1 -e MPLBACKEND=Agg -v "FTDT\\outputs:/app/outputs" fdtd-te
```
Open a frame:
```powershell
ii outputs\hz_step_00000.png
```

### Native run (optional)
```powershell
python -m pip install --upgrade pip numpy matplotlib

# Headless (frames only)
$env:HEADLESS='1'; $env:MPLBACKEND='Agg'
python examples\run_simulation.py

# GUI (live animation; requires a desktop session)
$env:HEADLESS='0'
python examples\run_simulation.py
```

### Parameters (edit `examples/run_simulation.py`)
- Grid: `nx, ny` (e.g., 200 x 200)
- Spacing: `dx, dy` (meters)
- Time step: computed via CFL (`compute_cfl_dt`), controlled by `cfl`
- CPML: `CPMLParams(thickness, order_m, R_target, kappa_max, alpha_max)`
- Source: `GaussianPulseSource(i0, j0, t0, spread, amplitude)`
- Duration: `total_steps`
- Frame interval (headless): `save_every` (default 25)

### Output
- Headless runs save PNG snapshots in `outputs/` as `hz_step_XXXXX.png`.
- GUI runs show a live colormap of `Hz`. MP4 export is disabled by default in Docker (no ffmpeg). To enable MP4 natively, install ffmpeg on your host and extend the example to set `save_path` in `animate_hz`.

### Tuning and stability
- Use smaller `cfl` (e.g., 0.6–0.8) if you see instability or overflow warnings
- Increase CPML `thickness` (≥ 10–20) for stronger absorption at higher frequencies
- Ensure ≥ 10 cells per wavelength for accuracy (decrease `dx, dy` for higher-frequency content)

### Troubleshooting
- Docker volume path on Windows must be absolute; example above uses an absolute host path
- If Docker build fails due to cache issues, try: `docker builder prune -af` and rebuild

### References
- A. Taflove and S. C. Hagness, Computational Electrodynamics (3rd ed.)
- J.-P. Berenger, “A perfectly matched layer for the absorption of electromagnetic waves,” JCP, 1994
- J. P. Berenger and CPML variants (Roden & Gedney), IEEE TAP, 2000




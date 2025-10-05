# Droplet Rebounds off a Fluid Bath: Kinematic‑Match Simulations and Experiments (MATLAB)

[![arXiv](https://img.shields.io/badge/View%20on%20arXiv-B31B1B?logo=arxiv&labelColor=gray)](https://arxiv.org/abs/2509.22826)




![Graphical Abstract](matlab/0_data/manual/GraphicalAbstract.jpg)

This repository accompanies the manuscript “Droplet rebounds off a fluid bath: kinematic match simulations and experiments.” It contains MATLAB code and precomputed data to simulate non‑coalescing droplet–bath impacts using the full kinematic‑match (KM) framework, now extended to include droplet deformation. The solver predicts the time‑evolving contact area, pressure distribution, and wave field on both bodies. Automation scripts generate parameterized folder trees, run sweeps (water/oil presets), and aid post‑processing and figure reproduction.

All inputs use CGS units (cm–g–s). Folder names encode parameters to keep runs organized and reproducible.

## Paper Context
- Scope: Non‑coalescing, axisymmetric droplet impacts on a deep bath at low Weber number, where surface tension and weak viscosity dominate.
- Model: Full KM with spherical‑harmonic deformation of the droplet, weakly viscous bath, and natural geometric/kinematic contact constraints. No explicit gas‑film resolution is used.
- Results: Detailed predictions (pressed area, pressures, wave field) that compare favorably with prior and new experiments for sub‑millimetric droplets at low impact speed.

## What’s Included
- Core solver routines (SH dynamics, contact/pressure solver, DtN coupling).
- Sweep scripts for water/oil setups with logging and optional email notifications.
- Precomputed matrices/templates to avoid recomputing the DtN operator for common domains/resolution.
- Post-processing utilities and figure scripts.

## Repository Layout
- `matlab/1_code/simulation_code/`
  - Core functions, e.g. `solve_motion.m`, `solveDD0.m`, `solvenDDCusp.m`, `solve_ODE_unkown.m`, `r_from_spherical.m`, `zs_from_spherical.m`, `theta_from_cylindrical.m`, `custom_project_amplitudes.m`, `zeta_generator.m`.
- `matlab/1_code/sweeper_water_2024.m`, `matlab/1_code/sweeper_oil_2024.m`
  - Parameter-sweep drivers that prepare folder trees and call `solve_motion`.
- `matlab/1_code/D*/`
  - Domain templates and cached data (e.g., `DomainMaker.m`, `ParRadDTNStops.m`, `DTN*.mat`). These serve as “safe” templates for generating new runs.
- `matlab/1_code/Figures/`
  - Figure scripts and extracted data for plotting (e.g., `impact_panel.m`).
- `matlab/0_data/`
  - `manual/` for figures/logs, `external/` for reference data. Logs typically go to `manual/Logger/`.
- `Read_me.txt`
  - Legacy notes describing the older workflow (still useful context).

## Requirements
- MATLAB (R2020b+ recommended). No special toolboxes are strictly required; Parallel Computing Toolbox is optional.
- Python 3 optional for email notifications (used by `sending_email.py`). On macOS, you may need to fix SSL certs (see `certificate_install.py`).
- Disk space: runs produce `.mat` files per step; DtN matrices can be large.

## Quickstart: Run a Water or Oil Sweep
1) Start MATLAB and set the working directory:
   - `cd matlab/1_code`
2) Open and edit either `sweeper_water_2024.m` or `sweeper_oil_2024.m` to choose:
   - Domain and resolution: `D`, `Quant`
   - Fluid: `rho`, `sigma`, `nu`, `muair`
   - Drop (solid) props: `RhoS`, `SigmaS`
   - Radius: `R` (cm; folder names reflect mm)
   - Impact setup: `Ang` (deg), `U` (cm/s), `modes` (SH count), `tol` (convergence)
3) Run the sweeper script in MATLAB. It will:
   - Create the parameterized folder tree.
   - Prepare any missing `.mat` inputs using template makers (`DomainMaker.m`, `BathMaker.m`, `DropFluidMaker.m`, `RoMaker.m`).
   - Call `simulation_code/solve_motion.m` for each parameter set.
   - Write logs under `../0_data/manual/Logger/` and results in the leaf run folders.

Re-running a sweeper will skip runs that already have results; set `force_sweep = false/true` inside the script as needed.

## Quickstart: Single Run (Manual)
If you prefer not to use a sweeper, you can set up and run a single configuration manually.

1) Create the domain and base data (inside a chosen `DXXQuantYY/`):
   - Run `DomainMaker.m` (generates `D.mat`, `nr.mat`, `dr.mat`, `Delta.mat`, `IntMat.mat`, etc.).
   - A precomputed DtN matrix `DTNnew345nr%D%drefp10.mat` is typically present. Avoid recomputing `ParRadDTNStops.m` unless necessary.
2) Create fluid property folder and files (inside `rhoXXXXsigmaYYYYnuZZZZmuairW/`):
   - Edit and run `BathMaker.m` to write `rho.mat`, `sigma.mat`, `nu.mat`, `muair.mat`, `g.mat`.
3) Create solid (drop) property folder and files (inside `RhoSXXXXSigmaSYYYY/`):
   - Edit and run `DropFluidMaker.m` to write `rhoS.mat`, `sigmaS.mat`, etc.
4) Create the radius folder `R0####mm/` and run `RoMaker.m` to write `Ro.mat`.
5) Create the run folder and call the solver:
   - Make `ImpDefCornerAng{Ang}U{U0}/N={N}tol={tol}/` and `cd` into it.
   - Run, e.g.: `solve_motion(U0, [], N, tol, pwd, false)`

Tip: `solve_motion.m` navigates up the folder tree to load all required `.mat` files. If the folder naming convention is followed, it “just works”.

## Folder Convention & Units (CGS)
The sweeper/manual setup expects the following folder pattern under a given domain:
```
D{D}Quant{Quant}/
  rho{1000*rho}sigma{round(100*sigma)}nu{round(10000*nu)}muair{muair}/
    RhoS{1000*rhoS}SigmaS{round(100*sigmaS)}/
      R0{Ro*10000 in mm}/
        ImpDefCornerAng{Ang}U{U0}/
          N={N}tol={tol}/
```
- Examples (water): `rho1000sigma7220nu98muair0`, `RhoS1000SigmaS7220`, `R0350mm`, `ImpDefCornerAng180U44.52`, `N=21tol=5.00e-05`.
- Inputs in CGS: `rho`~g/cm^3 (enter 1 for water), `sigma`~dyn/cm (enter 72.20), `nu`~cm^2/s (enter 0.00978), `U0`~cm/s, `R` in cm.

## Outputs
Each run folder contains:
- Time histories: `z.mat`, `vz.mat`, `etaOri.mat`, `tvec.mat`, `numl.mat`, `nlmax.mat`.
- Fields and amplitudes: `etas.mat`, `psMatPer.mat`, `oscillation_amplitudes.mat`, `pressure_amplitudes.mat`, `Rv.mat`.
- Summary: `ProblemConditions.mat` (nondimensional numbers, units, N, U0, Ang, dtb, and `PROBLEM_CONSTANTS`).
- On errors: `error_logU0=... .mat` and partial snapshots prefixed with `errored_`.

Post-processing helpers (examples):
- `SurfceSlicesAndDroplet.m` (surface + drop visualization)
- `PressViewer.m` (pressure visualization)
- `Figures/impact_panel.m` (figure reproduction)

## Tips & Notes
- DtN matrices are expensive to compute. Use the included `DTN*.mat` in the domain folders; only run `ParRadDTNStops.m` if you understand the cost and requirements.
- Sweeper scripts also support logging and optional email notifications when runs finish. To enable email, place a Gmail app-password at `matlab/0_data/credentials/lol.txt` and ensure Python 3 + SSL certs are set up. See `matlab/1_code/sending_email.py` and `matlab/1_code/certificate_install.py`.
- Many scripts assume the folder tree exists. If something is missing, re-run the relevant `*Maker.m` script in the right folder.
- Default contact angle is 180° (`Ang = 180`).

## Legacy Notes
`Read_me.txt` contains the original description of the folder system and an older entry-point (`VertPolarExact*.m`). The modern workflow routes runs through `simulation_code/solve_motion.m` and the 2024 sweepers.

## Cite This Work
If you use this code or reproduce results, please cite the companion paper. A BibTeX entry will be added upon publication. In the meantime, cite as “Agüero, Galeano‑Rios, Ragazzo, Gabbard, Harris, Milewski (2025), Droplet rebounds off a fluid bath: kinematic match simulations and experiments.”

## Contact
Questions or issues with running the code?
- Open an issue on this repository

# Molecular Communication Channel as a Physical Reservoir Computer

This repository contains the MATLAB simulation code accompanying the study of **molecular communication (MC) channel as physical reservoir computer**. 
---
The tasks used to benchmark the reservoir are standard nonlinear time-series benchmarks: **NARMA-2**, **NARMA-10**, and **Mackey–Glass** prediction, as well as information-theoretic capacity via **Integrated Polynomial Capacity (IPC)**.
---
## Dependency: Smoldyn (Required)
> **All stochastic simulation code in this repository requires [Smoldyn](https://www.smoldyn.org/) to be installed and accessible from the MATLAB working directory.**
Smoldyn is a particle-based spatial stochastic simulator that models individual molecule diffusion and reactions in 3D. MATLAB scripts generate Smoldyn configuration files (`.txt`) on the fly, invoke the Smoldyn executable, and read back the resulting molecule-count outputs.
**Installation:** Download and install Smoldyn from [https://www.smoldyn.org/download.html](https://www.smoldyn.org/download.html). Make sure the `smoldyn` binary is either on your system `PATH` or that the MATLAB scripts are updated to point to its location.

---

## Repository Structure

```
├── Point-Source(Monte Carlo)/          # Validates the point-source approximation
├── NARMA2_code/                        # NARMA-2 benchmark experiments
├── Masking/         # Virtual node / masking analysis
├── IPC_SWEEP/                      #   IPC sweep
├── MG/                    #   Mackey–Glass prediction
└── NARMA10/               #   NARMA-10 prediction
```

## Directory Descriptions

### `Point-Source(Monte Carlo)/`

Validates that the analytical **point-source / point-receiver** model accurately represents the full stochastic Smoldyn simulation. Three scenarios are compared side-by-side:

| Scenario | Transmitter | Receiver |
|---|---|---|
| Deterministic analytical | Point | Point |
| Stochastic (Smoldyn) | Point | Finite sphere |
| Stochastic (Smoldyn) | Finite sphere | Finite sphere |

Key files:

- [MonteCarloConfirmation.m](Point-Source\(Monte%20Carlo\)/MonteCarloConfirmation.m) — Top-level script; runs all three scenarios over 20 Monte Carlo trials and plots mean ± std traces.
- [deterministicPointModel.m](Point-Source\(Monte%20Carlo\)/deterministicPointModel.m) — Solves the 3D diffusion Green's function and the receptor binding ODE analytically.
- [generateSimConfig_MCConfirmation.m](Point-Source\(Monte%20Carlo\)/generateSimConfig_MCConfirmation.m) — Generates Smoldyn `.txt` config files for point or finite-radius molecule release.
- [runOneCase.m](Point-Source\(Monte%20Carlo\)/runOneCase.m) — Orchestrates one full simulation case (train + test), averaging over multiple Smoldyn runs.
- [computeTraceMetrics.m](Point-Source\(Monte%20Carlo\)/computeTraceMetrics.m) — Computes correlation, NRMSE, and peak-time error between analytical and stochastic traces.
- [randomPointsInSphere.m](Point-Source\(Monte%20Carlo\)/randomPointsInSphere.m) — Uniformly samples points inside a sphere (used for finite transmitter initialisation).

**Default physical parameters:** `D = 1×10⁻¹¹ m²/s`, `distance = 10 µm`, `N = 500` receptors, `k_on = 1×10⁻¹⁸ m³/s`, `k_off = 1 s⁻¹`, symbol period `T = 1 s`.

---

### `NARMA2_code/`

Applies the MC reservoir to the **NARMA-2** benchmark — a second-order nonlinear autoregressive moving-average task defined by:

```
q(n+1) = 0.4·q(n) + 0.4·q(n)·q(n-1) + 0.6·u(n)³ + 0.1
```

Key files:

- [NARMA2_Generate_Smoldyn_Data.m](NARMA2_code/NARMA2_Generate_Smoldyn_Data.m) — Converts the NARMA-2 input sequence to molecule counts, runs Smoldyn, and saves the resulting receptor occupancy traces.
- [NARMA2_Stochastic_Sweep.m](NARMA2_code/NARMA2_Stochastic_Sweep.m) — Parametric sweep over moving-average window size and virtual-node count; generates 2D NRMSE heatmaps comparing the deterministic ODE model to the Smoldyn stochastic model.
- [NARMA2_fig_generation_code.m](NARMA2_code/NARMA2_fig_generation_code.m) — Produces the publication-quality prediction figure (target vs. stochastic estimate over 300 test symbols).

---

### `Masking/`

Characterises the **information richness** of the virtual-node representation produced by temporal masking.

- [MaskingValidation.m](Masking/MaskingValidation.m) — Loads a saved receptor occupancy trace (100 virtual nodes per symbol) and computes:
  - Pairwise correlation matrix between nodes
  - Average off-diagonal correlation
  - Participation ratio (effective rank of the node covariance)
  - Effective rank via Shannon entropy of the eigenvalue spectrum
  - Eigenvalue spectrum of the correlation matrix

These metrics quantify the effective dimensionality of the reservoir state space and confirm that temporal masking produces sufficiently diverse, low-correlated features.

---

### `IPC_SWEEP/` — Information Processing Capacity

Measures the **Integrated Polynomial Capacity (IPC)** of the MC reservoir using Legendre polynomial basis functions. IPC decomposes total capacity into contributions from different orders of nonlinear computation.

- [IPC_sweep.m](Final%20Code%20Repo%20for%20Original%20Paper/IPC_SWEEP/IPC_sweep.m) — Sweeps `k_on` × `k_off` parameter space and produces IPC heatmaps plus mean receptor occupancy.
- [computeIPC_legendreNoZeroDup.m](Final%20Code%20Repo%20for%20Original%20Paper/IPC_SWEEP/computeIPC_legendreNoZeroDup.m) — Core IPC calculation: enumerates factor sets of Legendre polynomials up to a target degree and measures the squared correlation of each from a linear readout.
- [computeIPC_legendreNoZeroDup_sweepPrint.m](Final%20Code%20Repo%20for%20Original%20Paper/IPC_SWEEP/computeIPC_legendreNoZeroDup_sweepPrint.m) — Wrapper that reports swept parameters during a batch run.

**Parameter ranges swept:** `k_on ∈ [5×10⁻²⁰, 2×10⁻¹⁷] m³/s`, `k_off ∈ [0.1, 10] s⁻¹`.

### `MG/` — Mackey–Glass Time-Series Prediction

Benchmarks the reservoir on **Mackey–Glass** chaotic time-series prediction, a standard hard reservoir computing task:

```
dx/dt = β·x(t−τ) / (1 + x(t−τ)ⁿ) − γ·x(t)     (τ=17, β=0.2, γ=0.1, n=10)
```

- [createMGseries.m](Final%20Code%20Repo%20for%20Original%20Paper/MG/createMGseries.m) — Integrates the Mackey–Glass DDE with RK4 and saves a 5000-point normalised series.
- [NewMG_Generate_Smoldyn_Data.m](Final%20Code%20Repo%20for%20Original%20Paper/MG/NewMG_Generate_Smoldyn_Data.m) — Maps the MG series to molecule counts and runs the Smoldyn stochastic simulation.
- [NewMG_Stochastic_Sweep.m](Final%20Code%20Repo%20for%20Original%20Paper/MG/NewMG_Stochastic_Sweep.m) — 2D sweep over moving-average window and prediction horizon (`predictlength`); outputs NRMSE heatmaps.
- [NewMG_Only_Numerical_Sweep.m](Final%20Code%20Repo%20for%20Original%20Paper/MG/NewMG_Only_Numerical_Sweep.m) — Deterministic (ODE-only) version of the sweep for comparison.

### `NARMA10/` — NARMA-10 Time-Series Prediction

The most memory-demanding benchmark: **NARMA-10**, a 10th-order nonlinear task:

```
q(n+1) = 0.3·q(n) + 0.05·q(n)·Σq(n−9:n) + 1.5·u(n−9)·u(n) + 0.1
```

- [NARMA10_Generate_Smoldyn_Data.m](Final%20Code%20Repo%20for%20Original%20Paper/NARMA10/NARMA10_Generate_Smoldyn_Data.m) — Generates Smoldyn data for the NARMA-10 input sequence.
- [NARMA_Stochastic_Sweep.m](Final%20Code%20Repo%20for%20Original%20Paper/NARMA10/NARMA_Stochastic_Sweep.m) — Full stochastic parameter sweep with Smoldyn comparison.
- [NARMA10_Only_Numerical_Sweep.m](Final%20Code%20Repo%20for%20Original%20Paper/NARMA10/NARMA10_Only_Numerical_Sweep.m) — Numerical-only sweep that additionally computes and plots the theoretical receptor timescale `τ = 1/(k_on·⟨c⟩ + k_off)` and its autocorrelation alongside NRMSE heatmaps.

---

## Requirements

| Requirement | Notes |
|---|---|
| **MATLAB** | R2019b or later recommended |
| **Smoldyn** | **Required for all stochastic simulations** — [smoldyn.org](https://www.smoldyn.org/) |
| MATLAB Signal Processing Toolbox | Used for `movmean`, `xcorr` |
| MATLAB Statistics & Machine Learning Toolbox | Used for ridge regression utilities |

---

## Citation

If you use this code, please cite the associated paper (details to be added upon publication).

---

## License

See [LICENSE](LICENSE) for details.

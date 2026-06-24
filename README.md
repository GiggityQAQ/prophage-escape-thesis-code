# Prophage Escape Thesis Code

This repository contains the Julia source code used for the thesis:

**A Modeling Study of Prophage Escape under Periodic Antibiotic Pressure**

Author: **Zichang Li**
Repository: `https://github.com/GiggityQAQ/prophage-escape-thesis-code`

## Overview

This project implements deterministic ordinary differential equation models of bacteria, lysogens, and free phages under periodic antibiotic pressure.

The thesis compares two model variants:

* **Model A: Baseline model**
  All prophages induce at the same basal rate. Antibiotic exposure affects sensitive cells through antibiotic-mediated death, but does not increase prophage induction.

* **Model B: Escape model**
  ARG-negative prophages increase their induction rate during antibiotic exposure, while ARG-positive prophages keep the basal induction rate.

For each antibiotic schedule, four initial-condition tests are used:

1. `Sresident`
2. `Cresident`
3. `President`
4. `DirectPC`

Resident-based tests first run the resident community until convergence. After invaders are introduced, post-invasion dynamics are simulated for a fixed window of 200 complete antibiotic cycles. Final outcomes are classified using the final five complete cycles.

## Repository structure

```text
prophage-escape-thesis-code/
├── README.md
├── Project.toml
├── Manifest.toml
├── scripts/
│   ├── modelA_phase_diagram.jl
│   ├── modelB_phase_diagram.jl
│   ├── modelB_timeseries.jl
│   └── modelB_dominant_lysogen_check.jl
├── figures/
│   ├── 3.1.pdf
│   ├── 3.2.pdf
│   ├── 3.3.pdf
│   ├── 3.4.pdf
│   └── 3.5.pdf
└── results/
```

`Project.toml` and `Manifest.toml` are recommended for reproducibility. If they are not included, the required packages can be installed manually as described below.

## Scripts and thesis figures

| Thesis figure | Figure file       | Source script                              | Main generated output                                   | Description                                                                                                                                                                     |
| ------------- | ----------------- | ------------------------------------------ | ------------------------------------------------------- | ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| Figure 3.1    | `figures/3.1.pdf` | `scripts/modelA_phase_diagram.jl`          | `modelA_phase_SPC_dominant_square.png`                  | Generates the Model A baseline phase diagram. This figure tests whether ARG location alone produces a consistent long-term outcome in the absence of antibiotic-induced escape. |
| Figure 3.2    | `figures/3.2.pdf` | `scripts/modelB_phase_diagram.jl`          | `modelB_phase_lysogen_presence_by_test_square.png`      | Generates the Model B lysogen persistence-combination phase diagram. This figure records which lysogen types remain present in the final analysis window.                       |
| Figure 3.3    | `figures/3.3.pdf` | `scripts/modelB_timeseries.jl`             | `modelB_timeseries_Ton_600_Toff_100.png`                | Generates the representative Model B time-series simulation for `T_on = 600` and `T_off = 100` under the `DirectPC` initial condition.                                          |
| Figure 3.4    | `figures/3.4.pdf` | `scripts/modelB_timeseries.jl`             | `modelB_timeseries_Ton_600_Toff_900.png`                | Generates the representative Model B time-series simulation for `T_on = 600` and `T_off = 900` under the `DirectPC` initial condition.                                          |
| Figure 3.5    | `figures/3.5.pdf` | `scripts/modelB_dominant_lysogen_check.jl` | `modelB_phase_dominant_lysogen_type_by_test_square.png` | Generates the simplified Model B dominant-lysogen global check phase diagram. This figure records only the single dominant lysogen type for each run.                           |

## Requirements

The code was written in Julia.

Main Julia packages:

* `DifferentialEquations.jl`
* `Plots.jl`
* `Colors.jl`

The scripts also use Julia standard libraries such as:

* `Statistics`
* `Printf`

## Installation

Clone the repository:

```bash
git clone https://github.com/GiggityQAQ/prophage-escape-thesis-code.git
cd prophage-escape-thesis-code
```

Start Julia using the project environment:

```bash
julia --project=.
```

If `Project.toml` and `Manifest.toml` are included, instantiate the environment inside Julia:

```julia
using Pkg
Pkg.instantiate()
```

If the environment has not been created yet, install the required packages manually:

```julia
using Pkg
Pkg.add(["DifferentialEquations", "Plots", "Colors"])
```

## Running the scripts

All commands should be run from the root directory of the repository.

Run the Model A baseline phase diagram:

```bash
julia --project=. scripts/modelA_phase_diagram.jl
```

Run the Model B lysogen persistence-combination phase diagram:

```bash
julia --project=. scripts/modelB_phase_diagram.jl
```

Run the Model B representative time-series simulations:

```bash
julia --project=. scripts/modelB_timeseries.jl
```

Run the Model B dominant-lysogen global check:

```bash
julia --project=. scripts/modelB_dominant_lysogen_check.jl
```

## Expected outputs

The scripts generate the following main output files:

| Output file                                             | Generated by                               | Thesis figure |
| ------------------------------------------------------- | ------------------------------------------ | ------------- |
| `modelA_phase_SPC_dominant_square.png`                  | `scripts/modelA_phase_diagram.jl`          | Figure 3.1    |
| `modelB_phase_lysogen_presence_by_test_square.png`      | `scripts/modelB_phase_diagram.jl`          | Figure 3.2    |
| `modelB_timeseries_Ton_600_Toff_100.png`                | `scripts/modelB_timeseries.jl`             | Figure 3.3    |
| `modelB_timeseries_Ton_600_Toff_900.png`                | `scripts/modelB_timeseries.jl`             | Figure 3.4    |
| `modelB_phase_dominant_lysogen_type_by_test_square.png` | `scripts/modelB_dominant_lysogen_check.jl` | Figure 3.5    |

By default, the scripts save generated figures to the working directory from which the commands are run. The final thesis-ready versions of the figures can be stored in the `figures/` folder.

## Parameter scan settings

The phase-diagram scripts use the following default parameter grid:

```julia
T_on_values  = 10.0:80.0:1200.0
T_off_values = 10.0:80.0:1200.0
```

This scans antibiotic-on duration and antibiotic-off duration from 10 to 1200 with a step size of 80.

Some scripts include commented settings for higher-resolution scans, for example:

```julia
T_on_values  = 10.0:10.0:1200.0
T_off_values = 10.0:10.0:1200.0
```

The higher-resolution setting is computationally more expensive and is not used as the default.

## Model notes

The state variables are:

* `S_S`: non-lysogenic ARG-negative host
* `S_R`: non-lysogenic ARG-positive host
* `L_SS`: ARG-negative host carrying an ARG-negative prophage
* `L_SR`: ARG-negative host carrying an ARG-positive prophage
* `L_RS`: ARG-positive host carrying an ARG-negative prophage
* `L_RR`: ARG-positive host carrying an ARG-positive prophage
* `V_S`: ARG-negative free phage
* `V_R`: ARG-positive free phage
* `A`: antibiotic state, where `A = 0` indicates antibiotic absence and `A = 1` indicates antibiotic presence

In Model A, all prophages have the same basal induction rate.

In Model B, antibiotic-induced extra induction applies only to ARG-negative prophages:

```julia
eta_S = eta + epsilon * A
eta_R = eta
```

## Thesis

Zichang Li. *A Modeling Study of Prophage Escape under Periodic Antibiotic Pressure*. University of Auckland, 2026.

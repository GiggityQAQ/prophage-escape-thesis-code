# Prophage Escape Thesis Code

This repository contains the Julia source code used for the thesis:

**A Modeling Study of Prophage Escape under Periodic Antibiotic Pressure**

Author: **Zichang Li**
Repository: `https://github.com/GiggityQAQ/prophage-escape-thesis-code`

## Overview

This project implements deterministic ordinary differential equation (ODE)
models of temperate phage and bacterial host population dynamics under periodic
antibiotic exposure.

The thesis compares two model variants:

* **Model A: Baseline model**
  All prophages have the same basal induction rate. Antibiotic exposure affects
  sensitive cells through antibiotic-mediated death, but does not increase
  prophage induction.

* **Model B: Escape model**
  ARG-negative prophages increase their induction rate during antibiotic
  exposure, while ARG-positive prophages keep the basal induction rate.

The scripts reproduce the five main computational figures reported in Chapter 3
of the thesis.

## Repository structure

```text
prophage-escape-thesis-code/
├── README.md
├── run_all.sh
└── ODE/
    ├── modelA_phase_diagram.jl
    ├── modelB_phase_diagram.jl
    ├── modelB_timeseries.jl
    └── modelB_dominant_lysogen_global_check.jl
```

## Scripts and thesis figures

| Thesis figure | Script                                        | Main output                                             | Description                                                                                                               |
| ------------- | --------------------------------------------- | ------------------------------------------------------- | ------------------------------------------------------------------------------------------------------------------------- |
| Figure 3.1    | `ODE/modelA_phase_diagram.jl`                 | `modelA_phase_SPC_dominant_square.png`                  | Generates the Model A baseline phase diagram using the SPC dominant-attractor classification.                             |
| Figure 3.2    | `ODE/modelB_phase_diagram.jl`                 | `modelB_phase_lysogen_presence_by_test_square.png`      | Generates the Model B lysogen-presence phase diagram, recording which lysogen types persist in the final analysis window. |
| Figure 3.3    | `ODE/modelB_timeseries.jl`                    | `modelB_timeseries_Ton_600_Toff_100.png`                | Generates the Model B representative time-series simulation for `T_on = 600` and `T_off = 100`.                           |
| Figure 3.4    | `ODE/modelB_timeseries.jl`                    | `modelB_timeseries_Ton_600_Toff_900.png`                | Generates the Model B representative time-series simulation for `T_on = 600` and `T_off = 900`.                           |
| Figure 3.5    | `ODE/modelB_dominant_lysogen_global_check.jl` | `modelB_phase_dominant_lysogen_type_by_test_square.png` | Generates the simplified Model B dominant-lysogen global check phase diagram.                                             |

## Requirements

The code was written in Julia. Julia 1.9 or newer is recommended. You can
download it from <https://julialang.org/downloads/> and confirm the
installation by running `julia --version` in a terminal.

Main Julia packages:

* `DifferentialEquations.jl`
* `Plots.jl`
* `Colors.jl`
* `Measures.jl`

The scripts also use Julia standard libraries that ship with Julia and do not
need to be installed separately:

* `Statistics`
* `Printf`

## Installation

### 1. Clone the repository

Run the following in a terminal (not inside Julia):

```bash
git clone https://github.com/GiggityQAQ/prophage-escape-thesis-code.git
cd prophage-escape-thesis-code
```

### 2. Install the required packages

The commands below are **Julia commands, not terminal commands**. They must be
typed inside the Julia REPL, not directly into your shell.

First start the Julia REPL by typing `julia` in a terminal and pressing Enter:

```bash
julia
```

You should now see the Julia prompt, which looks like this:

```text
julia>
```

At this `julia>` prompt, run the following two lines:

```julia
using Pkg
Pkg.add(["DifferentialEquations", "Plots", "Colors", "Measures"])
```

Package installation may take several minutes the first time. After it
finishes, leave the REPL with `exit()` (or press `Ctrl+D`) to return to your
terminal:

```julia
exit()
```

You only need to do this installation step once per machine.

## Running the scripts

All commands in this section are terminal commands and should be run from the
root directory of the repository.

Run the Model A baseline phase diagram (Figure 3.1):

```bash
julia ODE/modelA_phase_diagram.jl
```

Run the Model B lysogen-presence phase diagram (Figure 3.2):

```bash
julia ODE/modelB_phase_diagram.jl
```

Run the Model B representative time-series simulations (Figures 3.3 and 3.4):

```bash
julia ODE/modelB_timeseries.jl
```

Run the Model B dominant-lysogen global check (Figure 3.5):

```bash
julia ODE/modelB_dominant_lysogen_global_check.jl
```

### Running everything at once

A helper script is provided to run all four scripts in sequence and generate
all five figures. From the repository root:

```bash
chmod +x run_all.sh
./run_all.sh
```

## Expected outputs

The scripts generate the following main figure files in the directory from
which they are run:

| Output file                                             | Generated by                                  | Thesis figure |
| ------------------------------------------------------- | --------------------------------------------- | ------------- |
| `modelA_phase_SPC_dominant_square.png`                  | `ODE/modelA_phase_diagram.jl`                 | Figure 3.1    |
| `modelB_phase_lysogen_presence_by_test_square.png`      | `ODE/modelB_phase_diagram.jl`                 | Figure 3.2    |
| `modelB_timeseries_Ton_600_Toff_100.png`                | `ODE/modelB_timeseries.jl`                    | Figure 3.3    |
| `modelB_timeseries_Ton_600_Toff_900.png`                | `ODE/modelB_timeseries.jl`                    | Figure 3.4    |
| `modelB_phase_dominant_lysogen_type_by_test_square.png` | `ODE/modelB_dominant_lysogen_global_check.jl` | Figure 3.5    |

### Approximate runtime

* The time-series script (`modelB_timeseries.jl`) is fast and usually finishes
  within a few minutes.
* The three phase-diagram scripts scan a two-dimensional parameter grid and are
  much slower. At the default grid resolution they may take from tens of
  minutes up to a few hours, depending on your machine.

## Parameter scan settings

The phase-diagram scripts use the following default parameter grid:

```julia
T_on_values  = 10.0:80.0:1200.0
T_off_values = 10.0:80.0:1200.0
```

This scans the antibiotic-on duration and the antibiotic-off duration from 10
to 1200 with a step size of 80.

Some scripts include commented settings for a higher-resolution scan:

```julia
T_on_values  = 10.0:10.0:1200.0
T_off_values = 10.0:10.0:1200.0
```

The higher-resolution scan is computationally more expensive and is not used as
the default setting.

## Model notes

The ODE state variables are:

* `S_S`: non-lysogenic ARG-negative host
* `S_R`: non-lysogenic ARG-positive host
* `L_SS`: ARG-negative host carrying an ARG-negative prophage
* `L_SR`: ARG-negative host carrying an ARG-positive prophage
* `L_RS`: ARG-positive host carrying an ARG-negative prophage
* `L_RR`: ARG-positive host carrying an ARG-positive prophage
* `V_S`: ARG-negative free phage
* `V_R`: ARG-positive free phage
* `A`: antibiotic state, where `A = 0` indicates antibiotic absence and
  `A = 1` indicates antibiotic presence

The models assume logistic growth with a carrying capacity normalized to 1, so
all population sizes lie between 0 and 1.

In Model A, all prophages have the same basal induction rate. In Model B,
antibiotic-induced extra induction applies only to ARG-negative prophages.

## Troubleshooting

* **`command not found: julia`** — Julia is not installed or not on your `PATH`.
  Install it from <https://julialang.org/downloads/> and reopen your terminal.
* **`UndefVarError` for a package function (e.g. `mean`)** — a required package
  was not installed. Re-run the installation step above inside the Julia REPL.
* **`Pkg.add` fails** — make sure you are running it at the `julia>` prompt, not
  in your terminal shell.

## Author

Zichang Li

## Thesis

Zichang Li. *A Modeling Study of Prophage Escape under Periodic Antibiotic
Pressure*. University of Auckland, 2026.

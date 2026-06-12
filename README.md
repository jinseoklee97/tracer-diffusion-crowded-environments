# Accessible pore geometry governs tracer diffusion in crowded environments

This repository contains analysis code, processed datasets, and simulation scripts associated with the manuscript:

**Accessible pore geometry governs tracer diffusion in crowded environments**
Jinseok Lee, Tong Lin, Mengyang Gu, and Yimin Luo

The repository provides workflows for analyzing tracer diffusion in crowded soft-particle and intracellular environments. The main analyses include multiple particle tracking (MPT), AIUQ-based MSD estimation, minimal Brownian-dynamics simulations, and PPGP-based prediction and inverse analysis of tracer mean-squared displacements (MSDs).

## Overview

Tracer diffusion in crowded environments is governed by the accessible pore space available to the tracer. This repository includes codes used to:

1. Analyze experimental tracer videos using MPT and AIUQ.
2. Simulate tracer diffusion in polydisperse soft-particle matrices.
3. Train and evaluate PPGP models that map matrix geometric variables to MSD curves.
4. Perform inverse analysis to identify matrix geometries that reproduce target MSDs.
5. Estimate tracer-defined percolation behavior from predicted MSDs.
6. Analyze pore-size distributions, dynamic pore statistics, van Hove functions, non-Gaussian parameters, and velocity autocorrelation functions.

The main geometric variables used in the simulation and PPGP analyses are:

* Matrix area fraction, `Phi`
* Mean matrix particle radius, `R_bar`
* Matrix particle polydispersity, `p`
* Matrix mobility factor, when dynamic matrix simulations are used

---

## Repository structure

```text
.
├── Experiment
│   ├── AIUQ
│   │   ├── Functions
│   │   │   ├── functions_DDM.R
│   │   │   └── functions_nonparametric.R
│   │   ├── Results
│   │   │   └── .gitkeep
│   │   ├── Videos
│   │   │   └── .gitkeep
│   │   └── AIUQ_nonparametric.R
│   │
│   ├── Example videos
│   │   └── .gitkeep
│   │
│   └── MPT
│       ├── Results
│       │   └── .gitkeep
│       ├── Videos
│       │   └── .gitkeep
│       ├── MPT_Main.m
│       ├── MPT_perform.m
│       ├── bpass.m
│       ├── feature2D.m
│       ├── fracshift.m
│       ├── localmax.m
│       ├── luberize.m
│       ├── rsqd.m
│       ├── thetarr.m
│       ├── trackmem.m
│       └── unq.m
│
├── PPGP
│   ├── Code
│   │   ├── functions
│   │   ├── Inverse_total.R
│   │   ├── forward_3dscatter.R
│   │   └── forward_percolation.R
│   │
│   ├── Data
│   │   ├── Inplane
│   │   │   ├── inplane_testing_MSDs.csv
│   │   │   └── inplane_testing_parameters.csv
│   │   ├── Input
│   │   │   └── 0.72_sim.csv
│   │   ├── Testing
│   │   │   ├── testing_MSD.csv
│   │   │   └── testing_input.csv
│   │   └── Training
│   │       ├── delta_t.csv
│   │       ├── training_MSD.csv
│   │       └── training_input.csv
│   │
│   └── Output
│       └── .gitkeep
│
├── Simulation
│   └── Code
│       ├── Simulation_setup.m
│       └── run_simulation.m
│
└── README.md
```

---

## Experiment

The `Experiment/` folder contains codes for analyzing experimental microscopy videos.

---

### MPT analysis

Folder:

```text
Experiment/MPT/
```

This folder contains the MATLAB implementation of multiple particle tracking.

Place TIFF videos in:

```text
Experiment/MPT/Videos/
```

Then edit the following user settings in `MPT_Main.m`:

```matlab
filename = 'SP 0.72 AF_2000frame';
maxdisp = 10;
nt = 2000;
dt = 0.5;
pxsz = 0.293;
```

Run in MATLAB:

```matlab
cd Experiment/MPT
MPT_Main
```

Main files:

* `MPT_Main.m`: main runner for MPT, MSD, NGP, and VACF analysis
* `MPT_perform.m`: performs image filtering, feature detection, trajectory linking, MSD calculation, van Hove displacement collection, NGP calculation, and VACF calculation
* `bpass.m`, `feature2D.m`, `trackmem.m`, and related files: particle detection and tracking functions

Default tracking parameters inside `MPT_perform.m` include:

```matlab
memory = 2;
Imin = 150;
rad = 3;
```

Outputs are saved in:

```text
Experiment/MPT/Results/
```

Typical output files include:

```text
<filename>_MSD.mat
<filename>_stDev.mat
<filename>_tau.mat
<filename>_lub.mat
<filename>_NGP.mat
<filename>_MSD_NGP_noEB_runner.csv
<filename>_VACF.mat
<filename>_VACF_runner.csv
<filename>_MPT_MSD_NGP_VACF_noEB_all.mat
<filename>_Trajectories.tiff
```

The CSV output `<filename>_MSD_NGP_noEB_runner.csv` contains MSD and non-Gaussian parameter values. The CSV output `<filename>_VACF_runner.csv` contains the velocity autocorrelation function.

---

### AIUQ analysis

Folder:

```text
Experiment/AIUQ/
```

This folder contains the R implementation of nonparametric AIUQ-based MSD estimation.

Place TIFF videos in:

```text
Experiment/AIUQ/Videos/
```

Then edit the user settings in `AIUQ_nonparametric.R`:

```r
file_name = "SP 0.72 AF_500frame"
mindt = 2
pxsz = 0.293
M = 100
q_thr = 0.999
```

Run in R:

```r
setwd("Experiment/AIUQ")
source("AIUQ_nonparametric.R")
```

or from the terminal:

```bash
cd Experiment/AIUQ
Rscript AIUQ_nonparametric.R
```

Main files:

* `AIUQ_nonparametric.R`: main script for nonparametric AIUQ-based MSD estimation
* `Functions/functions_DDM.R`: helper functions for DDM/AIUQ analysis
* `Functions/functions_nonparametric.R`: helper functions for nonparametric MSD estimation

The script reads:

```text
Videos/<file_name>.tif
```

and writes:

```text
results/<file_name>_AIUQ_nonparametric.csv
```

The output CSV contains:

```text
d_input
n_par_AIUQ
lower
upper
```

where `d_input` is the lag time, `n_par_AIUQ` is the estimated MSD, and `lower`/`upper` are uncertainty bounds.

---

## Simulation

Folder:

```text
Simulation/Code/
```

This folder contains the minimal simulation model used to generate tracer trajectories and MSDs in crowded soft-particle environments.

Main files:

```text
Simulation_setup.m
run_simulation.m
```

Run in MATLAB:

```matlab
cd Simulation/Code
Simulation_setup
```

`Simulation_setup.m` is the main runner. It defines the parameter sweep, output options, and base simulation settings, then calls `run_simulation.m`.

---

### Main simulation settings

The main sweep variables are defined near the top of `Simulation_setup.m`:

```matlab
phi_list = 0.717;
R_list   = 4.56;
p_list   = 0.286;
mobility_list = [0];
nRep     = 1;
```

These correspond to:

* `phi_list`: matrix area fraction
* `R_list`: mean matrix particle radius
* `p_list`: matrix particle polydispersity
* `mobility_list`: matrix mobility factor
* `nRep`: number of repeated simulations per condition

The base simulation settings include:

```matlab
cfg.L           = 150;
cfg.R_tr        = 1.00;
cfg.eta         = 1e-3;
cfg.dt_tr       = 0.5;
cfg.N_tr_target = 200;
cfg.N_steps     = 2000;
```

Hydrodynamic and vertical fluctuation parameters include:

```matlab
cfg.h0        = 1.25 * cfg.R_tr;
cfg.tau_h     = 2.0;
cfg.sigma_h   = 0.04 * cfg.R_tr;
cfg.blend_wall = 0.4;
cfg.alpha_hyd  = 0.6;
cfg.g0         = 0.30 * cfg.R_tr;
cfg.f_min      = 0.35;
```

---

### Simulation outputs

The main output folder is defined by:

```matlab
sweep_outdir = fullfile(pwd, '0.72');
```

Each simulation case is saved in:

```text
Simulation/Code/0.72/cases/<case_name>/
```

Typical per-case outputs include:

```text
info.csv
radii.csv
radii.mat
seed_holes.csv
emsd.csv
result.mat
traj.png
traj.fig
msd.png
msd.fig
sim.mp4
packing.mp4
```

Depending on the output switches, the simulation can also export:

```text
alpha2.csv
EB.csv
vacf.csv
vh.png
vh.fig
vh_1d_<lag>.csv
vh_radial_<lag>.csv
mat_msd.csv
mat_msd.png
mat_msd.fig
```

Pore-geometry outputs may include:

```text
preview_pass_size_map_<case_name>.png
FIG_pores_metric_1panel_INSET_<case_name>.png
FIG_pores_metric_1panel_NOINSET_<case_name>.png
holes_geometry_<case_name>.csv
holes_passable_<case_name>.csv
```

Dynamic pore outputs are saved under:

```text
pore_dyn/
```

and may include:

```text
summary.csv
exact.csv
pdf.csv
counts.csv
pore_data.mat
stats.png
area.png
count.png
pdf_heat.png
```

---

### Sweep-level outputs

After all cases are complete, `Simulation_setup.m` saves sweep-level files in the main sweep folder, including:

```text
summary.csv
all.mat
sweep_tau.csv
sweep_input.csv
sweep_MSD.csv
sweep_metrics.csv
sweep_input_labeled.csv
sweep_MSD_labeled.csv
sweep_metrics_labeled.csv
sweep_MSD_mean.csv
sweep_input_mean.csv
sweep_metrics_mean.csv
```

When experimental MSD comparison is enabled with:

```matlab
use_exp_msd = true;
```

the script can also generate:

```text
best_fit.png
best_fit.fig
rmse_heat_m*.png
rmse_heat_m*.fig
```

---

## PPGP

Folder:

```text
PPGP/
```

This folder contains the R code and data used for PPGP-based prediction, inverse analysis, percolation analysis, and visualization.

---

### Data files

Training data:

```text
PPGP/Data/Training/training_input.csv
PPGP/Data/Training/training_MSD.csv
PPGP/Data/Training/delta_t.csv
```

Testing data:

```text
PPGP/Data/Testing/testing_input.csv
PPGP/Data/Testing/testing_MSD.csv
```

In-plane testing data:

```text
PPGP/Data/Inplane/inplane_testing_parameters.csv
PPGP/Data/Inplane/inplane_testing_MSDs.csv
```

Target MSD for inverse analysis:

```text
PPGP/Data/Input/0.72_sim.csv
```

The input columns are assumed to be ordered as:

```text
column 1 = Phi
column 2 = p
column 3 = R_bar
```

The target MSD file for inverse analysis should contain:

```text
dt
msd
```

Optional uncertainty columns can also be included:

```text
lower
upper
```

---

## PPGP scripts

### `forward_3dscatter.R`

This script trains a PPGP model using simulation-generated MSDs and visualizes the predicted state space.

Run from the repository root:

```bash
Rscript PPGP/Code/forward_3dscatter.R
```

or run from R:

```r
setwd("PPGP/Code")
source("forward_3dscatter.R")
```

The script automatically detects paths when it is run from `PPGP/Code/`. Paths can also be overridden:

```bash
Rscript PPGP/Code/forward_3dscatter.R --root="PPGP" --data="PPGP/Data" --out="PPGP/Output"
```

This script:

1. Loads training, testing, and in-plane testing data.
2. Fits a PPGP model using `RobustGaSP::ppgasp`.
3. Trains on `log10(MSD + eps)`.
4. Computes the logarithmic MSD slope, `alpha`.
5. Generates a 3D training-data scatter plot.
6. Generates a 2D `Phi-p` slice at fixed `R_bar`.
7. Generates a 3D state-space surface visualization.
8. Validates predicted `alpha` values against held-out testing data.
9. Exports predicted MSDs and 95% predictive intervals for four user-selected conditions.

Important user settings include:

```r
R_slice = 4
n_grid_2d = 90
pred_cond_mat = rbind(
  c(0.10, 0.398, 3.31),
  c(0.31, 0.373, 3.86),
  c(0.549, 0.323, 3.73),
  c(0.717, 0.286, 4.56)
)
```

Typical outputs in `PPGP/Output/` include:

```text
3D_alpha_training.png
2D_alpha_slice_phi_p_R4.000.png
3D_state_space.png
validation_alphaPred_vs_alphaTrue_TEST.png
validation_alpha_true_vs_pred_TEST.csv
PredMSD_4conditions_with95PI.png
PredMSD_4conditions_with95PI_long.csv
PredMSD_4conditions_with95PI_wide.csv
PredMSD_<condition>_phi<Phi>_p<p>_R<R>.csv
```

---

### `forward_percolation.R`

This script uses the trained PPGP model to perform a dense `Phi` sweep at fixed `p` and `R_bar`, then identifies the tracer-defined percolation threshold using a knee-detection procedure.

Run:

```bash
Rscript PPGP/Code/forward_percolation.R
```

or in R:

```r
setwd("PPGP/Code")
source("forward_percolation.R")
```

Path override example:

```bash
Rscript PPGP/Code/forward_percolation.R --root="PPGP" --data="PPGP/Data" --out="PPGP/Output"
```

Important user settings include:

```r
p_fixed = 0.5
R_fixed = 3.76

N_PHI_DENSE = 400
EARLY_MAX = 10
LATE_MIN = 500

USE_LOGBIN_LATE = TRUE
N_BINS_LATE = 18
```

This script:

1. Loads simulation training data.
2. Fits a PPGP model on `log10(MSD + eps)`.
3. Predicts MSDs over a dense `Phi` grid at fixed `p` and `R_bar`.
4. Estimates early-time and late-time diffusion coefficients.
5. Computes the ratio `D_late / D_early`.
6. Uses a kneedle-style algorithm to identify the knee position.
7. Saves CSV summaries and diagnostic plots.

Typical outputs in `PPGP/Output/` include:

```text
phi_sweep_Dearly_Dlate_ratio_p<p>_R<R>.csv
phi_knee_summary_p<p>_R<R>.csv
phi_sweep_ratio_p<p>_R<R>.png
phi_sweep_metrics_p<p>_R<R>.png
phi_knee_kneedle_diagnostic_p<p>_R<R>.png
```

The resulting `phi_knee` is interpreted as a tracer-defined percolation threshold, where tracer motion begins to experience strong confinement due to disconnected accessible pore spaces.

---

### `Inverse_total.R`

This script performs inverse analysis by searching for matrix geometries that reproduce a target MSD.

Run:

```bash
Rscript PPGP/Code/Inverse_total.R
```

or in R:

```r
setwd("PPGP/Code")
source("Inverse_total.R")
```

Path override example:

```bash
Rscript PPGP/Code/Inverse_total.R --root="PPGP" --exp="0.72_sim.csv"
```

The script expects:

```text
PPGP/Data/Training/training_input.csv
PPGP/Data/Training/training_MSD.csv
PPGP/Data/Training/delta_t.csv
PPGP/Data/Testing/testing_input.csv
PPGP/Data/Testing/testing_MSD.csv
PPGP/Data/Input/0.72_sim.csv
PPGP/Code/functions/rcalibration_no_discrepancy_local.R
```

The target MSD file should contain:

```text
dt
msd
```

Optional columns:

```text
lower
upper
```

When `lower` and `upper` are provided, the inverse analysis uses uncertainty-based weights in log-MSD space. Otherwise, uniform weights are used.

The script:

1. Loads training and testing simulation data.
2. Loads a target MSD from `PPGP/Data/Input/`.
3. Fits a PPGP surrogate model.
4. Computes weighted RMSE between predicted and target MSDs in log10 space.
5. Samples the 3D parameter space using Latin hypercube sampling.
6. Identifies low-RMSE regions and representative local minima.
7. Generates a 3D inverse-map figure showing RMSE and selected local minima.

Default global search ranges are:

```r
phi_range = c(0.0, 0.75)
p_range   = c(0.1, 1.0)
R_range   = c(2.0, 6.0)
N_global  = 12000
```

Typical outputs in `PPGP/Output/` include:

```text
global_LHS_12000_samples_RMSE.csv
Fig4_best50_withLocalMinima_far3.png
```

---

## Requirements

### MATLAB

The MPT and simulation workflows require MATLAB.

Recommended:

* MATLAB R2021a or later
* Image Processing Toolbox, recommended for TIFF image-stack handling and microscopy-image analysis

---

### R

The AIUQ and PPGP workflows require R.

Core packages used by the PPGP scripts include:

```r
install.packages(c(
  "data.table",
  "RobustGaSP",
  "RobustCalibration",
  "lhs",
  "MASS",
  "mvtnorm",
  "scatterplot3d",
  "ggplot2",
  "dplyr",
  "tidyr",
  "readr",
  "pracma"
))
```

Additional packages used by the AIUQ workflow include:

```r
install.packages(c(
  "fftwtools",
  "geometry",
  "plot3D",
  "SuperGauss",
  "ijtiff",
  "here"
))
```

The AIUQ workflow also uses:

```r
AIUQ
bioimagetools
EBImage
```

Depending on the R environment, some of these packages may need to be installed from Bioconductor, GitHub, or a local source.

---

## Typical workflow

### 1. Analyze experimental videos with MPT

```matlab
cd Experiment/MPT
MPT_Main
```

This produces tracer trajectories, MSDs, NGP, and VACF outputs.

---

### 2. Analyze experimental videos with AIUQ

```bash
cd Experiment/AIUQ
Rscript AIUQ_nonparametric.R
```

This produces an AIUQ-estimated MSD curve with uncertainty bounds.

---

### 3. Run simulations

```matlab
cd Simulation/Code
Simulation_setup
```

This generates simulated tracer trajectories, MSDs, pore statistics, van Hove functions, NGP, VACF, and optional videos/figures depending on the selected output switches.

---

### 4. Train and visualize the PPGP state space

```bash
Rscript PPGP/Code/forward_3dscatter.R
```

This produces 3D and 2D visualizations of the predicted state space and validates the surrogate model against testing data.

---

### 5. Estimate tracer-defined percolation behavior

```bash
Rscript PPGP/Code/forward_percolation.R
```

This performs a dense `Phi` sweep at fixed `p` and `R_bar`, then identifies the knee in `D_late / D_early`.

---

### 6. Perform inverse MSD analysis

```bash
Rscript PPGP/Code/Inverse_total.R
```

This searches the 3D parameter space for geometries that reproduce a target MSD.

---

## Data notes

Large raw microscopy videos are not tracked directly in this repository.

Empty folders are retained using `.gitkeep` files so that the expected directory structure is preserved after cloning. Users should place their own microscopy videos in the appropriate `Videos/` folder before running MPT or AIUQ analysis.

---

## Reproducibility notes

The scripts are organized to run relative to the repository structure shown above. If folders are renamed or moved, file paths in the corresponding scripts may need to be updated.

For reproducible simulation results, check or set the random seed in the MATLAB simulation scripts before running large parameter sweeps. The simulation runner uses the random seed stored in the configuration structure.

The PPGP scripts support path overrides through command-line arguments such as:

```bash
--root
--data
--out
--exp
--function
```

They also support environment variables such as:

```text
SURROGATE_ROOT_DIR
SURROGATE_DATA_DIR
SURROGATE_OUT_DIR
SURROGATE_EXP_FILE
SURROGATE_FUNCTION_FILE
```

---

## Citation

If you use this repository, please cite the associated manuscript:

```text
Lee, J., Lin, T., Gu, M., and Luo, Y.
Accessible pore geometry governs tracer diffusion in crowded environments.
In preparation.
```

The final journal citation and DOI will be added after publication.

---

## Contact

For questions about the code or data, please contact:

**Jinseok Lee**
Department of Mechanical Engineering and Materials Science
Yale University
Email: [jinseok.lee@yale.edu](mailto:jinseok.lee@yale.edu)

**Yimin Luo**
Department of Mechanical Engineering
Yale University
Email: [yimin.luo@yale.edu](mailto:yimin.luo@yale.edu)


# Data and code for:
A framework for modeling and inferring tracer diffusion in crowded environments

Authors  
Jinseok Lee  
Tong Lin  
Mengyang Gu  
Yimin Luo  

---

## Overview
This repository contains experimental data, simulation outputs, and analysis code used in the manuscript.

The framework combines:
- tracer trajectory experiments
- minimal Monte Carlo simulations
- hydrodynamic corrections
- PPGP surrogate modeling
- inverse inference of crowded environments

---

## Repository structure

data/
    experimental trajectories, MSDs, and van Hove distributions

simulation/
    simulated MSD curves used for surrogate model training

cell_data/
    intracellular tracer MSD measurements

code/
    analysis, simulation, surrogate modeling, and inverse inference

figure_generation/
    scripts to reproduce manuscript figures

---

## Reproducing figures

Example workflow:

Figure 2:
    code/experiment_analysis/compute_MSD.m

Figure 3:
    code/simulation/run_simulation.m

Figure 4:
    code/surrogate_model/train_surrogate.R

---

## Software requirements

MATLAB R2023a  
R 4.3  

R packages:
RobustGaSP  
sensitivity  

---

## Data availability

All data and code required to reproduce the results of this study are available in this repository.

A permanent archived version will be released on Zenodo upon publication.

# Academic Workforce Dynamics in the Mathematical Sciences

This repository contains the data and computational code used to develop and analyze a discrete time dynamical model of the academic pipeline in the mathematical sciences. The study examines how degree production, hiring capacity, and faculty turnover jointly determine long term workforce structure.

## Overview

The model represents transitions across undergraduate, graduate, postdoctoral, and faculty stages using a discrete time compartmental framework. Observed national degree counts are used to reconstruct upstream population levels, while downstream dynamics are governed by completion rates, exit processes, and vacancy limited hiring.

Faculty hiring is constrained by turnover generated openings, and transition from postdoctoral positions depends on competition intensity. This structure produces nonlinear interactions between candidate supply and institutional capacity.

The computational analysis explores how variation in degree inflow and hiring parameters affects faculty growth, postdoctoral accumulation, and system level congestion.

## Repository Contents

* **NCESdata.csv**
  National Center for Education Statistics data used to reconstruct undergraduate and graduate degree flows in the mathematical sciences.

* **simulation.ipynb**
  Python notebook implementing the discrete time model, running simulations, and generating figures presented in the manuscript.

* **Sensitivity_global.m**
  MATLAB script used to perform global sensitivity analysis based on parameter sampling and correlation measures.

## Requirements

### Python

* numpy
* pandas
* matplotlib

### MATLAB

* MATLAB environment for running global sensitivity analysis script

## Usage

1. Run `simulation.ipynb` to reproduce model simulations and figures.
2. Use `NCESdata.csv` as input data for degree flow reconstruction.
3. Run `Sensitivity_global.m` in MATLAB to reproduce global sensitivity analysis.

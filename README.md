# Thesis: Robust Control Chart for Individual Observations in Phase I for High-Dimensional Processes
This repository contains the implementation of robust control charts for high-dimensional processes developed as part of a research thesis.

## Project Structure

The project is organized into two main folders, each focusing on different statistical distributions:

### ChartControlT2MRCDNormal/
This folder contains the implementation for normally distributed data.

**Main file:** `T2MRCDNormal.R`

Support files:
- `AlgorithmEBADI.R` - EBADI algorithm
- `AlgorithmRoMDP.R` - RoMDP algorithm
- `MethodUCLKernel.R` - UCL Kernel method
- `SignalProbability.R` - Signal probability calculation
- `SimulationT2ChartMedia.R` - T2 chart simulations

### ChartControlT2MRCDGamma/
This folder contains the implementation for gamma-distributed data.

**Main file:** `T2MRCDGamma.R`

Supporting files:
- `AlgorithmEBADI.R` - EBADI algorithm
- `AlgorithmRoMDP.R` - RoMDP algorithm
- `MethodUCLKernel.R` - UCL Kernel method
- `SignalProbability.R` - Signal probability calculation
- `SimulationT2ChartGammaMedia.R` - T² chart simulations for gamma distribution

## Usage

To run each project:

1. **For normal distribution:** Run the file `ChartControlT2MRCDNormal/T2MRCDNormal.R`
2. **For gamma distribution:** Run the file `ChartControlT2MRCDGamma/T2MRCDGamma.R`

## Description

This work implements robust T² control charts using the MRCD (Minimum Regularized Covariance Determinant) estimator for high-dimensional multivariate processes. The algorithms include robust outlier detection and parameter estimation methods to improve the performance of control charts in the presence of contaminated data.

## Requirements

- R (recommended version: 4.0 or higher)
- Required R packages (specified in each main file)

## Author

Developed as part of a research thesis on robust control charts for high-dimensional processes.
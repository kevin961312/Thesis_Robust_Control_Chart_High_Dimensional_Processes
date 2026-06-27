# Thesis: Robust Control Chart for Individual Observations in Phase I for High-Dimensional Processes

This repository contains Monte Carlo simulation notebooks for robust T² control charts applied to high-dimensional processes. Each notebook runs simulations comparing chart detection methods and saves results as `.RData` files. The `PlotResults.R` script loads those files and produces signal probability comparison plots saved as SVG.

## Project Structure

### ChartControlT2MRCDNormal/

Simulations for **normally distributed** multivariate data.

| File | What it does | Output `.RData` |
|------|-------------|-----------------|
| `Simulation(MRCD,TMOD,RMDP).ipynb` | Monte Carlo comparison of four methods — T²-MRCD, T²MOD, EBADI-Ui, and EBADI-Zi — across varying mean-shift magnitudes (delta) and outlier percentages. Parallel execution via `foreach`/`doParallel`. | `SinalprobabilityNormal{n}x{p}x{ε}x10000.RData` |
| `SimulationMRCDBRP.ipynb` | Same four-method comparison but uses a bootstrap procedure (K-MRCD / BRP) to automatically select the MRCD regularization parameter α for each shift level. | `SinalprobabilityNormalMC{n}x{p}x{ε}.RData` |
| `SimulationSVTDD.ipynb` | Monte Carlo evaluation of the SVTDD (Support Vector Tensor Data Description) chart. Observations are reshaped into tensors; `σ` and `C` are tuned via leave-one-out cross-validation. | `SinalprobabilitySVTDDNormalMC{n}x{p}x{ε}.RData` |

`rrcov.zip` — patched local build of the `rrcov` package (version 1.7-7.9000) required by all notebooks.

### ChartControlT2MRCDGamma/

Simulations for **gamma-distributed** multivariate data. Observations are generated via a Gaussian copula — each marginal follows a Gamma(shape, rate) distribution — allowing flexible dependence control through an AR(1) correlation matrix.

| File | What it does | Output `.RData` |
|------|-------------|-----------------|
| `SimulationGamma(MRCD,TMOD,RMDP).ipynb` | Monte Carlo comparison of four methods — T²-MRCD, T²MOD, EBADI-Ui, and EBADI-Zi — for gamma data. The shift is applied to the shape parameter (`rho + delta`), which changes the mean of each marginal. Default scenario: Gamma(2, 5). | `SinalprobabilityGamma(2,5){n}x{p}x{ε}x10000.RData` |
| `SimulationGammaMRCDBRP.ipynb` | Same MRCD simulation but uses the BRP bootstrap (K-MRCD) to automatically select the regularization parameter α per shift level. Works with gamma data via the same copula construction. | `SinalprobabilityGammaMRCDBRP(rho,beta){n}x{p}x{ε}.RData` |
| `SimulationGammaSTVDD.ipynb` | SVTDD chart for gamma data. Observations are reshaped into tensors (e.g., 25×10 for p=250); σ and C are tuned via LOO cross-validation. Uses `rmvgamma` for data generation. | `SinalprobabilitySVTDDGamma(2,5){n}x{p}x{ε}.RData` |

`rrcov.zip` — same patched package, copy kept per-folder for Colab uploads.

## Key variables inside each `.RData`

After loading an `.RData` file the following matrices are available:

| Variable | Column 1 | Column 2 | Cols 3–4 | Source notebook |
|----------|----------|----------|----------|-----------------|
| `MatrixDeltaMRCD` | delta | signal prob (UCL) | UCL-max, UCL-kernel | `Simulation…` / `SimulationMRCDBRP` |
| `MatrixDeltaT2MOD` | delta | signal prob (UCL) | UCL-max, UCL-kernel | `Simulation…` |
| `MatrixDeltaEBADIZI` | delta | signal prob (UCL) | UCL-max, UCL-kernel | `Simulation…` |
| `MatrixDelta` | row index | signal prob | — | `SimulationSVTDD` / `SimulationGammaSTVDD` |

> `MatrixDelta` from the SVTDD notebooks has a spurious first row (artifact of `rbind` initialization) — `PlotResults.R` removes it with `SVTDD$MatrixDelta <- SVTDD$MatrixDelta[-1, ]`.

### Google Drive output folders

| Notebook | Drive folder |
|----------|-------------|
| `Simulation(MRCD,TMOD,RMDP).ipynb` | `Colab10000simulacion/` |
| `SimulationMRCDBRP.ipynb` | `Colab10000BRPNormal/` |
| `SimulationSVTDD.ipynb` | `Colab10000SVTDDNormal/` |
| `SimulationGamma(MRCD,TMOD,RMDP).ipynb` | `Colab10000Gamma/` |
| `SimulationGammaMRCDBRP.ipynb` | `Colab10000BRPGamma/` |
| `SimulationGammaSTVDD.ipynb` | `Colab10000SVTDDGamma/` |

## Visualizing results

`PlotResults.R` loads the three `.RData` files for a given scenario and saves an SVG signal-probability-vs-delta comparison plot.

1. Open `PlotResults.R` and edit the configuration block at the top (paths and `dimension` string).
2. Run from the terminal:

```bash
Rscript PlotResults.R
```

The plot overlays five curves:

| Color | Method |
|-------|--------|
| Red | T²-MRCD (fixed α) |
| Purple | T²-MRCD (BRP α selection) |
| Green | T²MOD |
| Blue | EBADI-Zi |
| Orange | SVTDD |

The SVG is saved to the images directory specified in the script.

## How to run the simulations

The notebooks are designed for **Google Colab** — they authenticate with Google Drive and upload results automatically.

1. Upload the notebook and `rrcov.zip` to a Colab session.
2. Run all cells in order. The first cell authenticates with Google Drive.
3. Adjust `Observation`, `NumberVariable`, and `Percentoutliers` in the simulation cell as needed.
4. The `.RData` result is saved locally and uploaded to the Drive folder specified in the last cell.

## Requirements

- R ≥ 4.0
- R packages: `MASS`, `rrcov` (patched — install from `rrcov.zip`), `EnvStats`, `foreach`, `doParallel`, `Rfast`, `psych`, `expm`, `KernSmooth`, `StableMCD`, `quadprog`, `randcorr`, `mvtnorm`, `googledrive`, `svglite`
- Google Colab environment with multiple CPU cores (notebooks call `detectCores()`)

## Author

Developed as part of a research thesis on robust control charts for high-dimensional processes.

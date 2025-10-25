rm(list = ls())
library(here)
library(MASS)
library(clusterGeneration)
library(rrcov)
library(EnvStats)
library(randcorr)
library(KernSmooth)
library(psych)
library(expm)
library(Rfast)
library(copula)
library(lcmix)
library(foreach)
library(doParallel)

source(here::here("./MethodUCLKernel.R"))
source(here::here("./SignalProbability.R"))
source(here::here("./SimulationT2ChartGammaMedia.R"))
set.seed(123)

# Parameters Simulación
registerDoParallel(5)

Observation <- 100
NumberVariable <- 250
NumSimulation <- 1000
rho <- rep(1:5, times = NumberVariable / length(1:5))
beta <- rep(c(1, 3, 0.5, 2, 5), times = NumberVariable / length(1:5))
AlphaMRCD <- 0.75
AlphaPFA <- 0.05

Percentoutliers <- 0.2
NumSimulationDelta <- 1000

fun <- function(i, j) (0.5)^(abs(i - j))
rows <- 1:NumberVariable
cols <- 1:NumberVariable

SigmaCorr <- outer(rows, cols, FUN = fun)

# Parte 1 Calculo de los Limites de Probabilidad MRCD

valuesT2Total <- foreach(i = c("MRCD", "T2MOD", "EBADIUI", "EBADIZI"), .combine = "cbind") %dopar% {
  SimulationT2Chart(
    Observation,
    NumberVariable,
    NumSimulation,
    SigmaCorr,
    rho,
    AlphaMRCD,
    i, beta
  )
}

# MRCD
ValuesMRCDT2Matrix <- valuesT2Total[1, 1]
ValuesMRCDT2Total <- valuesT2Total[2, 1]
ValuesMRCDT2Max <- valuesT2Total[3, 1]
KernelMethodMRCD <- MethodUCLKernel(ValuesMRCDT2Total[[1]], AlphaPFA)
UCLMRCD <- qemp(p = 1 - AlphaPFA, obs = ValuesMRCDT2Total[[1]])
UCLMaxMRCD <- qemp(p = 1 - AlphaPFA, obs = ValuesMRCDT2Max[[1]])
UCLKernelMRCD <- KernelMethodMRCD$UCL
SignalProbabilityT2MRCD <- SignalProbability(ValuesMRCDT2Matrix[[1]], UCLMRCD, UCLMaxMRCD, UCLKernelMRCD)

# T2MOD
ValuesT2MODT2Matrix <- valuesT2Total[1, 2]
ValuesT2MODT2Total <- valuesT2Total[2, 2]
ValuesT2MODT2Max <- valuesT2Total[3, 2]
KernelMethodT2MOD <- MethodUCLKernel(ValuesT2MODT2Total[[1]], AlphaPFA)
UCLT2MOD <- qemp(p = 1 - AlphaPFA, obs = ValuesT2MODT2Total[[1]])
UCLMaxT2MOD <- qemp(p = 1 - AlphaPFA, obs = ValuesT2MODT2Max[[1]])
UCLKernelT2MOD <- KernelMethodT2MOD$UCL
SignalProbabilityT2MOD <- SignalProbability(ValuesT2MODT2Matrix[[1]], UCLT2MOD, UCLMaxT2MOD, UCLKernelT2MOD)

# EBADIUI
ValuesEBADIUIT2Matrix <- valuesT2Total[1, 3]
ValuesEBADIUIT2Total <- valuesT2Total[2, 3]
ValuesEBADIUIT2Max <- valuesT2Total[3, 3]
KernelMethodEBADIUI <- MethodUCLKernel(ValuesEBADIUIT2Total[[1]], AlphaPFA)
UCLEBADIUI <- qemp(p = 1 - AlphaPFA, obs = ValuesEBADIUIT2Total[[1]])
UCLMaxEBADIUI <- qemp(p = 1 - AlphaPFA, obs = ValuesEBADIUIT2Max[[1]])
UCLKernelEBADIUI <- KernelMethodEBADIUI$UCL
SignalProbabilityT2EBADIUI <- SignalProbability(ValuesEBADIUIT2Matrix[[1]], UCLEBADIUI, UCLMaxEBADIUI, UCLKernelEBADIUI)

# EBADIZI
ValuesEBADIZIT2Matrix <- valuesT2Total[1, 4]
ValuesEBADIZIT2Total <- valuesT2Total[2, 4]
ValuesEBADIZIT2Max <- valuesT2Total[3, 4]
KernelMethodEBADIZI <- MethodUCLKernel(ValuesEBADIZIT2Total[[1]], AlphaPFA)
UCLEBADIZI <- qemp(p = 1 - AlphaPFA, obs = ValuesEBADIZIT2Total[[1]])
UCLMaxEBADIZI <- qemp(p = 1 - AlphaPFA, obs = ValuesEBADIZIT2Max[[1]])
UCLKernelEBADIZI <- KernelMethodEBADIZI$UCL
SignalProbabilityT2EBADIZI <- SignalProbability(ValuesEBADIZIT2Matrix[[1]], UCLEBADIZI, UCLMaxEBADIZI, UCLKernelEBADIZI)

# Parte 2 Calculo del parametro de no centralidad
ArrayRhoShift <- c(
  0,
  1 / 100,
  5 / 100,
  10 / 100,
  15 / 100,
  25 / 100,
  50 / 100,
  60 / 100,
  75 / 100,
  90 / 100,
  1
)

Observation <- 10000
Data <- rmvgamma(Observation, rho, beta, SigmaCorr)
CovMatriz <- cov(Data)
Inverse <- solve(CovMatriz)


MatrixDeltaMRCD <- foreach(rhoShift = ArrayRhoShift, .combine = "rbind") %dopar% {
  DeltaNCP <- sqrt(t((rho + rhoShift) / beta - rho / beta) %*% Inverse %*% ((rho + rhoShift) / beta - rho / beta))
  rhosum <- rho + rhoShift
  SimulationT2ChartOutliers(
    Observation,
    NumberVariable,
    NumSimulationDelta,
    SigmaCorr,
    rho,
    rhosum,
    Percentoutliers,
    AlphaMRCD,
    "MRCD",
    beta,
    DeltaNCP[1, 1],
    UCLMRCD,
    UCLMaxMRCD,
    UCLKernelMRCD,
    rhoShift
  )
}

MatrixDeltaT2MOD <- foreach(rhoShift = ArrayRhoShift, .combine = "rbind") %dopar% {
  DeltaNCP <- sqrt(t((rho + rhoShift) / beta - rho / beta) %*% Inverse %*% ((rho + rhoShift) / beta - rho / beta))
  SimulationT2ChartOutliers(
    Observation,
    NumberVariable,
    NumSimulationDelta,
    SigmaCorr,
    rho,
    rho + rhoShift,
    Percentoutliers,
    AlphaMRCD,
    "T2MOD",
    beta,
    DeltaNCP[1, 1],
    UCLT2MOD,
    UCLMaxT2MOD,
    UCLKernelT2MOD,
    rhoShift
  )
}
MatrixDeltaEBADIUI <- foreach(rhoShift = ArrayRhoShift, .combine = "rbind") %dopar% {
  DeltaNCP <- sqrt(t((rho + rhoShift) / beta - rho / beta) %*% Inverse %*% ((rho + rhoShift) / beta - rho / beta))
  SimulationT2ChartOutliers(
    Observation,
    NumberVariable,
    NumSimulationDelta,
    SigmaCorr,
    rho,
    rho + rhoShift,
    Percentoutliers,
    AlphaMRCD,
    "EBADIUI",
    beta,
    DeltaNCP[1, 1],
    UCLEBADIUI,
    UCLMaxEBADIUI,
    UCLKernelEBADIUI,
    rhoShift
  )
}

MatrixDeltaEBADIZI <- foreach(rhoShift = ArrayRhoShift, .combine = "rbind") %dopar% {
  DeltaNCP <- sqrt(t((rho + rhoShift) / beta - rho / beta) %*% Inverse %*% ((rho + rhoShift) / beta - rho / beta))
  SimulationT2ChartOutliers(
    Observation,
    NumberVariable,
    NumSimulationDelta,
    SigmaCorr,
    rho,
    rho + rhoShift,
    Percentoutliers,
    AlphaMRCD,
    "EBADIZI",
    beta,
    DeltaNCP[1, 1],
    UCLEBADIZI,
    UCLMaxEBADIZI,
    UCLKernelEBADIZI,
    rhoShift
  )
}

p <- plot(MatrixDeltaMRCD[, 1], MatrixDeltaMRCD[, 2], ylim = c(0, 1), col = "red", pch = 16, type = "b", xlab = "Delta Values", ylab = "Signal Probability", lty = 1)
lines(MatrixDeltaT2MOD[, 1], MatrixDeltaT2MOD[, 2], col = "green", pch = 17, type = "b", lty = 4)
lines(MatrixDeltaEBADIUI[, 1], MatrixDeltaEBADIUI[, 2], col = "blue", pch = 18, type = "b", lty = 5)
lines(MatrixDeltaEBADIZI[, 1], MatrixDeltaEBADIZI[, 2], col = "orange", pch = 15, type = "b", lty = 6)


rm(list = ls())
library(svglite)
urlplot <- "/Users/kevin.pineda/Downloads/Tesis/TesisCode/Imagenes TG 2 KMRCD/Gamma(2,5)/ImagenGamma(2,5)50x250x0.2.svg"
load("~/Downloads/Tesis/TesisCode/KMRCD/DatosKMRCD/Parametro(2,5)/SinalprobabilityGamma(2,5)50x250x0.2.RData")
for (i in 1:length(ArrayRhoShift)) {
  rhoShift <- ArrayRhoShift[[i]]
  DeltaNCP <- sqrt(t(rep((rho + rhoShift) / beta, NumberVariable) - rep(rho / beta, NumberVariable)) %*% Inverse %*% (rep((rho + rhoShift) / beta, NumberVariable) - rep(rho / beta, NumberVariable)))
  MatrixDeltaMRCD[i, 1] <- DeltaNCP
  MatrixDeltaT2MOD[i, 1] <- DeltaNCP
  MatrixDeltaEBADIUI[i, 1] <- DeltaNCP
  MatrixDeltaEBADIZI[i, 1] <- DeltaNCP
}

par(mar = c(10, 5, 5, 5))
p <- plot(MatrixDeltaMRCD[, 1], MatrixDeltaMRCD[, 2], ylim = c(0, 1), col = "red", pch = 16, type = "b", xlab = "Delta Values", ylab = "Signal Probability", lty = 1)
lines(MatrixDeltaT2MOD[, 1], MatrixDeltaT2MOD[, 2], col = "green", pch = 17, type = "b", lty = 4)
lines(MatrixDeltaEBADIUI[, 1], MatrixDeltaEBADIUI[, 2], col = "blue", pch = 18, type = "b", lty = 5)
lines(MatrixDeltaEBADIZI[, 1], MatrixDeltaEBADIZI[, 2], col = "orange", pch = 15, type = "b", lty = 6)


svglite(urlplot, width = 8, height = 8)
p <- plot(MatrixDeltaMRCD[, 1], MatrixDeltaMRCD[, 2], ylim = c(0, 1), col = "red", pch = 16, type = "b", xlab = "Delta Values", ylab = "Signal Probability", lty = 1)
lines(MatrixDeltaT2MOD[, 1], MatrixDeltaT2MOD[, 2], col = "green", pch = 17, type = "b", lty = 4)
lines(MatrixDeltaEBADIUI[, 1], MatrixDeltaEBADIUI[, 2], col = "blue", pch = 18, type = "b", lty = 5)
lines(MatrixDeltaEBADIZI[, 1], MatrixDeltaEBADIZI[, 2], col = "orange", pch = 15, type = "b", lty = 6)

dev.off()

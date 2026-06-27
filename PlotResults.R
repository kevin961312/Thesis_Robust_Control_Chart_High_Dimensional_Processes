library(svglite)

# -- Configure these for each scenario ----------------------------------------
dimension  <- "(2,5)150x250x0.05"   # shape + obs x vars x outlier%

base_dir   <- "/Users/KEPINEDA/Downloads/Gamma organizado"
dir_MRCD   <- file.path(base_dir, "Colab10000Gamma")
dir_MRCDH  <- file.path(base_dir, "GammaBRP/Gamma(2,5)")
dir_SVTDD  <- file.path(base_dir, "Colab10000SVTDDGamma")
dir_images <- file.path(base_dir, "Gamma(2,5)Images")
# -----------------------------------------------------------------------------

dir.create(dir_images, showWarnings = FALSE, recursive = TRUE)

MRCD  <- new.env()
MRCDH <- new.env()
SVTDD <- new.env()

load(file.path(dir_MRCD,  paste0("SinalprobabilityGamma",        dimension, "x10000.RData")), envir = MRCD)
load(file.path(dir_MRCDH, paste0("SinalprobabilityGammaMRCDBRP", dimension, ".RData")),       envir = MRCDH)
load(file.path(dir_SVTDD, paste0("SinalprobabilitySVTDDGamma",   dimension, ".RData")),       envir = SVTDD)

SVTDD$MatrixDelta <- SVTDD$MatrixDelta[-1, ]

svg_path <- file.path(dir_images, paste0("Gamma", dimension, ".svg"))
svglite(svg_path, width = 8, height = 6)

par(mar = c(10, 5, 5, 5))

plot(
  MRCD$MatrixDeltaEBADIZI[, 1], MRCD$MatrixDeltaT2MOD[, 2],
  ylim = c(0, 1), col = "green", pch = 17, type = "b",
  xlab = "Delta Values", ylab = "Signal Probability", lty = 4
)
lines(MRCD$MatrixDeltaEBADIZI[, 1], MRCD$MatrixDeltaEBADIZI[, 2], col = "blue",   pch = 18, type = "b", lty = 5)
lines(MRCD$MatrixDeltaEBADIZI[, 1], MRCD$MatrixDeltaMRCD[, 2],    col = "red",    pch = 16, type = "b", lty = 1)
lines(MRCD$MatrixDeltaEBADIZI[, 1], SVTDD$MatrixDelta[, 2],       col = "orange", pch = 15, type = "b", lty = 6)
lines(MRCD$MatrixDeltaEBADIZI[, 1], MRCDH$MatrixDeltaMRCD[, 2],   col = "purple", pch = 16, type = "b", lty = 1)

legend(
  "bottomright",
  legend = c("T²MOD", "EBADI-Zi", "T²-MRCD", "SVTDD", "T²-MRCD (BRP)"),
  col    = c("green",      "blue",      "red",          "orange", "purple"),
  pch    = c(17,           18,           16,             15,       16),
  lty    = c(4,            5,            1,              6,        1),
  bty    = "n"
)

dev.off()
cat("Saved:", svg_path, "\n")

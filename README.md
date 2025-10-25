# Tesis: Carta de Control Robusta para Observaciones Individuales en Fase I para Procesos de Alta Dimensionalidad
Este repositorio contiene la implementación de cartas de control robustas para procesos de alta dimensionalidad desarrollada como parte de una tesis de investigación.

## Estructura del Proyecto

El proyecto está organizado en dos carpetas principales, cada una enfocada en diferentes distribuciones estadísticas:

### ChartControlT2MRCDNormal/
Esta carpeta contiene la implementación para datos con distribución normal.

**Archivo principal:** `T2MRCDNormal.R`

Archivos de soporte:
- `AlgorithmEBADI.R` - Algoritmo EBADI
- `AlgorithmRoMDP.R` - Algoritmo RoMDP
- `MethodUCLKernel.R` - Método UCL Kernel
- `SignalProbability.R` - Cálculo de probabilidades de señal
- `SimulationT2ChartMedia.R` - Simulaciones de carta T²

### ChartControlT2MRCDGamma/
Esta carpeta contiene la implementación para datos con distribución gamma.

**Archivo principal:** `T2MRCDGamma.R`

Archivos de soporte:
- `AlgorithmEBADI.R` - Algoritmo EBADI
- `AlgorithmRoMDP.R` - Algoritmo RoMDP
- `MethodUCLKernel.R` - Método UCL Kernel
- `SignalProbability.R` - Cálculo de probabilidades de señal
- `SimulationT2ChartGammaMedia.R` - Simulaciones de carta T² para distribución gamma

## Uso

Para ejecutar cada proyecto:

1. **Para distribución normal:** Ejecute el archivo `ChartControlT2MRCDNormal/T2MRCDNormal.R`
2. **Para distribución gamma:** Ejecute el archivo `ChartControlT2MRCDGamma/T2MRCDGamma.R`

## Descripción

Este trabajo implementa cartas de control T² robustas utilizando el estimador MRCD (Minimum Regularized Covariance Determinant) para procesos multivariados de alta dimensionalidad. Los algoritmos incluyen métodos robustos de detección de puntos atípicos y estimación de parámetros para mejorar el rendimiento de las cartas de control en presencia de datos contaminados.

## Requisitos

- R (versión recomendada: 4.0 o superior)
- Paquetes de R necesarios (especificados en cada archivo principal)

## Autor

Desarrollado como parte de una tesis de investigación sobre cartas de control robustas para procesos de alta dimensionalidad.
# GameStop Volatility Analysis

Financial time series analysis of GameStop stock using GARCH models to understand volatility patterns, particularly during the 2021 retail trading events.

## Overview

This project analyzes GameStop (GME) stock volatility from 2002 to 2023 using various GARCH models. The analysis focuses on volatility clustering, leverage effects, and parameter stability across different time periods.

## Repository Structure

```
gamestop_analysis/
├── README.md                           # Project documentation
├── code/
│   ├── gamestop.R                     # Main analysis script
│   ├── TSA-Predict-Student-Functions.R # Helper functions for time series analysis
│   └── TSA-Finance-Functions.R         # Finance-specific functions
└── report/
    └── gamestop.pdf                   # Analysis report (Italian)
```

## How to Run

1. **Install required R packages:**
```r
install.packages(c("tseries", "sandwich", "lmtest", "urca", 
                   "rugarch", "FinTS", "car", "forecast", 
                   "xts", "quantmod"))
```

2. **Download and set up:**
   - Clone this repository
   - Set your working directory to the `code/` folder
   - Ensure all three R files are in the same directory

3. **Run the analysis:**
```r
source("gamestop.R")
```

**Note:** Data is automatically downloaded from Yahoo Finance - no manual data files needed.

## Models and Methods

### Time Series Models
- **ARMA(1,1)** for mean equation modeling
- **sGARCH(1,1)** - Simple GARCH for volatility
- **GJR-GARCH(1,1)** - Asymmetric GARCH (best performing)
- **T-GARCH(1,1)** - Threshold GARCH

### Statistical Tests
- **Ljung-Box tests** for serial correlation
- **ARCH-LM tests** for heteroskedasticity  
- **Nyblom tests** for parameter stability
- **Sign bias tests** for leverage effects

### Key Analysis Features
- Volatility clustering identification
- News Impact Curve analysis
- Model comparison using information criteria
- Out-of-sample forecasting evaluation
- Garman-Klass volatility benchmark

## Key Findings

| Finding | Description |
|---------|-------------|
| **Volatility Clustering** | Strong evidence in GME returns, especially during 2021 |
| **Best Model** | GJR-GARCH outperforms other specifications |
| **Leverage Effects** | Significant asymmetric volatility response to shocks |
| **Parameter Stability** | Instability detected in full sample, improved in subsample (2008-2023) |
| **Distribution** | Student-t distribution fits better than normal |

## Data

- **Source:** Yahoo Finance (via `quantmod` package)
- **Symbol:** GME (GameStop Corp.)
- **Period:** January 3, 2002 - November 22, 2023
- **Frequency:** Daily
- **Variables:** Open, High, Low, Close, Adjusted Close, Volume

## Model Performance

Based on subsample analysis (2008-2023):

| Model | AIC | BIC |
|-------|-----|-----|
| sGARCH | 4.9886 | 4.9946 |
| **GJR-GARCH** | **4.9825** | **4.9897** |
| T-GARCH | 4.9831 | 4.9903 |

*Lower values indicate better fit*

## Technical Implementation

- **Language:** R
- **Key Libraries:** rugarch, quantmod, xts, forecast
- **Estimation Method:** Maximum Likelihood (SOLNP solver)
- **Distribution:** Student-t for innovation errors
- **Validation:** Extensive diagnostic testing and residual analysis

## Usage Notes

- All file paths are relative - no need to modify source paths
- Internet connection required for data download
- Analysis takes several minutes to complete
- All plots and results are displayed automatically

<!-- badges: start -->
[![R-CMD-check](https://github.com/mschemmel/qpcR/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/mschemmel/qpcR/actions/workflows/R-CMD-check.yaml)
[![codecov](https://codecov.io/gh/mschemmel/qpcR/graph/badge.svg?token=1PPYBCNCU7)](https://codecov.io/gh/mschemmel/qpcR)
<!-- badges: end -->


# qpcR

Calculate relative gene expression values of raw qPCR data. Base R. No dependencies. Detailed documentation will follow soon.

## Quick start

```r
library(devtools)
devtools::install_github("mschemmel/qpcR")

# load example data
qpcRdata <- read.table(system.file("extdata", "example2.tsv", package = "qpcR"), sep = "\t", head = TRUE)

# get mean relative expression
qpcR(qpcRdata, hkg = "HKG", groups = "dpi")
```

## Input data
The minimal structure of the input data has to contain following columns: `gene`, `treatment`, `cq`. Primer `efficiency` values are optional.

| Column | Description | Note |
|--------|-------------|------|
| `gene` | investigated genes | |
| `treatment` | variable to compare| |
| `cq` | cq value measured by qPCR machine | |
| `efficiency` | primer efficiency values (%)| (optional) assumed to be 100 % if not provided |

If `cq` contains NA, `brep`, `trep` or both columns are needed to properly exclude samples.

| Column | Description |
|--------|-------------|
| `brep` | number of biological replicate |
| `trep` | number of technical replicate |


Example datasets can be found in the `inst/extdata` folder.

## Output
qpcR outputs a data frame with following columns:

| Column | Description |
| ------ | ----------- |
| treatment | treatments used in your experiment |
| gene | all investigated genes except housekeeping gene(s) |
| dpi | your group(s) variable (if input data was grouped) |
| rexpr.mean | mean expression value |
| rexpr.sd | standard deviation of expression values |
| rexpr.se | standard error of expression values |
| rexpr.n | number observations per gene x group(s) x treatment variables |


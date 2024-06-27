<!-- badges: start -->
[![R-CMD-check](https://github.com/mschemmel/qpcR/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/mschemmel/qpcR/actions/workflows/R-CMD-check.yaml)
<!-- badges: end -->


# qpcR

Calculate relative gene expression values of raw qPCR data. Base R. No dependencies.

## Quick start
Clone the repo using `git clone REPO`.

```r
library(devtools)
devtools::load_all(PATH/TO/REPO)

# load example data
qpcRdata <- read.table("./data/example2.tsv", sep = "\t", head = TRUE)

# get mean relative expression
qpcR(qpcRdata, hkg = "HKG", groups = "dpi")
```
In order to get raw values set `aggregate = FALSE`.

## Input data
The minimal structure of the input data has to contain following columns: `gene`, `treatment`, `cq`. Primer `efficiency` values are optional.

| Column | Description | Note |
|--------|-------------|------|
| `gene` | investigated genes | |
| `treatment` | variable to compare|
| `cq` | cq value measured by qPCR machine | |
| `efficiency` | primer efficiency values (%)| (optional) assumed to be 100 % if not provided |

If `cq` contains NA, `brep`, `trep` columns are needed to properly exclude samples.

| Column | Description |
|--------|-------------|
| `brep` | number of biological replicate |
| `trep` | number of technical replicate |


An example dataset can be found in the `data` folder.



## Output
qpcR outputs a data frame with following columns:

| Column | Description |
| ------ | ----------- |
| treatment | treatments used in your experiment |
| gene | all investigated genes except housekeeping gene(s) |
| dpi | your group(s) variable |
| rexpr.mean | mean expression value |
| rexpr.sd | standard deviation of expression values |
| rexpr.se | standard error of expression values |
| rexpr.n | number observations per gene x group(s) x treatment variables |


qpcR
---

Calculate relative gene expression of qPCR data

Minimal requirement of input data
---
The minimal structure of the input data has to contain following columns: `gene`, `treatment`, `cq`.

| Column | Description | Note |
|--------|-------------|------|
| gene | investigated gene names |
| treatment | variable(s) to compare |
| cq | cq value measured by qPCR machine |
| efficiency | primer efficiency values (%)| (optional) assumed to be 100 % if not provided

If `cq` contains NA, `brep`, `trep` columns are needed to properly exclude samples.
| Column | Description |
|--------|-------------|
| brep | number of biological replicate |
| trep | number of technical replicate |


An example dataset can be found in the `data` folder.


Quick start
---
Clone the repo using `git clone REPO`. 

```r
devtools::load_all(PATH/TO/REPO)
qpcRdata <- read.table("./data/example2.tsv", sep = "\t", head = TRUE)

# get mean relative expression
qpcR(qpcRdata, hkg = c("HKG"), groups = "dpi")

# in order to get raw values set aggregate = FALSE
qpcR(qpcRdata, hkg = c("HKG"), groups = "dpi", aggregate = FALSE)
```

Testing your changes
---
In order to run all tests, please run:
```
testthat::test_dir("tests/testthat")
```

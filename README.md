qpcR
---

Calculate relative gene expression of qPCR data

Minimal requirement of input data
---
The minimal structure of the input data has to contain following information: `gene`, `treatment`, `cq`, `brep` and `trep`.

| Column | Description |
|--------|-------------|
| gene   | all investigated gene names |
| treatment | variable to compare | 
| cq | cq value measured by qPCR machine |
| brep | number of biological replication |
| trep | number of technical replication |
| efficiency | primer efficiency values |

An example dataset can be found in the `data` folder.


Quick start
---
Clone the repo using `git clone REPO`. 

```r
devtools::load_all(PATH/TO/REPO)
qpcRdata <- read.table("./data/example2.tsv", sep = "\t", head = TRUE)

# get mean relative expression
qpcR(qpcRdata, hkg = c("HKG"), groups = "dpi")
```

Testing your changes
---
In order to run all tests, please run:
```
testthat::test_dir("tests/testthat")
```

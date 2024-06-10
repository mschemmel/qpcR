qpcR
---

Calculate relative gene expression values of qPCR data

Minimal requirement of input data
---
The minimal structure of the input data has to contain following columns: `gene`, `treatment`, `cq`.

| Column | Description | Note |
|--------|-------------|------|
| `gene` | investigated gene names | |
| `treatment` | variable to compare to | default: 'control' |
| `cq` | cq value measured by qPCR machine | |
| `efficiency` | primer efficiency values (%)| (optional) assumed to be 100 % if not provided |

If `cq` contains NA, `brep`, `trep` columns are needed to properly exclude samples.
| Column | Description |
|--------|-------------|
| `brep` | number of biological replicate |
| `trep` | number of technical replicate |


An example dataset can be found in the `data` folder.


Quick start
---
Clone the repo using `git clone REPO`. 

```r
library(devtools)
devtools::load_all(PATH/TO/REPO)

# load example data
qpcRdata <- read.table("./data/example2.tsv", sep = "\t", head = TRUE)

# get mean relative expression
qpcR(qpcRdata, hkg = c("HKG"), groups = "dpi")
```
In order to get raw values set `aggregate = FALSE`.

Output
```
 treatment    gene   dpi  rexpr.mean   rexpr.sd   rexpr.n
   control    gene1   0   1.00000000  0.77378830     4
   control    gene1   3   1.00000000  1.24126171     5
  infected    gene1   0   1.36180622  0.60010957     5
  infected    gene1   3   0.06305240  0.07964586     5
inoculated    gene1   0   2.11208713  2.46449999     6
inoculated    gene1   3   0.07137122  0.08878563     5
```

Testing your changes
---
In order to run all tests, please run:
```
testthat::test_dir("tests/testthat")
```

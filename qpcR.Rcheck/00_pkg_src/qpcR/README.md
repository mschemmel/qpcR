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
treatment    gene   dpi rexpr.mean   rexpr.sd     rexpr.n
control      gene1  0   1.00000000   0.77378830   4
control      gene2  0   1.00000000   0.67041945   4
control      HKG2   0   1.00000000   0.72223068   4
control      gene1  3   1.00000000   1.24126171   5
control      gene2  3   1.00000000   1.52509029   5
control      HKG2   3   1.00000000   1.38556563   5
infected     gene1  0   1.36180622   0.60010957   5
infected     gene2  0   3.73508048   1.03314936   5
infected     HKG2   0   0.07212605   0.04191470   5
infected     gene1  3   0.06305240   0.07964586   5
infected     gene2  3   2.87064646   4.12356984   5
infected     HKG2   3   8.91031466   19.62060012  5
inoculated   gene1  0   2.11208713   2.46449999   6
inoculated   gene2  0   0.11087194   0.15580777   6
inoculated   HKG2   0   0.04133188   0.03066765   6
inoculated   gene1  3   0.07137122   0.08878563   5
inoculated   gene2  3   4.88018151   4.21235147   5
inoculated   HKG2   3   0.51919342   0.39638676   5

```

Testing your changes
---
In order to run all tests, please run:
```
testthat::test_dir("tests/testthat")
```

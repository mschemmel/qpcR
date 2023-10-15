qpcR
---

Calculate relative gene expression of qPCR data

Quick start
---

```r
library(qpcR)

# import data
# setwd() first
qpcRdata <- read.table("./data/example2.tsv", sep = "\t", head = TRUE)

# set primer efficiency
eff <- c("HKG" = 77.5, "gene1" = 79.8)

# get mean relative expression
print("######## Example ########")
print(qpcR(qpcRdata, hkg = c("HKG"), efficiency = eff))
```
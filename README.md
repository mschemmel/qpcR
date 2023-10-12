qpcR
---

Calculate relative gene expression of qPCR data

Quick start
---

```r
library(qpcR)

# import data
# setwd() first
qpcRdata1 <- read.table("./data/example1.tsv", sep = "\t", head = TRUE)
qpcRdata2 <- read.table("./data/example2.tsv", sep = "\t", head = TRUE)

# set primer efficiency
eff1 <- c("HKG" = 100)
eff2 <- c("HKG" = 77.5, "gene1" = 79.8)

# get mean relative expression
print("######## Example 1 ########")
print(qpcR(qpcRdata1, hkg = c("HKG"), efficiency = eff1))

print("######## Example 2 ########")
print(qpcR(qpcRdata2, hkg = c("HKG"), efficiency = eff2))
```
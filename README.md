# qpcR
---
Calculate relative gene expression of qPCR data

# Quick start
---
```r
library(qpcR)

# import data
data(qpcRdata1)
data(qpcRdata2)

# set primer efficiency
eff1 <- c("HKG" = 100)
eff2 <- c("HKG" = 77.5, "gene1" = 79.8)

# get mean relative expression
print("######## Example 1 ########")
print(qpcR(qpcRdata1, hkg = c("HKG"), efficiency = eff1))

print("######## Example 2 ########")
print(qpcR(qpcRdata2, hkg = c("HKG"), efficiency = eff2))
```
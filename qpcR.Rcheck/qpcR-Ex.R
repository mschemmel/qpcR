pkgname <- "qpcR"
source(file.path(R.home("share"), "R", "examples-header.R"))
options(warn = 1)
library('qpcR')

base::assign(".oldSearch", base::search(), pos = 'CheckExEnv')
base::assign(".old_wd", base::getwd(), pos = 'CheckExEnv')
cleanEx()
nameEx("delta_cq")
### * delta_cq

flush(stderr()); flush(stdout())

### Name: delta_cq
### Title: Calculate delta cq values per comparison
### Aliases: delta_cq

### ** Examples

delta_cq(df, contr_mean)



cleanEx()
nameEx("drop_columns")
### * drop_columns

flush(stderr()); flush(stdout())

### Name: drop_columns
### Title: Drop columns based on character vector
### Aliases: drop_columns

### ** Examples

drop_columns(df, cols)



cleanEx()
nameEx("geometric_mean")
### * geometric_mean

flush(stderr()); flush(stdout())

### Name: geometric_mean
### Title: Calculate geometric mean if multiple housekeeping genes are
###   given
### Aliases: geometric_mean

### ** Examples

geometric_mean(c(0.23,0.53,0.12))



cleanEx()
nameEx("get_e")
### * get_e

flush(stderr()); flush(stdout())

### Name: get_e
### Title: Calculate E value of standard samples
### Aliases: get_e

### ** Examples

get_e(named_list)



cleanEx()
nameEx("get_reference_group")
### * get_reference_group

flush(stderr()); flush(stdout())

### Name: get_reference_group
### Title: Select control group of data
### Aliases: get_reference_group

### ** Examples

get_control_group(df)



cleanEx()
nameEx("make_groups")
### * make_groups

flush(stderr()); flush(stdout())

### Name: make_groups
### Title: Create chunks of groups for calculation
### Aliases: make_groups

### ** Examples

make_groups(df, groups = c("plate", "diet", "timepoint"))



cleanEx()
nameEx("pair_wise")
### * pair_wise

flush(stderr()); flush(stdout())

### Name: pair_wise
### Title: Apply calculation for every hkg vs. target pair
### Aliases: pair_wise

### ** Examples

pair_wise(pair, e_values, hkg)



cleanEx()
nameEx("prepare")
### * prepare

flush(stderr()); flush(stdout())

### Name: prepare
### Title: Prepare input data
### Aliases: prepare

### ** Examples

prepare(df)



cleanEx()
nameEx("qpcR")
### * qpcR

flush(stderr()); flush(stdout())

### Name: qpcR
### Title: Main function to perform calculation of relative expression of
###   qpcr data
### Aliases: qpcR

### ** Examples

qpcr(df, hkg = c("HKG"))



cleanEx()
nameEx("ratio_by_mean_ratio")
### * ratio_by_mean_ratio

flush(stderr()); flush(stdout())

### Name: ratio_by_mean_ratio
### Title: Calculate the ratio compared to the mean ratio per gene
### Aliases: ratio_by_mean_ratio

### ** Examples

ratio_by_mean_ratio(df, d_cq, e_val, hkg_, treatm)



cleanEx()
nameEx("sanitize")
### * sanitize

flush(stderr()); flush(stdout())

### Name: sanitize
### Title: Remove NA sample wise as foundation of valid mean relative
###   expression calculation
### Aliases: sanitize

### ** Examples

sanitize(list_of_groups, khg_ = c("HKG"))



### * <FOOTER>
###
cleanEx()
options(digits = 7L)
base::cat("Time elapsed: ", proc.time() - base::get("ptime", pos = 'CheckExEnv'),"\n")
grDevices::dev.off()
###
### Local variables: ***
### mode: outline-minor ***
### outline-regexp: "\\(> \\)?### [*]+" ***
### End: ***
quit('no')

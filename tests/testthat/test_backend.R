dat <- read.table("example2.tsv", sep = "\t", head = TRUE)
test_that("make_groups() works", {
    expect_type(make_groups(dat, groups = NULL), "list")
    expect_type(make_groups(dat, groups = "dpi"), "list")
})

test_that("delta_cq() works", {
    expect_type(delta_cq(dat), "list")
})
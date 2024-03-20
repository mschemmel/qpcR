dat <- read.table("example2.tsv", sep = "\t", head = TRUE)
test_that("make_groups() works", {
    expect_type(make_groups(dat, hkg = NULL, groups = NULL), "list")
    expect_type(make_groups(dat, hkg = NULL, groups = "dpi"), "list")
    expect_error(make_groups(dat, hkg = NULL, groups = "dose"))
})

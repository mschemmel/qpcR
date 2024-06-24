dat <- read.table("example2.tsv", sep = "\t", head = TRUE)
test_that("make_groups() works", {
    expect_type(make_groups(dat, hkg = NULL, groups = NULL), "list")
    expect_type(make_groups(dat, hkg = NULL, groups = "dpi"), "list")
    expect_error(make_groups(dat, hkg = NULL, groups = "dose"))
    expect_error(make_groups(dat[dat$treatment == "control", ], groups = "dpi"))
})

test_that("comp_gene_pair() works", {
    expect_error(comp_gene_pair(dat[dat$treatment == "control", ], hkg = "HKG"))
    expect_type(comp_gene_pair(dat, hkg = "HKG"), "list")
})

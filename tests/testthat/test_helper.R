test_that("set_efficiency() works", {
    expect_equal(set_efficiency(c("gene1", "gene2", "HKG"), list("gene1" = 98, "gene2" = 95)), list("gene1" = 98, "gene2" = 95, "HKG" = 100))
    expect_equal(set_efficiency(c("gene1", "gene2", "HKG"), list("gene1" = 98)), list("gene1" = 98, "gene2" = 100, "HKG" = 100))
})

test_that("get_e() works", {
    expect_equal(get_e(100), 2)
    expect_equal(get_e(0), 1)
})

test_that("drop_columns() works", {
    dat <- read.table("example2.tsv", sep = "\t", head = TRUE)
    expect_equal(names(drop_columns(dat)), c("gene", "treatment", "cq", "brep", "trep"))
    expect_equal(names(drop_columns(dat, "treatment")), c("gene", "cq", "brep", "trep"))
    expect_equal(names(drop_columns(dat, c("treatment", "gene"))), c("cq", "brep", "trep"))
})
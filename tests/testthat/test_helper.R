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

test_that("geometric_mean() works", {
    expect_equal(geometric_mean(c(1:3)), 1.8171206)
    expect_equal(geometric_mean(c(1:3, NA)), 1.8171206)
})

test_that("prepare() works", {
    dat <- data.frame("gene" = rep("gene1", 3), cq = runif(3))
    expect_equal(unique(prepare(dat)$efficiency), 100)
})

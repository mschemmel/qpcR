dat <- read.table("example2.tsv", sep = "\t", head = TRUE)
test_that("get_e() works", {
    expect_equal(get_e(100), 2)
    expect_equal(get_e(0), 1)
    expect_true(all(!(is.na(get_e(c(56, NA))))))
    expect_type(get_e(1), "double")
    expect_type(get_e(NA), "double")
})

test_that("drop_columns() works", {
    expect_equal(names(drop_columns(dat)), c("gene", "treatment", "cq", "brep", "trep", "efficiency", "dpi"))
    expect_equal(names(drop_columns(dat, "treatment")), c("gene", "cq", "brep", "trep", "efficiency", "dpi"))
    expect_equal(names(drop_columns(dat, c("treatment", "gene"))), c("cq", "brep", "trep", "efficiency", "dpi"))
})

test_that("geometric_mean() works", {
    expect_equal(geometric_mean(c(1:3)), 1.8171206)
    expect_equal(geometric_mean(c(1:3, NA)), 1.8171206)
    expect_type(geometric_mean(c(1:3)), "double")
    expect_type(geometric_mean(c(1:3, NA)), "double")
})

test_that("prepare() works", {
    expect_equal(unique(qpcR(drop_columns(dat, "efficiency"), hkg = "HKG", reference = "control", groups = "dpi", aggregate = FALSE)$efficiency), 100)
    expect_error(qpcR(drop_columns(dat, "treatment"), hkg = "HKG", reference = "control", groups = "dpi", aggregate = FALSE)$efficiency)
})

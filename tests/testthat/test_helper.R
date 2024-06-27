dat <- read.table("example2.tsv", sep = "\t", head = TRUE)
assign("reference_group", "control", qenv)

test_that("get_e() works", {
    expect_equal(get_E(100), 2)
    expect_equal(get_E(0), 1)
    expect_type(get_E(1), "double")
})

test_that("drop_columns() works", {
    expect_equal(names(drop_columns(dat)), c("gene", "treatment", "cq", "efficiency", "dpi"))
    expect_equal(names(drop_columns(dat, "treatment")), c("gene", "cq", "brep", "trep", "efficiency", "dpi"))
    expect_equal(names(drop_columns(dat, c("treatment", "gene"))), c("cq", "brep", "trep", "efficiency", "dpi"))
})

test_that("geometric_mean() works", {
    expect_equal(geometric_mean(c(1:3)), 1.8171206)
    expect_equal(geometric_mean(c(1:3, NA)), 1.8171206)
    expect_type(geometric_mean(c(1:3)), "double")
    expect_type(geometric_mean(c(1:3, NA)), "double")
})

test_that("outlier_method() works", {
    num_range <- c(1:4, 100)
    expect_error(outlier_method(c(num_range, "test"), method = "interquartile"))
    expect_error(outlier_method(num_range, method = "test"))
    expect_equal(outlier_method(num_range, "interquartile"), c(FALSE, FALSE, FALSE, FALSE, TRUE))
    expect_equal(outlier_method(num_range, "z-score"), c(FALSE, FALSE, FALSE, FALSE, FALSE))
    expect_equal(outlier_method(num_range, "hampel"), c(FALSE, FALSE, FALSE, FALSE, TRUE))
    expect_type(outlier_method(num_range, "hampel"), "logical")
})

test_that("reference detection reacts properly", {
    expect_equal(detect_reference(c("control", "infected", "inoculated"), "control"), "control")
    expect_equal(detect_reference(c("control", "mock", "infected"), "mock"), "mock")
    expect_error(detect_reference(c("control", "mock", "infected"), NULL))
    expect_error(detect_reference(c("inoculated", "infected"), "test"))
    expect_error(detect_reference(c("inoculated", "infected"), NULL))
})

test_that("assert() works", {
    expect_no_error(assert(2 == 2))
    expect_no_error(assert(c(1 == 1, 2 == 2)))
    expect_error(assert(c(1 == 1, 1 == 2)))
    expect_error(assert())
    expect_error(assert(1 == 2))
    expect_error(assert(1 == 2, "error"))
    expect_error(assert(c(2 == 2, 2 == 3), "error"), "error")
})

test_that("in_range() works", {
    expect_equal(in_range(50, 1, 100), TRUE)
    expect_equal(in_range(1, 1, 100), TRUE)
    expect_equal(in_range(100, 1, 100), TRUE)
    expect_equal(in_range(0, 1, 100), FALSE)
    expect_equal(in_range(101, 1, 100), FALSE)
})
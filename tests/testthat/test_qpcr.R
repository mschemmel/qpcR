dat <- read.table("example2.tsv", sep = "\t", head = TRUE)

test_that("qpcR() works", {
    run <- get_reference_group(qpcR(dat, hkg = "HKG", reference = "control", groups = "dpi", aggregate = FALSE))
    control_groups <- aggregate(rexpr ~ gene + dpi, data = run, FUN = mean)
    expect_equal(unique(control_groups$rexpr), c(1, 1))

    run2 <- get_reference_group(qpcR(dat, hkg = c("HKG", "HKG2"), reference = "control", groups = "dpi", aggregate = FALSE))
    control_groups2 <- aggregate(rexpr ~ gene + dpi, data = run2, FUN = mean)
    expect_equal(unique(control_groups2$rexpr), c(1, 1))

    expect_no_error(qpcR(dat, hkg = "HKG", groups = "dpi", outlier.method = "interquartile"))
    expect_no_error(qpcR(dat, hkg = "HKG", groups = "dpi", outlier.method = "z-score"))
    expect_no_error(qpcR(dat, hkg = "HKG", groups = "dpi", outlier.method = "hampel"))

    expect_no_error(qpcR(dat, hkg = c("HKG", "HKG2"), groups = "dpi", outlier.method = "interquartile"))
    expect_no_error(qpcR(dat, hkg = c("HKG", "HKG2"), groups = "dpi", outlier.method = "z-score"))
    expect_no_error(qpcR(dat, hkg = c("HKG", "HKG2"), groups = "dpi", outlier.method = "hampel"))
})

test_that("output equals data.frame", {
    expect_type(qpcR(dat, hkg = "HKG", groups = "dpi"), "list")
    expect_type(qpcR(dat, hkg = c("HKG", "HKG2"), groups = "dpi"), "list")
})

test_that("all output columns present", {
    expect_equal(colnames(qpcR(dat, hkg = "HKG", groups = "dpi")), c("treatment", "gene", "dpi", "rexpr.mean", "rexpr.sd", "rexpr.se", "rexpr.n"))
})

test_that("valid input parameter present", {
    expect_error(qpcR(dat, groups = "dpi"))
    expect_error(qpcR(dat, hkg = "HKG", treatment = "Control", groups = "dpi"))
    expect_error(qpcR(dat, hkg = "HKG_", groups = "dpi"))
    expect_error(qpcR(dat, hkg = "HKG", groups = "dpI"))
})

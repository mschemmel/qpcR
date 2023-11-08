test_that("qpcR() works", {
    dat <- read.table("example2.tsv", sep = "\t", head = TRUE)
    run <- get_reference_group(qpcR(dat, hkg = "HKG", reference = "control", groups = "dpi"))
    control_groups <- aggregate(rexpr ~ gene + dpi, data = run, FUN = mean)
    expect_equal(unique(control_groups$rexpr), c(1, 1))

    run2 <- get_reference_group(qpcR(dat, hkg = "HKG", reference = "control"))
    control_groups2 <- aggregate(rexpr ~ gene, data = run2, FUN = mean)
    expect_equal(unique(control_groups2$rexpr), c(1, 1))
})
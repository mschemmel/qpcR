test_that("qpcR() works", {
    dat <- read.table("example2.tsv", sep = "\t", head = TRUE)
    run <- get_reference_group(qpcR(dat, hkg = "HKG", reference = "control", groups = "dpi", aggregate = FALSE))
    control_groups <- aggregate(rexpr ~ gene + dpi, data = run, FUN = mean)
    expect_equal(unique(control_groups$rexpr), c(1, 1))

    run2 <- get_reference_group(qpcR(dat, hkg = c("HKG", "HKG2"), reference = "control", groups = "dpi", aggregate = FALSE))
    control_groups2 <- aggregate(rexpr ~ gene + dpi, data = run2, FUN = mean)
    expect_equal(unique(control_groups2$rexpr), 1)

    run3 <- get_reference_group(qpcR(dat, hkg = "HKG", reference = "control", aggregate = FALSE))
    control_groups3 <- aggregate(rexpr ~ gene, data = run3, FUN = mean)
    expect_equal(unique(control_groups3$rexpr), c(1, 1))
})
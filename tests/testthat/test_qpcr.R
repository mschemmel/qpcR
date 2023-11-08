test_that("qpcR() works", {
    dat <- read.table("example2.tsv", sep = "\t", head = TRUE)
    run <- get_reference_group(qpcR(dat, hkg = "HKG", reference = "control", groups = "dpi"))
    control_groups <- as.numeric(tapply(run$rexpr, list(run$gene, run$dpi), FUN = mean))
    expect_equal(unique(control_groups), c(1, 1))

    run2 <- get_reference_group(qpcR(dat, hkg = "HKG", reference = "control", groups = "dpi"))
    control_groups2 <- as.numeric(tapply(run$rexpr, list(run$gene), FUN = mean))
    expect_equal(unique(control_groups2), c(1, 1))
})
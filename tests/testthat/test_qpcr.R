test_that("qpcR() works", {
    dat <- read.table("example2.tsv", sep = "\t", head = TRUE)
    run <- qpcR(dat, hkg = "HKG", reference = "control", groups = "dpi")
    run <- run[run$treatment == "control", ]
    control_groups <- as.numeric(tapply(run$rexpr, list(run$gene, run$dpi), FUN = mean))
    expect_equal(unique(control_groups), c(1, 1))
})
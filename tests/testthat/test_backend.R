dat <- read.table("example2.tsv", sep = "\t", head = TRUE)
test_that("make_groups() works", {
    expect_type(make_groups(dat, groups = NULL), "list")
    expect_type(make_groups(dat, groups = "dpi"), "list")
    expect_error(make_groups(dat, "dose"))
})

test_that("ratio() works", {
    dat <- dat[dat$dpi == 0, ]
    expect_type(ratio(delta_cq(prepare(dat))), "list")
})
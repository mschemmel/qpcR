dat <- read.table("example2.tsv", sep = "\t", head = TRUE)
test_that("make_groups() works", {
    expect_type(make_groups(dat, hkg = NULL, groups = NULL), "list")
    expect_type(make_groups(dat, hkg = NULL, groups = "dpi"), "list")
    expect_error(make_groups(dat, hkg = NULL, groups = "dose"))
    expect_error(make_groups(dat[dat$treatment == "control", ], groups = "dpi"))
})

test_that("comp_gene_pair() works", {
    expect_error(comp_gene_pair(dat[dat$treatment == "control", ], hkg = "HKG"))
    expect_type(comp_gene_pair(dat, hkg = "HKG"), "list")
})

test_that("prepare() works", {
    df <- data.frame(treatment = rep(c("Control", "Infected"), 2, each = 6),
                     gene = rep(c("Kinase","Actin"), 4, each = 3),
                     dpi = c(rep("2dpi", 12), rep("4dpi", 12)),
                     cq = sample(24))

    expect_no_error(prepare(df))
    expect_equal(unique(prepare(df)$efficiency), 100)
    expect_type(prepare(df)$cq, "double")
    expect_type(prepare(df)$efficiency, "double")
    expect_type(prepare(df)$E, "double")

    df$efficiency <- 0.4
    expect_equal(unique(prepare(df)$E), 1.4)

    df$efficiency <- 40
    expect_equal(unique(prepare(df)$E), 1.4)
})

test_that("plot_cpgs_coverage works with valid inputs", {
    skip_on_cran()

    # Create synthetic CpG data with ICR annotation
    set.seed(456)
    n_cpgs <- 100

    df_ICR_cpg <- data.frame(
        probeID = paste0("cg", sprintf("%08d", 1:n_cpgs)),
        cstart = seq(1000000, by = 50, length.out = n_cpgs),
        ICR = sample(c("ICR1", "ICR2", "ICR3"), n_cpgs, replace = TRUE),
        start = 1000000,
        end = 1005000
    )

    result <- plot_cpgs_coverage(
        df_ICR_cpg = df_ICR_cpg,
        bedmeth = "v1"
    )

    expect_type(result, "list")
    expect_named(result, c("plot_counts", "plot_percentage", "data"))
    expect_s3_class(result$plot_counts, "gg")
    expect_s3_class(result$plot_percentage, "gg")
    expect_s3_class(result$data, "data.frame")
})

test_that("plot_cpgs_coverage returns correct data structure", {
    skip_on_cran()

    set.seed(789)
    n_cpgs <- 50

    df_ICR_cpg <- data.frame(
        probeID = paste0("cg", sprintf("%08d", 1:n_cpgs)),
        cstart = seq(1000000, by = 100, length.out = n_cpgs),
        ICR = sample(c("ICR1", "ICR2"), n_cpgs, replace = TRUE),
        start = 1000000,
        end = 1005000
    )

    result <- plot_cpgs_coverage(df_ICR_cpg = df_ICR_cpg, bedmeth = "v1")

    # Check data structure
    expect_true("ICR" %in% colnames(result$data))
    expect_true("n_cpgs" %in% colnames(result$data))
    expect_true("Total_cov" %in% colnames(result$data))
    expect_true("Percentage_cov" %in% colnames(result$data))
})

test_that("plot_cpgs_coverage works with bedmeth v2", {
    skip_on_cran()

    set.seed(123)
    n_cpgs <- 75

    df_ICR_cpg <- data.frame(
        probeID = paste0("cg", sprintf("%08d", 1:n_cpgs)),
        cstart = seq(1000000, by = 50, length.out = n_cpgs),
        ICR = sample(c("ICR1", "ICR2", "ICR3"), n_cpgs, replace = TRUE),
        start = 1000000,
        end = 1005000
    )

    result <- plot_cpgs_coverage(df_ICR_cpg = df_ICR_cpg, bedmeth = "v2")

    expect_type(result, "list")
    expect_s3_class(result$plot_counts, "gg")
    expect_s3_class(result$plot_percentage, "gg")
})

test_that("plot_cpgs_coverage works with 450k array", {
    skip_on_cran()

    set.seed(321)
    n_cpgs <- 60

    df_ICR_cpg <- data.frame(
        probeID = paste0("cg", sprintf("%08d", 1:n_cpgs)),
        cstart = seq(1000000, by = 75, length.out = n_cpgs),
        ICR = sample(c("ICR1", "ICR2"), n_cpgs, replace = TRUE),
        start = 1000000,
        end = 1005000
    )

    result <- plot_cpgs_coverage(df_ICR_cpg = df_ICR_cpg, bedmeth = "450k")

    expect_type(result, "list")
    expect_s3_class(result$plot_counts, "gg")
})

test_that("plot_cpgs_coverage validates bedmeth parameter", {
    df_ICR_cpg <- data.frame(
        probeID = paste0("cg", sprintf("%08d", 1:10)),
        cstart = seq(1000000, by = 100, length.out = 10),
        ICR = rep("ICR1", 10),
        start = 1000000,
        end = 1005000
    )

    expect_error(
        plot_cpgs_coverage(df_ICR_cpg = df_ICR_cpg, bedmeth = "invalid"),
        "Invalid bedmeth version"
    )
})

test_that("plot_cpgs_coverage handles empty data frame", {
    skip_on_cran()

    df_ICR_cpg <- data.frame(
        probeID = character(0),
        cstart = numeric(0),
        ICR = character(0),
        start = numeric(0),
        end = numeric(0)
    )

    result <- plot_cpgs_coverage(df_ICR_cpg = df_ICR_cpg, bedmeth = "v1")

    expect_type(result, "list")
    expect_s3_class(result$data, "data.frame")
})

test_that("plot_cpgs_coverage calculates percentage correctly", {
    skip_on_cran()

    set.seed(111)
    n_cpgs <- 40

    df_ICR_cpg <- data.frame(
        probeID = paste0("cg", sprintf("%08d", 1:n_cpgs)),
        cstart = seq(1000000, by = 50, length.out = n_cpgs),
        ICR = sample(c("ICR1", "ICR2"), n_cpgs, replace = TRUE),
        start = 1000000,
        end = 1005000
    )

    result <- plot_cpgs_coverage(df_ICR_cpg = df_ICR_cpg, bedmeth = "v1")

    # Percentage should be between 0 and 100
    expect_true(all(result$data$Percentage_cov >= 0))
    expect_true(all(result$data$Percentage_cov <= 100))
})

test_that("plot_line_ICR works with valid inputs (static)", {
    skip_on_cran()

    # Create synthetic data
    set.seed(789)
    n_cpgs <- 20
    n_samples <- 6

    # Create significant DMPs data
    significantDMPs <- data.frame(
        cstart = seq(1000000, by = 100, length.out = n_cpgs),
        ICR = rep("KCNQ1OT1:TSS-DMR", n_cpgs),
        start = 1000000,
        end = 1002000,
        P.Value = runif(n_cpgs, 0.001, 0.049),
        stringsAsFactors = FALSE
    )

    # Create ICRcpg data
    meth_values <- matrix(
        runif(n_cpgs * n_samples, 0.3, 0.7),
        nrow = n_cpgs, ncol = n_samples
    )
    colnames(meth_values) <- paste0("Sample_", 1:n_samples)

    ICRcpg <- data.frame(
        meth_values,
        cstart = seq(1000000, by = 100, length.out = n_cpgs),
        ICR = rep("KCNQ1OT1:TSS-DMR", n_cpgs),
        start = 1000000,
        end = 1002000,
        stringsAsFactors = FALSE
    )

    sampleInfo <- c(rep("Control", 3), rep("Case", 3))

    result <- plot_line_ICR(
        significantDMPs = significantDMPs,
        ICRcpg = ICRcpg,
        ICR = "KCNQ1OT1:TSS-DMR",
        sampleInfo = sampleInfo,
        interactive = FALSE
    )

    expect_s3_class(result, "gg")
})

test_that("plot_line_ICR works with interactive plot", {
    skip_on_cran()

    set.seed(123)
    n_cpgs <- 15
    n_samples <- 4

    significantDMPs <- data.frame(
        cstart = seq(1000000, by = 100, length.out = n_cpgs),
        ICR = rep("TestICR", n_cpgs),
        start = 1000000,
        end = 1001500,
        P.Value = runif(n_cpgs, 0.001, 0.049),
        stringsAsFactors = FALSE
    )

    meth_values <- matrix(
        runif(n_cpgs * n_samples, 0.3, 0.7),
        nrow = n_cpgs, ncol = n_samples
    )
    colnames(meth_values) <- paste0("Sample_", 1:n_samples)

    ICRcpg <- data.frame(
        meth_values,
        cstart = seq(1000000, by = 100, length.out = n_cpgs),
        ICR = rep("TestICR", n_cpgs),
        start = 1000000,
        end = 1001500,
        stringsAsFactors = FALSE
    )

    sampleInfo <- c(rep("Control", 2), rep("Case", 2))

    result <- plot_line_ICR(
        significantDMPs = significantDMPs,
        ICRcpg = ICRcpg,
        ICR = "TestICR",
        sampleInfo = sampleInfo,
        interactive = TRUE
    )

    # plotly objects have class "plotly" and "htmlwidget"
    expect_true(any(class(result) %in% c("plotly", "htmlwidget")))
})

test_that("plot_line_ICR validates significantDMPs input", {
    # Empty data frame
    expect_error(
        plot_line_ICR(
            significantDMPs = data.frame(),
            ICRcpg = data.frame(),
            ICR = "test",
            sampleInfo = c("Control", "Case")
        ),
        "significantDMPs must be a non-empty data frame"
    )

    # Not a data frame
    expect_error(
        plot_line_ICR(
            significantDMPs = matrix(1:10, ncol = 2),
            ICRcpg = data.frame(x = 1),
            ICR = "test",
            sampleInfo = c("Control", "Case")
        ),
        "significantDMPs must be a non-empty data frame"
    )
})

test_that("plot_line_ICR validates ICRcpg input", {
    significantDMPs <- data.frame(
        cstart = 1:10,
        ICR = rep("test", 10),
        start = 1000,
        end = 2000
    )

    # Empty ICRcpg
    expect_error(
        plot_line_ICR(
            significantDMPs = significantDMPs,
            ICRcpg = data.frame(),
            ICR = "test",
            sampleInfo = c("Control", "Case")
        ),
        "ICRcpg must be a non-empty data frame"
    )
})

test_that("plot_line_ICR validates ICR parameter", {
    significantDMPs <- data.frame(
        cstart = 1:10,
        ICR = rep("ValidICR", 10),
        start = 1000,
        end = 2000
    )

    ICRcpg <- data.frame(
        Sample_1 = runif(10),
        cstart = 1:10,
        ICR = rep("ValidICR", 10),
        start = 1000,
        end = 2000
    )

    expect_error(
        plot_line_ICR(
            significantDMPs = significantDMPs,
            ICRcpg = ICRcpg,
            ICR = "InvalidICR",
            sampleInfo = c("Control")
        ),
        "ICR.*not found in significantDMPs"
    )
})

test_that("plot_line_ICR validates sampleInfo length", {
    n_cpgs <- 10
    n_samples <- 4

    significantDMPs <- data.frame(
        cstart = 1:n_cpgs,
        ICR = rep("TestICR", n_cpgs),
        start = 1000,
        end = 2000
    )

    meth_values <- matrix(runif(n_cpgs * n_samples), nrow = n_cpgs, ncol = n_samples)
    colnames(meth_values) <- paste0("Sample_", 1:n_samples)

    ICRcpg <- data.frame(
        meth_values,
        cstart = 1:n_cpgs,
        ICR = rep("TestICR", n_cpgs),
        start = 1000,
        end = 2000
    )

    # Wrong sampleInfo length
    expect_error(
        plot_line_ICR(
            significantDMPs = significantDMPs,
            ICRcpg = ICRcpg,
            ICR = "TestICR",
            sampleInfo = c("Control", "Case") # Only 2, need 4
        ),
        "sampleInfo length.*doesn't match"
    )
})

test_that("plot_line_ICR handles missing required columns", {
    significantDMPs <- data.frame(
        cstart = 1:10,
        ICR = rep("TestICR", 10),
        start = 1000,
        end = 2000
    )

    # Missing cstart column
    ICRcpg <- data.frame(
        Sample_1 = runif(10),
        ICR = rep("TestICR", 10),
        start = 1000,
        end = 2000
    )

    expect_error(
        plot_line_ICR(
            significantDMPs = significantDMPs,
            ICRcpg = ICRcpg,
            ICR = "TestICR",
            sampleInfo = c("Control")
        ),
        "Missing required columns in ICRcpg"
    )
})

test_that("plot_line_ICR handles data with NAs", {
    skip_on_cran()

    set.seed(456)
    n_cpgs <- 20
    n_samples <- 4

    significantDMPs <- data.frame(
        cstart = seq(1000000, by = 100, length.out = n_cpgs),
        ICR = rep("TestICR", n_cpgs),
        start = 1000000,
        end = 1002000,
        P.Value = runif(n_cpgs, 0.001, 0.049),
        stringsAsFactors = FALSE
    )

    meth_values <- matrix(
        runif(n_cpgs * n_samples, 0.3, 0.7),
        nrow = n_cpgs, ncol = n_samples
    )
    # Add some NAs
    meth_values[sample(length(meth_values), 5)] <- NA
    colnames(meth_values) <- paste0("Sample_", 1:n_samples)

    ICRcpg <- data.frame(
        meth_values,
        cstart = seq(1000000, by = 100, length.out = n_cpgs),
        ICR = rep("TestICR", n_cpgs),
        start = 1000000,
        end = 1002000,
        stringsAsFactors = FALSE
    )

    sampleInfo <- c(rep("Control", 2), rep("Case", 2))

    result <- plot_line_ICR(
        significantDMPs = significantDMPs,
        ICRcpg = ICRcpg,
        ICR = "TestICR",
        sampleInfo = sampleInfo,
        interactive = FALSE
    )

    expect_s3_class(result, "gg")
})

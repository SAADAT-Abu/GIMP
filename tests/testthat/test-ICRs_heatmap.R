test_that("ICRs_heatmap works with valid inputs", {
    skip_on_cran()

    # Use real ICR names from DMRs.hg19
    data(DMRs.hg19)
    real_icr_names <- head(DMRs.hg19$ICR, 10)

    set.seed(42)
    n_icrs <- length(real_icr_names)
    n_samples <- 8

    df_ICR <- matrix(runif(n_icrs * n_samples, 0.2, 0.8),
        nrow = n_icrs, ncol = n_samples
    )
    rownames(df_ICR) <- real_icr_names
    colnames(df_ICR) <- paste0("Sample_", 1:n_samples)
    df_ICR <- as.data.frame(df_ICR)

    sampleInfo <- c(rep("Control", 4), rep("Case", 4))

    result <- ICRs_heatmap(
        df_ICR = df_ICR,
        sampleInfo = sampleInfo,
        plot_type = "beta"
    )

    expect_s3_class(result, "gg")
})

test_that("ICRs_heatmap handles delta plot type", {
    skip_on_cran()

    # Use real ICR names from DMRs.hg19
    data(DMRs.hg19)
    real_icr_names <- head(DMRs.hg19$ICR, 8)

    set.seed(123)
    n_icrs <- length(real_icr_names)
    n_samples <- 6

    df_ICR <- matrix(runif(n_icrs * n_samples, 0.3, 0.7),
        nrow = n_icrs, ncol = n_samples
    )
    rownames(df_ICR) <- real_icr_names
    colnames(df_ICR) <- paste0("Sample_", 1:n_samples)
    df_ICR <- as.data.frame(df_ICR)

    sampleInfo <- c(rep("Control", 3), rep("Case", 3))

    result <- ICRs_heatmap(
        df_ICR = df_ICR,
        sampleInfo = sampleInfo,
        plot_type = "delta"
    )

    expect_s3_class(result, "gg")
})

test_that("ICRs_heatmap handles defect plot type", {
    skip_on_cran()

    # Use real ICR names from DMRs.hg19
    data(DMRs.hg19)
    real_icr_names <- head(DMRs.hg19$ICR, 6)

    set.seed(456)
    n_icrs <- length(real_icr_names)
    n_samples <- 8

    df_ICR <- matrix(runif(n_icrs * n_samples, 0.2, 0.8),
        nrow = n_icrs, ncol = n_samples
    )
    rownames(df_ICR) <- real_icr_names
    colnames(df_ICR) <- paste0("Sample_", 1:n_samples)
    df_ICR <- as.data.frame(df_ICR)

    sampleInfo <- c(rep("Control", 4), rep("Case", 4))

    result <- ICRs_heatmap(
        df_ICR = df_ICR,
        sampleInfo = sampleInfo,
        plot_type = "defect",
        sd_threshold = 3
    )

    expect_s3_class(result, "gg")
})

test_that("ICRs_heatmap validates sampleInfo length", {
    df_ICR <- matrix(runif(5 * 4, 0.3, 0.7), nrow = 5, ncol = 4)
    rownames(df_ICR) <- paste0("ICR_", 1:5)
    colnames(df_ICR) <- paste0("Sample_", 1:4)
    df_ICR <- as.data.frame(df_ICR)

    # Wrong length sampleInfo
    sampleInfo <- c(rep("Control", 2), rep("Case", 1)) # Only 3, need 4

    expect_error(
        ICRs_heatmap(df_ICR = df_ICR, sampleInfo = sampleInfo),
        "Length of 'sampleInfo' must match"
    )
})

test_that("ICRs_heatmap validates plot_type parameter", {
    df_ICR <- matrix(runif(5 * 4, 0.3, 0.7), nrow = 5, ncol = 4)
    rownames(df_ICR) <- paste0("ICR_", 1:5)
    colnames(df_ICR) <- paste0("Sample_", 1:4)
    df_ICR <- as.data.frame(df_ICR)

    sampleInfo <- c(rep("Control", 2), rep("Case", 2))

    expect_error(
        ICRs_heatmap(df_ICR = df_ICR, sampleInfo = sampleInfo, plot_type = "invalid"),
        "Invalid plot_type"
    )
})

test_that("ICRs_heatmap validates order_by parameter", {
    df_ICR <- matrix(runif(5 * 4, 0.3, 0.7), nrow = 5, ncol = 4)
    rownames(df_ICR) <- paste0("ICR_", 1:5)
    colnames(df_ICR) <- paste0("Sample_", 1:4)
    df_ICR <- as.data.frame(df_ICR)

    sampleInfo <- c(rep("Control", 2), rep("Case", 2))

    expect_error(
        ICRs_heatmap(
            df_ICR = df_ICR, sampleInfo = sampleInfo,
            order_by = "invalid"
        ),
        "Invalid order_by"
    )
})

test_that("ICRs_heatmap validates bedmeth parameter", {
    df_ICR <- matrix(runif(5 * 4, 0.3, 0.7), nrow = 5, ncol = 4)
    rownames(df_ICR) <- paste0("ICR_", 1:5)
    colnames(df_ICR) <- paste0("Sample_", 1:4)
    df_ICR <- as.data.frame(df_ICR)

    sampleInfo <- c(rep("Control", 2), rep("Case", 2))

    expect_error(
        ICRs_heatmap(df_ICR = df_ICR, sampleInfo = sampleInfo, bedmeth = "invalid"),
        "Invalid bedmeth version"
    )
})

test_that("ICRs_heatmap handles empty data frame", {
    df_ICR <- data.frame()
    sampleInfo <- character(0)

    expect_error(
        ICRs_heatmap(df_ICR = df_ICR, sampleInfo = sampleInfo),
        "df_ICR is empty"
    )
})

test_that("ICRs_heatmap works with bedmeth v2", {
    skip_on_cran()

    # Use real ICR names from DMRs.hg38 for v2
    data(DMRs.hg38)
    real_icr_names <- head(DMRs.hg38$ICR, 5)

    set.seed(789)
    n_icrs <- length(real_icr_names)
    n_samples <- 4

    df_ICR <- matrix(runif(n_icrs * n_samples, 0.3, 0.7),
        nrow = n_icrs, ncol = n_samples
    )
    rownames(df_ICR) <- real_icr_names
    colnames(df_ICR) <- paste0("Sample_", 1:n_samples)
    df_ICR <- as.data.frame(df_ICR)

    sampleInfo <- c(rep("Control", 2), rep("Case", 2))

    result <- ICRs_heatmap(
        df_ICR = df_ICR,
        sampleInfo = sampleInfo,
        bedmeth = "v2"
    )

    expect_s3_class(result, "gg")
})

test_that("ICRs_heatmap works with custom labels", {
    skip_on_cran()

    # Use real ICR names from DMRs.hg19
    data(DMRs.hg19)
    real_icr_names <- head(DMRs.hg19$ICR, 5)

    set.seed(321)
    n_icrs <- length(real_icr_names)
    n_samples <- 6

    df_ICR <- matrix(runif(n_icrs * n_samples, 0.3, 0.7),
        nrow = n_icrs, ncol = n_samples
    )
    rownames(df_ICR) <- real_icr_names
    colnames(df_ICR) <- paste0("Sample_", 1:n_samples)
    df_ICR <- as.data.frame(df_ICR)

    sampleInfo <- c(rep("Healthy", 3), rep("Disease", 3))

    result <- ICRs_heatmap(
        df_ICR = df_ICR,
        sampleInfo = sampleInfo,
        control_label = "Healthy",
        case_label = "Disease"
    )

    expect_s3_class(result, "gg")
})

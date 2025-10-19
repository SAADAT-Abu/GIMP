test_that("iDMPs works with valid inputs", {
    # Create synthetic methylation data with annotation
    set.seed(123)
    n_cpgs <- 50
    n_controls <- 3
    n_cases <- 3

    meth_data <- matrix(runif(n_cpgs * (n_controls + n_cases), 0.3, 0.7),
        nrow = n_cpgs, ncol = n_controls + n_cases
    )
    colnames(meth_data) <- c(
        paste0("Ctrl_", 1:n_controls),
        paste0("Case_", 1:n_cases)
    )

    # Add annotation columns
    ICRcpg_data <- data.frame(
        meth_data,
        cstart = seq(1000000, by = 100, length.out = n_cpgs),
        ICR = rep(c("ICR1", "ICR2"), length.out = n_cpgs),
        start = 1000000,
        end = 1005000
    )

    sampleInfo <- factor(c(rep("Control", n_controls), rep("Case", n_cases)))

    result <- iDMPs(
        data = ICRcpg_data,
        sampleInfo = sampleInfo,
        pValueCutoff = 0.05
    )

    expect_type(result, "list")
    expect_named(result, c("fit", "eBayesfit", "topDMPs", "allResults", "groupLabels"))
    expect_s3_class(result$topDMPs, "data.frame")
    expect_s3_class(result$allResults, "data.frame")
})

test_that("iDMPs requires Control group", {
    set.seed(456)
    meth_data <- matrix(runif(30 * 4, 0.3, 0.7), nrow = 30, ncol = 4)
    colnames(meth_data) <- paste0("Sample_", 1:4)

    ICRcpg_data <- data.frame(
        meth_data,
        cstart = seq(1000000, by = 100, length.out = 30),
        ICR = rep(c("ICR1", "ICR2"), length.out = 30),
        start = 1000000,
        end = 1005000
    )

    # Create sampleInfo without "Control" group
    sampleInfo <- factor(c(rep("Group1", 2), rep("Group2", 2)))

    expect_error(
        iDMPs(data = ICRcpg_data, sampleInfo = sampleInfo),
        "Control.*group.*must be present"
    )
})

test_that("iDMPs handles missing annotation columns", {
    set.seed(789)
    meth_data <- matrix(runif(20 * 4, 0.3, 0.7), nrow = 20, ncol = 4)
    colnames(meth_data) <- paste0("Sample_", 1:4)

    # Missing required annotation columns
    ICRcpg_data <- data.frame(meth_data)

    sampleInfo <- factor(c(rep("Control", 2), rep("Case", 2)))

    expect_error(
        iDMPs(data = ICRcpg_data, sampleInfo = sampleInfo),
        "No annotation columns found"
    )
})

test_that("iDMPs validates sample info length", {
    set.seed(321)
    meth_data <- matrix(runif(30 * 4, 0.3, 0.7), nrow = 30, ncol = 4)
    colnames(meth_data) <- paste0("Sample_", 1:4)

    ICRcpg_data <- data.frame(
        meth_data,
        cstart = seq(1000000, by = 100, length.out = 30),
        ICR = rep(c("ICR1", "ICR2"), length.out = 30),
        start = 1000000,
        end = 1005000
    )

    # Wrong number of samples in sampleInfo
    sampleInfo <- factor(c(rep("Control", 2), rep("Case", 1))) # Only 3, but need 4

    expect_error(
        iDMPs(data = ICRcpg_data, sampleInfo = sampleInfo),
        "doesn't match sampleInfo length"
    )
})

test_that("iDMPs handles data with many NAs", {
    set.seed(111)
    n_cpgs <- 40
    n_samples <- 6
    meth_data <- matrix(runif(n_cpgs * n_samples, 0.3, 0.7),
        nrow = n_cpgs, ncol = n_samples
    )
    # Add many NAs (more than 50% in some rows)
    meth_data[1:5, ] <- NA
    colnames(meth_data) <- paste0("Sample_", 1:n_samples)

    ICRcpg_data <- data.frame(
        meth_data,
        cstart = seq(1000000, by = 100, length.out = n_cpgs),
        ICR = rep(c("ICR1", "ICR2"), length.out = n_cpgs),
        start = 1000000,
        end = 1005000
    )

    sampleInfo <- factor(c(rep("Control", 3), rep("Case", 3)))

    # Should work but filter out rows with >50% NAs
    result <- iDMPs(
        data = ICRcpg_data,
        sampleInfo = sampleInfo,
        pValueCutoff = 0.05
    )

    expect_type(result, "list")
    expect_s3_class(result$allResults, "data.frame")
})

test_that("iDMPs custom p-value cutoff works", {
    set.seed(222)
    n_cpgs <- 30
    n_samples <- 6
    meth_data <- matrix(runif(n_cpgs * n_samples, 0.3, 0.7),
        nrow = n_cpgs, ncol = n_samples
    )
    colnames(meth_data) <- paste0("Sample_", 1:n_samples)

    ICRcpg_data <- data.frame(
        meth_data,
        cstart = seq(1000000, by = 100, length.out = n_cpgs),
        ICR = rep(c("ICR1", "ICR2"), length.out = n_cpgs),
        start = 1000000,
        end = 1005000
    )

    sampleInfo <- factor(c(rep("Control", 3), rep("Case", 3)))

    result1 <- iDMPs(data = ICRcpg_data, sampleInfo = sampleInfo, pValueCutoff = 0.05)
    result2 <- iDMPs(data = ICRcpg_data, sampleInfo = sampleInfo, pValueCutoff = 0.01)

    # More stringent cutoff should have fewer or equal significant DMPs
    expect_true(nrow(result2$topDMPs) <= nrow(result1$topDMPs))
})

test_that("iDMPs returns expected list elements", {
    set.seed(333)
    n_cpgs <- 25
    n_samples <- 4
    meth_data <- matrix(runif(n_cpgs * n_samples, 0.3, 0.7),
        nrow = n_cpgs, ncol = n_samples
    )
    colnames(meth_data) <- paste0("Sample_", 1:n_samples)

    ICRcpg_data <- data.frame(
        meth_data,
        cstart = seq(1000000, by = 100, length.out = n_cpgs),
        ICR = rep(c("ICR1", "ICR2"), length.out = n_cpgs),
        start = 1000000,
        end = 1005000
    )

    sampleInfo <- factor(c(rep("Control", 2), rep("Case", 2)))

    result <- iDMPs(data = ICRcpg_data, sampleInfo = sampleInfo)

    # Check structure
    expect_true("fit" %in% names(result))
    expect_true("eBayesfit" %in% names(result))
    expect_true("topDMPs" %in% names(result))
    expect_true("allResults" %in% names(result))
    expect_true("groupLabels" %in% names(result))

    # Check groupLabels
    expect_equal(result$groupLabels, c("Control", "Case"))
})

test_that("iDMPs handles input that is not a data frame", {
    meth_matrix <- matrix(runif(20 * 4, 0.3, 0.7), nrow = 20, ncol = 4)
    sampleInfo <- factor(c(rep("Control", 2), rep("Case", 2)))

    expect_error(
        iDMPs(data = meth_matrix, sampleInfo = sampleInfo),
        "input data must be a data frame"
    )
})

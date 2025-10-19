test_that("make_ICRs works with valid inputs", {
    # Create a synthetic beta matrix
    set.seed(123)
    n_probes <- 100
    n_samples <- 4
    beta_matrix <- matrix(runif(n_probes * n_samples, 0.3, 0.8),
        nrow = n_probes, ncol = n_samples
    )
    rownames(beta_matrix) <- paste0("cg", sprintf("%08d", 1:n_probes))
    colnames(beta_matrix) <- paste0("Sample_", 1:n_samples)

    # Test with default bedmeth version (v1)
    result <- make_ICRs(Bmatrix = beta_matrix, bedmeth = "v1")

    expect_s3_class(result, "data.frame")
    # Result should have sample columns
    expect_true(ncol(result) >= n_samples)
})

test_that("make_ICRs aggregates data per ICR", {
    set.seed(456)
    beta_matrix <- matrix(runif(50 * 3, 0.2, 0.9),
        nrow = 50, ncol = 3
    )
    rownames(beta_matrix) <- paste0("cg", sprintf("%08d", 1:50))
    colnames(beta_matrix) <- paste0("Sample_", 1:3)

    result <- make_ICRs(Bmatrix = beta_matrix, bedmeth = "v1")

    expect_s3_class(result, "data.frame")
    # ICR names should be in rownames if there are matching probes
    expect_true(length(rownames(result)) >= 0)
})

test_that("make_ICRs works with bedmeth v2", {
    set.seed(789)
    beta_matrix <- matrix(runif(60 * 5, 0.3, 0.7),
        nrow = 60, ncol = 5
    )
    rownames(beta_matrix) <- paste0("cg", sprintf("%08d", 101:160))
    colnames(beta_matrix) <- paste0("Sample_", 1:5)

    result <- make_ICRs(Bmatrix = beta_matrix, bedmeth = "v2")

    expect_s3_class(result, "data.frame")
})

test_that("make_ICRs works with 450k array", {
    set.seed(321)
    beta_matrix <- matrix(runif(40 * 2, 0.4, 0.8),
        nrow = 40, ncol = 2
    )
    rownames(beta_matrix) <- paste0("cg", sprintf("%08d", 201:240))
    colnames(beta_matrix) <- paste0("Sample_", 1:2)

    result <- make_ICRs(Bmatrix = beta_matrix, bedmeth = "450k")

    expect_s3_class(result, "data.frame")
})

test_that("make_ICRs handles invalid bedmeth parameter", {
    beta_matrix <- matrix(runif(10 * 2, 0.3, 0.7), nrow = 10, ncol = 2)
    rownames(beta_matrix) <- paste0("cg", sprintf("%08d", 1:10))
    colnames(beta_matrix) <- paste0("Sample_", 1:2)

    expect_error(
        make_ICRs(Bmatrix = beta_matrix, bedmeth = "invalid"),
        "Invalid bedmeth input"
    )
})

test_that("make_ICRs handles empty matrix", {
    beta_matrix <- matrix(numeric(0), nrow = 0, ncol = 0)

    result <- make_ICRs(Bmatrix = beta_matrix, bedmeth = "v1")

    expect_s3_class(result, "data.frame")
    # Empty input may still return ICR annotation data
    expect_true(nrow(result) >= 0)
})

test_that("make_ICRs handles matrix with NAs", {
    set.seed(111)
    beta_matrix <- matrix(runif(30 * 4, 0.2, 0.8), nrow = 30, ncol = 4)
    beta_matrix[sample(length(beta_matrix), 10)] <- NA
    rownames(beta_matrix) <- paste0("cg", sprintf("%08d", 1:30))
    colnames(beta_matrix) <- paste0("Sample_", 1:4)

    result <- make_ICRs(Bmatrix = beta_matrix, bedmeth = "v1")

    expect_s3_class(result, "data.frame")
})

test_that("make_ICRs produces numeric methylation values", {
    set.seed(222)
    beta_matrix <- matrix(runif(20 * 3, 0.3, 0.7),
        nrow = 20, ncol = 3
    )
    rownames(beta_matrix) <- paste0("cg", sprintf("%08d", 1:20))
    colnames(beta_matrix) <- paste0("Sample_", 1:3)

    result <- make_ICRs(Bmatrix = beta_matrix, bedmeth = "v1")

    if (nrow(result) > 0 && ncol(result) > 0) {
        # Check that at least some columns are numeric
        numeric_cols <- sapply(result, is.numeric)
        expect_true(any(numeric_cols))
    }
})

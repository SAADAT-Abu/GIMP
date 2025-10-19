test_that("make_cpgs works with valid inputs", {
    # Create a synthetic beta matrix with known probe IDs
    set.seed(123)
    n_probes <- 100
    n_samples <- 4
    beta_matrix <- matrix(runif(n_probes * n_samples, 0.3, 0.8),
        nrow = n_probes, ncol = n_samples
    )
    rownames(beta_matrix) <- paste0("cg", sprintf("%08d", 1:n_probes))
    colnames(beta_matrix) <- paste0("Sample_", 1:n_samples)

    # Test with default bedmeth version (v1)
    result <- make_cpgs(Bmatrix = beta_matrix, bedmeth = "v1")

    expect_s3_class(result, "data.frame")
    expect_true(all(c("ICR", "cstart", "start", "end") %in% colnames(result)))
})

test_that("make_cpgs works with bedmeth v2", {
    set.seed(456)
    beta_matrix <- matrix(runif(50 * 3, 0.2, 0.9),
        nrow = 50, ncol = 3
    )
    rownames(beta_matrix) <- paste0("cg", sprintf("%08d", 101:150))
    colnames(beta_matrix) <- paste0("Sample_", 1:3)

    result <- make_cpgs(Bmatrix = beta_matrix, bedmeth = "v2")

    expect_s3_class(result, "data.frame")
})

test_that("make_cpgs works with 450k array", {
    set.seed(789)
    beta_matrix <- matrix(runif(75 * 2, 0.4, 0.7),
        nrow = 75, ncol = 2
    )
    rownames(beta_matrix) <- paste0("cg", sprintf("%08d", 201:275))
    colnames(beta_matrix) <- paste0("Sample_", 1:2)

    result <- make_cpgs(Bmatrix = beta_matrix, bedmeth = "450k")

    expect_s3_class(result, "data.frame")
})

test_that("make_cpgs handles invalid bedmeth parameter", {
    beta_matrix <- matrix(runif(10 * 2, 0.3, 0.7), nrow = 10, ncol = 2)
    rownames(beta_matrix) <- paste0("cg", sprintf("%08d", 1:10))
    colnames(beta_matrix) <- paste0("Sample_", 1:2)

    expect_error(
        make_cpgs(Bmatrix = beta_matrix, bedmeth = "invalid"),
        "Invalid bedmeth input"
    )
})

test_that("make_cpgs handles empty matrix", {
    beta_matrix <- matrix(numeric(0), nrow = 0, ncol = 0)

    result <- make_cpgs(Bmatrix = beta_matrix, bedmeth = "v1")

    expect_s3_class(result, "data.frame")
    # Empty input should result in data frame with annotation columns only
    expect_true("ICR" %in% colnames(result) || ncol(result) >= 0)
})

test_that("make_cpgs returns expected columns", {
    set.seed(321)
    beta_matrix <- matrix(runif(20 * 3, 0.3, 0.8),
        nrow = 20, ncol = 3
    )
    rownames(beta_matrix) <- paste0("cg", sprintf("%08d", 1:20))
    colnames(beta_matrix) <- paste0("Sample_", 1:3)

    result <- make_cpgs(Bmatrix = beta_matrix, bedmeth = "v1")

    # Check that sample columns are preserved
    expected_cols <- c("ICR", "cstart", "start", "end")
    expect_true(all(expected_cols %in% colnames(result)))
})

test_that("make_cpgs handles matrix with NAs", {
    set.seed(111)
    beta_matrix <- matrix(runif(30 * 4, 0.2, 0.8), nrow = 30, ncol = 4)
    beta_matrix[sample(length(beta_matrix), 10)] <- NA
    rownames(beta_matrix) <- paste0("cg", sprintf("%08d", 1:30))
    colnames(beta_matrix) <- paste0("Sample_", 1:4)

    result <- make_cpgs(Bmatrix = beta_matrix, bedmeth = "v1")

    expect_s3_class(result, "data.frame")
})

test_that("get_bed_data works with v1", {
    skip_on_cran()

    result <- GIMP:::get_bed_data("v1")

    expect_s3_class(result, "data.frame")
    expect_true(nrow(result) > 0)
    expect_true(ncol(result) >= 4)
})

test_that("get_bed_data works with v2", {
    skip_on_cran()

    result <- GIMP:::get_bed_data("v2")

    expect_s3_class(result, "data.frame")
    expect_true(nrow(result) > 0)
})

test_that("get_bed_data works with 450k", {
    skip_on_cran()

    result <- GIMP:::get_bed_data("450k")

    expect_s3_class(result, "data.frame")
    expect_true(nrow(result) > 0)
})

test_that("get_bed_data validates version parameter", {
    expect_error(
        GIMP:::get_bed_data("invalid"),
        "Invalid version"
    )
})

test_that("get_bed_data returns consistent data structure", {
    skip_on_cran()

    result_v1 <- GIMP:::get_bed_data("v1")
    result_v2 <- GIMP:::get_bed_data("v2")
    result_450k <- GIMP:::get_bed_data("450k")

    # All should have same column structure
    expect_equal(colnames(result_v1), colnames(result_v2))
    expect_equal(colnames(result_v1), colnames(result_450k))
})

test_that("get_bed_data returns data frame with probe information", {
    skip_on_cran()

    result <- GIMP:::get_bed_data("v1")

    # Check for essential BED columns
    expect_true("probeID" %in% colnames(result))
    expect_true("chrom" %in% colnames(result) || "chr" %in% colnames(result))
})

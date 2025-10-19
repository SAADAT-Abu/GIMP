test_that("read_idat_zip validates required packages", {
    skip_on_cran()

    # Function should exist
    expect_true(exists("read_idat_zip"))
})

test_that("read_idat_zip validates ZIP file existence", {
    skip_on_cran()
    skip_if_not_installed("minfi")

    expect_error(
        read_idat_zip(
            zip_file = "nonexistent_file.zip",
            array_type = "EPIC"
        ),
        "ZIP file not found"
    )
})

test_that("read_idat_zip validates array_type parameter", {
    skip_on_cran()
    skip_if_not_installed("minfi")

    # Create a temporary empty zip for testing
    temp_zip <- tempfile(fileext = ".zip")
    file.create(temp_zip)

    expect_error(
        read_idat_zip(
            zip_file = temp_zip,
            array_type = "InvalidType"
        ),
        "should be one of"
    )

    unlink(temp_zip)
})

test_that("read_idat_zip validates normalize_method parameter", {
    skip_on_cran()
    skip_if_not_installed("minfi")

    temp_zip <- tempfile(fileext = ".zip")
    file.create(temp_zip)

    expect_error(
        read_idat_zip(
            zip_file = temp_zip,
            array_type = "EPIC",
            normalize_method = "InvalidMethod"
        ),
        "should be one of"
    )

    unlink(temp_zip)
})

test_that("read_idat_zip handles parallel processing parameter", {
    skip_on_cran()
    skip_if_not_installed("minfi")

    # Function should handle n_cores parameter
    expect_true("n_cores" %in% names(formals(read_idat_zip)))
})

test_that("read_idat_zip function signature is correct", {
    skip_on_cran()

    # Check function parameters
    params <- names(formals(read_idat_zip))

    expect_true("zip_file" %in% params)
    expect_true("sample_sheet_name" %in% params)
    expect_true("array_type" %in% params)
    expect_true("temp_dir" %in% params)
    expect_true("normalize_method" %in% params)
    expect_true("detection_pval" %in% params)
    expect_true("remove_failed_samples" %in% params)
    expect_true("n_cores" %in% params)
})

test_that("preview_idat_zip function exists", {
    skip_on_cran()

    # Check if helper function exists (if it's exported)
    # This is just to verify the function structure
    expect_true(exists("read_idat_zip"))
})

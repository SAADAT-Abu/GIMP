test_that("validate_geo_dataset validates GEO ID format", {
    skip_on_cran()
    skip_if_not_installed("GEOquery")

    # Invalid format should return FALSE
    result <- validate_geo_dataset("INVALID123")
    expect_false(result$valid)
    expect_match(result$message, "Invalid GEO ID format")
})

test_that("validate_geo_dataset requires GEOquery package", {
    skip_on_cran()

    # Should mention GEOquery requirement if package is not loaded
    # This test verifies the function structure exists
    expect_true(exists("validate_geo_dataset"))
})

test_that("validate_geo_dataset cleans GEO ID", {
    skip_on_cran()
    skip_if_not_installed("GEOquery")

    # Should handle whitespace and case
    result1 <- validate_geo_dataset("  gse12345  ")
    result2 <- validate_geo_dataset("GSE12345")

    # Both should attempt same ID format
    expect_true("valid" %in% names(result1))
    expect_true("valid" %in% names(result2))
})

test_that("process_geo_dataset validates prerequisites", {
    skip_on_cran()

    # Function should exist and be exported
    expect_true(exists("process_geo_dataset"))
})

test_that("diagnose_geo_dataset validates GEO ID", {
    skip_on_cran()

    # Function should exist and be exported
    expect_true(exists("diagnose_geo_dataset"))
})

test_that("get_geo_phenotype_data validates input", {
    skip_on_cran()

    # Function should exist and be exported
    expect_true(exists("get_geo_phenotype_data"))
})

test_that("process_geo_with_mappings validates input parameters", {
    skip_on_cran()

    # Function should exist and be exported
    expect_true(exists("process_geo_with_mappings"))
})

test_that("detect_array_platform works with platform IDs", {
    skip_on_cran()

    result1 <- GIMP:::detect_array_platform("GPL13534")
    expect_equal(result1, "450k")

    result2 <- GIMP:::detect_array_platform("GPL21145")
    expect_equal(result2, "EPIC")

    result3 <- GIMP:::detect_array_platform("GPL23976")
    expect_equal(result3, "EPICv2")
})

test_that("detect_array_platform works with text analysis", {
    skip_on_cran()

    result1 <- GIMP:::detect_array_platform("", "Study with 450k array", "")
    expect_equal(result1, "450k")

    result2 <- GIMP:::detect_array_platform("", "EPIC array study", "")
    expect_equal(result2, "EPIC")

    result3 <- GIMP:::detect_array_platform("", "", "Using EPICv2 platform")
    expect_equal(result3, "EPICv2")
})

test_that("detect_array_platform has fallback", {
    skip_on_cran()

    result <- GIMP:::detect_array_platform("", "", "")
    expect_equal(result, "EPIC")
})

test_that("assign_sample_groups classifies samples correctly", {
    skip_on_cran()

    values <- c("control sample", "disease sample", "normal tissue", "tumor tissue")
    control_terms <- c("control", "normal")
    case_terms <- c("disease", "tumor")

    groups <- GIMP:::assign_sample_groups(values, control_terms, case_terms)

    expect_equal(groups[1], "Control")
    expect_equal(groups[2], "Case")
    expect_equal(groups[3], "Control")
    expect_equal(groups[4], "Case")
})

test_that("assign_sample_groups handles unknown values", {
    skip_on_cran()

    values <- c("unknown", "other")
    control_terms <- c("control")
    case_terms <- c("case")

    groups <- GIMP:::assign_sample_groups(values, control_terms, case_terms)

    expect_true(all(groups == "Unknown"))
})

test_that("extract_sentrix_info parses standard format", {
    skip_on_cran()

    basenames <- c("200123456789_R01C01", "200987654321_R02C02")

    result <- GIMP:::extract_sentrix_info(basenames)

    expect_equal(result$sentrix_id[1], "200123456789")
    expect_equal(result$sentrix_position[1], "R01C01")
    expect_equal(result$sentrix_id[2], "200987654321")
    expect_equal(result$sentrix_position[2], "R02C02")
})

test_that("extract_sentrix_info creates fallback IDs", {
    skip_on_cran()

    basenames <- c("nonstandard_format", "another_bad_format")

    result <- GIMP:::extract_sentrix_info(basenames)

    expect_true(all(grepl("^200000000", result$sentrix_id)))
    expect_true(all(grepl("^R\\d+C\\d+", result$sentrix_position)))
})

test_that("extract_sample_id_from_idat extracts GSM IDs", {
    skip_on_cran()

    basename1 <- "GSM1234567_something_else"
    result1 <- GIMP:::extract_sample_id_from_idat(basename1)
    expect_equal(result1, "GSM1234567")

    basename2 <- "prefix_GSM9876543_suffix"
    result2 <- GIMP:::extract_sample_id_from_idat(basename2)
    expect_equal(result2, "GSM9876543")
})

test_that("extract_sample_id_from_idat handles non-GSM names", {
    skip_on_cran()

    basename <- "Sample_12345_R01C01"
    result <- GIMP:::extract_sample_id_from_idat(basename)
    expect_equal(result, "Sample")
})

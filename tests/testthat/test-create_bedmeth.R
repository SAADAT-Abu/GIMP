test_that("create_bedmeth works with v1", {
    skip_on_cran()
    skip_if_not_installed("minfi")
    skip_if_not_installed("IlluminaHumanMethylationEPICanno.ilm10b4.hg19")

    result <- create_bedmeth(version = "v1")

    expect_s3_class(result, "data.frame")
    expect_named(result, c("chrom", "start", "end", "probeID"))
    expect_true(nrow(result) > 0)
    expect_true(all(c("chrom", "start", "end", "probeID") %in% colnames(result)))
})

test_that("create_bedmeth works with v2", {
    skip_on_cran()
    skip_if_not_installed("minfi")
    skip_if_not_installed("IlluminaHumanMethylationEPICanno.ilm10b4.hg19")

    result <- create_bedmeth(version = "v2")

    expect_s3_class(result, "data.frame")
    expect_named(result, c("chrom", "start", "end", "probeID"))
    expect_true(nrow(result) > 0)
})

test_that("create_bedmeth works with 450k", {
    skip_on_cran()
    skip_if_not_installed("minfi")
    skip_if_not_installed("IlluminaHumanMethylation450kanno.ilmn12.hg19")

    result <- create_bedmeth(version = "450k")

    expect_s3_class(result, "data.frame")
    expect_named(result, c("chrom", "start", "end", "probeID"))
    expect_true(nrow(result) > 0)
})

test_that("create_bedmeth validates version parameter", {
    expect_error(
        create_bedmeth(version = "invalid"),
        "Invalid version"
    )
})

test_that("create_bedmeth returns proper BED format", {
    skip_on_cran()
    skip_if_not_installed("minfi")
    skip_if_not_installed("IlluminaHumanMethylation450kanno.ilmn12.hg19")

    result <- create_bedmeth(version = "450k")

    # Check column names match BED format
    expect_equal(colnames(result), c("chrom", "start", "end", "probeID"))

    # Check data types
    expect_type(result$start, "integer")
    expect_type(result$end, "integer")
    expect_type(result$probeID, "character")
})

test_that("create_bedmeth end equals start", {
    skip_on_cran()
    skip_if_not_installed("minfi")
    skip_if_not_installed("IlluminaHumanMethylation450kanno.ilmn12.hg19")

    result <- create_bedmeth(version = "450k")

    # For single-base CpG sites, end should equal start
    expect_equal(result$start, result$end)
})

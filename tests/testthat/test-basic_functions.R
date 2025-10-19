# Basic tests that don't require special dependencies and will run in all environments

test_that("Package loads successfully", {
    expect_true(require(GIMP, quietly = TRUE))
})

test_that("Core exported functions exist", {
    expect_true(exists("make_cpgs"))
    expect_true(exists("make_ICRs"))
    expect_true(exists("iDMPs"))
    expect_true(exists("ICRs_heatmap"))
    expect_true(exists("plot_cpgs_coverage"))
    expect_true(exists("plot_line_ICR"))
    expect_true(exists("create_bedmeth"))
})

test_that("Data objects are available", {
    data(DMRs.hg19)
    expect_true(exists("DMRs.hg19"))
    expect_s3_class(DMRs.hg19, "data.frame")
    expect_true(nrow(DMRs.hg19) > 0)

    data(DMRs.hg38)
    expect_true(exists("DMRs.hg38"))
    expect_s3_class(DMRs.hg38, "data.frame")
    expect_true(nrow(DMRs.hg38) > 0)
})

test_that("DMRs.hg19 has expected structure", {
    data(DMRs.hg19)
    expect_true("ICR" %in% colnames(DMRs.hg19))
    expect_true("chrom" %in% colnames(DMRs.hg19))
    expect_true("start" %in% colnames(DMRs.hg19))
    expect_true("end" %in% colnames(DMRs.hg19))
})

test_that("DMRs.hg38 has expected structure", {
    data(DMRs.hg38)
    expect_true("ICR" %in% colnames(DMRs.hg38))
    expect_true("chrom" %in% colnames(DMRs.hg38))
    expect_true("start" %in% colnames(DMRs.hg38))
    expect_true("end" %in% colnames(DMRs.hg38))
})

test_that("BED data objects can be loaded", {
    # Load bed data
    data(bed450k)
    data(bedEPICv1)
    data(bedEPICv2)

    # Check they exist
    expect_true(exists("bed450k"))
    expect_true(exists("bedEPICv1"))
    expect_true(exists("bedEPICv2"))

    # Check structure
    expect_s3_class(bed450k, "data.frame")
    expect_s3_class(bedEPICv1, "data.frame")
    expect_s3_class(bedEPICv2, "data.frame")

    # Check they have data
    expect_true(nrow(bed450k) > 0)
    expect_true(nrow(bedEPICv1) > 0)
    expect_true(nrow(bedEPICv2) > 0)
})

test_that("bed450k has correct columns", {
    data(bed450k)
    expect_true("chrom" %in% colnames(bed450k))
    expect_true("probeID" %in% colnames(bed450k))
})

test_that("bedEPICv1 has correct columns", {
    data(bedEPICv1)
    expect_true("chrom" %in% colnames(bedEPICv1))
    expect_true("probeID" %in% colnames(bedEPICv1))
})

test_that("bedEPICv2 has correct columns", {
    data(bedEPICv2)
    expect_true("chrom" %in% colnames(bedEPICv2))
    expect_true("probeID" %in% colnames(bedEPICv2))
})

test_that("make_cpgs validates bedmeth parameter", {
    # Create minimal test data
    beta_matrix <- matrix(0.5, nrow = 5, ncol = 2)
    rownames(beta_matrix) <- paste0("cg", 1:5)
    colnames(beta_matrix) <- c("Sample1", "Sample2")

    expect_error(make_cpgs(beta_matrix, bedmeth = "invalid"), "Invalid bedmeth")
})

test_that("make_ICRs validates bedmeth parameter", {
    # Create minimal test data
    beta_matrix <- matrix(0.5, nrow = 5, ncol = 2)
    rownames(beta_matrix) <- paste0("cg", 1:5)
    colnames(beta_matrix) <- c("Sample1", "Sample2")

    expect_error(make_ICRs(beta_matrix, bedmeth = "invalid"), "Invalid bedmeth")
})

test_that("iDMPs validates Control group requirement", {
    # Create test data without annotation
    meth_data <- matrix(runif(20 * 4, 0.3, 0.7), nrow = 20, ncol = 4)
    colnames(meth_data) <- paste0("Sample_", 1:4)

    ICRcpg_data <- data.frame(
        meth_data,
        cstart = seq(1000000, by = 100, length.out = 20),
        ICR = rep("ICR1", 20),
        start = 1000000,
        end = 1005000
    )

    # Without Control group
    sampleInfo <- factor(c(rep("Group1", 2), rep("Group2", 2)))

    expect_error(
        iDMPs(data = ICRcpg_data, sampleInfo = sampleInfo),
        "Control.*group.*must be present"
    )
})

test_that("iDMPs requires data frame input", {
    meth_matrix <- matrix(runif(20 * 4, 0.3, 0.7), nrow = 20, ncol = 4)
    sampleInfo <- factor(c(rep("Control", 2), rep("Case", 2)))

    expect_error(
        iDMPs(data = meth_matrix, sampleInfo = sampleInfo),
        "input data must be a data frame"
    )
})

test_that("iDMPs validates sampleInfo length", {
    meth_data <- matrix(runif(20 * 4, 0.3, 0.7), nrow = 20, ncol = 4)
    colnames(meth_data) <- paste0("Sample_", 1:4)

    ICRcpg_data <- data.frame(
        meth_data,
        cstart = seq(1000000, by = 100, length.out = 20),
        ICR = rep("ICR1", 20),
        start = 1000000,
        end = 1005000
    )

    # Wrong sampleInfo length
    sampleInfo <- factor(c(rep("Control", 2), rep("Case", 1)))

    expect_error(
        iDMPs(data = ICRcpg_data, sampleInfo = sampleInfo),
        "doesn't match sampleInfo length"
    )
})

test_that("GEO functions exist and are exported", {
    expect_true(exists("validate_geo_dataset"))
    expect_true(exists("process_geo_dataset"))
    expect_true(exists("diagnose_geo_dataset"))
    expect_true(exists("get_geo_phenotype_data"))
    expect_true(exists("process_geo_with_mappings"))
})

test_that("IDAT functions exist and are exported", {
    expect_true(exists("read_idat_zip"))
})

test_that("Shiny app function exists", {
    expect_true(exists("GIMP_app"))
})

test_that("Helper functions exist", {
    # These are internal but should exist in namespace
    expect_true(exists("get_bed_data", where = asNamespace("GIMP"), inherits = FALSE))
    expect_true(exists("detect_array_platform", where = asNamespace("GIMP"), inherits = FALSE))
    expect_true(exists("assign_sample_groups", where = asNamespace("GIMP"), inherits = FALSE))
})

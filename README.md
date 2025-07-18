# GIMP: Genomic Imprinting Methylation Patterns

[![R](https://img.shields.io/badge/R-4.0+-blue.svg)](https://cran.r-project.org/)
[![Bioconductor](https://img.shields.io/badge/Bioconductor-3.14+-green.svg)](https://bioconductor.org/)
[![License](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

<div align="center">
  <img src="GIMP.logo.png" alt="GIMP Logo" width="200"/>
</div>

**GIMP** (Genomic Imprinting Methylation Patterns) is an R package designed for the comprehensive analysis of Imprinting Control Regions (ICRs) from methylation array data. It provides a complete pipeline for extracting imprinted CpGs (iCpGs), computing coverage, and analyzing ICRs at both probe and sample-specific levels.

** This Repo is a fork of [ngsFC/GIMP.git](https://github.com/ngsFC/GIMP.git)

**🆕 NEW in v0.2.0**: GIMP now has a shinyApp supports **direct processing of raw IDAT files** from ZIP archives, making it a complete solution from raw data to specialized genomic imprinting analysis!

## Features

### Core Capabilities
- **Dual Input Support**: Process both preprocessed methylation data and raw IDAT files
- **Specialized ICR Analysis**: Focus on imprinting control regions with curated coordinates
- **Interactive Visualizations**: Comprehensive Shiny app with plotly integration
- **Multiple Array Support**: 450k, EPIC v1, and EPIC v2 arrays
- **Quality Control**: Automatic QC for IDAT data with detailed reporting
- **Flexible Analysis**: Beta values, delta-beta, and defect matrix visualizations

### Unique Imprinting Features
- **ICR-specific coordinates**: Based on [Joshi et al. 2016](https://doi.org/10.1080/15592294.2016.1264561)
- **Defect matrix analysis**: SD-based detection of imprinting disorders
- **Specialized heatmaps**: Designed for imprinting pattern visualization
- **Interactive region explorer**: Detailed methylation profiles across ICRs
- **Clinical interpretation**: Tools for analyzing imprinting disorders

## Installation

### Standard Installation

```r
# Install devtools if you don't have it
if (!requireNamespace("devtools", quietly = TRUE)) {
  install.packages("devtools")
}

# Install GIMP from GitHub
devtools::install_github("ngsFC/GIMP")

# Load the package
library(GIMP)
```

### For IDAT Processing (Additional Requirements)

```r
# Install Bioconductor packages for IDAT support
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

# Install required annotation packages
BiocManager::install(c(
  "minfi",
  "IlluminaHumanMethylation450kanno.ilmn12.hg19",
  "IlluminaHumanMethylationEPICanno.ilm10b4.hg19",
  "IlluminaHumanMethylationEPICv2anno.20a1.hg38"
))

# Test IDAT functionality
test_idat_functionality()
```

## Quick Start

### Option 1: Using the Shiny App (Recommended for Beginners)

```r
library(GIMP)

# Launch interactive app
GIMP_app()

# For large IDAT files, increase upload limit
GIMP_app(max_upload_size_mb = 1000)  # 1GB limit
```

### Option 2: Command Line Usage

#### For Preprocessed Data
```r
library(GIMP)

# Load your beta value matrix
df <- readRDS("methylation_data.rds")

# Standard GIMP workflow
ICRcpg <- make_cpgs(Bmatrix = df, bedmeth = "v1")
df_ICR <- make_ICRs(Bmatrix = df, bedmeth = "v1")

# Generate sample information
sampleInfo <- c(rep("Control", 10), rep("Case", 10))

# Create heatmap
ICRs_heatmap(df_ICR, sampleInfo = sampleInfo, plot_type = "beta")

# Differential analysis
dmps <- iDMPs(data = ICRcpg, sampleInfo = sampleInfo)
```

#### For Raw IDAT Files
```r
library(GIMP)

# Process IDAT files from ZIP archive
idat_data <- read_idat_zip(
  zip_file = "methylation_data.zip",
  array_type = "EPIC",
  normalize_method = "quantile"
)

# Extract processed data
beta_matrix <- idat_data$beta_matrix
sample_info <- idat_data$sample_info

# Continue with standard GIMP workflow
ICRcpg <- make_cpgs(Bmatrix = beta_matrix, bedmeth = "v1")
# ... rest of analysis
```

## IDAT File Processing Guide

### Required File Structure

Your ZIP file must contain:

```
methylation_data.zip
├── SampleID1_ChipID_Position_Red.idat
├── SampleID1_ChipID_Position_Grn.idat
├── SampleID2_ChipID_Position_Red.idat
├── SampleID2_ChipID_Position_Grn.idat
├── ...
└── samplesheet.csv
```

### Sample Sheet Requirements

#### Required Columns
- **Sample_Name**: Unique identifier for each sample
- **Sentrix_ID**: Chip/slide identifier (e.g., "200123456789")
- **Sentrix_Position**: Array position on chip (e.g., "R01C01")

#### Optional Columns
- **Sample_Group**: For automatic group assignment ("Control", "Case")
- **Sample_Plate**: Plate information
- **Sample_Well**: Well position

#### Sample Sheet Example

```csv
Sample_Name,Sentrix_ID,Sentrix_Position,Sample_Group
Control_01,200123456789,R01C01,Control
Control_02,200123456789,R02C01,Control
Control_03,200123456789,R03C01,Control
Case_01,200123456790,R01C01,Case
Case_02,200123456790,R02C01,Case
Case_03,200123456790,R03C01,Case
```

### Creating a Sample Sheet Template

```r
# Generate template with your sample information
template <- create_sample_sheet_template(
  sample_names = c("Sample_1", "Sample_2", "Sample_3", "Sample_4"),
  sentrix_ids = c("200123456789", "200123456789", "200123456790", "200123456790"),
  sentrix_positions = c("R01C01", "R02C01", "R01C01", "R02C01"),
  groups = c("Control", "Control", "Case", "Case")
)

# Save template
write.csv(template, "samplesheet.csv", row.names = FALSE)
```

### Common IDAT Issues and Solutions

#### Issue: "Missing IDAT files for samples"
**Solution**: Check that IDAT file names match your sample sheet

```r
# Diagnose your ZIP file structure
diagnose_idat_structure("your_file.zip")

# Generate sample sheet from actual IDAT files
new_sheet <- generate_samplesheet_from_idats("your_file.zip")
```

#### Issue: "Maximum upload size exceeded"
**Solution**: Increase upload limit

```r
# In Shiny app
GIMP_app(max_upload_size_mb = 1000)  # 1GB

# Or in R console before launching
options(shiny.maxRequestSize = 1000*1024^2)
GIMP_app()
```

#### Issue: File size too large
**Solutions**:
1. **Compress your ZIP file** with maximum compression
2. **Remove unnecessary files** from the ZIP
3. **Process in smaller batches**
4. **Use command line instead of Shiny app**

### Supported Array Types

| Array Type | GIMP Parameter | Genome Build | Typical File Size |
|-----------|----------------|--------------|-------------------|
| 450k | `"450k"` | hg19 | 50-150MB (10-50 samples) |
| EPIC v1 | `"EPIC"` or `"v1"` | hg19 | 100-500MB (10-50 samples) |
| EPIC v2 | `"EPICv2"` or `"v2"` | hg38 | 200-800MB (10-50 samples) |

## Analysis Workflows

### Standard Workflow

1. **Data Upload**
   - Upload processed data (CSV/RDS) OR raw IDAT files (ZIP)
   - Assign sample groups (Control vs Case)

2. **CpG Coverage Analysis**
   - Visualize probe coverage across ICRs
   - Assess data quality and completeness

3. **ICR Heatmap Analysis**
   - Generate methylation heatmaps
   - Choose from beta, delta-beta, or defect matrix views

4. **Differential Methylation**
   - Identify significantly different positions
   - Generate volcano plots and summary statistics

5. **Region Explorer**
   - Detailed visualization of specific ICRs
   - Interactive plots with DMP highlighting

### Analysis Functions

#### Core Functions
- `make_cpgs()`: Extract ICR CpG sites
- `make_ICRs()`: Create ICR-level methylation matrix
- `ICRs_heatmap()`: Generate ICR heatmaps
- `iDMPs()`: Identify differentially methylated positions
- `plot_line_region()`: Visualize specific ICR regions

#### IDAT Functions
- `read_idat_zip()`: Process IDAT files from ZIP
- `diagnose_idat_structure()`: Analyze ZIP file contents
- `generate_samplesheet_from_idats()`: Auto-generate sample sheets
- `test_idat_functionality()`: Verify IDAT processing setup

#### Utility Functions
- `plot_cpgs_coverage()`: Visualize CpG coverage
- `create_sample_sheet_template()`: Generate sample sheet templates
- `GIMP_app()`: Launch Shiny application

## Advanced Usage

### Custom Analysis Parameters

```r
# Advanced heatmap with defect detection
ICRs_heatmap(
  df_ICR = icr_data,
  sampleInfo = sample_groups,
  plot_type = "defect",
  sd_threshold = 2.5,  # More sensitive detection
  order_by = "meth"    # Cluster by methylation patterns
)

# Sensitive DMP detection
dmps <- iDMPs(
  data = ICRcpg,
  sampleInfo = sample_groups,
  pValueCutoff = 0.01  # More stringent threshold
)
```

### Working with Large Datasets

```r
# For large IDAT files
idat_data <- read_idat_zip(
  zip_file = "large_dataset.zip",
  normalize_method = "funnorm",  # Faster normalization
  detection_pval = 0.05,         # Less stringent QC
  remove_failed_samples = TRUE
)

# Batch processing
process_batch <- function(zip_files) {
  results <- list()
  for (zip_file in zip_files) {
    results[[basename(zip_file)]] <- read_idat_zip(zip_file)
  }
  return(results)
}
```

## Troubleshooting

### Installation Issues

```r
# Check GIMP installation
library(GIMP)

# Test IDAT functionality
test_idat_functionality()

# Fix common issues
fix_minfi_installation()

# Check system requirements
check_minfi_functions()
```

### Memory Requirements

| Dataset Size | Recommended RAM | Processing Time |
|-------------|-----------------|-----------------|
| 10-20 samples | 4GB | 2-5 minutes |
| 50-100 samples | 8GB | 5-15 minutes |
| 100+ samples | 16GB+ | 15+ minutes |

### Common Error Solutions

#### "minfi functions not found"
```r
# Reinstall minfi and dependencies
BiocManager::install("minfi", force = TRUE)
.rs.restartR()  # Restart R session
```

#### "No ICRs found"
- Check that your array type matches your data
- Verify probe IDs are in the correct format
- Try different `bedmeth` parameters

#### "Dimension mismatch errors"
- Check that sample information matches data dimensions
- Verify no missing values in critical columns
- Use `diagnose_idat_structure()` for IDAT files

## Data Sources

### ICR Coordinates
GIMP uses curated ICR coordinates from:
- **Joshi et al. (2016)**: "Detailed annotation of human Imprinting Control Regions" ([DOI](https://doi.org/10.1080/15592294.2016.1264561))
- Coordinates available for both hg19 and hg38 genome builds

### Compatible Data Sources
- **Illumina methylation arrays**: 450k, EPIC v1, EPIC v2
- **GEO datasets**: Direct processing of downloaded IDAT files
- **Preprocessed data**: From other methylation analysis pipelines
- **Clinical samples**: Hospital/research institution data

## Citation

If you use GIMP in your research, please cite:

```
Cecere, F. (2024). GIMP: Genomic Imprinting Methylation Patterns. 
R package version 0.2.0. https://github.com/ngsFC/GIMP
```

## Contributing

We welcome contributions! Please:

1. **Report bugs**: Use GitHub issues with detailed error messages
2. **Suggest features**: Describe your use case and proposed functionality
3. **Submit code**: Follow R package development best practices
4. **Share datasets**: Help us test with diverse methylation data

## Support

### Getting Help
- **GitHub Issues**: For bug reports and feature requests
- **Documentation**: `?function_name` for detailed help
- **Vignettes**: Comprehensive usage examples

### Frequently Asked Questions

**Q: Can GIMP analyze other genomic regions besides ICRs?**
A: GIMP is specifically designed for ICR analysis, but the core functions can be adapted for other regions of interest.

**Q: What's the difference between GIMP and shinyepico?**
A: GIMP focuses specifically on imprinting analysis with specialized ICR coordinates and imprinting disorder detection, while shinyepico is a general-purpose methylation analysis tool.

**Q: Can I use GIMP for non-human samples?**
A: Currently, GIMP includes human ICR coordinates. For other species, you would need to provide custom ICR coordinates.

## Acknowledgments

The GIMP package was developed by Francesco Cecere. We gratefully acknowledge:

- **National Centre for HPC, Big Data and Quantum Computing** for computational resources
- **Bioconductor community** for methylation analysis infrastructure
- **minfi developers** for IDAT processing capabilities
- **shinyepico** for inspiration on user-friendly methylation analysis

## License

GIMP is released under the MIT License. See LICENSE file for details.

---

**Keywords**: DNA methylation, genomic imprinting, IDAT processing, Illumina arrays, ICR analysis, R package, Bioconductor, shiny application

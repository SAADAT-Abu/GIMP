# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Repository Overview

GIMP (Genomic Imprinting Methylation Patterns) is an R package for analyzing Imprinting Control Regions (ICRs) from DNA methylation array data. It supports both preprocessed methylation data and raw IDAT files from Illumina arrays (450k, EPIC v1, EPIC v2).

## Key Development Commands

### R Package Development
```r
# Install package dependencies
install.packages("devtools")
BiocManager::install(c("minfi", "IlluminaHumanMethylation450kanno.ilmn12.hg19", 
                       "IlluminaHumanMethylationEPICanno.ilm10b4.hg19", 
                       "IlluminaHumanMethylationEPICv2anno.20a1.hg38"))

# Build and install package from source
devtools::install()

# Load package for development
devtools::load_all()

# Generate documentation
devtools::document()

# Build package data (run from data-raw/)
source("preparedata.R")

# Check package
devtools::check()

# Build vignettes
devtools::build_vignettes()
```

### Testing the Package
```r
# Launch the Shiny app
GIMP_app()

# Test IDAT functionality
test_idat_functionality()

# Test core functions
library(GIMP)
data(bedEPICv1)  # Load example data
```

## Architecture Overview

### Core Components

**Data Processing Pipeline:**
- `read_idat_zip()`: Processes raw IDAT files from ZIP archives with parallel processing support
- `make_cpgs()`: Extracts ICR-specific CpG sites from methylation data
- `make_ICRs()`: Creates ICR-level methylation matrices
- `create_bedmeth()`: Generates bed-format methylation data

**Analysis Functions:**
- `iDMPs()`: Identifies differentially methylated positions
- `ICRs_heatmap()`: Generates specialized heatmaps for imprinting patterns
- `plot_line_ICR()`: Creates detailed ICR region visualizations
- `plot_cpgs_coverage()`: Visualizes CpG coverage across ICRs

**Interactive Interface:**
- `GIMP_app()`: Launches Shiny application (located in `inst/shiny/`)
- Supports file uploads, interactive visualizations, and result export

**GEO Data Integration:**
- `validate_geo_dataset()`: Validates GEO datasets for IDAT availability
- `process_geo_dataset()`: Downloads and processes GEO methylation data
- `diagnose_geo_dataset()`: Analyzes GEO dataset structure

### Package Structure

**Core R Files:**
- `R/GIMP_app.r`: Main Shiny app launcher
- `R/read_idat_zip.r`: IDAT file processing with parallel support
- `R/geo_functions.r`: GEO dataset integration functions
- `R/iDMPs.R`: Differential methylation analysis
- `R/ICRs_heatmap.R`: Specialized heatmap generation
- `R/minfi_helpers.R`: Minfi package integration utilities

**Data Assets:**
- `data/DMRs.hg19.rda` & `data/DMRs.hg38.rda`: ICR coordinates for genome builds
- `data/bed450k.rda`, `data/bedEPICv1.rda`, `data/bedEPICv2.rda`: Array-specific probe mappings

**Shiny Application:**
- `inst/shiny/app.R`: Main app entry point
- `inst/shiny/ui.R` & `inst/shiny/server.R`: UI and server logic
- `inst/shiny/help.md`: In-app help documentation

### Key Dependencies

**Core Bioconductor:**
- `minfi`: IDAT processing and normalization
- Array annotation packages for 450k, EPIC v1, EPIC v2

**Analysis & Visualization:**
- `tidyverse`: Data manipulation
- `ggplot2`, `pheatmap`, `plotly`: Visualization
- `limma`: Statistical analysis for DMPs

**Parallel Processing:**
- `doParallel`: Multi-core IDAT processing
- `parallel`: Core detection and management

## Development Notes

### IDAT Processing
- Supports parallel processing with `n_cores` parameter
- Handles multiple normalization methods (quantile, SWAN, funnorm, noob)
- Includes QC metrics and failed sample detection
- Requires large memory for processing (see README for recommendations)

### ICR-Specific Analysis
- Uses curated ICR coordinates from Joshi et al. 2016
- Supports both hg19 and hg38 genome builds
- Specialized for genomic imprinting disorder analysis
- Includes defect matrix analysis with SD-based thresholds

### Array Support
- 450k arrays: ~450K probes, hg19 coordinates
- EPIC v1: ~850K probes, hg19 coordinates  
- EPIC v2: ~930K probes, hg38 coordinates

### File Upload Limits
- Default Shiny upload limit: 500MB
- Can be increased with `GIMP_app(max_upload_size_mb = 1000)`
- For large datasets, recommend command-line processing over Shiny interface
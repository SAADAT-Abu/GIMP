# GIMP Shiny App User Guide

## Overview

The GIMP (Genomic Imprinting Methylation Patterns) Shiny app provides an interactive interface for analyzing methylation patterns at Imprinted Control Regions (ICRs) from both preprocessed methylation data and raw IDAT files. This comprehensive tool is designed for researchers studying genomic imprinting disorders and DNA methylation patterns.

## Getting Started

### System Requirements
- **R Version**: 4.0 or higher
- **Memory**: 4GB minimum (8GB+ recommended for large datasets)
- **Browser**: Modern web browser with JavaScript enabled
- **Network**: Required for downloading annotation packages

### Launching the App
```r
library(GIMP)

# Standard launch
GIMP_app()

# For large files (adjust upload limit)
GIMP_app(max_upload_size_mb = 1000)  # 1GB limit
```

## Data Input Options

GIMP supports two types of input data:

### Option 1: Processed Methylation Data
**Best for**: GEO datasets, preprocessed data, quick analysis

**Supported Formats**:
- CSV files (rows = CpGs, columns = samples)
- RDS files (R data objects)
- Excel files (.xlsx)
- Tab-delimited text files

**Data Requirements**:
- Beta values between 0 and 1
- Row names must be probe IDs (e.g., cg00000029)
- Column names must be sample identifiers
- No missing values in probe IDs

### Option 2: Raw IDAT Files (ZIP Archive)
**Best for**: Raw Illumina data, complete processing control

**ZIP File Structure**:
```
methylation_data.zip
├── Sample1_ChipID_Position_Red.idat
├── Sample1_ChipID_Position_Grn.idat
├── Sample2_ChipID_Position_Red.idat
├── Sample2_ChipID_Position_Grn.idat
├── ...
└── samplesheet.csv
```

## Sample Sheet Creation Guide

### Required Columns

| Column Name | Description | Example | Notes |
|-------------|-------------|---------|--------|
| `Sample_Name` | Unique sample identifier | `Control_01` | No spaces or special characters recommended |
| `Sentrix_ID` | Chip/slide identifier | `200123456789` | 12-digit number from Illumina |
| `Sentrix_Position` | Position on chip | `R01C01` | Format: R##C## |

### Optional Columns

| Column Name | Description | Example | Usage |
|-------------|-------------|---------|--------|
| `Sample_Group` | Group assignment | `Control`, `Case` | Enables automatic group assignment |
| `Sample_Plate` | Plate information | `Plate1` | For batch effect tracking |
| `Sample_Well` | Well position | `A01` | For plate layout tracking |
| `Sample_Type` | Sample type | `Blood`, `Tissue` | For clinical annotation |

### Sample Sheet Examples

#### Basic Sample Sheet
```csv
Sample_Name,Sentrix_ID,Sentrix_Position
Sample_001,200123456789,R01C01
Sample_002,200123456789,R02C01
Sample_003,200123456789,R03C01
Sample_004,200123456790,R01C01
```

#### Complete Sample Sheet with Groups
```csv
Sample_Name,Sentrix_ID,Sentrix_Position,Sample_Group,Sample_Type,Sample_Plate
Control_01,200123456789,R01C01,Control,Blood,Plate1
Control_02,200123456789,R02C01,Control,Blood,Plate1
Control_03,200123456789,R03C01,Control,Blood,Plate1
BWS_Patient_01,200123456790,R01C01,Case,Blood,Plate1
BWS_Patient_02,200123456790,R02C01,Case,Blood,Plate1
SRS_Patient_01,200123456790,R03C01,Case,Blood,Plate1
```

#### GEO Dataset Sample Sheet
```csv
Sample_Name,Sentrix_ID,Sentrix_Position,Sample_Group
GSM3861617,200989060193,R02C01,Control
GSM3861618,200989060193,R08C01,Control
GSM3861619,200991620059,R02C01,Case
GSM3861669,200992350051,R07C01,Case
```

### Creating Sample Sheets in GIMP

#### Method 1: Generate Template
```r
# Create template with your sample information
template <- create_sample_sheet_template(
  sample_names = c("Sample_1", "Sample_2", "Sample_3", "Sample_4"),
  sentrix_ids = c("200123456789", "200123456789", "200123456790", "200123456790"),
  sentrix_positions = c("R01C01", "R02C01", "R01C01", "R02C01"),
  groups = c("Control", "Control", "Case", "Case")
)

# Save template
write.csv(template, "samplesheet.csv", row.names = FALSE)
```

#### Method 2: Generate from IDAT Files
```r
# Auto-generate from actual IDAT files in ZIP
new_sheet <- generate_samplesheet_from_idats("your_file.zip")

# Review and edit the generated sheet
View(new_sheet)

# Save corrected version
write.csv(new_sheet, "corrected_samplesheet.csv", row.names = FALSE)
```

## Workflow Guide

### Step 1: Data Upload

#### For Processed Data
1. Select **"Processed Data (CSV/RDS/Excel)"**
2. Upload your methylation file
3. Select appropriate array type:
   - **EPIC v1 (hg19)**: Most common, ~850k probes
   - **EPIC v2 (hg38)**: Newest version, ~935k probes
   - **450k (hg19)**: Older arrays, ~485k probes
4. Assign sample groups (Control vs Case)
5. Click **"Process Data"**

#### For IDAT Files
1. Select **"Raw IDAT Files (ZIP)"**
2. Upload ZIP file containing IDAT files and sample sheet
3. **Preview ZIP contents** to verify structure
4. Configure processing options:
   - **Array Type**: EPIC, 450k, or EPICv2
   - **Normalization**: Quantile (recommended), SWAN, funnorm, or noob
   - **Detection P-value**: Threshold for probe quality (default: 0.01)
5. Click **"Process IDAT Files"**
6. Wait for processing (can take 5-15 minutes)

### Step 2: CpG Coverage Analysis

**Purpose**: Assess how well your array covers ICR regions

1. Navigate to **"CpG Analysis"** tab
2. Click **"Generate CpG Analysis"**
3. Review three outputs:
   - **Coverage Counts**: Number of probes per ICR
   - **Coverage Percentage**: Percentage of possible probes covered
   - **Coverage Data**: Detailed statistics table

**Interpretation**:
- High coverage (>80%): Excellent for analysis
- Medium coverage (50-80%): Good for analysis
- Low coverage (<50%): Limited analysis power

### Step 3: ICR Heatmap Analysis

**Purpose**: Visualize methylation patterns across all ICRs

#### Heatmap Options

**Plot Types**:
- **Beta Values**: Raw methylation values (0-1 scale)
  - Use for: General methylation pattern visualization
- **Delta Beta**: Difference from control group mean
  - Use for: Highlighting changes relative to controls
- **Defect Matrix**: Binary visualization of abnormal methylation
  - Use for: Clinical analysis, disorder detection
  - SD Threshold: 2-3 (standard), 1.5 (sensitive), 4+ (stringent)

**Ordering Options**:
- **Coordinates**: ICRs ordered by genomic position
- **Methylation Values**: ICRs clustered by methylation similarity

#### Clinical Interpretation Guide

**Normal Imprinting Patterns**:
- Clear distinction between maternal and paternal alleles
- Consistent patterns within control groups
- Expected 50% methylation for most ICRs

**Loss of Imprinting (LOI)**:
- Abnormal methylation patterns in case samples
- Hypermethylation or hypomethylation of ICRs
- Patient-specific or syndrome-specific patterns

**Common Syndromes**:
- **Beckwith-Wiedemann Syndrome**: ICR1/ICR2 alterations
- **Silver-Russell Syndrome**: ICR1 hypomethylation
- **Prader-Willi/Angelman**: SNRPN ICR alterations

### Step 4: Differential Methylation Analysis

**Purpose**: Identify statistically significant differences between groups

#### Parameters
- **P-value Cutoff**: 0.05 (standard), 0.01 (stringent), 0.1 (exploratory)

#### Outputs
1. **Significant DMPs Table**: 
   - logFC: Log fold change in methylation
   - P.Value: Statistical significance
   - adj.P.Val: Multiple testing corrected p-value
   - ICR: Associated imprinting control region

2. **Volcano Plot**: 
   - X-axis: Log fold change (effect size)
   - Y-axis: -log10(p-value) (significance)
   - Red points: Significant DMPs

3. **Summary Statistics**: 
   - Total CpGs tested
   - Number of significant DMPs
   - Top hits by significance

### Step 5: Region Explorer

**Purpose**: Detailed visualization of specific ICR regions

#### Usage
1. Select ICR from dropdown (only ICRs with significant DMPs shown)
2. Choose plot type:
   - **Interactive (plotly)**: Zoom, pan, hover for details
   - **Static (ggplot2)**: High-quality publication plots
3. Click **"Plot Region"**

#### Plot Elements
- **Lines**: Individual sample methylation profiles
- **Colors**: Blue (Control), Red (Case)
- **Vertical dashed lines**: Significant DMP positions
- **Red rug marks**: Significant DMPs
- **Large points**: CpGs that are significant DMPs

## File Upload Guidelines

### File Size Limits
- **Default**: 500MB upload limit
- **Large datasets**: Use `GIMP_app(max_upload_size_mb = 1000)` for 1GB
- **Very large datasets**: Consider command-line processing

### Typical File Sizes
| Array Type | Samples | Processed Data | IDAT ZIP |
|-----------|---------|----------------|----------|
| 450k | 10-20 | 50-100MB | 100-200MB |
| EPIC v1 | 10-20 | 100-200MB | 200-400MB |
| EPIC v2 | 10-20 | 150-300MB | 300-600MB |
| EPIC v1 | 50+ | 500MB+ | 1GB+ |

### Optimization Tips
- **Compress ZIP files** with maximum compression
- **Remove unnecessary files** from ZIP archives
- **Process in batches** for very large datasets
- **Use command line** for routine processing

## Quality Control

### IDAT Data QC Metrics

#### Sample-Level QC
- **Detection P-values**: Proportion of failed probes per sample
- **Signal Intensity**: Overall signal strength
- **Failed Samples**: Samples with >10% failed probes (automatically removed)

#### Probe-Level QC
- **Detection Rate**: Percentage of samples with good signal
- **Variance Filtering**: Removal of invariant probes
- **Cross-hybridization**: Removal of problematic probes

#### QC Visualizations
1. **Failed Probes per Sample**: Bar plot showing sample quality
2. **Beta Value Distribution**: Histogram of methylation values
3. **Sample Correlation**: Heatmap of inter-sample correlations
4. **Probe Intensity**: Distribution of signal intensities

### Troubleshooting QC Issues

#### High Failure Rates
**Symptoms**: Many samples with >10% failed probes
**Solutions**:
- Check IDAT file integrity
- Verify array type selection
- Consider less stringent detection p-value (0.05 instead of 0.01)

#### Low Signal Intensity
**Symptoms**: Warning about low intensity samples
**Solutions**:
- Check sample storage conditions
- Verify DNA quality and quantity
- Consider excluding problematic samples

#### Batch Effects
**Symptoms**: Samples clustering by technical factors
**Solutions**:
- Include batch information in sample sheet
- Use appropriate normalization method
- Consider batch correction in downstream analysis

## Advanced Features

### Normalization Methods

| Method | Speed | Use Case | Recommendation |
|--------|-------|----------|----------------|
| **Quantile** | Fast | General purpose | ✅ Default choice |
| **SWAN** | Medium | Corrects probe bias | Good for mixed sample types |
| **Functional** | Slow | Removes unwanted variation | Best for large studies |
| **Noob** | Fast | Background correction | Good for low-quality data |

### Custom Analysis Parameters

#### Sensitive Analysis (Disorder Detection)
- P-value cutoff: 0.01
- SD threshold: 2.0
- Remove failed samples: Yes
- Normalization: Quantile

#### Exploratory Analysis
- P-value cutoff: 0.1
- SD threshold: 3.0
- Remove failed samples: No
- Normalization: SWAN

#### Publication-Quality Analysis
- P-value cutoff: 0.05
- SD threshold: 2.5
- Remove failed samples: Yes
- Normalization: Functional

## Export and Reporting

### Download Options
- **Heatmaps**: PDF format, publication-ready
- **DMP Results**: CSV format with all statistics
- **QC Reports**: Comprehensive quality metrics

### Integration with Other Tools
- **Excel**: Direct import of CSV results
- **R/Bioconductor**: Export data for custom analysis
- **Clinical databases**: Structured output for LIMS integration

## Common Issues and Solutions

### Upload Issues

#### "Maximum upload size exceeded"
```r
# Increase upload limit
GIMP_app(max_upload_size_mb = 1000)

# Or set globally
options(shiny.maxRequestSize = 1000*1024^2)
```

#### "File format not supported"
- Ensure file extensions are correct (.csv, .rds, .xlsx, .zip)
- Check file is not corrupted
- Verify data structure (rows=probes, columns=samples)

### IDAT Processing Issues

#### "Missing IDAT files for samples"
```r
# Diagnose ZIP file structure
diagnose_idat_structure("your_file.zip")

# Generate new sample sheet from IDAT files
new_sheet <- generate_samplesheet_from_idats("your_file.zip")
```

#### "minfi functions not found"
```r
# Reinstall required packages
BiocManager::install("minfi", force = TRUE)
BiocManager::install("IlluminaHumanMethylationEPICanno.ilm10b4.hg19")

# Test installation
test_idat_functionality()
```

### Analysis Issues

#### "No ICRs found"
- Check array type selection matches your data
- Verify probe IDs are in standard format
- Try different bedmeth parameters

#### "No significant DMPs"
- Increase p-value cutoff (try 0.1)
- Check sample group assignments
- Verify sufficient sample size per group

#### "Memory issues"
- Close other applications
- Process smaller batches
- Use command line instead of Shiny app

## Best Practices

### Study Design
- **Minimum samples**: 6 per group (more is better)
- **Balanced groups**: Equal numbers when possible
- **Technical replicates**: Include when available
- **Batch information**: Record all technical factors

### Data Quality
- **DNA quality**: Use high-quality, intact DNA
- **Storage conditions**: Consistent sample handling
- **Array processing**: Follow Illumina protocols
- **Documentation**: Maintain detailed sample records

### Analysis Strategy
1. **Start with QC**: Always review quality metrics first
2. **Explore data**: Use heatmaps before statistical testing
3. **Multiple approaches**: Try different plot types and thresholds
4. **Validate results**: Confirm findings with additional methods
5. **Clinical correlation**: Link results to phenotype information

## Getting Help

### Documentation Resources
- **Function Help**: Use `?function_name` in R console for detailed documentation
- **Package Vignettes**: Comprehensive examples and workflows
- **GitHub Repository**: [https://github.com/ngsFC/GIMP](https://github.com/ngsFC/GIMP)
- **Issue Tracker**: Report bugs and request features on GitHub

### Common Error Messages

#### "Error in read.csv: no lines available in input"
**Cause**: Empty or corrupted sample sheet
**Solution**: Check sample sheet format and content

#### "Error: subscript out of bounds"
**Cause**: Dimension mismatch between data and sample information
**Solution**: Verify sample numbers match between data and sample sheet

#### "Error: could not find function 'getDetectionP'"
**Cause**: minfi package issues
**Solution**: 
```r
BiocManager::install("minfi", force = TRUE)
.rs.restartR()  # Restart R session
```

#### "Warning: No group information found"
**Cause**: Missing Sample_Group column in sample sheet
**Solution**: Add Sample_Group column or assign groups manually in the app

### Performance Optimization

#### For Large Datasets
1. **Use command line processing** instead of Shiny app
2. **Process in batches** if memory is limited
3. **Use functional normalization** for speed
4. **Consider subset analysis** for initial exploration

#### Memory Management
```r
# Check available memory
memory.limit()  # Windows
gc()  # Garbage collection

# Increase memory limit (Windows)
memory.limit(size = 8000)  # 8GB
```

### Imprinting Disorder Diagnosis

#### Beckwith-Wiedemann Syndrome (BWS)
**Key ICRs**: KCNQ1OT1 (ICR2), H19 (ICR1)
**Expected patterns**:
- ICR1: Hypomethylation (2-7% of BWS cases)
- ICR2: Hypermethylation (50-60% of BWS cases)
- Both: <1% of BWS cases

**GIMP Analysis**:
1. Focus on KCNQ1OT1 and H19 regions
2. Look for patient-specific hypermethylation/hypomethylation

#### Silver-Russell Syndrome (SRS)
**Key ICRs**: H19 (ICR1), KCNQ1OT1 (ICR2)
**Expected patterns**:
- ICR1: Hypomethylation (30-60% of SRS cases)
- ICR2: Hypomethylation (rare)

**GIMP Analysis**:
1. Use beta value heatmap for initial screening
2. Check H19 region specifically

#### Prader-Willi/Angelman Syndromes
**Key ICRs**: SNRPN
**Expected patterns**:
- Methylation status depends on parental origin
- Loss of methylation or gain of methylation

#### Temple Syndrome
**Key ICRs**: MEG3/DLK1
**Expected patterns**:
- Hypermethylation of MEG3 ICR

### Community Contributions

#### How to Contribute
1. **Report bugs**: Use GitHub issues with detailed descriptions
2. **Suggest features**: Describe use cases and requirements
3. **Submit code**: Follow R package development guidelines
4. **Share data**: Contribute example datasets
5. **Write documentation**: Improve user guides and tutorials

#### Development Roadmap
- **Core functionality**: Maintain and improve existing features
- **Research features**: Add advanced analysis methods
- **User experience**: Improve interface and workflows

## Conclusion

The GIMP Shiny app provides a comprehensive platform for genomic imprinting analysis, supporting both novice and expert users in analyzing methylation patterns at ICRs. By combining user-friendly interfaces with powerful analytical capabilities, GIMP enables researchers and clinicians to:

- **Process diverse data types**: From raw IDAT files to processed methylation data
- **Perform quality control**: Comprehensive QC metrics and visualizations
- **Identify patterns**: Detect imprinting disorders and research findings
- **Generate reports**: Publication-ready figures and clinical reports
- **Integrate workflows**: Connect with existing laboratory and research pipelines

For additional support, please consult the documentation, GitHub repository, or contact the development team.

---

**Last Updated**: July 2025  
**Version**: 0.2.0  
**Contact**: francesco.cecerengs@gmail.com  
**GitHub**: [https://github.com/ngsFC/GIMP](https://github.com/ngsFC/GIMP)  
**Contact**: saadatabu1996@gmail.com  
**GitHub**: [https://github.com/SAADAT-Abu](https://github.com/SAADAT-Abu)

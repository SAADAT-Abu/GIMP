% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/read_idat_zip.r
\name{read_idat_zip}
\alias{read_idat_zip}
\title{Read IDAT Files from ZIP Archive
This function extracts and processes IDAT files from a ZIP archive containing
IDAT files and a sample sheet, returning a beta value matrix ready for GIMP analysis.}
\usage{
read_idat_zip(
  zip_file,
  sample_sheet_name = "samplesheet.csv",
  array_type = c("EPIC", "450k", "EPICv2"),
  temp_dir = NULL,
  normalize_method = c("quantile", "SWAN", "funnorm", "noob"),
  detection_pval = 0.01,
  remove_failed_samples = TRUE,
  n_cores = NULL
)
}
\arguments{
\item{zip_file}{Path to ZIP file containing IDAT files and sample sheet}

\item{sample_sheet_name}{Name of the sample sheet file in the ZIP (default: "samplesheet.csv")}

\item{array_type}{Array type for annotation ("450k", "EPIC", "EPICv2")}

\item{temp_dir}{Temporary directory for extraction (default: creates temporary directory)}

\item{normalize_method}{Normalization method for minfi ("quantile", "SWAN", "funnorm", "noob")}

\item{detection_pval}{P-value threshold for detection (default: 0.01)}

\item{remove_failed_samples}{Remove samples with >10 percent failed probes (default: TRUE)}

\item{n_cores}{Number of CPU cores to use for parallel processing (default: NULL for sequential processing)}
}
\value{
A list containing:
\item{beta_matrix}{Beta value matrix ready for GIMP analysis}
\item{sample_info}{Sample information from the sample sheet}
\item{qc_metrics}{Quality control metrics}
\item{failed_samples}{Names of samples that failed QC}
}
\description{
Read IDAT Files from ZIP Archive
This function extracts and processes IDAT files from a ZIP archive containing
IDAT files and a sample sheet, returning a beta value matrix ready for GIMP analysis.
}
\examples{
# Read IDAT files from ZIP
idat_data <- read_idat_zip("my_methylation_data.zip", array_type = "EPIC")
beta_matrix <- idat_data$beta_matrix

# Use parallel processing with 4 cores
idat_data <- read_idat_zip("my_methylation_data.zip", array_type = "EPIC", n_cores = 4)

# Use with GIMP functions
ICRcpg <- make_cpgs(Bmatrix = beta_matrix, bedmeth = "v1")
}

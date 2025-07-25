% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ICRs_heatmap.R
\name{ICRs_heatmap}
\alias{ICRs_heatmap}
\title{Generate Heatmap of ICRs Methylation}
\usage{
ICRs_heatmap(
  df_ICR,
  sampleInfo,
  control_label = "Control",
  case_label = "Case",
  bedmeth = "v1",
  order_by = "cord",
  annotation_col = NULL,
  plot_type = "beta",
  sd_threshold = 3
)
}
\arguments{
\item{df_ICR}{A data frame or matrix containing methylation beta values for ICRs.}

\item{sampleInfo}{A vector indicating the group labels (e.g., "Control" and "Case") 
  for each sample in `df_ICR`.
Each element in `sampleInfo` should correspond to a sample in `df_ICR`.}

\item{control_label}{A character string specifying the label for the control 
group in `sampleInfo`. Default is `"Control"`.}

\item{case_label}{A character string specifying the label for the case 
group in `sampleInfo`. Default is `"Case"`.}

\item{bedmeth}{A character string specifying the BED data version for DMR 
coordinates. Options are `"v1"`, `"v2"`, or `"450k"`. Default is `"v1"`.}

\item{order_by}{A character string specifying the ordering rows in the heatmap. 
Options are `"cord"` for coordinates or `"meth"` for methylation values. 
Default is `"cord"`.}

\item{annotation_col}{A named list of colors for each unique value in 
`sampleInfo`. If `NULL`, default colors are assigned using the "viridis" 
palette. Default is `NULL`.}

\item{plot_type}{A character string specifying the type of heatmap to generate. 
Options are `"beta"` for beta values, `"delta"` for values normalized against 
controls, and `"defect"` for defect matrix based on standard deviations. 
Default is `"beta"`.}

\item{sd_threshold}{A numeric value specifying the standard deviation threshold 
for detecting defects in the defect matrix. Only used if `plot_type` is 
`"defect"`. Default is `3`.}
}
\value{
A heatmap plot visualizing methylation of ICRs.
}
\description{
This function generates a heatmap for visualizing methylation data of ICRs.
}
\examples{
# Example sampleInfo with "Case" and "Control" labels for each sample
sampleInfo <- c(rep("Case", 10), rep("Control", 10))
ICRs_heatmap(df_ICR = my_ICR_data, sampleInfo = sampleInfo, 
             annotation_col = list(Sample = c("darkgreen", "darkred")))
}

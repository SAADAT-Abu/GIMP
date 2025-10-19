#' BED EPICv2 probes
#'
#' This dataset contains the EPICv2 array probes coordinates.
#' Data is now hosted on AnnotationHub and can be accessed using the
#' internal helper function get_bed_data("v2").
#'
#' @docType data
#' @usage data(bedEPICv2)
#' @format BED file with chromosome, start, end, and probeID columns.
#' @keywords datasets
#' @note For Bioconductor submissions, this data is available via AnnotationHub.
#'   Use functions like make_cpgs(), make_ICRs(), and plot_cpgs_coverage()
#'   which automatically retrieve the data from AnnotationHub when needed.
#'
#' @examples
#' \dontrun{
#' # Data is now accessed via AnnotationHub
#' # Use GIMP functions directly - they handle data loading:
#' # ICRcpg <- make_cpgs(Bmatrix = your_beta_matrix, bedmeth = "v2")
#' }
"bedEPICv2"

#' Create the ICR Matrix
#'
#' @description
#' This function generates an ICR (Imprinted Control Region) matrix from a given beta matrix, using specified BED data for CpG mapping.
#' The ICR matrix provides data organized by CpG probes and samples.
#' The coordinates of the Human Imprinted regions are taken from \url{https://doi.org/10.1080/15592294.2016.1264561}
#' @param Bmatrix A data frame or matrix containing methylation beta values. Rows should represent CpG probes, and columns represent samples.
#' @param bedmeth A character string indicating the BED data version to use for CpG mapping. Options are `"v1"` (EPIC v1), `"v2"` (EPIC v2), or `"450k"` (450k array). Default is `"v1"`.
#' @return A data frame representing the ICR matrix, structured by CpG probes and samples.
#' @examples
#' \donttest{
#' # Create sample beta matrix for demonstration
#' set.seed(123)
#' n_probes <- 1000
#' n_samples <- 6
#' 
#' # Generate random probe IDs that might overlap with ICRs
#' sample_probes <- paste0("cg", sprintf("%08d", sample(1:50000000, n_probes)))
#' beta_matrix <- matrix(runif(n_probes * n_samples, 0.3, 0.8), 
#'                       nrow = n_probes, ncol = n_samples)
#' rownames(beta_matrix) <- sample_probes
#' colnames(beta_matrix) <- paste0("Sample_", 1:n_samples)
#' 
#' # Generate the ICR matrix with default BED version (EPIC v1)
#' ICRmatrix <- make_ICRs(Bmatrix = beta_matrix, bedmeth = "v1")
#' }
#' 
#' # Simple usage with your own data:
#' # ICRmatrix <- make_ICRs(Bmatrix = your_beta_matrix, bedmeth = "v1")
#' @export

make_ICRs <- function(Bmatrix, bedmeth = "v1") {
  
  # Load the appropriate bedmeth data based on the bedmeth input
  if (bedmeth == "v1") {
    message("Loading bedEPICv1...")
    data(bedEPICv1)
    bedmeth_data <- bedEPICv1
  } else if (bedmeth == "v2") {
    message("Loading bedEPICv2...")
    data(bedEPICv2)
    bedmeth_data <- bedEPICv2
  } else if (bedmeth == "450k") {
    message("Loading bed450k...")
    data(bed450k)
    bedmeth_data <- bed450k
  } else {
    stop("Invalid bedmeth input. Choose from 'v1', 'v2', or '450k'.")
  }
  
  # Load the appropriate ICRs data based on bedmeth input
  if (bedmeth == "v1" || bedmeth == "450k") {
    message("Loading DMRs.hg19...")
    data(DMRs.hg19)
    ICRs <- DMRs.hg19
  } else if (bedmeth == "v2") {
    message("Loading DMRs.hg38...")
    data(DMRs.hg38)
    ICRs <- DMRs.hg38
  }

  # Perform bed intersection between ICRs and bedmeth_data
  probeICR <- bed_intersect(ICRs, bedmeth_data) %>%
    mutate(chr = gsub("chr", "", chrom)) %>%
    mutate(chr = as.numeric(chr)) %>%
    dplyr::arrange(chr, start.x) %>%
    dplyr::select(probeID.y, start.y, ICR.x, start.x, end.x) %>%
    as.data.frame()

  colnames(probeICR) <- c("probeID", "cstart", "ICR", "start", "end")

  # Create df.ICR matrix
  df.ICR <- as.data.frame(Bmatrix) %>%
  rownames_to_column("probeID") %>%
  full_join(probeICR, by = "probeID") %>%
  filter(!is.na(ICR)) %>%
  filter(!is.na(.[[2]])) %>%  # Dynamically filter based on the first data column
  dplyr::select(-probeID,-start,-end) %>%
  group_by(ICR) %>%
  summarise_all(mean) %>%
  ungroup() %>%
  full_join(ICRs, by = "ICR") %>%
  mutate(chr = gsub("chr", "", chrom)) %>%
  mutate(chr = as.numeric(chr)) %>%
  dplyr::arrange(across(c(chr, start))) %>%
  column_to_rownames("ICR") %>%
  filter(!is.na(.[[1]])) %>%  # Dynamically filter based on the first column of the grouped data
  dplyr::select(-c("chrom", "start", "end", "chr", "cstart")) %>%
  as.data.frame()

  return(df.ICR)
}

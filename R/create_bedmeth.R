#' Create BED File Data from Methylation Array Annotations
#'
#' This function generates a BED-format data frame from Illumina Human Methylation annotation files.
#' The BED data includes chromosome, position, and probe ID information, and supports multiple annotation versions.
#'
#' @param version A character string specifying the annotation version to use.
#' Options include `"450k"` for 450k array, `"v1"` for the EPIC version1 and `"v2"` for EPIC version2.
#' Default is `"v1"`.
#'
#' @return A data frame in BED format containing columns:
#'   \item{chr}{Chromosome name.}
#'   \item{pos}{Position on the chromosome.}
#'   \item{probeID}{Unique identifier for each probe.}
#'   \item{end}{End position, which is the same as `pos` in this output.}
#' @examples
#' # Create BED-format data with the default version (EPIC v1)
#' bed_data <- create_bedmeth()
#' head(bed_data)  # View the first few rows
#' 
#' # Use a different annotation version if available
#' bed_data_v2 <- create_bedmeth(version = "v2")
#' @export

create_bedmeth <- function(version = "v1") {
  
  if (version == "v1") {
    anno <- getAnnotation("IlluminaHumanMethylationEPICanno.ilm10b4.hg19")
    bedmeth <- as.data.frame(anno) %>%
      rownames_to_column("probeID") %>%
      dplyr::select(chr, pos, probeID) %>%
      mutate(end = pos) %>%
      relocate(end, .after = pos)
    
    colnames(bedmeth) <- c("chrom", "start", "end", "probeID")
    
  } else if (version == "v2") {
    anno <- getAnnotation("IlluminaHumanMethylationEPICanno.ilm10b4.hg19")
    bedmeth <- as.data.frame(anno) %>%
      rownames_to_column("probeID") %>%
      dplyr::select(chr, pos, probeID) %>%
      mutate(end = pos) %>%
      relocate(end, .after = pos)
    
    colnames(bedmeth) <- c("chrom", "start", "end", "probeID")
  } else if (version == "450k") {
    anno <- getAnnotation("IlluminaHumanMethylation450kanno.ilmn12.hg19")
    bedmeth <- as.data.frame(anno) %>%
      rownames_to_column("probeID") %>%
      dplyr::select(chr, pos, probeID) %>%
      mutate(end = pos) %>%
      relocate(end, .after = pos)
    
    colnames(bedmeth) <- c("chrom", "start", "end", "probeID")
    } else {
    stop("Invalid version. Use '450k','v1' or 'v2'.")
  }
  
  return(bedmeth)
  }

bed450k <- create_bedmeth(version = "450k")
write.csv(bed450k, "bed450k.csv", row.names = FALSE)
  
bedEPICv1 <- create_bedmeth(version = "v1")
write.csv(bedEPICv1, "bedEPICv1.csv", row.names = FALSE)

bedEPICv2 <- create_bedmeth(version = "v2")
write.csv(bedEPICv2, "bedEPICv2.csv", row.names = FALSE)


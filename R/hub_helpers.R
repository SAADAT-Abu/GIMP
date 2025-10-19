#' Get BED Data from AnnotationHub
#'
#' @description
#' Internal helper function to retrieve BED annotation data from AnnotationHub.
#' Falls back to bundled data if AnnotationHub is not available.
#'
#' @param version Character string: "450k", "v1" (EPIC v1), or "v2" (EPIC v2)
#' @return Data frame with BED coordinates
#' @keywords internal
get_bed_data <- function(version) {
    # Try to use AnnotationHub first
    if (requireNamespace("AnnotationHub", quietly = TRUE)) {
        tryCatch({
            ah <- AnnotationHub::AnnotationHub()

            # Query for GIMP BED files
            bed_title <- switch(version,
                "450k" = "bed450k: Illumina 450k Array CpG Coordinates",
                "v1" = "bedEPICv1: Illumina EPIC v1 Array CpG Coordinates",
                "v2" = "bedEPICv2: Illumina EPIC v2 Array CpG Coordinates",
                stop("Invalid version: ", version)
            )

            query <- AnnotationHub::query(ah, c("GIMP", bed_title))

            if (length(query) > 0) {
                message("Loading BED data from AnnotationHub...")
                bed_data <- query[[1]]
                return(bed_data)
            }
        }, error = function(e) {
            message("AnnotationHub not available, using bundled data")
        })
    }

    # Fallback to bundled data
    message("Using bundled BED data (this may be removed in future versions)")
    bed_var <- switch(version,
        "450k" = "bed450k",
        "v1" = "bedEPICv1",
        "v2" = "bedEPICv2",
        stop("Invalid version: ", version)
    )

    # Load from package data
    data(list = bed_var, package = "GIMP", envir = environment())
    return(get(bed_var))
}

#' Generate Heatmap of Imprinted DMRs Methylation
#'
#' This function generates a heatmap for visualizing methylation data of Imprinted Differentially Methylated Regions (DMRs).
#' It allows custom color schemes for group annotations, with default colors from the "viridis" palette.
#'
#' @param df_ICR A data frame or matrix containing methylation beta values for Imprinted DMRs.
#' @param group_vector A vector indicating the group labels (e.g., "Control" and "Case") for each sample in `df_ICR`.
#' Each element in `group_vector` should correspond to a sample in `df_ICR`.
#' @param control_label A character string specifying the label for the control group in `group_vector`. Default is `"Control"`.
#' @param case_label A character string specifying the label for the case group in `group_vector`. Default is `"Case"`.
#' @param bedmeth A character string specifying the BED data version for DMR coordinates. Options are `"v1"`, `"v2"`, or `"450k"`. Default is `"v1"`.
#' @param cluster_by A character string specifying the clustering method for rows in the heatmap. Options are `"cord"` or `"meth"`. Default is `"cord"`.
#' @param annotation_col A named list of colors for each unique value in `group_vector`. If `NULL`, default colors are assigned using the "viridis" palette. Default is `NULL`.
#' @return A heatmap plot visualizing methylation of Imprinted DMRs.
#' @examples
#' # Example group_vector with "Case" and "Control" labels for each sample
#' group_vector <- c(rep("Case", 10), rep("Control", 10))
#' DMR_heatmap(df_ICR = my_ICR_data, group_vector = group_vector, annotation_col = list(Sample = c("darkgreen", "darkred")))
#' @export

iDMR_heatmap <- function(df_ICR, group_vector, control_label = "Control", case_label = "Case", bedmeth = "v1", cluster_by = "cord", annotation_col = NULL) {
  
  # Load required libraries
  library(pheatmap)
  library(viridisLite)
  
  # Load BED data based on bedmeth version
  if (bedmeth == "v1" || bedmeth == "450k") {
    data(DMRs.hg19)
    ICR_cord <- DMRs.hg19
    odr <- ICR_cord$ICR
  } else if (bedmeth == "v2") {
    data(DMRs.hg38)
    ICR_cord <- DMRs.hg38
    odr <- ICR_cord$ICR
  } else {
    stop("Invalid bedmeth version. Choose from 'v1', 'v2', or '450k'.")
  }
  
  # Ensure that group_vector has the same length as the number of samples (columns) in df_ICR
  if (length(group_vector) != ncol(df_ICR)) {
    stop("Length of 'group_vector' must match the number of samples (columns) in 'df_ICR'.")
  }
  
  group_vector <- factor(group_vector, levels = c(control_label, case_label))
  annotation_name <- if (!is.null(annotation_col)) names(annotation_col)[1] else "Sample"
  
  # Create the annotation data
  mat_col <- data.frame(group_vector)
  colnames(mat_col) <- annotation_name
  rownames(mat_col) <- colnames(df_ICR)
  
  unique_groups <- levels(group_vector)
  if (is.null(annotation_col)) {
    default_colors <- viridisLite::viridis(length(unique_groups))
    annotation_col <- list(setNames(default_colors, unique_groups))
    names(annotation_col) <- annotation_name
  } else {
    if (!is.list(annotation_col) || length(annotation_col[[1]]) != length(unique_groups)) {
      stop("The 'annotation_col' list must have the same number of colors as the unique values in 'group_vector'.")
    }
    names(annotation_col[[1]]) <- unique_groups
  }
  
  # Define heatmap palette and breaks
  paletteLength <- 100
  myColor <- colorRampPalette(c("#785EF0", "white", "#9a031e"))(paletteLength)
  myBreaks <- c(seq(0, 0.5, length.out = ceiling(paletteLength / 2) + 1), 
                seq(0.500001, 1, length.out = floor(paletteLength / 2)))
  
  # Determine row clustering based on cluster_by parameter
  row_clust <- if (cluster_by == "cord") FALSE else if (cluster_by == "meth") TRUE else {
    stop("Please select a valid 'cluster_by' parameter: 'cord' to cluster by coordinates or 'meth' to cluster by methylation values.")
  }
  
  # Generate the heatmap
  pheatmap(
    mat = df_ICR[odr, ],
    annotation_col = mat_col,
    annotation_colors = annotation_col,
    color = myColor,
    breaks = myBreaks,
    border_color = "grey",
    main = "Methylation of Imprinted DMRs",
    annotation_legend = TRUE,
    annotation_names_col = FALSE,
    annotation_names_row = FALSE,
    drop_levels = FALSE,
    fontsize = 8,
    cluster_rows = row_clust,
    cluster_cols = TRUE,
    clustering_distance_rows = "euclidean",
    clustering_distance_cols = "euclidean", 
    clustering_method = "ward.D2"
  )
}


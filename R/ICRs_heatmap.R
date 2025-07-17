#' Generate Heatmap of ICRs Methylation
#'
#' This function generates a heatmap for visualizing methylation data of ICRs.
#'
#' @param df_ICR A data frame or matrix containing methylation beta values for ICRs.
#' @param sampleInfo A vector indicating the group labels (e.g., "Control" and "Case") 
#'   for each sample in `df_ICR`.
#' Each element in `sampleInfo` should correspond to a sample in `df_ICR`.
#' @param control_label A character string specifying the label for the control 
#'   group in `sampleInfo`. Default is `"Control"`.
#' @param case_label A character string specifying the label for the case 
#'   group in `sampleInfo`. Default is `"Case"`.
#' @param bedmeth A character string specifying the BED data version for DMR 
#'   coordinates. Options are `"v1"`, `"v2"`, or `"450k"`. Default is `"v1"`.
#' @param order_by A character string specifying the ordering rows in the heatmap. 
#'   Options are `"cord"` for coordinates or `"meth"` for methylation values. 
#'   Default is `"cord"`.
#' @param annotation_col A named list of colors for each unique value in 
#'   `sampleInfo`. If `NULL`, default colors are assigned using the "viridis" 
#'   palette. Default is `NULL`.
#' @param plot_type A character string specifying the type of heatmap to generate. 
#'   Options are `"beta"` for beta values, `"delta"` for values normalized against 
#'   controls, and `"defect"` for defect matrix based on standard deviations. 
#'   Default is `"beta"`.
#' @param sd_threshold A numeric value specifying the standard deviation threshold 
#'   for detecting defects in the defect matrix. Only used if `plot_type` is 
#'   `"defect"`. Default is `3`.
#' @return A heatmap plot visualizing methylation of ICRs.
#' @import pheatmap
#' @import viridisLite
#' @import ggplotify
#' @examples
#' # Example sampleInfo with "Case" and "Control" labels for each sample
#' sampleInfo <- c(rep("Case", 10), rep("Control", 10))
#' ICRs_heatmap(df_ICR = my_ICR_data, sampleInfo = sampleInfo, 
#'              annotation_col = list(Sample = c("darkgreen", "darkred")))
#' @export

ICRs_heatmap <- function(df_ICR, sampleInfo, control_label = "Control", 
                         case_label = "Case", bedmeth = "v1", 
                         order_by = "cord", annotation_col = NULL, 
                         plot_type = "beta", sd_threshold = 3) {
  
  # Debug info
  message("ICRs_heatmap Debug Info:\n")
  message("df_ICR dimensions:", dim(df_ICR), "\n")
  message("sampleInfo length:", length(sampleInfo), "\n")
  message("sampleInfo values:", paste(sampleInfo, collapse = ", "), "\n")
  message("bedmeth:", bedmeth, "\n")
  message("plot_type:", plot_type, "\n")
  
  # Check if df_ICR is empty
  if (nrow(df_ICR) == 0 || ncol(df_ICR) == 0) {
    stop("df_ICR is empty. Please check your data processing steps.")
  }
  
  
  cls_distance <- "euclidean" # setting default cluster distance
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
  
  # Ensure sampleInfo length matches number of samples in input data
  if (length(sampleInfo) != ncol(df_ICR)) {
    stop("Length of 'sampleInfo' must match the number of samples (columns) in 'df_ICR'.")
  }
  # Ensure plot_type parameter is valid
  if (!plot_type %in% c("beta", "delta", "defect")) {
    stop("Invalid plot_type. Choose from 'beta', 'delta', or 'defect'.")
  }
  # Ensure order_by parameter is valid
  if (!order_by %in% c("cord", "meth")) {
    stop("Invalid order_by. Choose from 'cord' or 'meth'.")
  }
  
  sampleInfo <- factor(sampleInfo, levels = c(control_label, case_label))
  annotation_name <- if (!is.null(annotation_col)) names(annotation_col)[1] else "Sample"
  
  # Preparing annotation data
  mat_col <- data.frame(sampleInfo)
  colnames(mat_col) <- annotation_name
  rownames(mat_col) <- colnames(df_ICR)
  
  unique_groups <- levels(sampleInfo)
  if (is.null(annotation_col)) {
    default_colors <- viridisLite::viridis(length(unique_groups))
    annotation_col <- list(setNames(default_colors, unique_groups))
    names(annotation_col) <- annotation_name
  } else {
    if (!is.list(annotation_col) || 
        length(annotation_col[[1]]) != length(unique_groups)) {
      stop("The 'annotation_col' list must have the same number of colors as ", 
           "the unique values in 'sampleInfo'.")
    }
    names(annotation_col[[1]]) <- unique_groups
  }
  
  available_icrs <- intersect(odr, rownames(df_ICR))
  df_ICR_filtered <- df_ICR[available_icrs, , drop = FALSE]
  
  # Remove rows with all NA values
  rows_with_data <- apply(df_ICR_filtered, 1, function(x) !all(is.na(x)))
  df_ICR_filtered <- df_ICR_filtered[rows_with_data, , drop = FALSE]
  
  # Check if we have any data left
  if (nrow(df_ICR_filtered) == 0) {
    stop("No ICRs with methylation data found. Please check your data.")
  }
  
  # Define heatmap palette and breaks
  paletteLength <- 100
  
  if (plot_type == "beta") {
    # Define palette and breaks for beta values
    colorPalette <- colorRampPalette(c("#785EF0", "white", 
                                      "#9a031e"))(paletteLength)
    breaks <- c(seq(0, 0.5, length.out = ceiling(paletteLength / 2) + 1), 
                seq(0.500001, 1, length.out = floor(paletteLength / 2)))
    heatmap_data <- df_ICR_filtered
    main_title <- "Methylation of Imprinted DMRs (Beta)"
    
  } else if (plot_type == "delta") {
    # Calculate DeltaBeta matrix
    control_indices <- which(sampleInfo == control_label)
    control_means <- rowMeans(df_ICR_filtered[, control_indices, drop = FALSE], na.rm = TRUE)
    heatmap_data <- sweep(df_ICR_filtered, 1, control_means)
    colorPalette <- colorRampPalette(c("blue", "white", "red"))(paletteLength)
    breaks <- seq(-0.5, 0.5, length.out = paletteLength + 1)
    main_title <- "Delta Beta of Imprinted DMRs"
    
  } else if (plot_type == "defect") {
    # Calculate binary defect matrix based on sd_threshold
    df_mvalues <- log2(df_ICR_filtered / (1 - df_ICR_filtered))
    control_indices <- which(sampleInfo == control_label)
    control_mvalues <- df_mvalues[, control_indices, drop = FALSE]
    control_means_m <- rowMeans(control_mvalues, na.rm = TRUE)
    control_sds_m <- apply(control_mvalues, 1, sd, na.rm = TRUE)
    upper_threshold <- control_means_m + sd_threshold * control_sds_m
    lower_threshold <- control_means_m - sd_threshold * control_sds_m
    heatmap_data <- (df_mvalues > upper_threshold | df_mvalues < lower_threshold) * 1
    # Convert NA to 0 for defect matrix
    heatmap_data[is.na(heatmap_data)] <- 0
    colorPalette <- c("white", "black")
    breaks <- NULL  # No breaks needed for binary data
    main_title <- "Defect Matrix of Imprinted DMRs"
    cls_distance <- "binary" # Updating cluster distance to binary
    
  } else {
    stop("Invalid plot_type. Choose from 'beta', 'delta', or 'defect'.")
  }
  
  # Ordering rows
  row_clust <- if (order_by == "cord") FALSE else if (order_by == "meth") TRUE else {
    stop("Please select a valid 'order_by' parameter: 'cord' to order by coordinates or 'meth' to order by methylation values.")
  }
  
  # Generate the heatmap
  pheatmap(
    mat = heatmap_data,
    annotation_col = mat_col,
    annotation_colors = annotation_col,
    color = colorPalette,
    breaks = breaks,
    border_color = "grey",
    main = main_title,
    annotation_legend = TRUE,
    annotation_names_col = FALSE,
    annotation_names_row = FALSE,
    drop_levels = FALSE,
    fontsize = 8,
    cluster_rows = row_clust,
    cluster_cols = TRUE,
    clustering_distance_rows = if (plot_type == "defect") "binary" else "euclidean",
    clustering_distance_cols = cls_distance, 
    clustering_method = "ward.D2"
  ) %>% as.ggplot()
}

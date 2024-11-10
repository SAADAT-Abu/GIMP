#' Define ICR defects
#'
#' This function generates heatmaps and delta matrices to analyze ICR (Imprinted Control Regions) defects based on the input data.
#'
#' @param df_ICR A data frame or matrix containing ICR data. Rows represent individual probes or features, and columns represent samples.
#' @param group_vector A vector indicating group labels (e.g., "Control" and "Case") for each sample in `df_ICR`.
#' @param control_label A character string specifying the label for the control group in `group_vector`. Default is "Control".
#' @return A list containing:
#'   \item{heatmap_ICR}{A heatmap displaying the ICR data for visualization.}
#'   \item{df_ICR_delta}{A delta matrix showing differences in beta values between groups.}
#'   \item{df_ICR_defect}{A binary matrix indicating the presence of ICR defects.}
#' @examples
#' # Assuming df_ICR is your matrix and you have a group vector like c("Control", "Case", "Control", "Case", ...)
#' result <- analyze_ICR(df_ICR, group_vector)
#' # View the heatmaps and matrices:
#' result$heatmap_ICR  # Heatmap of df_ICR
#' result$df_ICR_delta  # DeltaBeta matrix
#' result$df_ICR_defect  # Binary defect matrix
#' @export

# Function to create heatmaps and delta matrices
analyze_ICR <- function(df_ICR, group_vector, control_label = "Control", case_label = "Case") {

   # Load required library
  library(ggplot2)
  library(reshape2)
  
  # Ensure the group_vector is a factor
  group_vector <- factor(group_vector, levels = c(control_label, case_label))
  
  # Step 1: Plot a heatmap of df.ICR with Case/Control annotations, scale 0 to 1
  plot_heatmap <- function(df, group_vector, title = "Heatmap of df.ICR", scale_limits = c(0, 1), color_scheme = "beta") {
    df_melt <- melt(as.matrix(df))
    colnames(df_melt) <- c("ICR", "Sample", "BetaValue")
    
    # Add the group annotations (Case/Control)
    df_melt$Group <- group_vector[df_melt$Sample]
    
    # Adjust color scheme based on the heatmap type
    if (color_scheme == "delta") {
      # Use the same color scheme as Heatmap 1 (low: blue, mid: white, high: red)
      color_scale <- scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0, limits = scale_limits)
    } else {
      # Default for beta values (Heatmap 1)
      color_scale <- scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0.5, limits = scale_limits)
    }
    
    # Plot the heatmap
    ggplot(df_melt, aes(x = Sample, y = ICR, fill = BetaValue)) +
      geom_tile() +
      color_scale +
      labs(title = title, x = "Sample", y = "ICR", fill = "Beta Value") +
      theme_minimal() +
      theme(axis.text.x = element_text(angle = 90, hjust = 1))
  }
  
  # Step 2: Calculate DeltaBeta (df.ICR.delta) with scale -0.5 to 0.5
  calculate_delta_beta <- function(df, group_vector, control_label) {
    control_indices <- which(group_vector == control_label)
    
    # Calculate the mean for controls for each ICR (row-wise)
    control_means <- rowMeans(df[, control_indices, drop = FALSE])
    
    # Subtract control means from the entire matrix to get DeltaBeta
    df_delta <- sweep(df, 1, control_means)
    
    return(df_delta)
  }
  
  # Step 3: Create a binary defect matrix (df.ICR.defect) based on M-values
  calculate_defect_matrix <- function(df, group_vector, control_label) {
    control_indices <- which(group_vector == control_label)
    
    # Convert Beta values to M-values for controls and entire matrix
    df_mvalues <- log2(df/(1-df))
    control_mvalues <- df_mvalues[, control_indices, drop = FALSE]
    
    # Calculate the mean and 3*SD for controls (in M-value space)
    control_means <- rowMeans(control_mvalues)
    control_sds <- apply(control_mvalues, 1, sd)
    
    # Define thresholds for binary defect matrix
    upper_threshold <- control_means + 3 * control_sds
    lower_threshold <- control_means - 3 * control_sds
    
    # Create a binary matrix where BetaValue is outside the thresholds
    df_defect <- (df_mvalues > upper_threshold | df_mvalues < lower_threshold) * 1  # Binary matrix
    
    return(df_defect)
  }
  
  # Step 1: Heatmap of df.ICR with Case/Control annotations (scale 0 to 1)
  print("Plotting heatmap of df.ICR")
  heatmap_1 <- plot_heatmap(df_ICR, group_vector, title = "Heatmap of df.ICR", scale_limits = c(0, 1), color_scheme = "beta")
  
  # Step 2: Calculate DeltaBeta matrix and plot heatmap (scale -0.5 to 0.5) with the same color scheme as Heatmap 1
  df_ICR_delta <- calculate_delta_beta(df_ICR, group_vector, control_label)
  print("Plotting heatmap of DeltaBeta (df_ICR_delta)")
  heatmap_2 <- plot_heatmap(df_ICR_delta, group_vector, title = "Heatmap of DeltaBeta (df_ICR_delta)", scale_limits = c(-0.5, 0.5), color_scheme = "delta")
  
  # Step 3: Calculate binary defect matrix and plot heatmap
  df_ICR_defect <- calculate_defect_matrix(df_ICR, group_vector, control_label)
  print("Plotting heatmap of binary defect matrix (df_ICR_defect)")
  heatmap_3 <- plot_heatmap(df_ICR_defect, group_vector, title = "Heatmap of Defect Matrix (df_ICR_defect)", scale_limits = c(0, 1), color_scheme = "beta")
  
  # Return the three matrices and their heatmaps
  return(list(
    df_ICR_delta = df_ICR_delta,
    df_ICR_defect = df_ICR_defect,
    heatmap_ICR = heatmap_1,
    heatmap_delta = heatmap_2,
    heatmap_defect = heatmap_3
  ))
}



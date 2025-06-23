#' Plot ICR CpG Matrix with Counts and Percentage Coverage
#'
#' This function plots the CpG coverage for Imprinted Control Regions (ICRs) using the provided data frame of CpG counts.
#' It compares CpG counts in the specified BED data version for visual analysis and includes an additional plot for percentage coverage.
#'
#' @param df_ICR_cpg A data frame containing CpG counts for ICR regions. Each row represents a different CpG probe, and columns contain sample-related information.
#' @param bedmeth A character string specifying the BED data version to use for mapping CpG coverage. Options are `"v1"` (EPIC v1), `"v2"` (EPIC v2), or `"450k"` (450k array). Default is `"v1"`.
#' @return A list containing two plots (counts and percentage coverage) and the data frame with CpG counts and coverage information.
#' @examples
#' plot_cpgs_coverage(df_ICR_cpg_counts, bedmeth = "v1")
#' @export

plot_cpgs_coverage <- function(df_ICR_cpg, bedmeth = "v1") {

  # Load the bedmeth data based on the version specified
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
    stop("Invalid bedmeth version. Choose from 'v1', 'v2', or '450k'.")
  }

  # Count the number of CpGs per ICR in df_ICR_cpg (input data)
  df_ICR_cpg_counts <- df_ICR_cpg %>%
    group_by(ICR) %>%
    summarize(n_cpgs = n())

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

  bedmeth_counts <- bed_intersect(ICRs, bedmeth_data) %>%
    mutate(chr = gsub("chr", "", chrom)) %>%
    mutate(chr = as.numeric(chr)) %>%
    arrange(chr, start.x) %>%
    dplyr::select(probeID.y, start.y, ICR.x, start.x, end.x) %>%
    as.data.frame() %>%
    dplyr::count(ICR.x, sort = TRUE) %>%
    dplyr::rename("ICR" = "ICR.x",
                  "Total_cov" = "n")

  comparison_data <- full_join(df_ICR_cpg_counts, bedmeth_counts, by = "ICR")

  # Add a column for percentage coverage
  comparison_data <- comparison_data %>%
    mutate(Percentage_cov = (n_cpgs / Total_cov) * 100) %>%
    replace_na(list(n_cpgs = 0, Total_cov = 0, Percentage_cov = 0))

  # Reshape data for the count plot
  comparison_data_long <- comparison_data %>%
    pivot_longer(cols = c(n_cpgs, Total_cov), names_to = "Measure", values_to = "Count")

  # Plot 1: Number of CpGs Covered per ICR with dashed grid lines
  plot_counts <- ggplot(comparison_data_long, aes(x = reorder(ICR, Count), y = Count, fill = Measure)) +
    geom_bar(stat = "identity", position = position_dodge(width = 0.5)) +
    labs(title = "Number of CpGs Covered per ICR", x = "ICR", y = "Number of CpGs") +
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
    theme_bw() + 
    coord_flip() +
    scale_y_continuous(breaks = seq(0, max(comparison_data_long$Count, na.rm = TRUE), 10)) +
    theme(panel.grid.major.x = element_line(color = "grey", linetype = "dashed"))

  # Plot 2: Percentage of CpGs Covered per ICR
  plot_percentage <- ggplot(comparison_data, aes(x = reorder(ICR, Percentage_cov), y = Percentage_cov)) +
    geom_bar(stat = "identity", fill = "skyblue") +
    labs(title = "Percentage of CpGs Covered per ICR", x = "ICR", y = "Percentage of CpGs Covered (%)") +
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
    theme_bw() +
    coord_flip() +
    scale_y_continuous(breaks = seq(0, 100, 10)) +
    theme(panel.grid.major.x = element_line(color = "grey", linetype = "dashed"))

  # Return both plots and the data
  return(list(plot_counts = plot_counts, plot_percentage = plot_percentage, data = comparison_data))
}


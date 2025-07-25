#' Plot Line Plot for ICR Methylation
#'
#' This function generates a line plot to visualize methylation values across a specified ICR.
#' Users can choose between a static `ggplot2` plot or an interactive `plotly` plot.
#'
#' @param significantDMPs A data frame containing information about significant DMPs. Must include columns `ICR`, `start`, and `end`.
#' @param ICRcpg A data frame or matrix containing CpG methylation data. Includes CpG coordinates (`cstart`) and methylation values.
#' @param ICR A character string specifying the name of the ICR region to be plotted.
#' @param sampleInfo A character vector providing group labels (e.g., `"Control"` or `"Case"`) for each sample in the methylation data.
#' @param interactive A logical value indimessageing whether to return an interactive `plotly` plot (`TRUE`) or a static `ggplot2` plot (`FALSE`). Default is `TRUE`.
#' @return A plot representing the line plot of methylation values across the specified ICR region, highlighting significant DMPs. 
#' The plot is either a `ggplot2` object or a `plotly` object, depending on the value of `interactive`.
#' @examples
#' # Example data for significantDMPs
#' plot <- plot_line_ICR(significantDMPs, ICRcpg, ICR = "KCNQ1OT1:TSS-DMR", sampleInfo = sampleInfo, interactive = T)
#' print(plot)
#' @export
plot_line_ICR <- function(significantDMPs, ICRcpg, ICR, sampleInfo, interactive = TRUE) {

  # Input validation
  if (!is.data.frame(significantDMPs) || nrow(significantDMPs) == 0) {
    stop("significantDMPs must be a non-empty data frame.")
  }

  if (!is.data.frame(ICRcpg) || nrow(ICRcpg) == 0) {
    stop("ICRcpg must be a non-empty data frame.")
  }

  if (!ICR %in% significantDMPs$ICR) {
    stop(paste("ICR", ICR, "not found in significantDMPs. Available ICRs:",
               paste(unique(significantDMPs$ICR), collapse = ", ")))
  }

  # Filter DMPs for the specified ICR
  regionDMPs <- significantDMPs[significantDMPs$ICR == ICR, ]
  if (nrow(regionDMPs) == 0) {
    stop("No significant DMPs found for the specified ICR.")
  }

  message("Found ", nrow(regionDMPs), " DMPs for ICR ", ICR, "\n")

  # Get region boundaries
  regionStart <- min(regionDMPs$start, na.rm = TRUE)
  regionEnd <- max(regionDMPs$end, na.rm = TRUE)

  message("Region boundaries: ", regionStart, " to ", regionEnd, "\n")

  # Filter CpGs for the specified ICR region
  required_cols <- c("cstart", "ICR", "start", "end")
  missing_cols <- setdiff(required_cols, colnames(ICRcpg))
  if (length(missing_cols) > 0) {
    stop("Missing required columns in ICRcpg: ", paste(missing_cols, collapse = ", "))
  }

  # Filter by ICR name first, then by coordinates
  regionCpGs <- ICRcpg[ICRcpg$ICR == ICR, ]

  if (nrow(regionCpGs) == 0) {
    stop("No CpGs found for ICR ", ICR, " in ICRcpg data.")
  }

  message("Found ", nrow(regionCpGs), " CpGs for ICR ", ICR, "\n")

  # Identify methylation data columns (all except the last 4 annotation columns)
  annotation_cols <- c("cstart", "ICR", "start", "end")
  methylation_col_indices <- which(!colnames(regionCpGs) %in% annotation_cols)

  if (length(methylation_col_indices) == 0) {
    stop("No methylation data columns found in ICRcpg.")
  }

  # Check sampleInfo length matches methylation data
  if (length(sampleInfo) != length(methylation_col_indices)) {
    stop("sampleInfo length (", length(sampleInfo),
         ") doesn't match number of samples (", length(methylation_col_indices), ")")
  }

  # Extract methylation data and annotation
  methylationData <- regionCpGs[, methylation_col_indices, drop = FALSE]
  annotationData <- regionCpGs[, annotation_cols, drop = FALSE]

  # Remove rows with all NA values in methylation data
  valid_rows <- apply(methylationData, 1, function(x) !all(is.na(x)))
  if (sum(valid_rows) == 0) {
    stop("No valid methylation data found for the specified ICR region.")
  }

  methylationData <- methylationData[valid_rows, , drop = FALSE]
  annotationData <- annotationData[valid_rows, , drop = FALSE]

  message("After filtering, have ", nrow(methylationData), " valid CpGs\n")

  # Prepare data for plotting
  methylationData$cstart <- annotationData$cstart

  # Reshape data for plotting with better error handling
  tryCatch({
    methylationDataLong <- reshape2::melt(methylationData,
                                          id.vars = "cstart",
                                          variable.name = "Sample",
                                          value.name = "Methylation")
  }, error = function(e) {
    stop("Error reshaping data: ", e$message, call. = FALSE)
  })

  # Add sample type information - FIXED: proper mapping of sample to group
  # Create a mapping dataframe from sample names to groups
  sample_names <- colnames(methylationData)[colnames(methylationData) != "cstart"]
  sample_group_mapping <- data.frame(
    Sample = sample_names,
    SampleGroup = sampleInfo,
    stringsAsFactors = FALSE
  )

  # Merge to get correct group assignment for each sample
  methylationDataLong <- merge(methylationDataLong, sample_group_mapping, by = "Sample", all.x = TRUE)

  # Mark which CpGs are DMPs
  methylationDataLong$IsDMP <- methylationDataLong$cstart %in% regionDMPs$cstart

  # Remove rows with NA methylation values for plotting
  methylationDataLong <- methylationDataLong[!is.na(methylationDataLong$Methylation), ]

  if (nrow(methylationDataLong) == 0) {
    stop("No valid data points for plotting after removing NAs.")
  }

  message("Final data for plotting has ", nrow(methylationDataLong), " points\n")

  # Create the plot with DYNAMIC COLOR MAPPING
  tryCatch({

    # Extract unique labels from sampleInfo
    unique_labels <- unique(sampleInfo)

    # Validate that we have exactly 2 groups
    if (length(unique_labels) != 2) {
      stop(paste("Expected exactly 2 groups for plotting, but found", length(unique_labels),
                 "groups:", paste(unique_labels, collapse = ", ")))
    }

    # Determine which label should be control and which should be case
    if ("Control" %in% unique_labels) {
      control_label <- "Control"
      case_label <- setdiff(unique_labels, "Control")
    } else if ("control" %in% unique_labels) {
      control_label <- "control"
      case_label <- setdiff(unique_labels, "control")
    } else {
      sorted_labels <- sort(unique_labels)
      control_label <- sorted_labels[1]
      case_label <- sorted_labels[2]
    }

    # Create dynamic color mapping
    color_mapping <- c("yellow3", "purple4")
    names(color_mapping) <- c(control_label, case_label)

    # Create the plot - FIXED: group by Sample to get one line per sample
    plot <- ggplot(methylationDataLong, aes(x = cstart, y = Methylation, group = Sample, color = SampleGroup)) +
      geom_line(linetype = 1, alpha = 0.5, size = 0.5) +  # Now each sample gets its own line
      geom_point(data = methylationDataLong[!methylationDataLong$IsDMP, ],
                 size = 1, alpha = 0.7) +
      geom_point(data = methylationDataLong[methylationDataLong$IsDMP, ],
                 size = 1.5, alpha = 0.9) +
      geom_vline(data = data.frame(x = unique(regionDMPs$cstart)),
                 aes(xintercept = x), color = "lightgrey", linetype = "dashed", alpha = 0.7) +
      geom_rug(data = methylationDataLong[methylationDataLong$IsDMP, ],
               aes(x = cstart), color = "red", sides = "t", alpha = 0.8) +

      # DYNAMIC SCALES
      scale_color_manual("Group", values = color_mapping) +

      labs(
        title = paste("Methylation Profile:", ICR),
        subtitle = paste("Showing", sum(methylationDataLong$IsDMP), "significant DMPs |",
                         control_label, "(Yellow) vs", case_label, "(Purple)"),
        x = "CpG Coordinates (bp)",
        y = "Methylation Value (Beta)"
      ) +
      theme_minimal() +
      theme(
        legend.position = "bottom",
        plot.title = element_text(size = 14, face = "bold"),
        plot.subtitle = element_text(size = 12)
      ) +
      ylim(0, 1)  # Ensure y-axis shows full beta value range

    # Convert to interactive plot if requested - FIXED: corrected layout syntax
    if (interactive) {
      tryCatch({
        plot <- ggplotly(plot, tooltip = c("x", "y", "colour", "group"))
        # Apply layout modifications separately if needed
        # Note: ggplotly usually preserves the ggplot title and legend automatically
      }, error = function(e) {
        warning("Failed to create interactive plot: ", e$message, ". Returning static plot.")
        interactive <- FALSE  # Fall back to static
      })
    }

    return(plot)

  }, error = function(e) {
    stop("Error creating plot: ", e$message, call. = FALSE)
  })
}

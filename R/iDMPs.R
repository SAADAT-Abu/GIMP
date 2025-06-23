#' Identify Differentially Methylated Positions in ICRs
#'
#' This function identifies differentially methylated positions (DMPs) between
#' control and case groups using linear modeling and empirical Bayes methods.
#'
#' @param data A data frame containing CpG methylation data with annotation columns
#' @param sampleInfo A factor or character vector indicating sample groups
#' @param pValueCutoff P-value threshold for significance (default: 0.05)
#' @return A list containing:
#'   \item{fit}{Linear model fit object}
#'   \item{eBayesfit}{Empirical Bayes fit object}
#'   \item{topDMPs}{Data frame of significant DMPs}
#'   \item{allResults}{Data frame of all results}
#'   \item{groupLabels}{Group labels used in analysis}
#' @examples
#' # Run DMP analysis
#' dmps <- iDMPs(data = ICRcpg, sampleInfo = sample_groups)
#' significant_dmps <- dmps$topDMPs
#' @export
iDMPs <- function(data, sampleInfo, pValueCutoff = 0.05) {
  
  sampleInfo <- factor(sampleInfo)
  
  if (!"Control" %in% levels(sampleInfo)) {
    stop("The 'Control' group must be present in sampleInfo.")
  }
  sampleInfo <- relevel(sampleInfo, ref = "Control")
  
  if (!is.data.frame(data)) {
    stop("The input data must be a data frame.")
  }
  
  # Identify annotation columns more robustly ***
  # Look for the expected annotation columns
  expected_annotation_cols <- c("cstart", "ICR", "start", "end")
  annotation_col_indices <- which(colnames(data) %in% expected_annotation_cols)
  
  if (length(annotation_col_indices) == 0) {
    stop("No annotation columns found. Expected columns: cstart, ICR, start, end")
  }
  
  # Separate methylation data from annotation properly ***
  # All columns except the annotation columns are methylation data
  methylation_col_indices <- setdiff(1:ncol(data), annotation_col_indices)
  
  if (length(methylation_col_indices) == 0) {
    stop("No methylation data columns found.")
  }
  
  methylationData <- as.matrix(data[, methylation_col_indices])
  additionalColumns <- data[, annotation_col_indices, drop = FALSE]
  
  # Ensure we have the expected annotation columns
  missing_cols <- setdiff(expected_annotation_cols, colnames(additionalColumns))
  if (length(missing_cols) > 0) {
    stop(paste("Missing annotation columns:", paste(missing_cols, collapse = ", ")))
  }
  
  # Check dimensions match ***
  if (ncol(methylationData) != length(sampleInfo)) {
    stop(paste("Number of samples in methylation data (", ncol(methylationData), 
               ") doesn't match sampleInfo length (", length(sampleInfo), ")"))
  }
  
  # Improved beta to M-value conversion with error handling ***
  betaToM <- function(beta) {
    # Handle edge cases
    beta[beta <= 0] <- 0.001
    beta[beta >= 1] <- 0.999
    log2(beta / (1 - beta))
  }
  
  # Remove rows with too many NAs ***
  # Remove CpGs that have more than 50% missing values
  na_proportion <- apply(methylationData, 1, function(x) sum(is.na(x)) / length(x))
  valid_rows <- na_proportion < 0.5
  
  if (sum(valid_rows) == 0) {
    stop("No CpGs with sufficient data found.")
  }
  
  cat("Removing", sum(!valid_rows), "CpGs with >50% missing values\n")
  cat("Analyzing", sum(valid_rows), "CpGs\n")
  
  methylationData <- methylationData[valid_rows, ]
  additionalColumns <- additionalColumns[valid_rows, ]
  
  # Convert to M-values
  mValues <- betaToM(methylationData)
  
  # Handle remaining NAs in M-values ***
  # Impute remaining NAs with row means
  for (i in 1:nrow(mValues)) {
    if (any(is.na(mValues[i, ]))) {
      row_mean <- mean(mValues[i, ], na.rm = TRUE)
      if (!is.na(row_mean)) {
        mValues[i, is.na(mValues[i, ])] <- row_mean
      }
    }
  }
  
  # Create design matrix
  design <- model.matrix(~sampleInfo)
  
  # Add row names to ensure proper matching ***
  rownames(mValues) <- rownames(data)[valid_rows]
  
  # Fit linear model
  fit <- lmFit(mValues, design)
  eBayesfit <- eBayes(fit)
  
  groupLabels <- levels(sampleInfo)
  coefName <- paste0("sampleInfo", groupLabels[2])  # Select second group relative to "Control"
  
  if (!coefName %in% colnames(fit$coefficients)) {
    stop(paste("Coefficient", coefName, "not found in the fit model. Check group labels."))
  }
  
  # Get all results first, then filter ***
  allResults <- topTable(eBayesfit, number = Inf, coef = coefName, sort.by = "P")
  
  # Add annotation columns back to results
  # Match by row names
  matched_indices <- match(rownames(allResults), rownames(additionalColumns))
  topDMPs_with_annotation <- cbind(allResults, additionalColumns[matched_indices, ])
  
  # Filter by p-value
  significantDMPs <- topDMPs_with_annotation[topDMPs_with_annotation$P.Value <= pValueCutoff, ]
  
  cat("Found", nrow(significantDMPs), "significant DMPs at p <", pValueCutoff, "\n")
  
  return(list(
    fit = fit,
    eBayesfit = eBayesfit,
    topDMPs = significantDMPs,
    allResults = topDMPs_with_annotation,  # Include all results for reference
    groupLabels = groupLabels
  ))
}
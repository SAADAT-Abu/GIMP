#' Validate IDAT ZIP File Contents
#'
#' @description
#' Helper function to validate the contents of a ZIP file before processing
#' IDAT files. Checks for required files and proper structure.
#'
#' @param zip_file Path to ZIP file
#' @param sample_sheet_name Expected sample sheet filename
#' @return List with validation results
#' @keywords internal

validate_idat_zip <- function(zip_file, sample_sheet_name = "samplesheet.csv") {
  
  if (!file.exists(zip_file)) {
    return(list(valid = FALSE, message = "ZIP file not found"))
  }
  
  # Get ZIP contents without extracting
  zip_contents <- unzip(zip_file, list = TRUE)
  
  # Check for IDAT files
  idat_files <- grep("\\.idat$", zip_contents$Name, value = TRUE)
  
  if (length(idat_files) == 0) {
    return(list(valid = FALSE, message = "No IDAT files found in ZIP"))
  }
  
  # Check for sample sheet
  sample_sheet_found <- any(grepl(sample_sheet_name, zip_contents$Name, ignore.case = TRUE))
  
  if (!sample_sheet_found) {
    # Look for alternative sample sheet names
    alt_sheets <- c("SampleSheet.csv", "sample_sheet.csv", "samples.csv", "metadata.csv")
    alt_found <- vapply(alt_sheets, function(x) any(grepl(x, zip_contents$Name, ignore.case = TRUE)))
    
    if (any(alt_found)) {
      found_sheet <- names(alt_found)[which(alt_found)[1]]
      return(list(valid = TRUE, 
                 message = paste("Found alternative sample sheet:", found_sheet),
                 idat_count = length(idat_files),
                 sample_sheet = found_sheet))
    } else {
      return(list(valid = FALSE, 
                 message = paste("Sample sheet not found. Expected:", sample_sheet_name)))
    }
  }
  
  # Estimate number of samples (each sample should have Red and Grn IDAT)
  estimated_samples <- length(idat_files) / 2
  
  return(list(
    valid = TRUE,
    message = "ZIP file validation passed",
    idat_count = length(idat_files),
    estimated_samples = floor(estimated_samples),
    sample_sheet = sample_sheet_name
  ))
}

#' Create Sample Sheet Template
#'
#' @description
#' Creates a template sample sheet for IDAT files to help users format their data correctly.
#'
#' @param sample_names Vector of sample names
#' @param sentrix_ids Vector of Sentrix IDs (slide IDs)
#' @param sentrix_positions Vector of Sentrix positions
#' @param groups Optional vector of sample groups
#' @return Data frame with sample sheet structure
#' @examples
#' # Create a basic template
#' template <- create_sample_sheet_template(
#'   sample_names = c("Sample1", "Sample2", "Sample3"),
#'   sentrix_ids = c("200123456", "200123456", "200123457"),
#'   sentrix_positions = c("R01C01", "R02C01", "R01C01"),
#'   groups = c("Control", "Control", "Case")
#' )
#' print(template)
#' 
#' # Create template with default values (minimal input)
#' simple_template <- create_sample_sheet_template(
#'   sample_names = c("Ctrl_01", "Case_01")
#' )
#' head(simple_template)
#' @export

create_sample_sheet_template <- function(sample_names = NULL, 
                                       sentrix_ids = NULL, 
                                       sentrix_positions = NULL,
                                       groups = NULL) {
  
  # Validate input lengths
  if (!is.null(sentrix_ids) && length(sentrix_ids) != length(sample_names)) {
    stop("sentrix_ids length must match sample_names length")
  }
  
  if (!is.null(sentrix_positions) && length(sentrix_positions) != length(sample_names)) {
    stop("sentrix_positions length must match sample_names length")
  }
  
  # Create template
  template <- data.frame(
    Sample_Name = sample_names,
    Sentrix_ID = sentrix_ids %||% paste0("20012345678", seq_len(length(sample_names))),
    Sentrix_Position = sentrix_positions %||% paste0("R0", ((seq_len(length(sample_names)) - 1) %% 6) + 1, "C01"),
    stringsAsFactors = FALSE
  )
  
  # Add group information if provided
  if (!is.null(groups)) {
    template$Sample_Group <- groups
  }
  
  # Add additional useful columns
  template$Sample_Plate <- "Plate1"
  template$Sample_Well <- paste0(LETTERS[((seq_len(nrow(template)) - 1) %/% 12) + 1], 
                                sprintf("%02d", ((seq_len(nrow(template)) - 1) %% 12) + 1))
  
  return(template)
}

# Helper operator for NULL coalescing
`%||%` <- function(x, y) if (is.null(x)) y else x

#' Preview IDAT ZIP Contents
#'
#' @description
#' Preview the contents of an IDAT ZIP file without processing it.
#' Useful for checking file structure before full processing.
#'
#' @param zip_file Path to ZIP file
#' @return List with ZIP contents summary
#' @examples
#' \donttest{
#' # Preview IDAT ZIP file contents
#' # preview_info <- preview_idat_zip("methylation_data.zip")
#' # print(preview_info)
#' }
#' @export

preview_idat_zip <- function(zip_file) {
  
  validation <- validate_idat_zip(zip_file)
  
  if (!validation$valid) {
    return(validation)
  }
  
  # Get detailed ZIP contents
  zip_contents <- unzip(zip_file, list = TRUE)
  
  # Categorize files
  idat_files <- grep("\\.idat$", zip_contents$Name, value = TRUE)
  csv_files <- grep("\\.csv$", zip_contents$Name, value = TRUE)
  other_files <- setdiff(zip_contents$Name, c(idat_files, csv_files))
  
  # Analyze IDAT files
  red_files <- grep("_Red\\.idat$", idat_files, value = TRUE)
  grn_files <- grep("_Grn\\.idat$", idat_files, value = TRUE)
  
  # Extract sample identifiers from IDAT files
  sample_ids <- unique(gsub("_(Red|Grn)\\.idat$", "", basename(idat_files)))
  
  return(list(
    valid = TRUE,
    total_files = nrow(zip_contents),
    idat_files = length(idat_files),
    red_files = length(red_files),
    grn_files = length(grn_files),
    csv_files = csv_files,
    estimated_samples = length(sample_ids),
    sample_identifiers = head(sample_ids, 10),  # Show first 10
    other_files = other_files,
    total_size_mb = round(sum(zip_contents$Length) / (1024^2), 2)
  ))
}

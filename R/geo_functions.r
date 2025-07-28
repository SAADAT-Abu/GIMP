#' Check if GEO ID exists and has IDAT data
#'
#' @description
#' Validates a GEO accession ID and checks if it contains raw IDAT files
#' for methylation array analysis.
#'
#' @param geo_id Character string of GEO accession (e.g., "GSE12345")
#' @return List with validation results and metadata
#' @examples
#' \donttest{
#' # Check a GEO dataset
#' result <- validate_geo_dataset("GSE68777")
#' if (result$valid && result$has_idats) {
#'   message("Dataset is suitable for GIMP analysis")
#' }
#' }
#' @export
validate_geo_dataset <- function(geo_id) {

  if (!requireNamespace("GEOquery", quietly = TRUE)) {
    stop("GEOquery package is required. Install with: BiocManager::install('GEOquery')")
  }

  # Clean up GEO ID format
  geo_id <- toupper(trimws(geo_id))
  if (!grepl("^GSE\\d+$", geo_id)) {
    return(list(
      valid = FALSE,
      message = paste("Invalid GEO ID format. Expected format: GSExxxxx, got:", geo_id)
    ))
  }

  message("Validating GEO dataset: ", geo_id)

  tryCatch({
    # Get GSE metadata without downloading expression data
    gse <- GEOquery::getGEO(geo_id, GSEMatrix = FALSE)

    # Extract basic information
    title <- GEOquery::Meta(gse)$title
    summary <- GEOquery::Meta(gse)$summary
    sample_count <- length(GEOquery::GSMList(gse))

    message("Dataset found: ", title)
    message("Sample count: ", sample_count)

    # Check if it's a methylation study
    is_methylation <- any(grepl("methylation|methyl|450k|epic",
                               c(title, summary), ignore.case = TRUE))

    if (!is_methylation) {
      return(list(
        valid = TRUE,
        has_idats = FALSE,
        message = "Dataset found but doesn't appear to be a supported methylation study",
        title = title,
        sample_count = sample_count
      ))
    }

    # Get sample information to check for IDAT availability
    gsm_list <- GEOquery::GSMList(gse)
    sample_info <- lapply(gsm_list[1:min(3, length(gsm_list))], function(gsm) {
      meta <- GEOquery::Meta(gsm)
      list(
        title = meta$title %||% "Unknown",
        platform = meta$platform_id %||% "Unknown",
        supplementary = meta$supplementary_file %||% character(0)
      )
    })

    # Check for IDAT files in supplementary data
    has_idats <- any(sapply(sample_info, function(s) {
      any(grepl("\\.idat", s$supplementary, ignore.case = TRUE))
    }))

    # Detect array platform
    platforms <- unique(sapply(sample_info, function(s) s$platform))
    array_type <- detect_array_platform(platforms, title, summary)

    return(list(
      valid = TRUE,
      has_idats = has_idats,
      is_methylation = is_methylation,
      title = title,
      summary = if(nchar(summary) > 200) paste0(substr(summary, 1, 200), "...") else summary,
      sample_count = sample_count,
      array_type = array_type,
      platforms = platforms,
      message = if(has_idats) "Dataset suitable for GIMP analysis" else "Dataset not suitable for Imprinting analysis"
    ))

  }, error = function(e) {
    if (grepl("not found", e$message, ignore.case = TRUE)) {
      return(list(
        valid = FALSE,
        message = paste("GEO dataset", geo_id, "not found or not accessible")
      ))
    } else {
      return(list(
        valid = FALSE,
        message = paste("Error accessing GEO dataset:", e$message)
      ))
    }
  })
}

#' Detect array platform from GEO metadata
#'
#' @param platforms Vector of platform IDs
#' @param title Dataset title
#' @param summary Dataset summary
#' @return Detected array type
#' @keywords internal
detect_array_platform <- function(platforms, title = "", summary = "") {

  # Platform ID mappings
  platform_map <- list(
    "GPL13534" = "450k",    # HumanMethylation450 BeadChip
    "GPL8490" = "450k",     # HumanMethylation450 BeadChip
    "GPL21145" = "EPIC",    # HumanMethylationEPIC BeadChip
    "GPL23976" = "EPICv2"   # HumanMethylationEPIC v2.0 BeadChip
  )

  # Check platform IDs first
  for (platform in platforms) {
    if (platform %in% names(platform_map)) {
      return(platform_map[[platform]])
    }
  }

  # Fallback to text analysis
  combined_text <- paste(title, summary, collapse = " ")

  if (grepl("epicv2|epic.?v2|epic.?2", combined_text, ignore.case = TRUE)) {
    return("EPICv2")
  } else if (grepl("epic|850k", combined_text, ignore.case = TRUE)) {
    return("EPIC")
  } else if (grepl("450k|450", combined_text, ignore.case = TRUE)) {
    return("450k")
  }

  # Default to EPIC if unsure but looks like methylation
  return("EPIC")
}

#' Download and process GEO methylation dataset
#'
#' @description
#' Downloads IDAT files and phenotypic data from a GEO dataset and processes
#' them for GIMP analysis.
#'
#' @param geo_id Character string of GEO accession (e.g., "GSE12345")
#' @param download_dir Directory to download files (default: temporary directory)
#' @param group_column Column name for sample groups (auto-detect if NULL)
#' @param control_terms Terms that indicate control samples
#' @param case_terms Terms that indicate case/treatment samples
#' @param max_samples Maximum number of samples to process (NULL for all)
#' @param normalize_method Normalization method for minfi
#' @param n_cores Number of CPU cores for parallel processing
#'
#' @return List containing processed beta matrix and sample information
#' @examples
#' \donttest{
#' # Process a GEO dataset
#' geo_data <- process_geo_dataset(
#'   geo_id = "GSE68777",
#'   group_column = "characteristics_ch1.1",
#'   control_terms = c("normal", "control"),
#'   case_terms = c("tumor", "cancer")
#' )
#'
#' # Use with GIMP functions
#' ICRcpg <- make_cpgs(Bmatrix = geo_data$beta_matrix, bedmeth = "v1")
#' }
#' @export
process_geo_dataset <- function(geo_id,
                                download_dir = NULL,
                                group_column = NULL,
                                control_terms = c("control", "normal", "healthy", "ctrl"),
                                case_terms = c("case", "disease", "tumor", "cancer", "patient"),
                                max_samples = NULL,
                                normalize_method = "quantile",
                                n_cores = NULL) {

  # Validate prerequisites
  if (!requireNamespace("GEOquery", quietly = TRUE)) {
    stop("GEOquery package is required. Install with: BiocManager::install('GEOquery')")
  }

  # Validate dataset first
  validation <- validate_geo_dataset(geo_id)
  if (!validation$valid) {
    stop("Dataset validation failed: ", validation$message)
  }

  if (!validation$has_idats) {
    stop("Dataset does not contain IDAT files. GIMP requires raw IDAT data.")
  }

  message("=== PROCESSING GEO DATASET: ", geo_id, " ===")
  message("Title: ", validation$title)
  message("Array type: ", validation$array_type)
  message("Sample count: ", validation$sample_count)

  # Setup directories
  if (is.null(download_dir)) {
    download_dir <- file.path(tempdir(), paste0("geo_", geo_id))
  }

  if (!dir.exists(download_dir)) {
    dir.create(download_dir, recursive = TRUE)
  }

  idat_dir <- file.path(download_dir, "idats")
  if (!dir.exists(idat_dir)) {
    dir.create(idat_dir, recursive = TRUE)
  }

  # Download and extract data
  message("\n=== DOWNLOADING GEO DATA ===")

  tryCatch({
    # Get detailed GSE information
    gse <- GEOquery::getGEO(geo_id, destdir = download_dir, GSEMatrix = TRUE)

    # Extract phenotypic data from the first (usually only) series matrix
    if (length(gse) > 1) {
      message("Multiple series detected, using the first one")
    }

    eset <- gse[[1]]
    pheno_data <- Biobase::pData(eset)

    message("Phenotypic data extracted: ", nrow(pheno_data), " samples, ",
            ncol(pheno_data), " variables")

    # Download supplementary files (IDAT files)
    message("\n=== DOWNLOADING IDAT FILES ===")
    supp_files <- GEOquery::getGEOSuppFiles(geo_id, baseDir = download_dir,
                                           fetch_files = TRUE)

    # Extract any compressed files
    extract_compressed_files(file.path(download_dir, geo_id))

    # Find and organize IDAT files
    idat_files <- find_and_organize_idats(file.path(download_dir, geo_id), idat_dir)

    if (length(idat_files) == 0) {
      stop("No IDAT files found in downloaded supplementary data")
    }

    message("Found ", length(idat_files), " IDAT files")

    # Create sample sheet from phenotypic data
    message("\n=== CREATING SAMPLE SHEET ===")
    sample_sheet <- create_geo_sample_sheet(
      pheno_data = pheno_data,
      idat_files = idat_files,
      group_column = group_column,
      control_terms = control_terms,
      case_terms = case_terms,
      max_samples = max_samples
    )

    message("Sample sheet created with ", nrow(sample_sheet), " samples")

    # Save sample sheet
    sample_sheet_path <- file.path(idat_dir, "samplesheet.csv")
    write.csv(sample_sheet, sample_sheet_path, row.names = FALSE)

    # Create ZIP file for compatibility with existing GIMP functions
    message("\n=== CREATING ZIP FOR PROCESSING ===")
    zip_path <- create_geo_zip(idat_dir, download_dir, geo_id)

    # Process with existing GIMP IDAT pipeline
    message("\n=== PROCESSING WITH GIMP PIPELINE ===")

    # Map GEO array type to GIMP bedmeth parameter
    bedmeth_mapping <- list(
      "450k" = "450k",
      "EPIC" = "v1",
      "EPICv2" = "v2"
    )

    idat_results <- read_idat_zip(
      zip_file = zip_path,
      sample_sheet_name = "samplesheet.csv",
      array_type = validation$array_type,
      normalize_method = normalize_method,
      detection_pval = 0.01,
      remove_failed_samples = TRUE,
      n_cores = n_cores
    )

    # Add GEO-specific metadata
    idat_results$geo_metadata <- list(
      geo_id = geo_id,
      title = validation$title,
      array_type = validation$array_type,
      original_sample_count = validation$sample_count,
      processed_sample_count = ncol(idat_results$beta_matrix),
      download_dir = download_dir,
      group_detection = sample_sheet$group_detection_info[1]
    )

    message("\n=== GEO PROCESSING COMPLETED ===")
    message("✅ Successfully processed ", ncol(idat_results$beta_matrix), " samples")
    message("✅ Beta matrix dimensions: ", paste(dim(idat_results$beta_matrix), collapse = " x "))

    return(idat_results)

  }, error = function(e) {
    # Cleanup on error
    if (dir.exists(download_dir)) {
      unlink(download_dir, recursive = TRUE)
    }
    stop("GEO processing failed: ", e$message, call. = FALSE)
  })
}

#' Extract compressed files in a directory
#'
#' @param dir_path Directory to search for compressed files
#' @keywords internal
extract_compressed_files <- function(dir_path) {

  message("Searching for compressed files in: ", dir_path)

  # Find compressed files (including nested directories)
  tar_files <- list.files(dir_path, pattern = "\\.(tar|tar\\.gz|tgz)$",
                         full.names = TRUE, recursive = TRUE)
  zip_files <- list.files(dir_path, pattern = "\\.zip$",
                         full.names = TRUE, recursive = TRUE)
  gz_files <- list.files(dir_path, pattern = "\\.gz$",
                        full.names = TRUE, recursive = TRUE)

  message("Found compressed files: ", length(tar_files), " TAR, ",
          length(zip_files), " ZIP, ", length(gz_files), " GZ")

  # Extract tar files
  for (tar_file in tar_files) {
    message("Extracting TAR: ", basename(tar_file))
    tryCatch({
      utils::untar(tar_file, exdir = dirname(tar_file), verbose = TRUE)
      message("✅ Successfully extracted: ", basename(tar_file))
    }, error = function(e) {
      message("❌ Failed to extract ", basename(tar_file), ": ", e$message)
    })
  }

  # Extract zip files
  for (zip_file in zip_files) {
    message("Extracting ZIP: ", basename(zip_file))
    tryCatch({
      utils::unzip(zip_file, exdir = dirname(zip_file))
      message("✅ Successfully extracted: ", basename(zip_file))
    }, error = function(e) {
      message("❌ Failed to extract ", basename(zip_file), ": ", e$message)
    })
  }

  # Handle individual .gz files (might be compressed IDAT files)
  for (gz_file in gz_files) {
    # Skip if it's part of a .tar.gz file (already handled above)
    if (grepl("\\.tar\\.gz$", gz_file)) next

    message("Extracting GZ: ", basename(gz_file))
    tryCatch({
      # Use R.utils if available, otherwise system gunzip
      if (requireNamespace("R.utils", quietly = TRUE)) {
        R.utils::gunzip(gz_file, destname = gsub("\\.gz$", "", gz_file),
                       overwrite = TRUE, remove = FALSE)
      } else {
        # Fallback to system gunzip
        system(paste("gunzip -c", shQuote(gz_file), ">",
                    shQuote(gsub("\\.gz$", "", gz_file))))
      }
      message("✅ Successfully extracted: ", basename(gz_file))
    }, error = function(e) {
      message("❌ Failed to extract ", basename(gz_file), ": ", e$message)
    })
  }

  # After extraction, check what files we have
  all_files <- list.files(dir_path, recursive = TRUE, full.names = TRUE)
  idat_files <- grep("\\.idat$", all_files, value = TRUE, ignore.case = TRUE)
  message("After extraction, found ", length(idat_files), " IDAT files")

  if (length(idat_files) > 0) {
    message("Example IDAT files found:")
    for (i in seq_len(min(5, length(idat_files)))) {
      message("  ", basename(idat_files[i]))
    }
  }
}

#' Find and organize IDAT files
#'
#' @param source_dir Directory to search for IDAT files
#' @param target_dir Directory to copy organized IDAT files
#' @return Vector of IDAT file paths
#' @keywords internal
find_and_organize_idats <- function(source_dir, target_dir) {

  message("=== SEARCHING FOR IDAT FILES ===")
  message("Source directory: ", source_dir)
  message("Target directory: ", target_dir)

  # First, let's see what's in the source directory
  all_files <- list.files(source_dir, recursive = TRUE, full.names = TRUE)
  message("Total files found: ", length(all_files))

  # Find all IDAT files recursively (case insensitive)
  idat_files <- grep("\\.idat$", all_files, value = TRUE, ignore.case = TRUE)

  message("IDAT files found: ", length(idat_files))

  if (length(idat_files) == 0) {
    # Maybe they're compressed? Let's check for compressed IDAT files
    compressed_idats <- grep("\\.(idat\\.gz|idat\\.bz2)$", all_files, value = TRUE, ignore.case = TRUE)

    if (length(compressed_idats) > 0) {
      message("Found ", length(compressed_idats), " compressed IDAT files")

      # Decompress them
      for (compressed_file in compressed_idats) {
        message("Decompressing: ", basename(compressed_file))

        tryCatch({
          if (grepl("\\.gz$", compressed_file)) {
            # Handle .gz files
            if (requireNamespace("R.utils", quietly = TRUE)) {
              decompressed_file <- gsub("\\.gz$", "", compressed_file)
              R.utils::gunzip(compressed_file, destname = decompressed_file,
                             overwrite = TRUE, remove = FALSE)
              idat_files <- c(idat_files, decompressed_file)
            } else {
              warning("R.utils package not available for decompressing .gz files")
            }
          } else if (grepl("\\.bz2$", compressed_file)) {
            # Handle .bz2 files
            if (requireNamespace("R.utils", quietly = TRUE)) {
              decompressed_file <- gsub("\\.bz2$", "", compressed_file)
              R.utils::bunzip2(compressed_file, destname = decompressed_file,
                              overwrite = TRUE, remove = FALSE)
              idat_files <- c(idat_files, decompressed_file)
            } else {
              warning("R.utils package not available for decompressing .bz2 files")
            }
          }
        }, error = function(e) {
          message("Failed to decompress ", basename(compressed_file), ": ", e$message)
        })
      }
    }

    # Re-search for IDAT files after decompression
    all_files <- list.files(source_dir, recursive = TRUE, full.names = TRUE)
    idat_files <- grep("\\.idat$", all_files, value = TRUE, ignore.case = TRUE)
    message("After decompression, IDAT files found: ", length(idat_files))
  }

  if (length(idat_files) == 0) {
    message("❌ No IDAT files found")

    # Debug information
    message("\n=== DEBUG INFORMATION ===")
    message("Directory structure:")

    # Show directory structure (first 20 files)
    for (i in seq_len(min(20, length(all_files)))) {
      rel_path <- gsub(paste0("^", source_dir, "/?"), "", all_files[i])
      message("  ", rel_path)
    }

    if (length(all_files) > 20) {
      message("  ... and ", length(all_files) - 20, " more files")
    }

    # Check for files that might be IDAT-related
    possible_files <- grep("(idat|methylation|450|epic)", all_files, value = TRUE, ignore.case = TRUE)
    if (length(possible_files) > 0) {
      message("\nPossible methylation-related files:")
      for (file in head(possible_files, 10)) {
        rel_path <- gsub(paste0("^", source_dir, "/?"), "", file)
        message("  ", rel_path)
      }
    }

    return(character(0))
  }

  # Show examples of found IDAT files
  message("Found IDAT files (first 10):")
  for (i in seq_len(min(10, length(idat_files)))) {
    message("  ", basename(idat_files[i]))
  }

  # Verify we have both Red and Grn files
  red_files <- grep("_Red\\.idat$", idat_files, value = TRUE, ignore.case = TRUE)
  grn_files <- grep("_Grn\\.idat$", idat_files, value = TRUE, ignore.case = TRUE)

  message("Red channel files: ", length(red_files))
  message("Green channel files: ", length(grn_files))

  if (length(red_files) == 0 || length(grn_files) == 0) {
    warning("Missing Red or Green channel files. Found ", length(red_files),
            " Red and ", length(grn_files), " Green files.")
  }

  # Copy IDAT files to organized directory
  message("Copying IDAT files to target directory...")
  copied_files <- character(0)

  for (idat_file in idat_files) {
    target_file <- file.path(target_dir, basename(idat_file))

    tryCatch({
      file.copy(idat_file, target_file, overwrite = TRUE)
      copied_files <- c(copied_files, target_file)
    }, error = function(e) {
      message("Failed to copy ", basename(idat_file), ": ", e$message)
    })
  }

  message("Successfully copied ", length(copied_files), " IDAT files")

  # Return list of organized files
  return(copied_files)
}

#' Create sample sheet from GEO phenotypic data
#'
#' @param pheno_data Phenotypic data from GEO
#' @param idat_files Vector of IDAT file paths
#' @param group_column Column name for groups
#' @param control_terms Terms indicating controls
#' @param case_terms Terms indicating cases
#' @param max_samples Maximum samples to include
#' @return Data frame with sample sheet
#' @keywords internal
create_geo_sample_sheet <- function(pheno_data, idat_files, group_column,
                                   control_terms, case_terms, max_samples) {

  # Extract sample identifiers from IDAT files
  idat_basenames <- unique(gsub("_(Red|Grn)\\.idat$", "", basename(idat_files)))

  message("Creating sample sheet...")
  message("Phenotypic data samples: ", nrow(pheno_data))
  message("IDAT basenames found: ", length(idat_basenames))

  # Match samples between phenotypic data and IDAT files
  # Try multiple matching strategies
  sample_matches <- match_geo_samples(pheno_data, idat_basenames)

  if (nrow(sample_matches) == 0) {
    stop("Could not match any samples between phenotypic data and IDAT files")
  }

  message("Matched ", nrow(sample_matches), " samples between pheno data and IDAT files")

  # Limit samples if requested
  if (!is.null(max_samples) && nrow(sample_matches) > max_samples) {
    sample_matches <- sample_matches[seq_len(max_samples), ]
    message("Limited to first ", max_samples, " samples")
  }

  # Detect or use specified group column
  if (is.null(group_column)) {
    group_info <- auto_detect_groups(pheno_data, sample_matches$pheno_id,
                                   control_terms, case_terms)
    group_column <- group_info$column
    sample_groups <- group_info$groups
    group_detection_msg <- group_info$message
  } else {
    if (!group_column %in% colnames(pheno_data)) {
      stop("Specified group column '", group_column, "' not found in phenotypic data")
    }
    sample_groups <- assign_sample_groups(pheno_data[sample_matches$pheno_id, group_column],
                                        control_terms, case_terms)
    group_detection_msg <- paste("Used specified column:", group_column)
  }

  # Extract sentrix info with better error handling
  sentrix_info <- extract_sentrix_info(sample_matches$idat_id)

  # Validate sentrix_info - FIXED: Handle potential NA values
  if (any(is.na(sentrix_info$sentrix_id)) || any(is.na(sentrix_info$sentrix_position))) {
    warning("Some Sentrix IDs or positions could not be extracted. Using fallback values.")

    # Fill in missing values
    na_indices <- which(is.na(sentrix_info$sentrix_id))
    for (i in na_indices) {
      sentrix_info$sentrix_id[i] <- paste0("200000000", sprintf("%02d", i))
    }

    na_indices <- which(is.na(sentrix_info$sentrix_position))
    for (i in na_indices) {
      sentrix_info$sentrix_position[i] <- paste0("R", sprintf("%02d", ((i-1) %% 8) + 1), "C01")
    }
  }

  # Create sample sheet with proper error handling
  sample_sheet <- data.frame(
    Sample_Name = sample_matches$pheno_id,
    Sentrix_ID = sentrix_info$sentrix_id,
    Sentrix_Position = sentrix_info$sentrix_position,
    Sample_Group = sample_groups,
    GEO_Sample = sample_matches$pheno_id,
    IDAT_Basename = sample_matches$idat_id,
    group_detection_info = group_detection_msg,
    stringsAsFactors = FALSE
  )

  # Validate sample sheet more carefully - FIXED: Better NA handling
  valid_rows <- !is.na(sample_sheet$Sample_Name) &
                sample_sheet$Sample_Name != "" &
                !is.na(sample_sheet$Sentrix_ID) &
                sample_sheet$Sentrix_ID != "" &
                !is.na(sample_sheet$Sentrix_Position) &
                sample_sheet$Sentrix_Position != ""

  # Check if any rows are invalid
  if (any(!valid_rows)) {
    message("Warning: ", sum(!valid_rows), " samples have invalid data and will be removed")
    sample_sheet <- sample_sheet[valid_rows, ]
  }

  if (nrow(sample_sheet) == 0) {
    stop("No valid samples remain after filtering. Check sample matching and data quality.")
  }

  # Print group assignment summary
  group_table <- table(sample_sheet$Sample_Group)
  message("Group assignment: ", paste(names(group_table), "=", group_table, collapse = ", "))
  message("Group detection: ", group_detection_msg)

  return(sample_sheet)
}

#' Match samples between phenotypic data and IDAT files
#'
#' @param pheno_data Phenotypic data
#' @param idat_basenames IDAT file basenames
#' @return Data frame with matched samples
#' @keywords internal
match_geo_samples <- function(pheno_data, idat_basenames) {

  matches <- data.frame(
    pheno_id = character(0),
    idat_id = character(0),
    stringsAsFactors = FALSE
  )

  # Strategy 1: Direct GEO sample ID matching
  geo_samples <- rownames(pheno_data)
  for (geo_sample in geo_samples) {
    # Try exact match
    exact_match <- grep(paste0("^", geo_sample, "_"), idat_basenames, value = TRUE)
    if (length(exact_match) > 0) {
      matches <- rbind(matches, data.frame(
        pheno_id = geo_sample,
        idat_id = exact_match[1],
        stringsAsFactors = FALSE
      ))
      next
    }

    # Try partial match
    partial_match <- grep(geo_sample, idat_basenames, value = TRUE)
    if (length(partial_match) > 0) {
      matches <- rbind(matches, data.frame(
        pheno_id = geo_sample,
        idat_id = partial_match[1],
        stringsAsFactors = FALSE
      ))
    }
  }

  # Strategy 2: If few matches, try reverse matching
  if (nrow(matches) < length(geo_samples) * 0.5) {
    for (idat_base in idat_basenames) {
      # Extract potential sample ID from IDAT filename
      potential_id <- extract_sample_id_from_idat(idat_base)

      if (potential_id %in% geo_samples && !potential_id %in% matches$pheno_id) {
        matches <- rbind(matches, data.frame(
          pheno_id = potential_id,
          idat_id = idat_base,
          stringsAsFactors = FALSE
        ))
      }
    }
  }

  # Remove duplicates
  matches <- matches[!duplicated(matches$pheno_id), ]
  matches <- matches[!duplicated(matches$idat_id), ]

  return(matches)
}

#' Auto-detect sample groups from phenotypic data
#'
#' @param pheno_data Phenotypic data
#' @param sample_ids Sample IDs to analyze
#' @param control_terms Terms indicating controls
#' @param case_terms Terms indicating cases
#' @return List with group information
#' @keywords internal
auto_detect_groups <- function(pheno_data, sample_ids, control_terms, case_terms) {

  message("Auto-detecting sample groups...")

  # Look for promising columns
  potential_columns <- grep("group|condition|treatment|disease|status|type|class",
                           colnames(pheno_data), ignore.case = TRUE, value = TRUE)

  if (length(potential_columns) == 0) {
    # Fallback to characteristics columns
    potential_columns <- grep("characteristics", colnames(pheno_data),
                            ignore.case = TRUE, value = TRUE)
  }

  if (length(potential_columns) == 0) {
    # Use all character/factor columns as potential group columns
    potential_columns <- colnames(pheno_data)[sapply(pheno_data, function(x) {
      is.character(x) || is.factor(x)
    })]
  }

  message("Found ", length(potential_columns), " potential group columns")

  # Score each column for group potential
  best_column <- NULL
  best_score <- 0
  best_groups <- NULL

  for (col in potential_columns) {
    if (col %in% colnames(pheno_data)) {

      # FIXED: Better handling of sample subsetting
      tryCatch({
        # Get column data for the matched samples
        if (all(sample_ids %in% rownames(pheno_data))) {
          col_data <- pheno_data[sample_ids, col]
        } else {
          # Fallback: use all data if rowname matching fails
          col_data <- pheno_data[[col]]
          if (length(col_data) != length(sample_ids)) {
            col_data <- col_data[seq_len(min(length(col_data), length(sample_ids)))]
          }
        }

        # Convert to character and handle NAs
        col_data <- as.character(col_data)
        col_data[is.na(col_data)] <- "Unknown"

        groups <- assign_sample_groups(col_data, control_terms, case_terms)

        # Score based on balance and clarity
        group_table <- table(groups)
        if (length(group_table) == 2 && !all(names(group_table) == "Unknown")) {
          balance_score <- min(group_table) / max(group_table)
          clarity_score <- sum(groups != "Unknown") / length(groups)
          total_score <- balance_score * clarity_score

          message("Column '", col, "' score: ", round(total_score, 3),
                 " (", paste(names(group_table), "=", group_table, collapse = ", "), ")")

          if (total_score > best_score) {
            best_score <- total_score
            best_column <- col
            best_groups <- groups
          }
        }

      }, error = function(e) {
        message("Error processing column '", col, "': ", e$message)
      })
    }
  }

  if (is.null(best_column)) {
    # Fallback: assign all as Unknown
    best_groups <- rep("Unknown", length(sample_ids))
    message_text <- "Could not auto-detect sample groups. Manual assignment required."
  } else {
    message_text <- paste("Auto-detected groups from column:", best_column,
                         "(", paste(names(table(best_groups)), collapse = " vs "), ")")
  }

  message("Group detection result: ", message_text)

  return(list(
    column = best_column,
    groups = best_groups,
    message = message_text
  ))
}

#' Assign sample groups based on terms
#'
#' @param values Vector of values to classify
#' @param control_terms Terms indicating controls
#' @param case_terms Terms indicating cases
#' @return Vector of group assignments
#' @keywords internal
assign_sample_groups <- function(values, control_terms, case_terms) {

  values <- as.character(values)
  groups <- rep("Unknown", length(values))

  # Check for control terms
  for (term in control_terms) {
    control_matches <- grepl(term, values, ignore.case = TRUE)
    groups[control_matches] <- "Control"
  }

  # Check for case terms
  for (term in case_terms) {
    case_matches <- grepl(term, values, ignore.case = TRUE)
    groups[case_matches] <- "Case"
  }

  return(groups)
}

#' Extract Sentrix information from IDAT basename
#'
#' @param idat_basenames Vector of IDAT basenames
#' @return List with Sentrix ID and position
#' @keywords internal
extract_sentrix_info <- function(idat_basenames) {

  # Try to extract standard Sentrix format: SentrixID_SentrixPosition
  sentrix_pattern <- "^(\\d+)_([A-Z]\\d+[A-Z]\\d+)$"

  sentrix_ids <- character(length(idat_basenames))
  sentrix_positions <- character(length(idat_basenames))

  for (i in seq_along(idat_basenames)) {
    basename <- idat_basenames[i]

    # Try standard pattern
    if (grepl(sentrix_pattern, basename)) {
      matches <- regmatches(basename, regexec(sentrix_pattern, basename))[[1]]
      sentrix_ids[i] <- matches[2]
      sentrix_positions[i] <- matches[3]
    } else {
      # Fallback: create synthetic IDs
      sentrix_ids[i] <- paste0("200000000", sprintf("%02d", i))
      sentrix_positions[i] <- paste0("R", sprintf("%02d", ((i-1) %% 8) + 1), "C01")
    }
  }

  return(list(
    sentrix_id = sentrix_ids,
    sentrix_position = sentrix_positions
  ))
}

#' Extract sample ID from IDAT filename
#'
#' @param idat_basename IDAT file basename
#' @return Potential sample ID
#' @keywords internal
extract_sample_id_from_idat <- function(idat_basename) {

  # Remove common prefixes and suffixes
  clean_name <- gsub("^.*?(GSM\\d+).*$", "\\1", idat_basename)

  if (clean_name != idat_basename) {
    return(clean_name)
  }

  # Fallback: return first part before underscore
  parts <- strsplit(idat_basename, "_")[[1]]
  return(parts[1])
}

#' Create ZIP file for GIMP processing
#'
#' @param idat_dir Directory containing organized IDAT files
#' @param base_dir Base directory for ZIP creation
#' @param geo_id GEO ID for naming
#' @return Path to created ZIP file
#' @keywords internal
create_geo_zip <- function(idat_dir, base_dir, geo_id) {

  zip_path <- file.path(base_dir, paste0(geo_id, "_processed.zip"))

  # Get all files to include
  files_to_zip <- list.files(idat_dir, full.names = TRUE)

  # Create ZIP file
  oldwd <- getwd()
  setwd(idat_dir)

  tryCatch({
    utils::zip(zip_path, files = basename(files_to_zip), flags = "-r9Xq")
    setwd(oldwd)
  }, error = function(e) {
    setwd(oldwd)
    stop("Failed to create ZIP file: ", e$message)
  })

  if (!file.exists(zip_path)) {
    stop("ZIP file creation failed")
  }

  message("Created ZIP file: ", zip_path, " (", round(file.size(zip_path)/1024^2, 1), " MB)")

  return(zip_path)
}

#' Diagnose GEO dataset structure
#'
#' @description
#' Downloads and examines a GEO dataset to understand its file structure
#' without processing the data. Useful for troubleshooting.
#'
#' @param geo_id Character string of GEO accession (e.g., "GSE12345")
#' @param temp_dir Temporary directory for download (default: temporary directory)
#' @return List with diagnostic information
#' @examples
#' \donttest{
#' # Diagnose a problematic dataset
#' diag <- diagnose_geo_dataset("GSE289527")
#' print(diag$summary)
#' }
#' @export
diagnose_geo_dataset <- function(geo_id, temp_dir = NULL) {
  # Null coalescing operator if not already defined
  `%||%` <- function(a, b) if (!is.null(a)) a else b

  if (!requireNamespace("GEOquery", quietly = TRUE)) {
    stop("GEOquery package is required. Install with: BiocManager::install('GEOquery')")
  }

  geo_id <- toupper(trimws(geo_id))

  if (is.null(temp_dir)) {
    temp_dir <- file.path(tempdir(), paste0("geo_diag_", geo_id))
  }

  if (!dir.exists(temp_dir)) {
    dir.create(temp_dir, recursive = TRUE)
  }

  message("=== DIAGNOSING GEO DATASET: ", geo_id, " ===")

  tryCatch({
    # Download GSE information
    message("Downloading GSE metadata...")
    gse <- GEOquery::getGEO(geo_id, destdir = temp_dir, GSEMatrix = TRUE)

    # Get basic information
    eset <- gse[[1]]
    pheno_data <- Biobase::pData(eset)

    message("Downloading supplementary files...")
    supp_files <- GEOquery::getGEOSuppFiles(geo_id, baseDir = temp_dir,
                                            fetch_files = TRUE)

    # Analyze downloaded files
    geo_data_dir <- file.path(temp_dir, geo_id)

    message("Extracting compressed files...")
    extract_compressed_files(geo_data_dir)

    # Get comprehensive file listing
    all_files <- list.files(geo_data_dir, recursive = TRUE, full.names = TRUE)
    file_info <- file.info(all_files)
    file_info$filename <- basename(all_files)
    file_info$extension <- tools::file_ext(all_files)
    file_info$relative_path <- gsub(paste0("^", geo_data_dir, "/?"), "", all_files)

    # Categorize files
    idat_files <- grep("\\.idat$", all_files, value = TRUE, ignore.case = TRUE)
    compressed_files <- grep("\\.(gz|bz2|zip|tar)$", all_files, value = TRUE, ignore.case = TRUE)
    csv_files <- grep("\\.csv$", all_files, value = TRUE, ignore.case = TRUE)
    txt_files <- grep("\\.txt$", all_files, value = TRUE, ignore.case = TRUE)

    # Create summary
    summary_info <- list(
      geo_id = geo_id,
      title = GEOquery::Meta(gse)$title %||% "Unknown",
      samples = nrow(pheno_data),
      phenotype_columns = ncol(pheno_data),
      total_files = length(all_files),
      idat_files = length(idat_files),
      compressed_files = length(compressed_files),
      csv_files = length(csv_files),
      txt_files = length(txt_files),
      largest_files = head(file_info[order(-file_info$size), c("filename", "size")], 10),
      file_extensions = sort(table(file_info$extension), decreasing = TRUE)
    )

    # Print diagnostic report
    cat("\n=== DIAGNOSTIC REPORT ===\n")
    cat("GEO ID:", summary_info$geo_id, "\n")
    cat("Title:", summary_info$title, "\n")
    cat("Samples:", summary_info$samples, "\n")
    cat("Total files downloaded:", summary_info$total_files, "\n")
    cat("IDAT files found:", summary_info$idat_files, "\n")

    cat("\nFile types:\n")
    for (ext in names(summary_info$file_extensions)) {
      ext_label <- if (ext == "") "(no extension)" else paste0(".", ext)
      cat("  ", ext_label, ":", summary_info$file_extensions[ext], "\n")
    }

    cat("\nLargest files:\n")
    for (i in seq_len(min(5, nrow(summary_info$largest_files)))) {
      cat("  ", summary_info$largest_files$filename[i],
          " (", round(summary_info$largest_files$size[i] / 1024^2, 1), " MB)\n")
    }

    if (length(idat_files) > 0) {
      cat("\nIDAT files found:\n")
      red_files <- grep("_Red\\.idat$", idat_files, ignore.case = TRUE)
      grn_files <- grep("_Grn\\.idat$", idat_files, ignore.case = TRUE)

      cat("  Red channel:", length(red_files), "\n")
      cat("  Green channel:", length(grn_files), "\n")

      if (length(red_files) > 0) {
        cat("  Example Red files:\n")
        for (i in seq_len(min(3, length(red_files)))) {
          cat("    ", basename(idat_files[red_files[i]]), "\n")
        }
      }
    } else {
      cat("\n❌ NO IDAT FILES FOUND\n")

      # Suggest possible issues
      cat("\nPossible issues:\n")
      cat("1. Dataset contains only processed data (no raw IDAT files)\n")
      cat("2. IDAT files are in unexpected compressed format\n")
      cat("3. Files are stored in subdirectories we didn't check\n")
      cat("4. Dataset is restricted access\n")

      # Look for methylation-related files
      meth_files <- grep("(methylation|450|epic|raw)", all_files, value = TRUE, ignore.case = TRUE)
      if (length(meth_files) > 0) {
        cat("\nMethylation-related files found:\n")
        for (file in head(meth_files, 5)) {
          rel_path <- gsub(paste0("^", geo_data_dir, "/?"), "", file)
          cat("  ", rel_path, "\n")
        }
      }
    }

    # Check phenotypic data
    cat("\nPhenotypic data columns (first 10):\n")
    pheno_cols <- colnames(pheno_data)
    for (col in head(pheno_cols, 10)) {
      unique_vals <- length(unique(pheno_data[[col]]))
      cat("  ", col, " (", unique_vals, " unique values)\n")
    }

    # Cleanup
    unlink(temp_dir, recursive = TRUE)

    return(list(
      summary = summary_info,
      file_info = file_info,
      phenotype_data = pheno_data,
      idat_files = basename(idat_files),
      diagnosis = if (length(idat_files) > 0) "SUITABLE" else "NO_IDATS"
    ))

  }, error = function(e) {
    message("Error diagnosing GEO dataset: ", e$message)
    unlink(temp_dir, recursive = TRUE)
    return(list(
      summary = NULL,
      file_info = NULL,
      phenotype_data = NULL,
      idat_files = NULL,
      diagnosis = "ERROR",
      error_message = e$message
    ))
  })
}

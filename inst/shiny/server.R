library(shiny)
library(shinydashboard)
library(DT)
library(GIMP)
library(tidyverse)
library(plotly)
library(ggplot2)
library(reshape2)

server <- function(input, output, session) {

  options(shiny.maxRequestSize = 500*1024^2)  # 500MB limit
  # Reactive values to store data
  values <- reactiveValues(
    betaMatrix = NULL,
    sampleInfo = NULL,
    ICRcpg = NULL,
    df.ICR = NULL,
    cpgs_analysis = NULL,
    dmps_results = NULL,
    processed = FALSE,
    qc_metrics = NULL,
    sample_sheet = NULL,
    zip_preview = NULL,
    geo_validation = NULL,
    geo_metadata = NULL,
    # New GEO multi-step values
    geo_pheno_data = NULL,
    geo_step2_ready = FALSE,
    geo_step3_ready = FALSE,
    geo_step4_ready = FALSE,
    selected_group_column = NULL,
    group_mappings = NULL,
    unique_group_values = NULL
  )

  # Check if data is loaded
  output$dataLoaded <- reactive({
    !is.null(values$betaMatrix)
  })
  outputOptions(output, "dataLoaded", suspendWhenHidden = FALSE)

  # Check if DMPs analysis is done
  output$dmpsDone <- reactive({
    !is.null(values$dmps_results)
  })
  outputOptions(output, "dmpsDone", suspendWhenHidden = FALSE)

  # Check if ZIP preview is available
  output$zipPreviewAvailable <- reactive({
    !is.null(input$idatZip) && !is.null(values$zip_preview)
  })
  outputOptions(output, "zipPreviewAvailable", suspendWhenHidden = FALSE)

  # Check if GEO is validated
  output$geoValidated <- reactive({
    !is.null(values$geo_validation)
  })
  outputOptions(output, "geoValidated", suspendWhenHidden = FALSE)

  # Check if GEO has IDATs
  output$geoHasIdats <- reactive({
    !is.null(values$geo_validation) && values$geo_validation$valid && values$geo_validation$has_idats
  })
  outputOptions(output, "geoHasIdats", suspendWhenHidden = FALSE)

  # New GEO step reactive outputs
  output$geoStep2Ready <- reactive({
    values$geo_step2_ready
  })
  outputOptions(output, "geoStep2Ready", suspendWhenHidden = FALSE)

  output$geoStep3Ready <- reactive({
    values$geo_step3_ready
  })
  outputOptions(output, "geoStep3Ready", suspendWhenHidden = FALSE)

  output$geoStep4Ready <- reactive({
    values$geo_step4_ready
  })
  outputOptions(output, "geoStep4Ready", suspendWhenHidden = FALSE)

  # ============================================================================
  # DATA UPLOAD - PROCESSED FILES
  # ============================================================================

  # File upload and preview for processed data
  observeEvent(input$betaMatrix, {
    req(input$betaMatrix)

    file_ext <- tools::file_ext(input$betaMatrix$datapath)

    tryCatch({
      df <- switch(file_ext,
        "csv" = read.csv(input$betaMatrix$datapath, row.names = 1, check.names = FALSE),
        "txt" = read.table(input$betaMatrix$datapath, header = TRUE, row.names = 1, sep = "\t"),
        "rds" = readRDS(input$betaMatrix$datapath),
        "xlsx" = {
          library(readxl)
          read_excel(input$betaMatrix$datapath, sheet = 1) %>%
            column_to_rownames(var = names(.)[1])
        },
        stop("Unsupported file format")
      )

      # Validate data
      if (!all(df >= 0 & df <= 1, na.rm = TRUE)) {
        showNotification("Warning: Values outside 0-1 range detected. Please check your data.",
                         type = "warning", duration = 10)
      }

      values$betaMatrix <- as.matrix(df)
      values$qc_metrics <- NULL  # Clear IDAT QC metrics
      values$geo_metadata <- NULL  # Clear GEO metadata

      # Generate sample group UI for processed data
      output$sampleGroupUI <- renderUI({
        samples <- colnames(values$betaMatrix)

        tagList(
          selectInput("controlSamples", "Select Control Samples:",
                      choices = samples,
                      multiple = TRUE,
                      selected = NULL
          ),
          helpText("Samples not selected as controls will be assigned to the Case group.")
        )
      })

      updateDataSummary()

    }, error = function(e) {
      showNotification(paste("Error loading file:", e$message), type = "error", duration = 10)
    })
  })

  # Process data for processed files
  observeEvent(input$processData, {
    req(values$betaMatrix, input$controlSamples)

    showModal(modalDialog(
      title = "Processing Data",
      "Extracting CpGs and creating ICR matrices...",
      footer = NULL
    ))

    # Create sample info for processed data
    samples <- colnames(values$betaMatrix)
    values$sampleInfo <- ifelse(samples %in% input$controlSamples, "Control", "Case")

    # Process data with GIMP functions
    tryCatch({
      values$ICRcpg <- make_cpgs(Bmatrix = values$betaMatrix, bedmeth = input$arrayType)
      values$df.ICR <- make_ICRs(Bmatrix = values$betaMatrix, bedmeth = input$arrayType)
      values$processed <- TRUE

      removeModal()
      showNotification("Data processed successfully!", type = "message")

    }, error = function(e) {
      removeModal()
      showNotification(paste("Processing failed:", e$message), type = "error", duration = 10)
    })
  })

  # ============================================================================
  # DATA UPLOAD - IDAT FILES
  # ============================================================================

  # ZIP preview
  observeEvent(input$previewZip, {
    req(input$idatZip)

    tryCatch({
      preview_results <- preview_idat_zip(input$idatZip$datapath)
      values$zip_preview <- preview_results

      output$zipPreview <- renderPrint({
        if (preview_results$valid) {
          cat("ZIP File Preview\n")
          cat("===============\n\n")
          cat("File Structure:\n")
          cat("  Total files:", preview_results$total_files, "\n")
          cat("  IDAT files:", preview_results$idat_files, "\n")
          cat("  - Red channel:", preview_results$red_files, "\n")
          cat("  - Green channel:", preview_results$grn_files, "\n")
          cat("  CSV files:", length(preview_results$csv_files), "\n")
          cat("  Other files:", length(preview_results$other_files), "\n")
          cat("  Total size:", preview_results$total_size_mb, "MB\n\n")

          cat("Sample Information:\n")
          cat("  Estimated samples:", preview_results$estimated_samples, "\n")

          if (length(preview_results$sample_identifiers) > 0) {
            cat("  Sample IDs (first 10):\n")
            for (i in seq_len(min(10, length(preview_results$sample_identifiers)))) {
              cat("   ", i, ":", preview_results$sample_identifiers[i], "\n")
            }
          }

          if (length(preview_results$csv_files) > 0) {
            cat("\n  CSV files found:\n")
            for (csv in preview_results$csv_files) {
              cat("   -", csv, "\n")
            }
          }

        } else {
          cat("ZIP Validation Failed\n")
          cat("====================\n\n")
          cat("Error:", preview_results$message, "\n")
        }
      })

      showNotification("ZIP file preview generated", type = "message")

    }, error = function(e) {
      showNotification(paste("Preview failed:", e$message), type = "error", duration = 10)
    })
  })

  # IDAT file processing
  observeEvent(input$processIDAT, {
    req(input$idatZip)

    # Pre-validation
    if (is.null(values$zip_preview)) {
      showNotification("Please preview the ZIP file first to validate its contents.",
                       type = "warning", duration = 8)
      return()
    }

    if (!values$zip_preview$valid) {
      showNotification("ZIP file validation failed. Please check the file contents.",
                       type = "error", duration = 10)
      return()
    }

    # Estimate processing time
    estimated_time <- ceiling(values$zip_preview$estimated_samples * 0.5)

    showModal(modalDialog(
      title = "Processing IDAT Files",
      div(
        h4("Processing in progress..."),
        p(paste("Estimated samples:", values$zip_preview$estimated_samples)),
        p(paste("Estimated time:", estimated_time, "minutes")),
        br(),
        div(class = "progress progress-striped active",
            div(class = "progress-bar progress-bar-info",
                style = "width: 100%", "Processing...")
        ),
        br(),
        p("Steps:"),
        tags$ol(
          tags$li("Extracting ZIP file"),
          tags$li("Reading sample sheet"),
          tags$li("Loading IDAT files"),
          tags$li("Quality control"),
          tags$li("Normalization"),
          tags$li("Generating beta matrix")
        )
      ),
      footer = NULL,
      size = "m"
    ))

    tryCatch({
      # Map array type for GIMP bedmeth parameter
      bedmeth_mapping <- list(
        "450k" = "450k",
        "EPIC" = "v1",
        "EPICv2" = "v2"
      )

      # Process IDAT files
      idat_results <- read_idat_zip(
        zip_file = input$idatZip$datapath,
        sample_sheet_name = input$sampleSheetName,
        array_type = input$idatArrayType,
        normalize_method = input$normMethod,
        detection_pval = input$detectionPval,
        remove_failed_samples = input$removeFailedSamples,
        n_cores = if (input$enableParallel) input$nCores else NULL
      )

      # Store results
      values$betaMatrix <- idat_results$beta_matrix
      values$sample_sheet <- idat_results$sample_sheet
      values$qc_metrics <- idat_results$qc_metrics
      values$geo_metadata <- NULL  # Clear GEO metadata for IDAT uploads

      # Auto-assign sample groups based on sample sheet if available
      if ("Sample_Group" %in% colnames(idat_results$sample_sheet)) {
        values$sampleInfo <- idat_results$sample_sheet$Sample_Group
        names(values$sampleInfo) <- idat_results$sample_sheet$Sample_Name

        # Auto-process if groups are available
        bedmeth_type <- bedmeth_mapping[[input$idatArrayType]]
        values$ICRcpg <- make_cpgs(Bmatrix = values$betaMatrix, bedmeth = bedmeth_type)
        values$df.ICR <- make_ICRs(Bmatrix = values$betaMatrix, bedmeth = bedmeth_type)
        values$processed <- TRUE

      } else {
        # Create UI for manual group assignment
        output$sampleGroupUI <- renderUI({
          samples <- colnames(values$betaMatrix)

          tagList(
            div(class = "alert alert-warning",
                icon("exclamation-triangle"),
                " No group information found in sample sheet. Please assign groups manually."
            ),
            selectInput("controlSamples", "Select Control Samples:",
                        choices = samples,
                        multiple = TRUE,
                        selected = NULL
            ),
            helpText("Samples not selected as controls will be assigned to the Case group."),
            actionButton("processDataIDAT", "Process with Group Assignment",
                        class = "btn-primary", icon = icon("play"))
          )
        })
      }

      removeModal()
      updateDataSummary()

      success_msg <- paste("IDAT files processed successfully!",
                          ncol(values$betaMatrix), "samples,",
                          nrow(values$betaMatrix), "probes")
      showNotification(success_msg, type = "message", duration = 8)

      # Generate QC outputs
      generateQCOutputs()

    }, error = function(e) {
      removeModal()
      error_msg <- paste("IDAT processing failed:", e$message)
      showNotification(error_msg, type = "error", duration = 15)

      # Detailed error logging
      cat("=== IDAT PROCESSING ERROR ===\n")
      cat("Error:", e$message, "\n")
      cat("============================\n")
      print(traceback())
    })
  })

  # Handle manual group assignment for IDAT data
  observeEvent(input$processDataIDAT, {
    req(values$betaMatrix, input$controlSamples, values$qc_metrics)

    showModal(modalDialog(
      title = "Processing IDAT Data for GIMP Analysis",
      "Extracting ICRs and creating analysis matrices...",
      footer = NULL
    ))

    # Create sample info for IDAT data
    samples <- colnames(values$betaMatrix)
    values$sampleInfo <- ifelse(samples %in% input$controlSamples, "Control", "Case")

    # Process data with GIMP functions
    tryCatch({
      # Map IDAT array type to GIMP bedmeth
      bedmeth_type <- switch(values$qc_metrics$array_type,
                            "450k" = "450k",
                            "EPIC" = "v1",
                            "EPICv2" = "v2")

      values$ICRcpg <- make_cpgs(Bmatrix = values$betaMatrix, bedmeth = bedmeth_type)
      values$df.ICR <- make_ICRs(Bmatrix = values$betaMatrix, bedmeth = bedmeth_type)
      values$processed <- TRUE

      removeModal()

      success_msg <- paste("IDAT data processed for GIMP analysis!",
                          "Groups:", paste(table(values$sampleInfo), collapse = " vs "))
      showNotification(success_msg, type = "message")

      # Update QC summary to include group info
      generateQCOutputs()

    }, error = function(e) {
      removeModal()
      showNotification(paste("GIMP processing failed:", e$message), type = "error", duration = 10)
    })
  })

  # ============================================================================
  # GEO DATASET PROCESSING - MULTI-STEP WORKFLOW
  # ============================================================================

  # Step 1: GEO validation (improved)
  observeEvent(input$validateGEO, {
    req(input$geoId)

    if (nchar(trimws(input$geoId)) == 0) {
      showNotification("Please enter a GEO ID first.", type = "warning")
      return()
    }

    # Reset step state
    values$geo_step2_ready <- FALSE
    values$geo_step3_ready <- FALSE
    values$geo_step4_ready <- FALSE
    values$geo_pheno_data <- NULL

    showModal(modalDialog(
      title = "Validating GEO Dataset",
      div(
        h4("Checking GEO dataset..."),
        p("Validating:", input$geoId),
        br(),
        div(class = "progress progress-striped active",
            div(class = "progress-bar progress-bar-info",
                style = "width: 100%", "Validating...")
        )
      ),
      footer = NULL,
      size = "m"
    ))

    tryCatch({
      # First validate the dataset
      validation_result <- validate_geo_dataset(input$geoId)
      values$geo_validation <- validation_result

      if (validation_result$valid && validation_result$has_idats) {
        # If valid, immediately get phenotypic data for step 2
        pheno_result <- get_geo_phenotype_data(input$geoId)
        
        if (pheno_result$success) {
          values$geo_pheno_data <- pheno_result
          values$geo_step2_ready <- TRUE
          
          removeModal()
          showNotification("Dataset validated! Proceeding to phenotypic data review.", 
                          type = "message", duration = 5)
        } else {
          removeModal()
          showNotification(paste("Failed to get phenotypic data:", pheno_result$error_message), 
                          type = "error", duration = 10)
        }
      } else {
        removeModal()
        if (validation_result$valid) {
          showNotification("Dataset found but no IDAT files detected. GIMP requires raw IDAT data.",
                           type = "warning", duration = 10)
        } else {
          showNotification(paste("Validation failed:", validation_result$message),
                           type = "error", duration = 10)
        }
      }

    }, error = function(e) {
      removeModal()
      showNotification(paste("GEO validation error:", e$message), type = "error", duration = 15)
    })
  })

  # UI Output for validation results in step 1
  output$geoValidationUI <- renderUI({
    req(values$geo_validation)
    
    validation_result <- values$geo_validation
    
    if (validation_result$valid && validation_result$has_idats) {
      div(class = "success-box",
        h5(icon("check-circle"), " Dataset Validated Successfully"),
        p(strong("Title:"), validation_result$title),
        p(strong("Samples:"), validation_result$sample_count),
        p(strong("Array Type:"), validation_result$array_type),
        p("✅ Contains IDAT files - suitable for GIMP analysis")
      )
    } else if (validation_result$valid) {
      div(class = "warning-box",
        h5(icon("exclamation-triangle"), " Dataset Found - No IDAT Files"),
        p(strong("Title:"), validation_result$title),
        p("❌ No IDAT files detected - not suitable for GIMP analysis")
      )
    } else {
      div(class = "warning-box",
        h5(icon("times-circle"), " Validation Failed"),
        p(validation_result$message)
      )
    }
  })

  # Step 2: Phenotypic data table and column selection
  output$phenoDataTable <- DT::renderDataTable({
    req(values$geo_pheno_data)
    
    pheno_data <- values$geo_pheno_data$pheno_data
    
    # Show only first 100 rows and key columns for performance
    display_data <- pheno_data[1:min(100, nrow(pheno_data)), ]
    
    DT::datatable(display_data, 
                  options = list(
                    pageLength = 10,
                    scrollX = TRUE,
                    scrollY = "300px",
                    dom = 'tip'
                  ),
                  rownames = TRUE)
  })

  # Update group column choices when step 2 is ready
  observe({
    req(values$geo_pheno_data)
    
    pheno_data <- values$geo_pheno_data$pheno_data
    column_choices <- colnames(pheno_data)
    names(column_choices) <- column_choices
    
    # Prioritize likely group columns
    potential_cols <- values$geo_pheno_data$potential_group_columns
    if (length(potential_cols) > 0) {
      priority_choices <- column_choices[column_choices %in% potential_cols]
      other_choices <- column_choices[!column_choices %in% potential_cols]
      column_choices <- c(priority_choices, other_choices)
    }
    
    updateSelectInput(session, "groupColumn", 
                     choices = column_choices,
                     selected = NULL)
  })

  # Preview selected column
  output$selectedColumnPreview <- renderPrint({
    req(input$groupColumn, values$geo_pheno_data)
    
    pheno_data <- values$geo_pheno_data$pheno_data
    
    if (input$groupColumn %in% colnames(pheno_data)) {
      col_data <- pheno_data[[input$groupColumn]]
      unique_vals <- unique(as.character(col_data))
      unique_vals <- unique_vals[!is.na(unique_vals) & unique_vals != ""]
      
      cat("Unique values in '", input$groupColumn, "':\n")
      cat("Total unique values:", length(unique_vals), "\n\n")
      
      if (length(unique_vals) <= 20) {
        for (val in unique_vals) {
          count <- sum(col_data == val, na.rm = TRUE)
          cat("- ", val, " (", count, " samples)\n")
        }
      } else {
        cat("First 20 values:\n")
        for (val in unique_vals[1:20]) {
          count <- sum(col_data == val, na.rm = TRUE)
          cat("- ", val, " (", count, " samples)\n")
        }
        cat("... and", length(unique_vals) - 20, "more values\n")
      }
    }
  })

  # Step 3: Proceed to group mapping
  observeEvent(input$proceedToStep3, {
    req(input$groupColumn, values$geo_pheno_data)
    
    if (input$groupColumn == "") {
      showNotification("Please select a grouping column first.", type = "warning")
      return()
    }
    
    pheno_data <- values$geo_pheno_data$pheno_data
    
    if (!input$groupColumn %in% colnames(pheno_data)) {
      showNotification("Selected column not found in data.", type = "error")
      return()
    }
    
    # Get unique values for mapping
    col_data <- pheno_data[[input$groupColumn]]
    unique_vals <- unique(as.character(col_data))
    unique_vals <- unique_vals[!is.na(unique_vals) & unique_vals != ""]
    
    if (length(unique_vals) == 0) {
      showNotification("Selected column contains no valid values.", type = "error")
      return()
    }
    
    if (length(unique_vals) > 50) {
      showNotification("Selected column has too many unique values (>50). Please choose a different column.", type = "warning")
      return()
    }
    
    values$selected_group_column <- input$groupColumn
    values$unique_group_values <- unique_vals
    values$geo_step3_ready <- TRUE
    
    showNotification("Proceeding to group mapping step.", type = "message")
  })

  # Generate group mapping UI for step 3
  output$groupMappingUI <- renderUI({
    req(values$unique_group_values)
    
    unique_vals <- values$unique_group_values
    
    # Create mapping inputs for each unique value
    mapping_inputs <- lapply(unique_vals, function(val) {
      div(style = "margin-bottom: 10px;",
        fluidRow(
          column(4, strong(val)),
          column(8,
            selectInput(paste0("mapping_", make.names(val)),
                       label = NULL,
                       choices = list("Case" = "Case", "Control" = "Control", "Exclude" = "Exclude"),
                       selected = "Exclude",
                       width = "100%")
          )
        )
      )
    })
    
    div(
      p("Map each unique value to Case, Control, or Exclude:"),
      do.call(tagList, mapping_inputs)
    )
  })

  # Group mapping summary
  output$groupMappingSummary <- renderPrint({
    req(values$unique_group_values, values$geo_pheno_data)
    
    unique_vals <- values$unique_group_values
    pheno_data <- values$geo_pheno_data$pheno_data
    col_data <- pheno_data[[values$selected_group_column]]
    
    # Get current mappings from inputs
    case_count <- 0
    control_count <- 0
    exclude_count <- 0
    
    for (val in unique_vals) {
      input_id <- paste0("mapping_", make.names(val))
      mapping <- input[[input_id]]
      
      if (!is.null(mapping)) {
        sample_count <- sum(col_data == val, na.rm = TRUE)
        
        if (mapping == "Case") {
          case_count <- case_count + sample_count
        } else if (mapping == "Control") {
          control_count <- control_count + sample_count
        } else {
          exclude_count <- exclude_count + sample_count
        }
      }
    }
    
    cat("Mapping Summary:\n===============\n\n")
    cat("Case samples:", case_count, "\n")
    cat("Control samples:", control_count, "\n")
    cat("Excluded samples:", exclude_count, "\n")
    cat("Total included:", case_count + control_count, "\n\n")
    
    if (case_count + control_count == 0) {
      cat("Warning: No samples assigned to Case or Control!\n")
    } else if (case_count == 0) {
      cat("Warning: No Case samples assigned.\n")
    } else if (control_count == 0) {
      cat("Warning: No Control samples assigned.\n")
    } else {
      cat("Ready for processing!\n")
    }
  })

  # Step 4: Process with mappings
  observeEvent(input$proceedToProcessing, {
    req(values$unique_group_values, values$geo_pheno_data, input$geoId)
    
    # Collect group mappings from UI
    unique_vals <- values$unique_group_values
    mappings <- list()
    
    for (val in unique_vals) {
      input_id <- paste0("mapping_", make.names(val))
      mapping <- input[[input_id]]
      if (!is.null(mapping)) {
        mappings[[val]] <- mapping
      }
    }
    
    # Validate mappings
    pheno_data <- values$geo_pheno_data$pheno_data
    col_data <- pheno_data[[values$selected_group_column]]
    
    case_count <- 0
    control_count <- 0
    
    for (val in unique_vals) {
      if (mappings[[val]] %in% c("Case", "Control")) {
        sample_count <- sum(col_data == val, na.rm = TRUE)
        if (mappings[[val]] == "Case") {
          case_count <- case_count + sample_count
        } else {
          control_count <- control_count + sample_count
        }
      }
    }
    
    if (case_count + control_count == 0) {
      showNotification("Please assign at least some samples to Case or Control.", type = "warning")
      return()
    }
    
    values$group_mappings <- mappings
    values$geo_step4_ready <- TRUE
    
    showNotification("Starting GEO dataset processing...", type = "message")
  })

  # Process GEO dataset with user mappings
  observe({
    req(values$geo_step4_ready, values$group_mappings, input$geoId)
    
    if (values$geo_step4_ready && !is.null(values$group_mappings)) {
      showModal(modalDialog(
        title = "Processing GEO Dataset",
        div(
          h4("Downloading and processing GEO data..."),
          p("This may take 10-30 minutes. Please be patient."),
          div(class = "progress progress-striped active",
              div(class = "progress-bar progress-bar-success",
                  style = "width: 100%", "Processing...")
          )
        ),
        footer = NULL,
        size = "m"
      ))
      
      tryCatch({
        n_cores <- if (input$geoEnableParallel) input$geoNCores else NULL
        
        geo_result <- process_geo_with_mappings(
          geo_id = input$geoId,
          group_column = values$selected_group_column,
          group_mappings = values$group_mappings,
          max_samples = input$geoMaxSamples,
          normalize_method = input$geoNormMethod,
          n_cores = n_cores
        )
        
        # Store results
        values$betaMatrix <- geo_result$beta_matrix
        values$sampleInfo <- geo_result$sample_info  # sample_info is already the Sample_Group vector
        values$qc_metrics <- geo_result$qc_metrics
        values$geo_metadata <- geo_result$geo_metadata
        
        # Process for GIMP analysis
        bedmeth_type <- switch(geo_result$geo_metadata$array_type,
                              "450k" = "450k",
                              "EPIC" = "v1",
                              "EPICv2" = "v2")
        
        values$ICRcpg <- make_cpgs(Bmatrix = values$betaMatrix, bedmeth = bedmeth_type)
        values$df.ICR <- make_ICRs(Bmatrix = values$betaMatrix, bedmeth = bedmeth_type)
        values$processed <- TRUE
        
        removeModal()
        updateDataSummary()
        generateQCOutputs()
        generateGEOMetadata()
        
        success_msg <- paste("GEO dataset processed successfully!",
                             ncol(values$betaMatrix), "samples from", input$geoId)
        showNotification(success_msg, type = "message", duration = 10)
        
      }, error = function(e) {
        removeModal()
        showNotification(paste("GEO processing failed:", e$message), type = "error", duration = 20)
      })
      
      values$geo_step4_ready <- FALSE
    }
  })

  # Reset GEO process
  observeEvent(input$resetGeoProcess, {
    values$geo_validation <- NULL
    values$geo_pheno_data <- NULL
    values$geo_step2_ready <- FALSE
    values$geo_step3_ready <- FALSE
    values$geo_step4_ready <- FALSE
    values$selected_group_column <- NULL
    values$group_mappings <- NULL
    values$unique_group_values <- NULL
    
    updateTextInput(session, "geoId", value = "")
    showNotification("GEO analysis reset. You can start a new analysis.", type = "message")
  })

  # Processing log output
  output$geoProcessingLog <- renderPrint({
    if (values$geo_step4_ready) {
      cat("Processing steps:\n")
      cat("Dataset validated\n")
      cat("Phenotypic data reviewed\n")
      cat("Groups mapped\n")
      cat("Downloading and processing...\n")
    }
  })

  # ============================================================================
  # HELPER FUNCTIONS
  # ============================================================================

  # Update data summary function
  updateDataSummary <- function() {
    output$dataSummary <- renderPrint({
      if (!is.null(values$betaMatrix)) {
        cat("Data dimensions:\n")
        cat(sprintf("  CpGs: %d\n", nrow(values$betaMatrix)))
        cat(sprintf("  Samples: %d\n", ncol(values$betaMatrix)))
        cat("\nSample names:\n")
        cat(paste(colnames(values$betaMatrix), collapse = ", "))

        if (!is.null(values$geo_metadata)) {
          cat("\n\nData Source: GEO dataset")
          cat("\nGEO ID:", values$geo_metadata$geo_id)
          cat("\nArray Type:", values$geo_metadata$array_type)
        } else if (!is.null(values$qc_metrics)) {
          cat("\n\nData Source: IDAT files")
          cat("\nArray Type:", values$qc_metrics$array_type)
          cat("\nNormalization:", values$qc_metrics$normalization)
        } else {
          cat("\n\nData Source: Processed data")
        }
      } else {
        cat("No data loaded.")
      }
    })

    output$dataPreview <- DT::renderDataTable({
      if (!is.null(values$betaMatrix)) {
        DT::datatable(
          head(values$betaMatrix, 100),
          options = list(
            pageLength = 10,
            scrollX = TRUE
          )
        )
      }
    })
  }

  # Generate QC outputs for IDAT/GEO data
  generateQCOutputs <- function() {
    output$qcSummary <- renderPrint({
      if (!is.null(values$geo_metadata)) {
        cat("GEO Dataset Processing Report\n")
        cat("============================\n\n")
        cat("GEO Information:\n")
        cat("  GEO ID:", values$geo_metadata$geo_id, "\n")
        cat("  Title:", values$geo_metadata$title, "\n")
        cat("  Array type:", values$geo_metadata$array_type, "\n")
        cat("  Original samples:", values$geo_metadata$original_sample_count, "\n")
        cat("  Processed samples:", values$geo_metadata$processed_sample_count, "\n")
        if (!is.null(values$geo_metadata$group_detection)) {
          cat("  Group detection:", values$geo_metadata$group_detection, "\n")
        }
        cat("\n")
      }

      if (!is.null(values$qc_metrics)) {
        if (is.null(values$geo_metadata)) {
          cat("IDAT Processing Report\n")
          cat("=====================\n\n")
        }

        cat("Processing Details:\n")
        cat("  Array type:", values$qc_metrics$array_type, "\n")
        cat("  Normalization:", values$qc_metrics$normalization, "\n")
        cat("  Samples processed:", values$qc_metrics$total_samples, "\n")
        cat("  Samples passed QC:", ncol(values$betaMatrix), "\n")
        cat("  Probes retained:", values$qc_metrics$total_probes, "\n")
        cat("  Detection p-value threshold:",
            if (!is.null(input$detectionPval)) input$detectionPval else 0.01, "\n\n")

        if (length(values$qc_metrics$failed_samples) > 0) {
          cat("Quality Control Issues:\n")
          cat("  Failed samples (>10% failed probes):\n")
          for (failed in values$qc_metrics$failed_samples) {
            cat("   -", failed, "\n")
          }
          cat("\n")
        }
      }

      if (!is.null(values$sampleInfo)) {
        cat("Sample Groups:\n")
        group_table <- table(values$sampleInfo)
        for (i in seq_len(length(group_table))) {
          cat("  ", names(group_table)[i], ":", group_table[i], "\n")
        }

        if (values$processed) {
          cat("\nGIMP Analysis Status: Ready for ICR analysis\n")
        }
      } else {
        cat("Sample Groups: Not assigned (manual assignment required)\n")
      }
    })

    # Enhanced QC plots
    output$qcPlots <- renderPlot({
      par(mfrow = c(2, 2))

      # 1. Failed probe proportion per sample
      if (!is.null(values$qc_metrics$failed_sample_proportion)) {
        failed_prop <- values$qc_metrics$failed_sample_proportion
        barplot(failed_prop,
                main = "Failed Probes per Sample",
                ylab = "Proportion of Failed Probes",
                xlab = "Samples",
                las = 2,
                cex.names = 0.5,
                col = ifelse(failed_prop > 0.1, "red", "lightblue"))
        abline(h = 0.1, col = "red", lty = 2)
        legend("topright", legend = c("Pass", "Fail"),
               fill = c("lightblue", "red"), cex = 0.8)
      }

      # 2. Beta value distribution
      beta_sample <- values$betaMatrix[sample(nrow(values$betaMatrix),
                                            min(10000, nrow(values$betaMatrix))), ]
      # Remove non-finite values (NaN, Inf, NA)
      beta_sample_clean <- beta_sample[is.finite(beta_sample)]
      
      if (length(beta_sample_clean) > 0) {
        hist(beta_sample_clean,
             main = "Beta Value Distribution",
             xlab = "Beta Values",
             ylab = "Frequency",
             breaks = 50,
             col = "lightgreen",
             border = "darkgreen")
      } else {
        plot(1, type = "n", axes = FALSE, xlab = "", ylab = "", main = "Beta Value Distribution")
        text(1, 1, "No finite beta values found", cex = 1.2, col = "red")
      }

      # 3. Sample correlation heatmap (subset of samples if too many)
      if (ncol(values$betaMatrix) <= 20) {
        cor_matrix <- cor(values$betaMatrix, use = "complete.obs")
        heatmap(cor_matrix,
                main = "Sample Correlation",
                col = colorRampPalette(c("blue", "white", "red"))(50),
                cexRow = 0.7, cexCol = 0.7)
      } else {
        # Plot mean correlation per sample instead
        sample_subset <- sample(ncol(values$betaMatrix), 20)
        cor_matrix <- cor(values$betaMatrix[, sample_subset], use = "complete.obs")
        heatmap(cor_matrix,
                main = "Sample Correlation (20 random samples)",
                col = colorRampPalette(c("blue", "white", "red"))(50),
                cexRow = 0.7, cexCol = 0.7)
      }

      # 4. Probe intensity distribution
      if (ncol(values$betaMatrix) > 0) {
        probe_means <- rowMeans(values$betaMatrix, na.rm = TRUE)
        # Remove non-finite values
        probe_means_clean <- probe_means[is.finite(probe_means)]
        
        if (length(probe_means_clean) > 0) {
          hist(probe_means_clean,
               main = "Probe Mean Intensity",
               xlab = "Mean Beta Value",
               ylab = "Frequency",
               breaks = 50,
               col = "lightcoral")
        } else {
          plot(1, type = "n", axes = FALSE, xlab = "", ylab = "", main = "Probe Mean Intensity")
          text(1, 1, "No finite probe means found", cex = 1.2, col = "red")
        }
      }
    })
  }

  # Generate GEO metadata output
  generateGEOMetadata <- function() {
    output$geoMetadata <- renderPrint({
      if (!is.null(values$geo_metadata)) {
        cat("GEO Dataset Metadata\n")
        cat("===================\n\n")
        cat("Dataset Information:\n")
        cat("  GEO Accession:", values$geo_metadata$geo_id, "\n")
        cat("  Title:", values$geo_metadata$title, "\n")
        cat("  Array Platform:", values$geo_metadata$array_type, "\n")
        cat("  Original Sample Count:", values$geo_metadata$original_sample_count, "\n")
        cat("  Processed Sample Count:", values$geo_metadata$processed_sample_count, "\n")

        if (!is.null(values$geo_metadata$group_detection)) {
          cat("  Group Detection Method:", values$geo_metadata$group_detection, "\n")
        }

        cat("\nProcessing Information:\n")
        cat("  Download Directory:", values$geo_metadata$download_dir, "\n")
        cat("  Processing Date:", Sys.time(), "\n")

        if (!is.null(values$sampleInfo)) {
          cat("\nSample Group Distribution:\n")
          group_table <- table(values$sampleInfo)
          for (i in seq_len(length(group_table))) {
            cat("  ", names(group_table)[i], ":", group_table[i], "samples\n")
          }
        }

        cat("\nGIMP Analysis Status:")
        if (values$processed) {
          cat(" ✅ Ready for ICR analysis\n")
        } else {
          cat(" ⏳ Waiting for group assignment\n")
        }

      } else {
        cat("No GEO metadata available.\n")
        cat("This output is only available for datasets processed from GEO.")
      }
    })
  }

  # ============================================================================
  # CPG COVERAGE ANALYSIS
  # ============================================================================

  # CpG Coverage Analysis
  observeEvent(input$runCpGAnalysis, {
    req(values$ICRcpg)

    showModal(modalDialog(
      title = "Running Analysis",
      "Calculating CpG coverage...",
      footer = NULL
    ))

    tryCatch({
      # Determine bedmeth type
      if (!is.null(values$geo_metadata)) {
        bedmeth_type <- switch(values$geo_metadata$array_type,
                              "450k" = "450k",
                              "EPIC" = "v1",
                              "EPICv2" = "v2")
      } else if (!is.null(values$qc_metrics)) {
        bedmeth_type <- switch(values$qc_metrics$array_type,
                              "EPIC" = "v1",
                              "450k" = "450k",
                              "EPICv2" = "v2")
      } else {
        bedmeth_type <- input$arrayType
      }

      values$cpgs_analysis <- plot_cpgs_coverage(values$ICRcpg, bedmeth = bedmeth_type)

      removeModal()

      # Render plots
      output$cpgCountPlot <- renderPlotly({
        ggplotly(values$cpgs_analysis$plot_counts)
      })

      output$cpgPercentagePlot <- renderPlotly({
        ggplotly(values$cpgs_analysis$plot_percentage)
      })

      output$cpgCoverageTable <- DT::renderDataTable({
        DT::datatable(
          values$cpgs_analysis$data,
          options = list(pageLength = 25)
        )
      })

    }, error = function(e) {
      removeModal()
      showNotification(paste("Error:", e$message), type = "error", duration = 10)
    })
  })

  # ============================================================================
  # ICR HEATMAP ANALYSIS
  # ============================================================================

  # ICR Heatmap
  observeEvent(input$generateHeatmap, {
    req(values$df.ICR, values$sampleInfo)

    # Add validation
    if (is.null(values$df.ICR) || nrow(values$df.ICR) == 0) {
      showNotification("No ICR data available. Please process data first.", type = "error")
      return()
    }

    if (is.null(values$sampleInfo) || length(values$sampleInfo) != ncol(values$df.ICR)) {
      showNotification("Sample information doesn't match data dimensions.", type = "error")
      return()
    }

    showModal(modalDialog(
      title = "Generating Heatmap",
      "Creating ICR methylation heatmap...",
      footer = NULL
    ))

    tryCatch({
      # Determine bedmeth type
      if (!is.null(values$geo_metadata)) {
        bedmeth_type <- switch(values$geo_metadata$array_type,
                              "450k" = "450k",
                              "EPIC" = "v1",
                              "EPICv2" = "v2")
      } else if (!is.null(values$qc_metrics)) {
        bedmeth_type <- switch(values$qc_metrics$array_type,
                              "450k" = "450k",
                              "EPIC" = "v1",
                              "EPICv2" = "v2")
      } else {
        bedmeth_type <- input$arrayType
      }

      # Extract unique labels from sampleInfo dynamically
      unique_labels <- unique(values$sampleInfo)

      # Validate that we have exactly 2 groups
      if (length(unique_labels) != 2) {
        stop(paste("Expected exactly 2 groups, but found", length(unique_labels),
                 "groups:", paste(unique_labels, collapse = ", ")))
      }

      # Determine which label should be control and which should be case
      # Priority order: "Control" > "control" > first alphabetically
      if ("Control" %in% unique_labels) {
        control_label <- "Control"
        case_label <- setdiff(unique_labels, "Control")
      } else if ("control" %in% unique_labels) {
        control_label <- "control"
        case_label <- setdiff(unique_labels, "control")
      } else {
        # Sort alphabetically and assign first as control, second as case
        sorted_labels <- sort(unique_labels)
        control_label <- sorted_labels[1]
        case_label <- sorted_labels[2]

        # Notify user about the assignment
        showNotification(
          paste("Auto-assigned groups: Control =", control_label, ", Case =", case_label),
          type = "message", duration = 5
        )
      }

      # Generate the heatmap
      p <- ICRs_heatmap(
        df_ICR = values$df.ICR,
        sampleInfo = values$sampleInfo,
        control_label = control_label,
        case_label = case_label,
        bedmeth = bedmeth_type,
        order_by = input$orderBy,
        plot_type = input$plotType,
        sd_threshold = ifelse(input$plotType == "defect", input$sdThreshold, 3)
      )

      # Render the plot
      output$icrHeatmap <- renderPlot({
        print(p)
      }, height = 800)

      removeModal()
      showNotification("Heatmap generated successfully!", type = "message")

    }, error = function(e) {
      removeModal()

      # More detailed error reporting
      error_msg <- paste("Heatmap generation failed:", e$message)
      showNotification(error_msg, type = "error", duration = 15)

      # Print detailed error to console for debugging
      cat("=== HEATMAP ERROR DETAILS ===\n")
      cat("Error message:", e$message, "\n")
      cat("df_ICR dimensions:", dim(values$df.ICR), "\n")
      cat("sampleInfo:", paste(values$sampleInfo, collapse = ", "), "\n")
      cat("============================\n")
      print(traceback())
    })
  })

  # Download heatmap
  output$downloadHeatmap <- downloadHandler(
    filename = function() {
      paste0("ICR_heatmap_", input$plotType, "_", Sys.Date(), ".pdf")
    },
    content = function(file) {
      tryCatch({
        # Determine bedmeth type
        if (!is.null(values$geo_metadata)) {
          bedmeth_type <- switch(values$geo_metadata$array_type,
                                "450k" = "450k",
                                "EPIC" = "v1",
                                "EPICv2" = "v2")
        } else if (!is.null(values$qc_metrics)) {
          bedmeth_type <- switch(values$qc_metrics$array_type,
                                "450k" = "450k",
                                "EPIC" = "v1",
                                "EPICv2" = "v2")
        } else {
          bedmeth_type <- input$arrayType
        }

        unique_labels <- unique(values$sampleInfo)
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

        pdf(file, width = 12, height = 10)
        p <- ICRs_heatmap(
          df_ICR = values$df.ICR,
          sampleInfo = values$sampleInfo,
          control_label = control_label,
          case_label = case_label,
          bedmeth = bedmeth_type,
          order_by = input$orderBy,
          plot_type = input$plotType,
          sd_threshold = ifelse(input$plotType == "defect", input$sdThreshold, 3)
        )
        print(p)
        dev.off()
      }, error = function(e) {
        showNotification(paste("Download failed:", e$message), type = "error", duration = 10)
      })
    }
  )

  # ============================================================================
  # DIFFERENTIAL METHYLATION ANALYSIS
  # ============================================================================

  # DMP Analysis
  observeEvent(input$runDMPAnalysis, {
    req(values$ICRcpg, values$sampleInfo)

    showModal(modalDialog(
      title = "Running DMP Analysis",
      "Identifying differentially methylated positions...",
      footer = NULL
    ))

    tryCatch({
      values$dmps_results <- iDMPs(
        data = values$ICRcpg,
        sampleInfo = values$sampleInfo,
        pValueCutoff = input$pValueCutoff
      )

      removeModal()

      # Check if we have any significant results
      if (nrow(values$dmps_results$topDMPs) == 0) {
        showNotification(paste0("No significant DMPs found at p < ", input$pValueCutoff,
                               ". Try increasing the p-value cutoff."),
                         type = "warning", duration = 10)

        # Still show all results table for reference
        output$dmpsTable <- DT::renderDataTable({
          if (nrow(values$dmps_results$allResults) > 0) {
            DT::datatable(
              head(values$dmps_results$allResults, 100),
              options = list(
                pageLength = 25,
                scrollX = TRUE
              ),
              caption = "Top 100 results (none meet significance threshold)"
            ) %>%
              DT::formatRound(columns = c("logFC", "AveExpr", "t", "P.Value", "adj.P.Val"), digits = 4)
          } else {
            DT::datatable(data.frame(Message = "No results available"))
          }
        })

        # Clear ICR choices since no significant DMPs
        updateSelectInput(session, "selectedICR", choices = NULL)

      } else {
        # We have significant results
        # Update ICR choices for region explorer
        icr_choices <- unique(values$dmps_results$topDMPs$ICR)
        updateSelectInput(session, "selectedICR",
                          choices = icr_choices,
                          selected = icr_choices[1]
        )

        # Render results table
        output$dmpsTable <- DT::renderDataTable({
          DT::datatable(
            values$dmps_results$topDMPs,
            options = list(
              pageLength = 25,
              scrollX = TRUE
            )
          ) %>%
            DT::formatRound(columns = c("logFC", "AveExpr", "t", "P.Value", "adj.P.Val"), digits = 4)
        })
      }

      # Volcano plot - use all results
      output$volcanoPlot <- renderPlotly({
        if (nrow(values$dmps_results$allResults) > 0) {
          plot_data <- values$dmps_results$allResults
          plot_data$neg_log_p <- -log10(plot_data$P.Value)
          plot_data$significant <- plot_data$P.Value < input$pValueCutoff

          p <- ggplot(plot_data, aes(x = logFC, y = neg_log_p, color = significant,
                                     text = paste("ICR:", ICR, "<br>",
                                                  "CpGID:", rownames(plot_data), "<br>",
                                                  "logFC:", round(logFC, 3), "<br>",
                                                  "P.Value:", format(P.Value, scientific = TRUE)))) +
            geom_point(alpha = 0.6) +
            scale_color_manual(values = c("FALSE" = "grey", "TRUE" = "red")) +
            geom_hline(yintercept = -log10(input$pValueCutoff), linetype = "dashed") +
            geom_vline(xintercept = c(-0.5, 0.5), linetype = "dashed") +
            labs(x = "Log Fold Change", y = "-log10(p-value)", title = "Volcano Plot") +
            theme_bw()

          ggplotly(p, tooltip = "text")
        } else {
          # Empty plot if no results
          ggplot() +
            annotate("text", x = 0, y = 0, label = "No results available") +
            theme_void()
        }
      })

      # Summary
      output$dmpsSummary <- renderPrint({
        cat("DMP Analysis Summary\n")
        cat("==================\n\n")
        cat(sprintf("Total CpGs tested: %d\n", nrow(values$dmps_results$allResults)))
        cat(sprintf("Significant DMPs (p < %g): %d\n", input$pValueCutoff, nrow(values$dmps_results$topDMPs)))

        if (nrow(values$dmps_results$topDMPs) > 0) {
          cat(sprintf("\nTop 5 DMPs by p-value:\n"))
          top5 <- head(values$dmps_results$topDMPs[order(values$dmps_results$topDMPs$P.Value),
                                                   c("ICR", "logFC", "P.Value")], 5)
          print(top5)
        } else {
          cat("\nNo significant DMPs found. Consider:\n")
          cat("- Increasing the p-value cutoff\n")
          cat("- Checking sample group assignments\n")
          cat("- Verifying data quality\n")
        }
      })

    }, error = function(e) {
      removeModal()
      showNotification(paste("DMP Analysis Error:", e$message), type = "error", duration = 15)

      # Print debug info
      cat("=== DMP ANALYSIS ERROR ===\n")
      cat("Error:", e$message, "\n")
      cat("ICRcpg dimensions:", dim(values$ICRcpg), "\n")
      cat("ICRcpg columns:", colnames(values$ICRcpg), "\n")
      cat("sampleInfo:", paste(values$sampleInfo, collapse = ", "), "\n")
      cat("========================\n")
      print(traceback())
    })
  })

  # Download DMP results
  output$downloadDMPs <- downloadHandler(
    filename = function() {
      paste0("DMP_results_", Sys.Date(), ".csv")
    },
    content = function(file) {
      write.csv(values$dmps_results$topDMPs, file, row.names = FALSE)
    }
  )

  # ============================================================================
  # REGION EXPLORER
  # ============================================================================

  # Refresh ICR list
  observeEvent(input$refreshICRList, {
    if (!is.null(values$dmps_results) && nrow(values$dmps_results$topDMPs) > 0) {
      icr_choices <- unique(values$dmps_results$topDMPs$ICR)
      updateSelectInput(session, "selectedICR",
                        choices = icr_choices,
                        selected = icr_choices[1]
      )
      showNotification(paste("Found", length(icr_choices), "ICRs with significant DMPs"),
                       type = "message")
    } else {
      updateSelectInput(session, "selectedICR", choices = NULL)
      showNotification("No significant DMPs found. Run DMP analysis first.", type = "warning")
    }
  })

  # Region Explorer
  observeEvent(input$plotRegion, {
    req(input$selectedICR, values$dmps_results, values$ICRcpg)

    # Validation checks
    if (is.null(input$selectedICR) || input$selectedICR == "") {
      showNotification("Please select an ICR first.", type = "warning")
      return()
    }

    if (nrow(values$dmps_results$topDMPs) == 0) {
      showNotification("No significant DMPs available for plotting.", type = "warning")
      return()
    }

    showModal(modalDialog(
      title = "Generating Plot",
      "Creating region plot...",
      footer = NULL
    ))

    # Clear previous plots
    output$regionPlot <- renderPlot({ NULL })
    output$regionPlotUI <- renderUI({ NULL })

    tryCatch({
      # Check if selected ICR exists in DMPs
      if (!input$selectedICR %in% values$dmps_results$topDMPs$ICR) {
        stop(paste("Selected ICR", input$selectedICR, "not found in significant DMPs"))
      }

      # Create the plot
      plot_interactive <- as.logical(input$plotInteractive)

      plot <- plot_line_ICR(
        significantDMPs = values$dmps_results$topDMPs,
        ICRcpg = values$ICRcpg,
        ICR = input$selectedICR,
        sampleInfo = values$sampleInfo,
        interactive = plot_interactive
      )

      removeModal()

      # Render the appropriate UI based on plot type
      if (plot_interactive) {
        output$regionPlotUI <- renderUI({
          tagList(
            h4(paste("Interactive Plot:", input$selectedICR)),
            plotlyOutput("regionPlotInteractive", height = "600px"),
            br(),
            helpText("Use plotly controls to zoom, pan, and hover over data points.")
          )
        })

        output$regionPlotInteractive <- renderPlotly({
          plot
        })
      } else {
        output$regionPlotUI <- renderUI({
          tagList(
            h4(paste("Static Plot:", input$selectedICR)),
            plotOutput("regionPlotStatic", height = "600px")
          )
        })

        output$regionPlotStatic <- renderPlot({
          print(plot)
        }, height = 600)
      }

      showNotification("Region plot generated successfully!", type = "message")

    }, error = function(e) {
      removeModal()

      error_msg <- paste("Region plot failed:", e$message)
      showNotification(error_msg, type = "error", duration = 15)

      # Detailed error logging
      cat("=== REGION PLOT ERROR ===\n")
      cat("Error message:", e$message, "\n")
      cat("Selected ICR:", input$selectedICR, "\n")
      cat("Interactive setting:", input$plotInteractive, "\n")
      cat("DMPs available:", nrow(values$dmps_results$topDMPs), "\n")
      cat("Stack trace:\n")
      print(traceback())
      cat("========================\n")

      # Show fallback message
      output$regionPlotUI <- renderUI({
        div(
          style = "text-align: center; padding: 50px;",
          h4("Plot Generation Failed", style = "color: red;"),
          p("Error:", e$message),
          p("Please try selecting a different ICR or check the console for detailed error information.")
        )
      })
    })
  })

  # Add a reset when ICR selection changes
  observeEvent(input$selectedICR, {
    # Clear plots when ICR selection changes
    output$regionPlotUI <- renderUI({
      div(
        style = "text-align: center; padding: 50px;",
        h4("Select Plot Type and Click 'Plot Region'"),
        p("Choose an ICR and plot type, then click the 'Plot Region' button to generate the visualization.")
      )
    })
  })
}

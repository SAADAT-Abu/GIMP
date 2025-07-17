#' Launch GIMP Shiny Application
#'
#' @description 
#' Launches an interactive Shiny application for GIMP analysis. The app provides
#' a graphical user interface for analyzing methylation patterns at Imprinted 
#' Control Regions (ICRs) without requiring R programming knowledge.
#'
#' @param max_upload_size_mb Maximum file upload size in MB (default: 500)
#' 
#' @details 
#' The GIMP Shiny app includes the following features:
#' \itemize{
#'   \item Upload methylation beta matrices from CSV files or raw IDAT files
#'   \item Analyze CpG coverage at ICRs
#'   \item Generate methylation heatmaps (beta, delta, and defect plots)
#'   \item Identify differentially methylated positions (DMPs)
#'   \item Explore specific ICR regions with interactive visualizations
#'   \item Export results and plots
#' }
#'
#' @return 
#' Opens the Shiny application in the default web browser. The function returns
#' invisibly once the app is closed.
#'
#' @examples
#' # Launch the GIMP Shiny app
#' GIMP_app()
#' 
#' # Launch with larger upload limit
#' GIMP_app(max_upload_size_mb = 1000)
#'
#' @export
GIMP_app <- function(max_upload_size_mb = 500) {
  # Check if required packages are installed
  required_packages <- c("shiny", "shinydashboard", "DT", "plotly",
                         "ggplot2", "reshape2", "tidyverse")
  missing_packages <- required_packages[!(required_packages %in% 
                                         rownames(installed.packages()))]
  
  if (length(missing_packages) > 0) {
    stop(paste("The following required packages are not installed:", 
               paste(missing_packages, collapse = ", "), 
               "\nPlease install them using: install.packages(c(", 
               paste(paste0('"', missing_packages, '"'), collapse = ", "), 
               "))"),
         call. = FALSE)
  }
  
  # Set upload size limit
  options(shiny.maxRequestSize = max_upload_size_mb * 1024^2)
  
  # Find the app directory
  appDir <- system.file("shiny", package = "GIMP")
  
  if (appDir == "") {
    stop("Could not find the Shiny app directory. Try re-installing the GIMP package.", 
         call. = FALSE)
  }
  
  # Check if app files exist
  required_files <- c("app.R", "ui.R", "server.R")
  missing_files <- required_files[!file.exists(file.path(appDir, 
                                                        required_files))]
  
  if (length(missing_files) > 0) {
    stop(paste("The following required app files are missing:", 
               paste(missing_files, collapse = ", ")),
         call. = FALSE)
  }
  
  # Launch the app
  message("Launching GIMP Shiny application...")
  message("Upload size limit set to: ", max_upload_size_mb, " MB")
  shiny::runApp(appDir, display.mode = "normal", launch.browser = TRUE)
}
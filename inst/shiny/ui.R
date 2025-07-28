library(shiny)
library(shinydashboard)
library(DT)
library(plotly)

ui <- dashboardPage(
  dashboardHeader(title = "GIMP"),

  dashboardSidebar(
    sidebarMenu(
      menuItem("Data Upload", tabName = "upload", icon = icon("upload")),
      menuItem("CpG Analysis", tabName = "cpg", icon = icon("dna")),
      menuItem("ICR Analysis", tabName = "icr", icon = icon("chart-bar")),
      menuItem("Differential Analysis", tabName = "dmps", icon = icon("chart-line")),
      menuItem("Region Explorer", tabName = "explorer", icon = icon("search")),
      menuItem("Help", tabName = "help", icon = icon("question-circle"))
    )
  ),

  dashboardBody(
    tags$head(
      tags$style(HTML("
        .content-wrapper, .right-side {
          background-color: #f4f4f4;
        }
        .box-header {
          background-color: #3c8dbc;
          color: white;
        }
        .file-validation {
          margin-top: 5px;
          font-size: 12px;
        }
        .info-box {
          background-color: #e8f4fd;
          border: 1px solid #b8daff;
          border-radius: 5px;
          padding: 15px;
          margin: 10px 0;
        }
        .warning-box {
          background-color: #fff3cd;
          border: 1px solid #ffeaa7;
          border-radius: 5px;
          padding: 15px;
          margin: 10px 0;
        }
        .geo-box {
          background-color: #e8f5e8;
          border: 1px solid #c3e6c3;
          border-radius: 5px;
          padding: 15px;
          margin: 10px 0;
        }
        .success-box {
          background-color: #d4edda;
          border: 1px solid #c3e6cb;
          border-radius: 5px;
          padding: 15px;
          margin: 10px 0;
        }
      ")),

      tags$script(HTML("
        // Custom message handler for updating array type
        Shiny.addCustomMessageHandler('updateArrayType', function(value) {
          $('#arrayType').val(value).trigger('change');
        });

        // File validation feedback
        $(document).on('change', '#idatZip', function() {
          var file = this.files[0];
          if (file) {
            if (file.name.toLowerCase().endsWith('.zip')) {
              $('#zip-validation').html('<i class=\"fa fa-check text-success\"></i> ZIP file selected');
            } else {
              $('#zip-validation').html('<i class=\"fa fa-times text-danger\"></i> Please select a ZIP file');
            }
          }
        });

        // GEO ID validation
        $(document).on('input', '#geoId', function() {
          var geoId = $(this).val().toUpperCase();
          if (geoId.match(/^GSE\\d+$/)) {
            $('#geo-validation').html('<i class=\"fa fa-check text-success\"></i> Valid GEO ID format');
          } else if (geoId.length > 0) {
            $('#geo-validation').html('<i class=\"fa fa-times text-danger\"></i> Invalid format (use GSExxxxx)');
          } else {
            $('#geo-validation').html('');
          }
        });
      "))
    ),

    tabItems(
      # Data Upload Tab
      tabItem(tabName = "upload",
        fluidRow(
          box(
            title = "Upload Methylation Data",
            status = "primary",
            solidHeader = TRUE,
            width = 12,

            # Data source selection
            fluidRow(
              column(12,
                radioButtons("dataSource", "Select Data Source:",
                  choices = list(
                    "Processed Data (CSV/RDS/Excel)" = "processed",
                    "Raw IDAT Files (ZIP)" = "idat",
                    "GEO Dataset (Download from NCBI)" = "geo"
                  ),
                  selected = "processed",
                  inline = TRUE
                ),
                hr()
              )
            ),

            # Processed data upload
            conditionalPanel(
              condition = "input.dataSource == 'processed'",

              fluidRow(
                column(6,
                  fileInput("betaMatrix",
                    "Choose Beta Matrix File",
                    accept = c(".csv", ".txt", ".rds", ".xlsx"),
                    placeholder = "Rows: CpGs, Columns: Samples"
                  ),

                  helpText("Upload a file where rows are CpG probes and columns are samples.
                           Supports CSV, TXT, RDS, and Excel formats."),

                  radioButtons("arrayType", "Array Type:",
                    choices = list(
                      "450k (hg19)" = "450k",
                      "EPIC v1 (hg19)" = "v1",
                      "EPIC v2 (hg38)" = "v2"
                    ),
                    selected = "v1"
                  )
                ),

                column(6,
                  h4("Sample Group Assignment"),
                  uiOutput("sampleGroupUI"),

                  actionButton("processData", "Process Data",
                    class = "btn-primary",
                    icon = icon("play")
                  )
                )
              )
            ),

            # IDAT file upload
            conditionalPanel(
              condition = "input.dataSource == 'idat'",

              fluidRow(
                column(6,
                  fileInput("idatZip",
                    "Upload ZIP File with IDAT Files",
                    accept = c(".zip"),
                    placeholder = "ZIP containing IDAT files + sample sheet"
                  ),

                  div(id = "zip-validation", class = "file-validation"),

                  div(class = "info-box",
                    h5(icon("info-circle"), " ZIP File Requirements:"),
                    tags$ul(
                      tags$li("IDAT files (*.idat) - both Red and Grn for each sample"),
                      tags$li("Sample sheet CSV file with required columns:"),
                      tags$ul(
                        tags$li(strong("Sample_Name"), " - Unique sample identifiers"),
                        tags$li(strong("Sentrix_ID"), " - Slide/chip ID (e.g., 200123456789)"),
                        tags$li(strong("Sentrix_Position"), " - Array position (e.g., R01C01)")
                      ),
                      tags$li("Optional: Sample_Group column for automatic group assignment")
                    )
                  ),

                  textInput("sampleSheetName", "Sample Sheet Filename:",
                    value = "samplesheet.csv",
                    placeholder = "e.g., samplesheet.csv"
                  ),

                  helpText("If your sample sheet has a different name, specify it here."),

                  actionButton("previewZip", "Preview ZIP Contents",
                    class = "btn-info btn-sm",
                    icon = icon("eye")
                  )
                ),

                column(6,
                  selectInput("idatArrayType", "Array Type:",
                    choices = list(
                      "450k" = "450k",
                      "EPIC (850k)" = "EPIC",
                      "EPICv2 (935k)" = "EPICv2"
                    ),
                    selected = "EPIC"
                  ),

                  selectInput("normMethod", "Normalization Method:",
                    choices = list(
                      "Quantile (Recommended)" = "quantile",
                      "SWAN" = "SWAN",
                      "Functional Normalization" = "funnorm",
                      "Noob" = "noob"
                    ),
                    selected = "quantile"
                  ),

                  div(class = "warning-box",
                    h5(icon("exclamation-triangle"), " Processing Time:"),
                    p("IDAT processing can take 5-15 minutes depending on:"),
                    tags$ul(
                      tags$li("Number of samples"),
                      tags$li("Normalization method chosen"),
                      tags$li("Computer performance")
                    ),
                    p("Please be patient and don't close the browser.")
                  ),

                  hr(),

                  h5("Advanced Options:"),

                  numericInput("detectionPval", "Detection P-value Threshold:",
                    value = 0.01, min = 0.001, max = 0.05, step = 0.001
                  ),
                  helpText("Probes with detection p-value above this threshold are considered failed."),

                  checkboxInput("removeFailedSamples", "Remove Failed Samples", value = TRUE),
                  helpText("Remove samples with >10% failed probes."),

                  checkboxInput("enableParallel", "Enable Parallel Processing", value = FALSE),
                  helpText("Use multiple CPU cores to speed up processing (requires more memory)."),

                  conditionalPanel(
                    condition = "input.enableParallel == true",
                    numericInput("nCores", "Number of CPU Cores:",
                      value = 2, min = 1, max = 16, step = 1
                    ),
                    helpText("Number of CPU cores to use. System will automatically limit to available cores.")
                  ),

                  br(),
                  actionButton("processIDAT", "Process IDAT Files",
                    class = "btn-success btn-lg",
                    icon = icon("cogs"),
                    style = "width: 100%;"
                  )
                )
              ),

              # ZIP preview results
              conditionalPanel(
                condition = "output.zipPreviewAvailable",
                hr(),
                fluidRow(
                  column(12,
                    box(
                      title = "ZIP File Preview",
                      status = "info",
                      solidHeader = TRUE,
                      width = 12,
                      collapsible = TRUE,

                      verbatimTextOutput("zipPreview")
                    )
                  )
                )
              )
            ),

            # GEO dataset processing
            conditionalPanel(
              condition = "input.dataSource == 'geo'",

              fluidRow(
                column(6,
                  div(class = "geo-box",
                    h5(icon("database"), " Download from GEO (Gene Expression Omnibus)"),
                    p("Automatically download and process methylation datasets from NCBI GEO."),
                    p("GIMP will:"),
                    tags$ul(
                      tags$li("Validate the GEO ID and check for IDAT files"),
                      tags$li("Download raw IDAT files and phenotypic data"),
                      tags$li("Auto-detect sample groups when possible"),
                      tags$li("Process data through the standard GIMP pipeline")
                    )
                  ),

                  textInput("geoId",
                    "Enter GEO Accession ID:",
                    placeholder = "e.g., GSE68777",
                    value = ""
                  ),

                  div(id = "geo-validation", class = "file-validation"),

                  helpText("Enter a GEO Series (GSE) accession number. The dataset must contain raw IDAT files."),

                  div(
                    style = "margin-top: 15px;",
                    actionButton("validateGEO", "Validate GEO Dataset",
                      class = "btn-info",
                      icon = icon("search")
                    ),
                    helpText("Check if the dataset exists and contains IDAT files", style = "font-size: 11px; margin-top: 5px;")
                  )
                ),

                column(6,
                  h5("Processing Options:"),

                  textInput("geoGroupColumn", "Group Column (optional):",
                    placeholder = "e.g., characteristics_ch1.1",
                    value = ""
                  ),
                  helpText("Specify a phenotype column for sample groups, or leave blank for auto-detection."),

                  textInput("geoControlTerms", "Control Terms:",
                    value = "control,normal,healthy,ctrl",
                    placeholder = "comma-separated terms"
                  ),
                  helpText("Terms that identify control samples (comma-separated)."),

                  textInput("geoCaseTerms", "Case Terms:",
                    value = "case,disease,tumor,cancer,patient",
                    placeholder = "comma-separated terms"
                  ),
                  helpText("Terms that identify case/treatment samples (comma-separated)."),

                  numericInput("geoMaxSamples", "Maximum Samples:",
                    value = NULL, min = 10, max = 500, step = 10
                  ),
                  helpText("Limit the number of samples to process (leave blank for all)."),

                  selectInput("geoNormMethod", "Normalization Method:",
                    choices = list(
                      "Quantile (Recommended)" = "quantile",
                      "SWAN" = "SWAN",
                      "Functional Normalization" = "funnorm",
                      "Noob" = "noob"
                    ),
                    selected = "quantile"
                  ),

                  checkboxInput("geoEnableParallel", "Enable Parallel Processing", value = FALSE),

                  conditionalPanel(
                    condition = "input.geoEnableParallel == true",
                    numericInput("geoNCores", "Number of CPU Cores:",
                      value = 2, min = 1, max = 16, step = 1
                    )
                  )
                )
              ),

              # GEO validation results
              conditionalPanel(
                condition = "output.geoValidated",
                hr(),
                fluidRow(
                  column(12,
                    box(
                      title = "GEO Dataset Information",
                      status = "success",
                      solidHeader = TRUE,
                      width = 12,
                      collapsible = TRUE,

                      verbatimTextOutput("geoValidationResult"),

                      br(),

                      conditionalPanel(
                        condition = "output.geoHasIdats",
                        div(class = "success-box",
                          h5(icon("check-circle"), " Dataset Ready for Processing"),
                          p("This dataset contains IDAT files and is suitable for GIMP analysis."),

                          div(class = "warning-box",
                            h5(icon("clock"), " Expected Processing Time:"),
                            p("GEO download and processing can take 10-30 minutes depending on:"),
                            tags$ul(
                              tags$li("Dataset size and number of samples"),
                              tags$li("Internet connection speed"),
                              tags$li("GEO server response time"),
                              tags$li("Computer performance")
                            ),
                            p(strong("Please be patient and don't close the browser during processing."))
                          ),

                          actionButton("processGEO", "Download and Process GEO Dataset",
                            class = "btn-success btn-lg",
                            icon = icon("download"),
                            style = "width: 100%; margin-top: 10px;"
                          )
                        )
                      )
                    )
                  )
                )
              )
            )
          )
        ),

        # Data summary and preview
        fluidRow(
          box(
            title = "Data Summary",
            status = "info",
            solidHeader = TRUE,
            width = 12,

            tabsetPanel(
              tabPanel("Data Summary",
                verbatimTextOutput("dataSummary")
              ),
              tabPanel("Data Preview",
                DT::dataTableOutput("dataPreview")
              ),
              tabPanel("Quality Control",
                conditionalPanel(
                  condition = "input.dataSource == 'idat' || input.dataSource == 'geo'",
                  verbatimTextOutput("qcSummary"),
                  plotOutput("qcPlots", height = "400px")
                ),
                conditionalPanel(
                  condition = "input.dataSource == 'processed'",
                  p("Quality control metrics available for IDAT and GEO data only.")
                )
              ),
              tabPanel("GEO Metadata",
                conditionalPanel(
                  condition = "input.dataSource == 'geo' && output.dataLoaded",
                  verbatimTextOutput("geoMetadata")
                ),
                conditionalPanel(
                  condition = "input.dataSource != 'geo' || !output.dataLoaded",
                  p("GEO metadata available only for GEO datasets.")
                )
              )
            )
          )
        )
      ),

      # CpG Analysis Tab
      tabItem(tabName = "cpg",
        fluidRow(
          box(
            title = "CpG Coverage Analysis",
            status = "primary",
            solidHeader = TRUE,
            width = 12,

            conditionalPanel(
              condition = "output.dataLoaded",

              actionButton("runCpGAnalysis", "Analyse CpG Coverage",
                class = "btn-success",
                icon = icon("calculator")
              ),

              br(), br(),

              tabsetPanel(
                tabPanel("Coverage Counts",
                  plotlyOutput("cpgCountPlot", height = "600px")
                ),
                tabPanel("Coverage Percentage",
                  plotlyOutput("cpgPercentagePlot", height = "600px")
                ),
                tabPanel("Coverage Data",
                  DT::dataTableOutput("cpgCoverageTable")
                )
              )
            ),

            conditionalPanel(
              condition = "!output.dataLoaded",
              div(
                style = "text-align: center; padding: 50px;",
                icon("upload", style = "font-size: 48px; color: #3c8dbc;"),
                h4("No Data Loaded"),
                p("Please upload data first in the Data Upload tab.")
              )
            )
          )
        )
      ),

      # ICR Analysis Tab
      tabItem(tabName = "icr",
        fluidRow(
          box(
            title = "ICR Methylation Heatmap",
            status = "primary",
            solidHeader = TRUE,
            width = 12,

            conditionalPanel(
              condition = "output.dataLoaded",

              fluidRow(
                column(4,
                  selectInput("plotType", "Plot Type:",
                    choices = list(
                      "Beta Values" = "beta",
                      "Delta Beta (vs Control)" = "delta",
                      "Defect Matrix" = "defect"
                    ),
                    selected = "beta"
                  ),
                  helpText("Beta: Raw methylation (0-1), Delta: Difference from controls, Defect: Binary abnormal patterns")
                ),

                column(4,
                  selectInput("orderBy", "Order By:",
                    choices = list(
                      "Coordinates" = "cord",
                      "Methylation Values" = "meth"
                    ),
                    selected = "meth"
                  ),
                  helpText("Coordinates: Genomic position, Methylation: Cluster by similarity")
                ),

                column(4,
                  conditionalPanel(
                    condition = "input.plotType == 'defect'",
                    div(
                      h5("SD Threshold:", style = "margin-bottom: 5px;"),
                      sliderInput("sdThreshold",
                        label = NULL,
                        value = 3.0,
                        min = 1.0,
                        max = 5.0,
                        step = 0.1,
                        ticks = TRUE,
                        animate = FALSE
                      ),
                      helpText("Lower = more sensitive detection", style = "margin-top: -10px; font-size: 11px;"),
                      div(
                        style = "font-size: 10px; color: #666; margin-top: 5px;",
                        "1.5-2.0: Very sensitive | 2.5-3.0: Standard | 3.5+: Conservative"
                      )
                    )
                  )
                )
              ),

              actionButton("generateHeatmap", "Generate Heatmap",
                class = "btn-success",
                icon = icon("fire")
              ),

              br(), br(),

              plotOutput("icrHeatmap", height = "800px"),

              br(),

              downloadButton("downloadHeatmap", "Download Heatmap",
                class = "btn-info")
            ),

            conditionalPanel(
              condition = "!output.dataLoaded",
              div(
                style = "text-align: center; padding: 50px;",
                icon("upload", style = "font-size: 48px; color: #3c8dbc;"),
                h4("No Data Loaded"),
                p("Please upload data first in the Data Upload tab.")
              )
            )
          )
        )
      ),

      # Differential Analysis Tab
      tabItem(tabName = "dmps",
        fluidRow(
          box(
            title = "Differential Methylation Analysis",
            status = "primary",
            solidHeader = TRUE,
            width = 12,

            conditionalPanel(
              condition = "output.dataLoaded",

              fluidRow(
                column(6,
                  numericInput("pValueCutoff", "P-value Cutoff:",
                    value = 0.05, min = 0.001, max = 0.1, step = 0.01
                  )
                ),

                column(6,
                  actionButton("runDMPAnalysis", "Run DMP Analysis",
                    class = "btn-success",
                    icon = icon("calculator")
                  )
                )
              ),

              br(),

              tabsetPanel(
                tabPanel("Significant DMPs",
                  DT::dataTableOutput("dmpsTable")
                ),
                tabPanel("Volcano Plot",
                  plotlyOutput("volcanoPlot", height = "600px")
                ),
                tabPanel("Summary Statistics",
                  verbatimTextOutput("dmpsSummary")
                )
              ),

              br(),

              downloadButton("downloadDMPs", "Download DMP Results",
                class = "btn-info")
            ),

            conditionalPanel(
              condition = "!output.dataLoaded",
              div(
                style = "text-align: center; padding: 50px;",
                icon("upload", style = "font-size: 48px; color: #3c8dbc;"),
                h4("No Data Loaded"),
                p("Please upload data first in the Data Upload tab.")
              )
            )
          )
        )
      ),

      # Region Explorer Tab
      tabItem(tabName = "explorer",
        fluidRow(
          box(
            title = "ICR Region Explorer",
            status = "primary",
            solidHeader = TRUE,
            width = 12,

            conditionalPanel(
              condition = "output.dmpsDone",

              fluidRow(
                column(4,
                  selectInput("selectedICR", "Select ICR:",
                    choices = NULL,
                    selected = NULL
                  ),
                  helpText("Only ICRs with significant DMPs are shown.")
                ),

                column(4,
                  radioButtons("plotInteractive", "Plot Type:",
                    choices = list(
                      "Interactive (plotly)" = TRUE,
                      "Static (ggplot2)" = FALSE
                    ),
                    selected = TRUE
                  ),
                  helpText("Interactive plots allow zooming and hovering.")
                ),

                column(4,
                  br(),
                  actionButton("plotRegion", "Plot Region",
                    class = "btn-success",
                    icon = icon("chart-line"),
                    style = "margin-top: 5px;"
                  ),
                  br(), br(),
                  actionButton("refreshICRList", "Refresh ICR List",
                    class = "btn-info btn-sm",
                    icon = icon("refresh")
                  )
                )
              ),

              hr(),

              # Dynamic plot output area
              uiOutput("regionPlotUI"),

              # Add some helpful information
              br(),
              div(
                style = "background-color: #f8f9fa; padding: 15px; border-radius: 5px;",
                h5("Plot Information:"),
                tags$ul(
                  tags$li("Blue lines/points: Control samples"),
                  tags$li("Red lines/points: Case samples"),
                  tags$li("Vertical dashed lines: Significant DMP positions"),
                  tags$li("Red rug marks: Significant DMPs"),
                  tags$li("Larger points: CpGs that are significant DMPs")
                )
              )
            ),

            conditionalPanel(
              condition = "!output.dmpsDone",
              div(
                style = "text-align: center; padding: 50px;",
                icon("exclamation-triangle", style = "font-size: 48px; color: #f39c12;"),
                h4("DMP Analysis Required"),
                p("Please run DMP analysis first in the Differential Analysis tab."),
                p("The Region Explorer shows detailed methylation patterns for significant DMPs.")
              )
            )
          )
        )
      ),

      # Help Tab
      tabItem(tabName = "help",
        fluidRow(
          box(
            title = "GIMP Shiny App Help",
            status = "info",
            solidHeader = TRUE,
            width = 12,

            h3("GIMP: Genomic Imprinting Methylation Patterns"),
            p("GIMP is a specialized tool for analyzing methylation patterns at Imprinting Control Regions (ICRs)."),

            h4("Data Sources Supported:"),
            tags$ul(
              tags$li(strong("Processed Data:"), " Upload CSV/RDS/Excel files with beta values"),
              tags$li(strong("Raw IDAT Files:"), " Upload ZIP archives containing IDAT files and sample sheets"),
              tags$li(strong("GEO Datasets:"), " Automatically download and process data from NCBI GEO")
            ),

            h4("GEO Dataset Processing:"),
            p("GIMP can automatically download and process methylation datasets from the Gene Expression Omnibus (GEO):"),
            tags$ol(
              tags$li("Enter a GEO Series ID (e.g., GSE68777)"),
              tags$li("Validate the dataset to check for IDAT files"),
              tags$li("Configure processing options (group detection, normalization)"),
              tags$li("Download and process automatically")
            ),

            div(class = "info-box",
              h5(icon("info-circle"), " GEO Processing Notes:"),
              tags$ul(
                tags$li("Only datasets with raw IDAT files are supported"),
                tags$li("Group auto-detection looks for terms like 'control', 'case', 'tumor', etc."),
                tags$li("Processing time depends on dataset size (typically 10-30 minutes)"),
                tags$li("Requires internet connection and sufficient disk space")
              )
            ),

            h4("Analysis Workflow:"),
            tags$ol(
              tags$li(strong("Data Upload:"), " Choose your data source and upload/process"),
              tags$li(strong("CpG Coverage:"), " Analyze probe coverage across ICRs"),
              tags$li(strong("ICR Heatmaps:"), " Visualize methylation patterns with multiple plot types"),
              tags$li(strong("Differential Analysis:"), " Identify significantly different methylation positions"),
              tags$li(strong("Region Explorer:"), " Detailed visualization of specific ICRs")
            ),

            h4("Troubleshooting:"),
            div(class = "warning-box",
              h5("Common Issues:"),
              tags$ul(
                tags$li(strong("GEO validation fails:"), " Check internet connection and GEO ID format"),
                tags$li(strong("No IDAT files found:"), " Dataset may only have processed data"),
                tags$li(strong("Processing timeout:"), " Try smaller datasets or increase browser timeout"),
                tags$li(strong("Memory errors:"), " Close other browser tabs and try again")
              )
            ),

            h4("About ICRs:"),
            p("Imprinting Control Regions (ICRs) are genomic regions that regulate parent-of-origin-specific gene expression. GIMP uses curated ICR coordinates from published literature to focus analysis on these biologically important regions."),

            h4("Citations:"),
            p("If you use GIMP in your research, please cite:"),
            tags$ul(
              tags$li("Cecere, F. (2024). GIMP: Genomic Imprinting Methylation Patterns. R package."),
              tags$li("ICR coordinates from: Joshi et al. (2016) Epigenetics, DOI: 10.1080/15592294.2016.1264561")
            )
          )
        )
      )
    )
  )
)

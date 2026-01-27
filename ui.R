dashboardPage(
  dashboardHeader(title = "NeRNA - Negative RNA Generator"),
  
  dashboardSidebar(
    sidebarMenu(id = "sidebar",
      menuItem("About", tabName = "about", icon = icon("info-circle"), selected = TRUE),
      menuItem("Upload FASTA", tabName = "upload", icon = icon("upload")),
      menuItem("Overview", tabName = "overview", icon = icon("tachometer-alt")),
      menuItem("Secondary Structures", tabName = "structures", icon = icon("dna")),
      menuItem("Negative Generation", tabName = "negative", icon = icon("flask")),
      menuItem("Overview Comparison", tabName = "overview_comparison", icon = icon("chart-pie")),
      menuItem("Result Comparison", tabName = "comparison", icon = icon("chart-bar"))
    )
  ),
  
  dashboardBody(
    tabItems(
      # About Tab
      tabItem(tabName = "about",
        fluidRow(
          box(
            title = "Welcome to NeRNA - Negative RNA Generator Tool",
            status = "primary",
            solidHeader = TRUE,
            width = 12,
            tags$div(
              style = "padding: 15px;",
              tags$h3("About NeRNA"),
              tags$p("NeRNA is a comprehensive tool for generating negative RNA datasets for RNA sequence analysis and research. 
                     The tool provides three different methods to generate negative RNA sequences, allowing researchers to 
                     choose the most appropriate method for their specific research needs."),
              
              tags$h4("Key Features:"),
              tags$ul(
                tags$li(tags$strong("Multiple Negative Generation Methods:"), 
                       " NeRNA offers three different ways to generate negative RNA sequences:"),
                tags$ul(
                  tags$li(tags$strong("NeRNA Method:"), " A novel approach that shifts the sequence by a specified amount, preserving many structural features while providing effective negative sequences."),
                  tags$li(tags$strong("Dinucleotide Shuffling:"), " Preserves dinucleotide frequencies, suitable for studies where composition is important."),
                  tags$li(tags$strong("Random Shuffling:"), " Basic randomization method that completely disrupts the sequence order.")
                ),
                tags$li(tags$strong("Secondary Structure Analysis:"), " Visualize and analyze RNA secondary structures using RNAfold."),
                tags$li(tags$strong("Comprehensive Comparisons:"), " Compare different methods using metrics like MFE, sequence length, base pairing, and GC content."),
                tags$li(tags$strong("Exportable Results:"), " Download results in multiple formats (CSV, Excel, FASTA) for further analysis.")
              ),
              
              tags$h4("How to Use NeRNA:"),
              tags$ol(
                tags$li("Start by uploading a FASTA file containing your RNA sequences."),
                tags$li("Explore your dataset in the Overview tab to understand basic statistics."),
                tags$li("Generate secondary structures to analyze the original RNA sequences."),
                tags$li("Use the Negative Generation tab to create negative sequences with your preferred method."),
                tags$li("Compare the results in the Overview Comparison and Result Comparison tabs."),
                tags$li("Download your results for publication or further analysis.")
              ),
              
              tags$h4("Citation:"),
              tags$p("If you use NeRNA in your research, please cite: [Citation information]"),
              
              tags$div(
                style = "margin-top: 20px; text-align: center;",
                actionButton("getStartedBtn", "Get Started", 
                            class = "btn-lg btn-success", 
                            style = "margin-right: 10px;"),
                actionButton("githubBtn", "View on GitHub", 
                            class = "btn-lg btn-info", 
                            onclick = "window.open('https://github.com/yourusername/nerna', '_blank')")
              )
            )
          )
        ),
        fluidRow(
          box(
            title = "Navigation Guide",
            status = "info",
            solidHeader = TRUE,
            width = 12,
            tags$div(
              style = "display: flex; flex-wrap: wrap; justify-content: space-between;",
              
              tags$div(
                style = "width: 30%; min-width: 250px; margin-bottom: 20px;",
                tags$h4(icon("upload"), " Upload FASTA"),
                tags$p("Start here to upload your RNA sequence FASTA files. Choose the RNA type for optimal analysis.")
              ),
              
              tags$div(
                style = "width: 30%; min-width: 250px; margin-bottom: 20px;",
                tags$h4(icon("tachometer-alt"), " Overview"),
                tags$p("Get a quick overview of your uploaded sequences, including length distribution and GC content.")
              ),
              
              tags$div(
                style = "width: 30%; min-width: 250px; margin-bottom: 20px;",
                tags$h4(icon("dna"), " Secondary Structures"),
                tags$p("Generate and visualize RNA secondary structures using RNAfold algorithm.")
              ),
              
              tags$div(
                style = "width: 30%; min-width: 250px; margin-bottom: 20px;",
                tags$h4(icon("flask"), " Negative Generation"),
                tags$p("Generate negative RNA sequences using NeRNA, Dinucleotide Shuffling, or Random Shuffling methods.")
              ),
              
              tags$div(
                style = "width: 30%; min-width: 250px; margin-bottom: 20px;",
                tags$h4(icon("chart-pie"), " Overview Comparison"),
                tags$p("Compare the statistical properties of original and generated sequences across all methods.")
              ),
              
              tags$div(
                style = "width: 30%; min-width: 250px; margin-bottom: 20px;",
                tags$h4(icon("chart-bar"), " Result Comparison"),
                tags$p("Detailed comparison of individual sequences, including secondary structure visualization.")
              )
            )
          )
        )
      ),
      
      # Upload Tab
      tabItem(tabName = "upload",
        fluidRow(
          box(
            title = "Upload FASTA File",
            status = "primary",
            solidHeader = TRUE,
            width = 12,
            tags$div(
              style = "border-left: 4px solid #17a2b8; padding: 10px; margin-bottom: 15px; background-color: #f8f9fa;",
              tags$h4(style = "margin-top: 0;", "Getting Started"),
              tags$p("Upload your FASTA format RNA sequence file here. Your file should contain RNA sequences in standard FASTA format. Choose the appropriate RNA type for optimal analysis. The tool can process various RNA types including miRNA, circRNA, lncRNA and tRNA."),
              tags$div(
                style = "background-color: #fff3cd; border-left: 4px solid #ffc107; padding: 8px; margin-top: 10px;",
                tags$strong(icon("exclamation-triangle"), " Important Limits:"),
                tags$ul(style = "margin-bottom: 0; margin-top: 5px;",
                  tags$li(tags$strong("Maximum sequences per file: 250"), " - If your file contains more than 250 sequences, it will be automatically split into batches for download."),
                  tags$li(tags$strong("Maximum sequence length: 1000 nucleotides"), " - Longer sequences will be filtered out for system stability.")
                )
              )
            ),
            fileInput("fastaFile", "Choose FASTA File",
                      accept = c(".fasta", ".fa", ".txt")),
            selectInput("rnaType", "RNA Type:",
                        choices = list(
                          "miRNA" = "mirna",
                          "circRNA" = "circrna",
                          "lncRNA" = "lncrna",
                          "tRNA" = "trna",
                          "Other" = "other"
                        ),
                        selected = "mirna"),
            actionButton("loadBtn", "Load Sequences", 
                        class = "btn-primary"),
            tags$div(
              style = "margin-top: 15px;",
              conditionalPanel(
                condition = "output.batchesReady == true",
                tags$div(
                  style = "background-color: #f8d7da; border: 2px solid #f5c6cb; border-radius: 5px; padding: 15px; margin-bottom: 10px;",
                  tags$h4(style = "color: #721c24; margin-top: 0;", 
                          icon("exclamation-circle"), " Too Many Sequences Detected!"),
                  tags$p(style = "color: #721c24;", 
                         "Your file contains more than 250 sequences. The sequences have been automatically split into batches."),
                  tags$p(style = "color: #721c24; font-weight: bold;", 
                         "Click the button below to download all batches as a ZIP file:")
                ),
                downloadButton("downloadBatchesZip", 
                              "Download Batches (ZIP)",
                              class = "btn-danger btn-lg",
                              style = "width: 100%; margin-bottom: 10px;"),
                tags$div(
                  style = "background-color: #d1ecf1; border-left: 4px solid #0c5460; padding: 10px; margin-top: 10px;",
                  tags$strong(icon("info-circle"), " Instructions:"),
                  tags$ol(style = "margin-bottom: 0; margin-top: 5px;",
                    tags$li("Download the ZIP file"),
                    tags$li("Extract all FASTA files"),
                    tags$li("Upload each batch separately"),
                    tags$li("Process each batch individually"),
                    tags$li("Combine results for final analysis")
                  )
                )
              )
            )
          ),
          box(
            title = "Parallel Processing Settings",
            status = "info",
            solidHeader = TRUE,
            width = 6,
            collapsible = TRUE,
            collapsed = TRUE,
            tags$div(
              style = "border-left: 4px solid #17a2b8; padding: 10px; margin-bottom: 15px; background-color: #f8f9fa;",
              tags$h4(style = "margin-top: 0;", icon("microchip"), " Multi-Core Processing"),
              tags$p("Enable parallel processing to speed up computations. Automatically uses 75% of available CPU cores.")
            ),
            checkboxInput("enableParallel", 
                         "Enable Parallel Processing (for 10+ sequences)", 
                         value = TRUE),
            uiOutput("parallelInfoUI"),
            actionButton("refreshParallelInfo", "Refresh CPU Info", 
                        class = "btn-info btn-sm", 
                        style = "width: 100%;")
          ),
          box(
            title = "Session Management",
            status = "success",
            solidHeader = TRUE,
            width = 6,
            collapsible = TRUE,
            collapsed = TRUE,
            tags$div(
              style = "border-left: 4px solid #28a745; padding: 10px; margin-bottom: 15px; background-color: #f8f9fa;",
              tags$h4(style = "margin-top: 0;", icon("save"), " Save & Load Sessions"),
              tags$p("Save your current analysis session to continue later. All sequences, structures, and negative results will be preserved.")
            ),
            downloadButton("saveSession", "Save Current Session (.rds)", 
                          class = "btn-success", icon = icon("save"),
                          style = "width: 100%; margin-bottom: 10px;"),
            fileInput("loadSessionFile", "Load Session File (.rds)",
                     accept = c(".rds", ".RDS")),
            actionButton("loadSessionBtn", "Load Session", 
                        class = "btn-success", icon = icon("folder-open"),
                        style = "width: 100%;")
          ),
          box(
            title = "Batch Processing",
            status = "warning",
            solidHeader = TRUE,
            width = 6,
            collapsible = TRUE,
            collapsed = TRUE,
            tags$div(
              style = "border-left: 4px solid #ffc107; padding: 10px; margin-bottom: 15px; background-color: #f8f9fa;",
              tags$h4(style = "margin-top: 0;", icon("layer-group"), " Batch FASTA Upload"),
              tags$p("Upload multiple FASTA files at once for batch processing. All files will be combined and processed together.")
            ),
            fileInput("batchFastaFiles", "Choose Multiple FASTA Files",
                     accept = c(".fasta", ".fa", ".txt"),
                     multiple = TRUE),
            actionButton("loadBatchBtn", "Load Batch Files", 
                        class = "btn-warning", icon = icon("upload"),
                        style = "width: 100%;")
          )
        )
      ),

      # Overview Tab
      tabItem(tabName = "overview",
        conditionalPanel(
          condition = "output.dataLoaded != true",
          fluidRow(
            box(
              title = "Data Required",
              status = "warning",
              solidHeader = TRUE,
              width = 12,
              tags$div(
                style = "text-align: center; padding: 20px;",
                tags$i(class = "fa fa-exclamation-circle fa-5x", style = "color: #f39c12; margin-bottom: 20px;"),
                tags$h3("No Data Loaded"),
                tags$p("Please upload a FASTA file and click 'Load Sequences' in the Upload tab to proceed."),
                tags$br(),
                actionButton("goToUploadBtn", "Go to Upload", class = "btn-primary")
              )
            )
          )
        ),
        conditionalPanel(
          condition = "output.dataLoaded == true",
          fluidRow(
            box(
              title = "Overview Information",
              status = "success",
              solidHeader = TRUE,
              width = 12,
              tags$div(
                style = "border-left: 4px solid #28a745; padding: 10px; margin-bottom: 15px; background-color: #f8f9fa;",
                tags$h4(style = "margin-top: 0;", "Data Overview"),
                tags$p("This tab provides a summary of your uploaded RNA sequences. You can view sequence details, length distribution, GC content, and other basic statistics. Sequences longer than 1,000 nucleotides are filtered out automatically for system stability and are listed separately below.")
              ),
              collapsible = TRUE
            )
          ),
          fluidRow(
            box(
              title = "Sequence Overview",
              status = "info",
              solidHeader = TRUE,
              width = 12,
              collapsible = TRUE,
              DTOutput("sequenceTable")
            )
          ),
          fluidRow(
            box(
              title = "Filtered Sequences (>1,000 nucleotides)",
              status = "warning",
              solidHeader = TRUE,
              width = 12,
              collapsible = TRUE,
              collapsed = FALSE,
              uiOutput("filteredSequencesInfo")
            )
          ),
          fluidRow(
            valueBoxOutput("seqCountBox", width = 4),
            valueBoxOutput("avgLengthBox", width = 4),
            valueBoxOutput("gcContentBox", width = 4)
          )
        )
      ),
      
      # Secondary Structures Tab
      tabItem(tabName = "structures",
        conditionalPanel(
          condition = "output.dataLoaded != true",
          fluidRow(
            box(
              title = "Data Required",
              status = "warning",
              solidHeader = TRUE,
              width = 12,
              tags$div(
                style = "text-align: center; padding: 20px;",
                tags$i(class = "fa fa-exclamation-circle fa-5x", style = "color: #f39c12; margin-bottom: 20px;"),
                tags$h3("No Data Loaded"),
                tags$p("Please upload a FASTA file and click 'Load Sequences' in the Upload tab to proceed."),
                tags$br(),
                actionButton("goToUploadBtn2", "Go to Upload", class = "btn-primary")
              )
            )
          )
        ),
        conditionalPanel(
          condition = "output.dataLoaded == true",
          fluidRow(
            box(
              title = "RNA Secondary Structures",
              status = "warning",
              solidHeader = TRUE,
              width = 12,
              tags$div(
                style = "border-left: 4px solid #ffc107; padding: 10px; margin-bottom: 15px; background-color: #f8f9fa;",
                tags$h4(style = "margin-top: 0;", "Secondary Structure Analysis"),
                tags$p("This tool uses RNAfold to predict the secondary structure of your RNA sequences. The results include the minimum free energy (MFE) value and the predicted secondary structure in dot-bracket notation. Click the 'Generate Secondary Structures' button to begin processing. Results can be downloaded in CSV, Excel, or FASTA formats.")
              ),
              actionButton("runRNAfoldBtn", "Generate Secondary Structures",
                          class = "btn-warning"),
              br(), br(),
              # Add download buttons for CSV and Excel
              fluidRow(
                column(3, 
                  downloadButton("downloadCSV", "Download CSV", 
                                class = "btn-info"),
                  tags$style(type = "text/css", "#downloadCSV {margin-right: 15px;}")
                ),
                column(3,
                  downloadButton("downloadExcel", "Download Excel", 
                                class = "btn-success")
                ),
                column(3,
                  downloadButton("downloadFASTA", "Download FASTA", 
                                class = "btn-primary")
                )
              ),
              br(),
              DTOutput("rnafoldTable")
            )
          ),
          fluidRow(
            box(
              title = "Real-time Structure Viewer",
              status = "info",
              solidHeader = TRUE,
              width = 12,
              collapsible = TRUE,
              collapsed = FALSE,
              tags$div(
                style = "border-left: 4px solid #17a2b8; padding: 10px; margin-bottom: 15px; background-color: #f8f9fa;",
                tags$h4(style = "margin-top: 0;", icon("eye"), " Interactive Structure Viewer"),
                tags$p("Select a sequence from the dropdown to view its secondary structure in detail. The viewer shows the sequence, structure notation, base pairs, and structural statistics.")
              ),
              selectInput("structureViewerSeq", "Select Sequence to View:",
                         choices = NULL,
                         width = "100%"),
              hr(),
              tabsetPanel(
                tabPanel("3D Arc Diagram",
                  tags$div(style = "margin-top: 15px;"),
                  plotlyOutput("structure3DPlot", height = "500px"),
                  tags$div(
                    style = "background-color: #e7f3ff; padding: 10px; border-radius: 5px; margin-top: 15px;",
                    tags$h5(icon("info-circle"), " Arc Diagram Guide:"),
                    tags$p(tags$strong("🔴 Red:"), " Adenine (A)"),
                    tags$p(tags$strong("🔵 Blue:"), " Uracil (U)"),
                    tags$p(tags$strong("🟢 Green:"), " Guanine (G)"),
                    tags$p(tags$strong("🟠 Orange:"), " Cytosine (C)"),
                    tags$p(tags$strong("Blue arcs:"), " Base pairing (hydrogen bonds)"),
                    tags$p("Hover over bases or arcs to see details. Use mouse to zoom and pan.")
                  )
                ),
                tabPanel("Text View",
                  tags$div(style = "margin-top: 15px;"),
                  verbatimTextOutput("structureViewer", placeholder = TRUE),
                  tags$div(
                    style = "background-color: #e7f3ff; padding: 10px; border-radius: 5px; margin-top: 15px;",
                    tags$h5(icon("info-circle"), " Text Legend:"),
                    tags$p(tags$strong("( and )"), " = Paired bases (forming stems)"),
                    tags$p(tags$strong("."), " = Unpaired bases (loops, bulges)"),
                    tags$p(tags$strong("MFE"), " = Minimum Free Energy (lower = more stable)")
                  )
                )
              )
            )
          )
        )
      ),
      
      # Negative Generation Tab (Combined)
      tabItem(tabName = "negative",
        conditionalPanel(
          condition = "output.dataLoaded != true",
          fluidRow(
            box(
              title = "Data Required",
              status = "warning",
              solidHeader = TRUE,
              width = 12,
              tags$div(
                style = "text-align: center; padding: 20px;",
                tags$i(class = "fa fa-exclamation-circle fa-5x", style = "color: #f39c12; margin-bottom: 20px;"),
                tags$h3("No Data Loaded"),
                tags$p("Please upload a FASTA file and click 'Load Sequences' in the Upload tab to proceed."),
                tags$br(),
                actionButton("goToUploadBtn3", "Go to Upload", class = "btn-primary")
              )
            )
          )
        ),
        conditionalPanel(
          condition = "output.dataLoaded == true",
          fluidRow(
            column(width = 3,
              box(
                title = "Methods",
                status = "primary",
                solidHeader = TRUE,
                width = 12,
                tags$div(
                  style = "border-left: 4px solid #007bff; padding: 10px; margin-bottom: 15px; background-color: #f8f9fa;",
                  tags$h4(style = "margin-top: 0;", "Negative RNA Generation"),
                  tags$p("Generate negative RNA sequences using three different methods. Each method has different characteristics and may be suitable for different research questions.")
                ),
                radioButtons("negativeMethod", "Select Method:",
                            choices = list(
                              "NeRNA" = "nerna",
                              "Dinucleotide Shuffling" = "dinucleotide",
                              "Random Shuffling" = "random"
                            ),
                            selected = "nerna")
              ),
              box(
                title = "Batch Operations",
                status = "info",
                solidHeader = TRUE,
                width = 12,
                actionButton("generateAllBtn", "Generate All Methods", 
                            class = "btn-warning", 
                            style = "width: 100%; margin-bottom: 10px;"),
                downloadButton("downloadAllCSV", "Download All (CSV)", 
                            class = "btn-info", 
                            style = "width: 100%; margin-bottom: 10px;"),
                downloadButton("downloadAllExcel", "Download All (Excel)", 
                            class = "btn-success", 
                            style = "width: 100%; margin-bottom: 10px;"),
                downloadButton("downloadBatchFASTA", "Download All (FASTA)", 
                            class = "btn-primary", 
                            style = "width: 100%; margin-bottom: 10px;"),
                actionButton("showOverviewComparisonBtn", "Result Comparison", 
                            class = "btn-primary", 
                            style = "width: 100%;")
              )
            ),
            column(width = 9,
              # NeRNA Panel - Show when NeRNA is selected
              conditionalPanel(
                condition = "input.negativeMethod == 'nerna'",
                box(
                  title = "NeRNA Method",
                  status = "success",
                  solidHeader = TRUE,
                  width = 12,
                  tags$div(
                    style = "border-left: 4px solid #28a745; padding: 10px; margin-bottom: 15px; background-color: #f8f9fa;",
                    tags$h4(style = "margin-top: 0;", "NeRNA Algorithm"),
                    tags$p("Structural transformation and octal shifting method for generating negative RNA sequences. Uses secondary structure information to preserve biological relevance while creating effective negative controls."),
                    tags$div(
                      style = "background-color: #fff3cd; border: 1px solid #ffeaa7; padding: 10px; border-radius: 4px; margin-top: 10px;",
                      tags$strong("⚠️ Important: "),
                      "You must run secondary structure analysis first (in the 'Secondary Structures' tab) before using NeRNA method."
                    )
                  ),
                  fluidRow(
                    column(width = 6,
                      tags$label("Shift Amount:"),
                      sliderInput("nernaShiftingSize", 
                                label = NULL,
                                min = 1, 
                                max = 72, 
                                value = 1,
                                width = "100%"),
                      tags$p("Number of bits to shift in the octal representation", 
                            style = "font-size: 90%; color: #6c757d;")
                    ),
                    column(width = 6,
                      numericInput("nernaSeed", "Random Seed:", 
                                  value = 123, min = 1, step = 1),
                      tags$p("For reproducible results", 
                            style = "font-size: 90%; color: #6c757d;")
                    )
                  ),
                  br(),
                  actionButton("runNeRNABtn", "Generate with NeRNA",
                              class = "btn-success"),
                  br(), br(),
                  # Add download buttons for NeRNA
                  fluidRow(
                    column(3, 
                      downloadButton("downloadNeRNACSV", "Download CSV", 
                                  class = "btn-info"),
                      tags$style(type = "text/css", "#downloadNeRNACSV {margin-right: 15px;}")
                    ),
                    column(3,
                      downloadButton("downloadNeRNAExcel", "Download Excel", 
                                  class = "btn-success")
                    ),
                    column(3,
                      downloadButton("downloadNeRNAFASTA", "Download FASTA", 
                                  class = "btn-primary")
                    )
                  ),
                  br(),
                  DTOutput("nernaOutput")
                )
              ),
              
              # Dinucleotide Shuffling Panel - Show when Dinucleotide is selected
              conditionalPanel(
                condition = "input.negativeMethod == 'dinucleotide'",
                box(
                  title = "Dinucleotide Shuffling Method",
                  status = "primary",
                  solidHeader = TRUE,
                  width = 12,
                  fluidRow(
                    column(width = 6,
                      numericInput("dinucleotideSeed", "Random Seed:", 
                                  value = 123, min = 1, step = 1)
                    )
                  ),
                  br(),
                  actionButton("runDinucleotideBtn", "Generate with Dinucleotide Shuffling",
                              class = "btn-primary"),
                  br(), br(),
                  # Add download buttons for Dinucleotide
                  fluidRow(
                    column(3, 
                      downloadButton("downloadDinucleotideCSV", "Download CSV", 
                                  class = "btn-info"),
                      tags$style(type = "text/css", "#downloadDinucleotideCSV {margin-right: 15px;}")
                    ),
                    column(3,
                      downloadButton("downloadDinucleotideExcel", "Download Excel", 
                                  class = "btn-success")
                    ),
                    column(3,
                      downloadButton("downloadDinucleotideFASTA", "Download FASTA", 
                                  class = "btn-primary")
                    )
                  ),
                  br(),
                  DTOutput("dinucleotideOutput")
                )
              ),
              
              # Random Shuffling Panel - Show when Random is selected
              conditionalPanel(
                condition = "input.negativeMethod == 'random'",
                box(
                  title = "Random Shuffling Method",
                  status = "danger",
                  solidHeader = TRUE,
                  width = 12,
                  fluidRow(
                    column(width = 6,
                      numericInput("randomSeed", "Random Seed:", 
                                  value = 123, min = 1, step = 1)
                    )
                  ),
                  br(),
                  actionButton("runRandomBtn", "Generate with Random Shuffling",
                              class = "btn-danger"),
                  br(), br(),
                  # Add download buttons for Random
                  fluidRow(
                    column(3, 
                      downloadButton("downloadRandomCSV", "Download CSV", 
                                  class = "btn-info"),
                      tags$style(type = "text/css", "#downloadRandomCSV {margin-right: 15px;}")
                    ),
                    column(3,
                      downloadButton("downloadRandomExcel", "Download Excel", 
                                  class = "btn-success")
                    ),
                    column(3,
                      downloadButton("downloadRandomFASTA", "Download FASTA", 
                                  class = "btn-primary")
                    )
                  ),
                  br(),
                  DTOutput("randomOutput")
                )
              )
            )
          )
        )
      ),
      
      # Comparison Tab
      tabItem(tabName = "comparison",
        conditionalPanel(
          condition = "output.dataLoaded != true",
          fluidRow(
            box(
              title = "Data Required",
              status = "warning",
              solidHeader = TRUE,
              width = 12,
              tags$div(
                style = "text-align: center; padding: 20px;",
                tags$i(class = "fa fa-exclamation-circle fa-5x", style = "color: #f39c12; margin-bottom: 20px;"),
                tags$h3("No Data Loaded"),
                tags$p("Please upload a FASTA file and click 'Load Sequences' in the Upload tab to proceed."),
                tags$br(),
                actionButton("goToUploadBtn4", "Go to Upload", class = "btn-primary")
              )
            )
          )
        ),
        conditionalPanel(
          condition = "output.dataLoaded == true",
          fluidRow(
            box(
              title = "Sequence Selection",
              status = "primary",
              solidHeader = TRUE,
              width = 12,
              tags$div(
                style = "border-left: 4px solid #17a2b8; padding: 10px; margin-bottom: 15px; background-color: #f8f9fa;",
                tags$h4(style = "margin-top: 0;", "Sequence-Specific Comparison"),
                tags$p("This page allows you to compare a specific sequence across all negative generation methods. Select a sequence from the dropdown menu and click 'Compare' to visualize and analyze the differences between the original and generated sequences.")
              ),
              fluidRow(
                column(width = 5,
                  selectInput("comparisonSequence", "Select Sequence for Detailed Comparison:",
                              choices = NULL)
                ),
                column(width = 3, 
                  actionButton("compareBtn", "Compare Selected Sequence", 
                              class = "btn-primary",
                              style = "margin-top: 25px; width: 100%;")
                ),
                column(width = 4,
                  downloadButton("downloadSequenceSpecificExcel", "Download Sequence-Specific Excel", 
                              class = "btn-success",
                              style = "margin-top: 25px; width: 100%;")
                )
              )
            )
          ),
          
          fluidRow(
            box(
              title = "Sequence-Specific Comparison",
              status = "success",
              solidHeader = TRUE,
              width = 12,
              uiOutput("selectedSequenceComparison")
            )
          ),
          
          fluidRow(
            box(
              title = "3D Structure Comparison: Original vs All Methods",
              status = "primary",
              solidHeader = TRUE,
              width = 12,
              collapsible = TRUE,
              collapsed = FALSE,
              tags$div(
                style = "border-left: 4px solid #007bff; padding: 10px; margin-bottom: 15px; background-color: #f8f9fa;",
                tags$h4(style = "margin-top: 0;", icon("cubes"), " 4-Way Structure Comparison"),
                tags$p("Compare the secondary structures of original sequence with all three negative generation methods. This visualization helps you understand how each method affects RNA folding patterns and thermodynamic stability.")
              ),
              fluidRow(
                column(6,
                  tags$div(
                    style = "background-color: #e3f2fd; padding: 10px; border-radius: 5px; margin-bottom: 10px;",
                    tags$h5(style = "margin: 0; color: #1976d2;", icon("dna"), " Original Sequence")
                  ),
                  plotlyOutput("comp4OriginalPlot", height = "350px"),
                  verbatimTextOutput("comp4OriginalText")
                ),
                column(6,
                  tags$div(
                    style = "background-color: #fff3e0; padding: 10px; border-radius: 5px; margin-bottom: 10px;",
                    tags$h5(style = "margin: 0; color: #f57c00;", icon("sync"), " NeRNA Method")
                  ),
                  plotlyOutput("comp4NeRNAPlot", height = "350px"),
                  verbatimTextOutput("comp4NeRNAText")
                )
              ),
              hr(),
              fluidRow(
                column(6,
                  tags$div(
                    style = "background-color: #e8f5e9; padding: 10px; border-radius: 5px; margin-bottom: 10px;",
                    tags$h5(style = "margin: 0; color: #388e3c;", icon("random"), " Dinucleotide Shuffling")
                  ),
                  plotlyOutput("comp4DinucPlot", height = "350px"),
                  verbatimTextOutput("comp4DinucText")
                ),
                column(6,
                  tags$div(
                    style = "background-color: #fce4ec; padding: 10px; border-radius: 5px; margin-bottom: 10px;",
                    tags$h5(style = "margin: 0; color: #c2185b;", icon("random"), " Random Shuffling")
                  ),
                  plotlyOutput("comp4RandomPlot", height = "350px"),
                  verbatimTextOutput("comp4RandomText")
                )
              ),
              hr(),
              tags$div(
                style = "background-color: #f3e5f5; padding: 15px; border-radius: 5px;",
                tags$h5(icon("chart-bar"), " Structure Comparison Summary"),
                tableOutput("structureComparisonTable")
              )
            )
          ),
          
          fluidRow(
            box(
              title = "Minimum Free Energy (MFE) Comparison",
              status = "warning",
              solidHeader = TRUE,
              width = 6,
              plotOutput("mfeVisualization", height = "300px")
            ),
            box(
              title = "Sequence Length Comparison",
              status = "primary",
              solidHeader = TRUE,
              width = 6,
              plotOutput("lengthVisualization", height = "300px")
            )
          ),
          
          fluidRow(
            box(
              title = "Base Pairing Comparison",
              status = "info",
              solidHeader = TRUE,
              width = 6,
              plotOutput("pairingVisualization", height = "300px")
            ),
            box(
              title = "GC Content Comparison",
              status = "success",
              solidHeader = TRUE,
              width = 6,
              plotOutput("gcVisualization", height = "300px")
            )
          ),
          
          fluidRow(
            box(
              title = "Sequence Feature Overview",
              status = "danger",
              solidHeader = TRUE,
              width = 6,
              plotOutput("radarVisualization", height = "400px")
            ),
            box(
              title = "Spider Chart Comparison",
              status = "primary",
              solidHeader = TRUE,
              width = 6,
              plotOutput("spiderVisualization", height = "400px")
            )
          ),
          
          # Export Section
          fluidRow(
            box(
              title = "Export Analysis Results",
              status = "primary",
              solidHeader = TRUE,
              width = 12,
              tags$div(
                style = "border-left: 4px solid #007bff; padding: 10px; margin-bottom: 15px; background-color: #f8f9fa;",
                tags$h4(style = "margin-top: 0;", icon("download"), " Download Analysis Data"),
                tags$p("Export your analysis results in various formats for further analysis or publication.")
              ),
              
              # Excel Export
              fluidRow(
                column(width = 6,
                  tags$div(
                    style = "background-color: #e8f5e9; padding: 15px; border-radius: 5px; margin-bottom: 15px;",
                    tags$h5(style = "margin-top: 0; color: #2e7d32;", icon("file-excel"), " Comprehensive Excel Report"),
                    tags$p("Complete analysis with 12 detailed sheets including all sequences, metrics, and statistics."),
                    downloadButton("downloadComprehensiveExcel", "Download Excel Report", 
                                class = "btn-success btn-lg", icon = icon("file-excel"), 
                                style = "width: 100%; padding: 10px;")
                  )
                ),
                column(width = 6,
                  tags$div(
                    style = "background-color: #fff3e0; padding: 15px; border-radius: 5px; margin-bottom: 15px;",
                    tags$h5(style = "margin-top: 0; color: #f57c00;", icon("file-csv"), " Summary CSV"),
                    tags$p("Quick summary table with key metrics for all methods."),
                    downloadButton("downloadComparisonCSV", "Download CSV Summary", 
                                class = "btn-warning btn-lg", icon = icon("file-csv"), 
                                style = "width: 100%; padding: 10px;")
                  )
                )
              ),
              
              # FASTA Downloads
              fluidRow(
                column(width = 12,
                  tags$div(
                    style = "background-color: #f3e5f5; padding: 15px; border-radius: 5px;",
                    tags$h5(style = "margin-top: 0; color: #7b1fa2;", icon("dna"), " FASTA Sequence Downloads"),
                    tags$p("Download generated sequences in FASTA format for each method."),
                    fluidRow(
                      column(width = 3,
                        downloadButton("downloadNeRNAFASTA", "NeRNA Sequences", 
                                    class = "btn-info", icon = icon("download"), 
                                    style = "width: 100%; margin-bottom: 10px;")
                      ),
                      column(width = 3,
                        downloadButton("downloadDinucFASTA", "Dinucleotide Sequences", 
                                    class = "btn-info", icon = icon("download"), 
                                    style = "width: 100%; margin-bottom: 10px;")
                      ),
                      column(width = 3,
                        downloadButton("downloadRandomFASTA", "Random Sequences", 
                                    class = "btn-info", icon = icon("download"), 
                                    style = "width: 100%; margin-bottom: 10px;")
                      ),
                      column(width = 3,
                        downloadButton("downloadAllFASTA", "All Sequences", 
                                    class = "btn-primary", icon = icon("download"), 
                                    style = "width: 100%; margin-bottom: 10px;")
                      )
                    )
                  )
                )
              )
            )
          )
        )
      ),
      
      # Overview Comparison Tab
      tabItem(tabName = "overview_comparison",
        conditionalPanel(
          condition = "output.dataLoaded != true",
          fluidRow(
            box(
              title = "Data Required",
              status = "warning",
              solidHeader = TRUE,
              width = 12,
              tags$div(
                style = "text-align: center; padding: 20px;",
                tags$i(class = "fa fa-exclamation-circle fa-5x", style = "color: #f39c12; margin-bottom: 20px;"),
                tags$h3("No Data Loaded"),
                tags$p("Please upload a FASTA file and click 'Load Sequences' in the Upload tab to proceed."),
                tags$br(),
                actionButton("goToUploadBtn5", "Go to Upload", class = "btn-primary")
              )
            )
          )
        ),
        conditionalPanel(
          condition = "output.dataLoaded == true",
          fluidRow(
            box(
              title = "Methods Comparison Summary",
              status = "primary",
              solidHeader = TRUE,
              width = 12,
              tags$div(
                style = "border-left: 4px solid #6f42c1; padding: 10px; margin-bottom: 15px; background-color: #f8f9fa;",
                tags$h4(style = "margin-top: 0;", "Methods Comparison"),
                tags$p("This page provides a statistical comparison of all negative RNA generation methods. The table shows average values for key metrics like sequence length, GC content, MFE, and base pairing. Use this overview to determine which method best preserves or alters the properties you're interested in.")
              ),
              p("This table shows a comparison of average statistics for all RNA sequence generation methods."),
              br(),
              tags$div(
                style = "background-color: #e7f3ff; border-left: 4px solid #2196F3; padding: 15px; margin-bottom: 20px;",
                tags$h4(style = "color: #1976D2; margin-top: 0;", icon("dna"), " Secondary Structure Calculation"),
                tags$p(style = "color: #1565C0; margin-bottom: 10px;",
                  "Click the button below to calculate secondary structures (MFE) for all original and negative sequences. ",
                  "This is required for accurate MFE values and structure comparisons."
                ),
                actionButton("runAllRNAfoldBtn", 
                            HTML("<i class='fa fa-calculator'></i> Run RNAfold for All Sequences"),
                            class = "btn-info btn-lg",
                            style = "width: 100%; font-weight: bold;"),
                tags$div(
                  id = "rnafoldProgress",
                  style = "margin-top: 15px; display: none;",
                  tags$div(class = "progress",
                    tags$div(class = "progress-bar progress-bar-striped active",
                            role = "progressbar",
                            style = "width: 100%",
                            "Calculating structures...")
                  )
                ),
                uiOutput("rnafoldAllStatus")
              ),
              fluidRow(
                column(6,
                  downloadButton("downloadComparisonExcel", "Download Summary (Excel)", class = "btn-success"),
                  style = "margin-bottom: 20px;"
                ),
                column(6,
                  actionButton("gotoSequenceComparisonBtn", "Sequence Specific Comparison", 
                              class = "btn-primary",
                              style = "width: 100%;"),
                  style = "margin-bottom: 20px;"
                )
              ),
              conditionalPanel(
                condition = "output.allStructuresCalculated == true",
                DTOutput("comparisonSummaryTable")
              ),
              conditionalPanel(
                condition = "output.allStructuresCalculated != true",
                tags$div(
                  style = "background-color: #fff3cd; border-left: 4px solid #ffc107; padding: 20px; margin-top: 15px; text-align: center;",
                  tags$h4(style = "color: #856404; margin-top: 0;", icon("exclamation-triangle"), " Structure Calculation Required"),
                  tags$p(style = "color: #856404; font-size: 16px;",
                    "Secondary structures (MFE values) have not been calculated yet for all sequences."
                  ),
                  tags$p(style = "color: #856404;",
                    "Please click the ", tags$strong("'Run RNAfold for All Sequences'"), " button above to calculate structures and display accurate MFE values."
                  )
                )
              )
            )
          ),
          conditionalPanel(
            condition = "output.allStructuresCalculated == true",
            fluidRow(
              box(
                title = "Method Effectiveness Score",
                status = "success",
                solidHeader = TRUE,
                width = 12,
                p("This visualization compares the effectiveness of each method based on various metrics."),
                plotOutput("methodComparisonPlot", height = "400px")
              )
            )
          )
        )
      )
    )
  )
) 
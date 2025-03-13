library(shiny)
library(shinythemes)

ui <- fluidPage(
  tags$head(
    tags$meta(charset="UTF-8"),
    tags$style(HTML("
      .method-panel {
        border: 1px solid #ddd;
        border-radius: 5px;
        padding: 15px;
        margin-bottom: 20px;
        background-color: #f9f9f9;
      }
      .method-title {
        font-weight: bold;
        margin-bottom: 10px;
        color: #2c3e50;
      }
      .download-section {
        margin-top: 15px;
      }
      .results-panel {
        margin-top: 20px;
        border-top: 1px solid #eee;
        padding-top: 15px;
      }
      .nav-tabs {
        margin-bottom: 20px;
      }
    "))
  ),
  
  theme = shinytheme("flatly"),
  
  titlePanel("RNA Structure Prediction & Negative Data Generation"),
  
  sidebarLayout(
    sidebarPanel(
      fileInput("fasta_file", "Upload FASTA File",
                accept = c(".fasta", ".fa", ".txt")),
      
      conditionalPanel(
        condition = "input.mainTabset == 'overview'",
        hr(),
        h4("RNA Type Selection"),
        selectInput("rnaType", "Select RNA Type:", 
                   choices = c("tRNA", "circRNA", "miRNA", "lncRNA", "Other")),
        
        hr(),
        h4("Structure Predictions"),
        tableOutput("structure_summary"),
        
        actionButton("goToNegative", "Generate Negative Data", 
                    class = "btn-primary btn-block", 
                    style = "margin-top: 20px;")
      ),
      
      conditionalPanel(
        condition = "input.mainTabset == 'negative'",
        hr(),
        h4("Negative Data Methods"),
        tabsetPanel(id = "negativeMethodTabs",
          tabPanel("NeRNA", 
            sliderInput("nernaShift", "Shift Amount:", 
                       min = 1, max = 10, value = 1, step = 1),
            p("Shift amount determines how many positions to shift the sequence."),
            actionButton("generateNeRNA", "Generate NeRNA", class = "btn-primary"),
            hr(),
            downloadButton("downloadNeRNA", "Download NeRNA Sequences")
          ),
          
          tabPanel("Random Shuffling",
            numericInput("shuffleSeed", "Random Seed:", value = 123),
            p("Using the same seed will produce the same shuffled sequences."),
            actionButton("generateShuffle", "Generate Shuffled", class = "btn-primary"),
            hr(),
            downloadButton("downloadShuffle", "Download Shuffled Sequences")
          ),
          
          tabPanel("Dinucleotide",
            numericInput("dinucSeed", "Random Seed:", value = 456),
            p("Preserves dinucleotide frequencies while shuffling."),
            actionButton("generateDinucleotide", "Generate Dinucleotide", class = "btn-primary"),
            hr(),
            downloadButton("downloadDinucleotide", "Download Dinucleotide Sequences")
          )
        )
      ),
      
      conditionalPanel(
        condition = "input.mainTabset == 'results'",
        hr(),
        h4("Export Results"),
        downloadButton("downloadExcel", "Download Full Results (Excel)"),
        downloadButton("downloadCSV", "Download Full Results (CSV)")
      )
    ),
    
    mainPanel(
      tabsetPanel(id = "mainTabset",
        tabPanel("Overview", value = "overview",
          h4("File Information"),
          textOutput("sequence_info"),
          
          hr(),
          
          h4("Structure Predictions"),
          tableOutput("structure_results")
        ),
        
        tabPanel("Negative Data Generation", value = "negative",
          h4("Generated Negative Sequences"),
          
          tabsetPanel(id = "negativeResultTabs",
            tabPanel("NeRNA Results", 
              tableOutput("nernaTable"),
              plotOutput("nernaPlot")
            ),
            
            tabPanel("Random Shuffling Results", 
              tableOutput("shuffleTable"),
              plotOutput("shufflePlot")
            ),
            
            tabPanel("Dinucleotide Results", 
              tableOutput("dinucleotideTable"),
              plotOutput("dinucleotidePlot")
            )
          )
        ),
        
        tabPanel("Results Comparison", value = "results",
          h4("Comparison of Methods"),
          plotOutput("comparisonPlot", height = "400px"),
          
          hr(),
          
          h4("Summary Statistics"),
          tableOutput("comparisonTable")
        ),
        
        tabPanel("Instructions", value = "instructions",
          h4("How to Use This Tool"),
          tags$ol(
            tags$li("Upload your RNA sequence file in FASTA format."),
            tags$li("Select the RNA type from the dropdown menu."),
            tags$li("Review the structure predictions in the Overview tab."),
            tags$li("Click 'Generate Negative Data' to proceed to the negative data generation."),
            tags$li("In the Negative Data Generation tab:"),
            tags$ul(
              tags$li(strong("NeRNA:"), " Set the shift amount and generate shifted sequences."),
              tags$li(strong("Random Shuffling:"), " Set a random seed and generate shuffled sequences."),
              tags$li(strong("Dinucleotide Shuffling:"), " Set a random seed and generate sequences with preserved dinucleotide frequencies.")
            ),
            tags$li("Compare results across methods in the Results Comparison tab."),
            tags$li("Download results using the download buttons.")
          )
        )
      )
    )
  )
) 
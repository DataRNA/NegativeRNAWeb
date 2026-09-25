function(input, output, session) {


  
  # Reactive values to store data
  values <- reactiveValues(
    fastaData = NULL,
    rnafoldResults = NULL,
    nernaResults = NULL,
    dinucleotideResults = NULL,
    randomResults = NULL,
    filteredSeqs = NULL,
    rnaType = NULL,
    dataLoaded = FALSE,
    sequenceBatches = NULL,
    batchesReady = FALSE,
    parallelEnabled = TRUE,
    parallelInfo = NULL
  )
  
  # ============================================================
  # PARALLEL PROCESSING INITIALIZATION
  # ============================================================
  
  # Initialize parallel processing on startup
  observe({
    tryCatch({
      n_cores <- setup_parallel_processing()
      values$parallelInfo <- get_parallel_info()
      message(sprintf("Parallel processing initialized with %d workers", n_cores))
    }, error = function(e) {
      message("Failed to initialize parallel processing: ", e$message)
      values$parallelEnabled <- FALSE
    })
  })
  
  # Update parallel info display
  output$parallelInfoUI <- renderUI({
    if (is.null(values$parallelInfo)) {
      return(tags$p("Loading parallel processing info..."))
    }
    
    info <- values$parallelInfo
    
    tags$div(
      style = "background-color: #e7f3ff; padding: 10px; border-radius: 5px; margin-top: 10px;",
      tags$p(
        tags$strong("Current Plan: "), 
        tags$span(info$plan, style = "color: #0066cc;")
      ),
      tags$p(
        tags$strong("Workers Available: "), 
        tags$span(info$workers, style = "color: #0066cc;")
      ),
      tags$p(
        tags$strong("Status: "), 
        if (info$is_parallel) {
          tags$span("✓ Parallel Processing Active", style = "color: #28a745; font-weight: bold;")
        } else {
          tags$span("Sequential Processing", style = "color: #ffc107;")
        }
      )
    )
  })
  
  # Refresh parallel info button
  observeEvent(input$refreshParallelInfo, {
    tryCatch({
      if (input$enableParallel) {
        n_cores <- setup_parallel_processing()
        showNotification(sprintf("Parallel processing enabled with %d cores", n_cores), 
                        type = "success")
      } else {
        disable_parallel_processing()
        showNotification("Parallel processing disabled", type = "info")
      }
      values$parallelInfo <- get_parallel_info()
    }, error = function(e) {
      showNotification(paste("Error updating parallel settings:", e$message), 
                      type = "error")
    })
  })
  
  # Output variable to control UI visibility
  output$dataLoaded <- reactive({
    return(values$dataLoaded)
  })
  outputOptions(output, "dataLoaded", suspendWhenHidden = FALSE)
  
  # Output variable to show batch download panel
  output$batchesReady <- reactive({
    return(values$batchesReady)
  })
  outputOptions(output, "batchesReady", suspendWhenHidden = FALSE)

  # Clear results that belong to a previously loaded dataset.
  reset_analysis_state <- function() {
    values$rnafoldResults <- NULL
    values$nernaResults <- NULL
    values$dinucleotideResults <- NULL
    values$randomResults <- NULL
    values$nernaStructures <- NULL
    values$dinucStructures <- NULL
    values$randomStructures <- NULL
    values$selectedComparison <- NULL
    values$comparisonSummary <- NULL
    values$sequenceBatches <- NULL
    values$batchesReady <- FALSE
  }
  
  # Output variable to check if RNAfold has been run for all sequences
  output$allStructuresCalculated <- reactive({
    cache_is_complete <- function(cache, sequences) {
      if (is.null(cache) || length(cache) != length(sequences)) return(FALSE)

      all(vapply(seq_along(sequences), function(i) {
        result <- cache[[i]]
        has_valid_structure(result$structure, nchar(sequences[i])) &&
          has_valid_mfe(result$mfe)
      }, logical(1)))
    }

    if (is.null(values$rnafoldResults) || is.null(values$fastaData) ||
        nrow(values$rnafoldResults) != nrow(values$fastaData)) {
      return(FALSE)
    }

    original_ready <- all(vapply(seq_len(nrow(values$rnafoldResults)), function(i) {
      has_valid_structure(
        values$rnafoldResults$SecondaryStructure[i],
        nchar(values$rnafoldResults$Sequence[i])
      ) && has_valid_mfe(values$rnafoldResults$MFE[i])
    }, logical(1)))

    if (!original_ready) return(FALSE)

    if (!is.null(values$nernaResults) &&
        !cache_is_complete(values$nernaStructures, values$nernaResults$NegativeSequence)) {
      return(FALSE)
    }

    if (!is.null(values$dinucleotideResults) &&
        !cache_is_complete(values$dinucStructures, values$dinucleotideResults$NegativeSequence)) {
      return(FALSE)
    }

    if (!is.null(values$randomResults) &&
        !cache_is_complete(values$randomStructures, values$randomResults$NegativeSequence)) {
      return(FALSE)
    }

    TRUE
  })
  outputOptions(output, "allStructuresCalculated", suspendWhenHidden = FALSE)
  
  # Load FASTA file when button is clicked
  observeEvent(input$loadBtn, {
    req(input$fastaFile)
    
    # Read FASTA file and convert to table (automatically filters sequences > 1000 nt, and cleans non-standard bases for tRNA)
    fasta_data <- readFastaToTable(input$fastaFile$datapath, rnaType = input$rnaType)
    
    # Check if any sequences remain after filtering
    if (is.null(fasta_data) || nrow(fasta_data) == 0) {
      showNotification(
        "All sequences were filtered out (longer than 1000 nucleotides). Please provide shorter sequences.",
        type = "error",
        duration = 10
      )
      return(NULL)
    }
    
    # Check if sequences exceed maximum allowed (250)
    if (attr(fasta_data, "exceeds_limit")) {
      total_seqs <- attr(fasta_data, "total_sequences")
      max_seqs <- attr(fasta_data, "max_sequences")
      num_batches <- ceiling(total_seqs / max_seqs)
      
      # Check if any sequences were filtered by length
      original_count <- attr(fasta_data, "original_count")
      filtered_count <- if (!is.null(original_count)) original_count - total_seqs else 0
      
      # Build notification message
      notification_msg <- sprintf(
        "<strong>Too Many Sequences!</strong><br/>
        Your file contains <strong>%d valid sequences</strong> (after filtering)",
        total_seqs
      )
      
      if (filtered_count > 0) {
        notification_msg <- paste0(
          notification_msg,
          sprintf("<br/><small><i>Note: %d sequences were filtered out (>1000 nucleotides)</i></small>",
                  filtered_count)
        )
      }
      
      notification_msg <- paste0(
        notification_msg,
        sprintf("<br/><br/>
        This exceeds the maximum allowed (<strong>%d sequences</strong>).<br/>
        <br/>
        The sequences have been split into <strong>%d batches</strong> of up to 250 sequences each.<br/>
        Please click the <strong>'Download Batches (ZIP)'</strong> button below to download all batches.<br/>
        <br/>
        Then upload and process each batch separately.",
        max_seqs, num_batches)
      )
      
      # Show warning notification
      showNotification(
        HTML(notification_msg),
        type = "error",
        duration = NULL,  # Keep notification visible
        closeButton = TRUE
      )
      
      # Split sequences into batches
      batches <- split_sequences_to_batches(fasta_data, batch_size = max_seqs)
      
      # Store batches for download
      reset_analysis_state()
      values$fastaData <- NULL
      values$filteredSeqs <- NULL
      values$rnaType <- input$rnaType
      values$sequenceBatches <- batches
      values$batchesReady <- TRUE
      values$dataLoaded <- FALSE  # Don't allow analysis
      
      return(NULL)
    }
    
    # If within limit, proceed normally
    reset_analysis_state()
    values$fastaData <- fasta_data
    
    # Store filtered sequences for display in Overview
    filtered_sequences <- attr(fasta_data, "filtered_sequences")
    if (!is.null(filtered_sequences) && nrow(filtered_sequences) > 0) {
      values$filteredSeqs <- filtered_sequences
    } else {
      values$filteredSeqs <- NULL
    }
    
    # Show success notification
    filtered_count <- attr(fasta_data, "filtered_count")
    if (!is.null(filtered_count) && filtered_count > 0) {
      showNotification(
        paste("Loaded", nrow(values$fastaData), "sequences.", filtered_count, "sequences were filtered (>1000 nt)."),
        type = "warning",
        duration = 8
      )
    } else {
      showNotification(
        paste("Loaded", nrow(values$fastaData), "sequences."),
        type = "message",
        duration = 5
      )
    }
    
    # Store selected RNA type
    values$rnaType <- input$rnaType
    
    # Set data as loaded to show other menu items
    values$dataLoaded <- TRUE
    
    # Check if the loaded data is DNA
    isDNA <- attr(values$fastaData, "isDNA")
    
    if (isDNA) {
      # Show warning notification
      showNotification(
        "Warning: DNA sequence detected. Please upload RNA sequences for better results.", 
        type = "warning",
        duration = 10
      )
    } else {
      # Show success notification
      showNotification("RNA sequence loaded successfully", type = "message")
    }
    
    # Automatically switch to overview tab
    updateTabItems(session, "sidebar", "overview")
  })

  # Load a small bundled dataset for a reproducible first-use example.
  load_example_data <- function(selected_type) {
    example_files <- c(
      mirna = file.path("cases", "miRNA.fa"),
      circrna = file.path("cases", "circRNA_1000.fa"),
      lncrna = file.path("cases", "lncRNA_1000.fa"),
      trna = file.path("cases", "tRNA.fasta")
    )

    example_file <- unname(example_files[selected_type])

    if (length(example_file) != 1 || is.na(example_file)) {
      showNotification(
        "Example data are available for miRNA, circRNA, lncRNA, and tRNA. Please select one of these RNA types.",
        type = "warning",
        duration = 8
      )
      return(NULL)
    }

    if (!file.exists(example_file)) {
      showNotification(
        paste("Bundled example file was not found:", example_file),
        type = "error",
        duration = 10
      )
      return(NULL)
    }

    fasta_data <- readFastaToTable(example_file, rnaType = selected_type)
    if (is.null(fasta_data) || nrow(fasta_data) == 0) {
      showNotification("The bundled example file contains no usable sequences.", type = "error")
      return(NULL)
    }

    fasta_data <- fasta_data[seq_len(min(10L, nrow(fasta_data))), , drop = FALSE]
    attr(fasta_data, "exceeds_limit") <- FALSE
    attr(fasta_data, "total_sequences") <- nrow(fasta_data)
    attr(fasta_data, "max_sequences") <- 250L
    attr(fasta_data, "original_count") <- nrow(fasta_data)
    attr(fasta_data, "filtered_count") <- 0L
    attr(fasta_data, "filtered_sequences") <- fasta_data[0, , drop = FALSE]

    reset_analysis_state()
    values$fastaData <- fasta_data
    values$filteredSeqs <- NULL
    values$rnaType <- selected_type
    values$dataLoaded <- TRUE
    updateSelectInput(session, "rnaType", selected = selected_type)

    showNotification(
      paste("Loaded", nrow(fasta_data), selected_type, "example sequences."),
      type = "message",
      duration = 6
    )
    updateTabItems(session, "sidebar", "overview")
  }

  observeEvent(input$loadExampleBtn, {
    load_example_data(input$rnaType)
  })

  observeEvent(input$loadHomeExampleBtn, {
    load_example_data(input$homeExampleType)
  })
  
  # Display sequence table
  output$sequenceTable <- renderDT({
    req(values$fastaData)
    datatable(values$fastaData, 
              options = list(scrollX = TRUE, pageLength = 5),
              rownames = FALSE)
  })
  
  # Display information about filtered sequences
  output$filteredSequencesInfo <- renderUI({
    if (!is.null(values$filteredSeqs) && nrow(values$filteredSeqs) > 0) {
      # Create a summary data frame for display
      filtered_summary <- data.frame(
        `Seq#` = 1:nrow(values$filteredSeqs),
        Accession = values$filteredSeqs$Accession,
        `Length (nt)` = sapply(values$filteredSeqs$Sequence, nchar),
        Status = rep("Filtered (>1000 nt)", nrow(values$filteredSeqs)),
        stringsAsFactors = FALSE,
        check.names = FALSE
      )
      
      tagList(
        tags$div(
          style = "background-color: #fff3cd; border-left: 4px solid #ffc107; padding: 15px; margin-bottom: 15px;",
          tags$h4(style = "color: #856404; margin-top: 0;", icon("exclamation-triangle"), " Filtered Sequences"),
          tags$p(style = "color: #856404; margin-bottom: 5px;",
            tags$strong(paste0(nrow(values$filteredSeqs), " sequences")),
            " were automatically filtered out because they exceed the maximum allowed length of 1,000 nucleotides."
          ),
          tags$p(style = "color: #856404; margin-bottom: 0; font-size: 14px;",
            icon("info-circle"),
            " These sequences are not included in the analysis. If you need to analyze them, please split them into smaller segments."
          )
        ),
        tags$div(
          style = "margin-top: 15px;",
          DTOutput("filteredSequencesTable")
        ),
        tags$div(
          style = "margin-top: 15px; padding: 10px; background-color: #d1ecf1; border-left: 4px solid #0c5460;",
          tags$p(style = "margin-bottom: 5px; color: #0c5460;",
            tags$strong("Note:"), " If your file contained more than 250 sequences (including filtered ones),",
            " you would have received batch files with a complete list of filtered sequences in the ZIP archive."
          )
        )
      )
    } else {
      tags$div(
        style = "background-color: #d4edda; border-left: 4px solid #28a745; padding: 15px;",
        tags$p(style = "color: #155724; margin-bottom: 0;",
          icon("check-circle"), 
          tags$strong(" All sequences passed the length filter!"),
          " No sequences were filtered out."
        )
      )
    }
  })
  
  # Render filtered sequences table
  output$filteredSequencesTable <- renderDT({
    req(values$filteredSeqs)
    
    filtered_summary <- data.frame(
      `Seq#` = 1:nrow(values$filteredSeqs),
      Accession = values$filteredSeqs$Accession,
      `Length (nt)` = sapply(values$filteredSeqs$Sequence, nchar),
      Status = rep("Filtered", nrow(values$filteredSeqs)),
      stringsAsFactors = FALSE,
      check.names = FALSE
    )
    
    datatable(
      filtered_summary,
      options = list(
        pageLength = 10,
        scrollX = TRUE,
        dom = 'Bfrtip',
        buttons = c('copy', 'csv', 'excel'),
        columnDefs = list(
          list(className = 'dt-center', targets = c(0, 2, 3))
        )
      ),
      rownames = FALSE,
      class = 'cell-border stripe',
      extensions = 'Buttons'
    ) %>%
      formatStyle(
        'Length (nt)',
        backgroundColor = styleInterval(c(1000), c('#fff3cd', '#f8d7da')),
        fontWeight = 'bold'
      ) %>%
      formatStyle(
        'Status',
        color = '#721c24',
        backgroundColor = '#f8d7da',
        fontWeight = 'bold'
      )
  })
  
  # Calculate sequence statistics
  output$seqCountBox <- renderValueBox({
    req(values$fastaData)
    valueBox(
      nrow(values$fastaData),
      "Number of Sequences",
      icon = icon("list"),
      color = "blue"
    )
  })
  
  output$avgLengthBox <- renderValueBox({
    req(values$fastaData)
    avgLength <- round(mean(nchar(values$fastaData$Sequence)), 1)
    valueBox(
      avgLength,
      "Average Sequence Length",
      icon = icon("ruler"),
      color = "green"
    )
  })
  
  output$gcContentBox <- renderValueBox({
    req(values$fastaData)
    # Calculate GC content
    gcCount <- sum(sapply(values$fastaData$Sequence, function(seq) {
      sum(gregexpr("[GC]", seq)[[1]] > 0)
    }))
    totalBases <- sum(nchar(values$fastaData$Sequence))
    gcContent <- round(gcCount / totalBases * 100, 1)
    
    valueBox(
      paste0(gcContent, "%"),
      "GC Content",
      icon = icon("percent"),
      color = "yellow"
    )
  })
  
  # Run RNAfold when button is clicked
  observeEvent(input$runRNAfoldBtn, {
    req(values$fastaData)
    
    # Create a temporary FASTA file with only the filtered (valid) sequences
    temp_fasta <- tempfile(fileext = ".fasta")
    
    tryCatch({
      # Write only the valid sequences (those in values$fastaData) to temp file
      fasta_content <- character(nrow(values$fastaData) * 2)
      for (i in 1:nrow(values$fastaData)) {
        fasta_content[(i - 1) * 2 + 1] <- paste0(">", values$fastaData$Accession[i])
        fasta_content[(i - 1) * 2 + 2] <- values$fastaData$Sequence[i]
      }
      writeLines(fasta_content, temp_fasta)
      
      # Show processing notification with sequence count
      filtered_count <- attr(values$fastaData, "filtered_count")
      if (!is.null(filtered_count) && filtered_count > 0) {
        notification_msg <- sprintf("Running RNAfold on %d sequences (%d sequences were filtered out)...", 
                                   nrow(values$fastaData), filtered_count)
      } else {
        notification_msg <- sprintf("Running RNAfold on %d sequences...", nrow(values$fastaData))
      }
      notification <- showNotification(notification_msg, type = "message", duration = NULL)
      
      # Run RNAfold with the filtered sequences (with parallel processing support)
      values$rnafoldResults <- runRNAfold(temp_fasta, 
                                          rnaType = values$rnaType,
                                          use_parallel = input$enableParallel)
      
      # Clean up temp file
      unlink(temp_fasta)
      
      # Remove notification
      removeNotification(notification)
      
      if (!is.null(values$rnafoldResults)) {
        # Show success notification
        showNotification(
          sprintf("Secondary structures generated successfully for %d sequences", 
                  nrow(values$rnafoldResults)), 
          type = "message",
          duration = 5
        )
      } else {
        # Show error notification
        showNotification("Error generating secondary structures", type = "error", duration = 10)
      }
    }, error = function(e) {
      # Clean up temp file in case of error
      if (file.exists(temp_fasta)) {
        unlink(temp_fasta)
      }
      showNotification(paste("Error:", e$message), type = "error", duration = 10)
    })
  })
  
  # Display RNAfold results as a table
  output$rnafoldTable <- renderDT({
    req(values$rnafoldResults)
    datatable(values$rnafoldResults, 
              options = list(scrollX = TRUE, pageLength = 10),
              rownames = FALSE)
  })
  
  # Download handlers for CSV and Excel
  output$downloadCSV <- downloadHandler(
    filename = function() {
      paste("RNAfold-Results-", Sys.Date(), ".csv", sep = "")
    },
    content = function(file) {
      if (!is.null(values$rnafoldResults)) {
        write.csv(values$rnafoldResults, file, row.names = FALSE)
      } else {
        # Create empty dataframe with message
        df <- data.frame(Message = "No RNAfold results available. Please generate secondary structures first.")
        write.csv(df, file, row.names = FALSE)
      }
    }
  )
  
  output$downloadExcel <- downloadHandler(
    filename = function() {
      paste("RNAfold-Results-", Sys.Date(), ".xlsx", sep = "")
    },
    content = function(file) {
      if (!is.null(values$rnafoldResults)) {
        # Create workbook
        wb <- createWorkbook()
        
        # Add a worksheet
        addWorksheet(wb, "RNAfold Results")
        
        # Write data
        writeData(wb, "RNAfold Results", values$rnafoldResults, rowNames = FALSE)
        
        # Style headers
        headerStyle <- createStyle(
          fontColour = "#FFFFFF", fgFill = "#4F81BD",
          halign = "center", valign = "center",
          textDecoration = "bold"
        )
        
        # Apply style to header row
        addStyle(wb, "RNAfold Results", headerStyle, rows = 1, cols = 1:ncol(values$rnafoldResults), gridExpand = TRUE)
        
        # Save workbook
        saveWorkbook(wb, file, overwrite = TRUE)
      } else {
        # Create empty workbook with message
        wb <- createWorkbook()
        addWorksheet(wb, "No Data")
        writeData(wb, "No Data", "No RNAfold results available. Please generate secondary structures first.")
        saveWorkbook(wb, file, overwrite = TRUE)
      }
    }
  )
  
  # Download FASTA format
  output$downloadFASTA <- downloadHandler(
    filename = function() {
      paste("RNAfold-Sequences-", Sys.Date(), ".fasta", sep = "")
    },
    content = function(file) {
      if (!is.null(values$rnafoldResults)) {
        # Create a connection to the file
        con <- file(file, "w")
        
        # Write each sequence in FASTA format
        for (i in 1:nrow(values$rnafoldResults)) {
          # Header line with sequence ID and MFE info
          writeLines(paste0(">", values$rnafoldResults$Accession[i], " | MFE: ", values$rnafoldResults$MFE[i]), con)
          
          # Sequence (adding line breaks every 60 characters for readability)
          seq <- values$rnafoldResults$Sequence[i]
          seq_wrapped <- gsub("(.{60})", "\\1\n", seq)
          writeLines(seq_wrapped, con)
        }
        
        # Close the connection
        close(con)
      } else {
        # Create empty file with message
        writeLines("No RNAfold results available. Please generate secondary structures first.", file)
      }
    }
  )
  
  # Run NeRNA method when button is clicked
  observeEvent(input$runNeRNABtn, {
    req(values$fastaData)
    
    # Show processing notification
    notification <- showNotification("Generating negative RNA sequences with Advanced NeRNA...", type = "message", duration = NULL)
    
    # Get parameters from input
    shifting_size <- input$nernaShiftingSize
    seed <- input$nernaSeed
    
    # Run NeRNA with correct algorithm (with parallel processing support)
    generatedSeqs <- generateNeRNA(
      values$fastaData$Sequence, 
      shifting_size = shifting_size,
      seed = seed,
      rnafold_results = values$rnafoldResults,
      use_parallel = input$enableParallel
    )
    
    # Create result data frame using helper function
    values$nernaResults <- create_results_dataframe(
      values$fastaData$Sequence, 
      generatedSeqs, 
      values$fastaData$Accession,
      additional_cols = list(ShiftAmount = shifting_size)
    )
    
    # Remove notification
    removeNotification(notification)
    
    # Show success notification with method details
    showNotification(paste("NeRNA sequences generated successfully with shift =", shifting_size), type = "message")
  })
  
  # Display NeRNA results
  output$nernaOutput <- renderDT({
    req(values$nernaResults)
    datatable(values$nernaResults, 
              options = list(scrollX = TRUE, pageLength = 5),
              rownames = FALSE)
  })
  
  # Download handlers for NeRNA
  output$downloadNeRNACSV <- downloadHandler(
    filename = function() {
      paste("NeRNA-Results-", Sys.Date(), ".csv", sep = "")
    },
    content = function(file) {
      if (!is.null(values$nernaResults)) {
        write.csv(values$nernaResults, file, row.names = FALSE)
      } else {
        # Create empty dataframe with message
        df <- data.frame(Message = "No NeRNA results available. Please generate with NeRNA first.")
        write.csv(df, file, row.names = FALSE)
      }
    }
  )
  
  output$downloadNeRNAExcel <- downloadHandler(
    filename = function() {
      paste("NeRNA-Results-", Sys.Date(), ".xlsx", sep = "")
    },
    content = function(file) {
      if (!is.null(values$nernaResults)) {
        # Create workbook
        wb <- createWorkbook()
        
        # Add a worksheet
        addWorksheet(wb, "NeRNA Results")
        
        # Write data
        writeData(wb, "NeRNA Results", values$nernaResults, rowNames = FALSE)
        
        # Style headers
        headerStyle <- createStyle(
          fontColour = "#FFFFFF", fgFill = "#4F81BD",
          halign = "center", valign = "center",
          textDecoration = "bold"
        )
        
        # Apply style to header row
        addStyle(wb, "NeRNA Results", headerStyle, rows = 1, cols = 1:ncol(values$nernaResults), gridExpand = TRUE)
        
        # Save workbook
        saveWorkbook(wb, file, overwrite = TRUE)
      } else {
        # Create empty workbook with message
        wb <- createWorkbook()
        addWorksheet(wb, "No Data")
        writeData(wb, "No Data", "No NeRNA results available. Please generate with NeRNA first.")
        saveWorkbook(wb, file, overwrite = TRUE)
      }
    }
  )
  
  # Download NeRNA results in FASTA format
  output$downloadNeRNAFASTA <- downloadHandler(
    filename = function() {
      paste("NeRNA-Sequences-", Sys.Date(), ".fasta", sep = "")
    },
    content = function(file) {
      if (!is.null(values$nernaResults)) {
        # Create a connection to the file
        con <- file(file, "w")
        
        # Write each sequence in FASTA format
        for (i in 1:nrow(values$nernaResults)) {
          # Header line with sequence ID and negative type
          method_info <- if("Method" %in% colnames(values$nernaResults)) {
            paste0("Method: ", values$nernaResults$Method[i], " | Shift: ", values$nernaResults$ShiftAmount[i])
          } else {
            paste0("Method: NeRNA | Shift: ", input$nernaShiftingSize)
          }
          writeLines(paste0(">", values$nernaResults$OriginalAccession[i], " | ", method_info), con)
          
          # Sequence (adding line breaks every 60 characters for readability)
          seq <- values$nernaResults$NegativeSequence[i]
          seq_wrapped <- gsub("(.{60})", "\\1\n", seq)
          writeLines(seq_wrapped, con)
        }
        
        # Close the connection
        close(con)
      } else {
        # Create empty file with message
        writeLines("No NeRNA results available. Please generate with NeRNA first.", file)
      }
    }
  )
  
  # Run Dinucleotide Shuffling method when button is clicked
  observeEvent(input$runDinucleotideBtn, {
    req(values$fastaData)
    
    # Show processing notification
    notification <- showNotification("Generating negative RNA sequences with Dinucleotide Shuffling...", type = "message", duration = NULL)
    
    # Get the seed from input
    seed <- input$dinucleotideSeed
    
    # Run Dinucleotide Shuffling with seed (with parallel processing support)
    generatedSeqs <- dinucleotideShuffle(values$fastaData$Sequence, 
                                         seed = seed,
                                         use_parallel = input$enableParallel)
    
    # Create result data frame using helper function
    values$dinucleotideResults <- create_results_dataframe(
      values$fastaData$Sequence, 
      generatedSeqs, 
      values$fastaData$Accession
    )
    
    # Remove notification
    removeNotification(notification)
    
    # Show success notification
    showNotification("Dinucleotide shuffled sequences generated successfully", type = "message")
  })
  
  # Display Dinucleotide Shuffling results
  output$dinucleotideOutput <- renderDT({
    req(values$dinucleotideResults)
    datatable(values$dinucleotideResults, 
              options = list(scrollX = TRUE, pageLength = 5),
              rownames = FALSE)
  })
  
  # Download handlers for Dinucleotide Shuffling
  output$downloadDinucleotideCSV <- downloadHandler(
    filename = function() {
      paste("Dinucleotide-Results-", Sys.Date(), ".csv", sep = "")
    },
    content = function(file) {
      if (!is.null(values$dinucleotideResults)) {
        write.csv(values$dinucleotideResults, file, row.names = FALSE)
      } else {
        # Create empty dataframe with message
        df <- data.frame(Message = "No Dinucleotide Shuffling results available. Please generate with Dinucleotide Shuffling first.")
        write.csv(df, file, row.names = FALSE)
      }
    }
  )
  
  output$downloadDinucleotideExcel <- downloadHandler(
    filename = function() {
      paste("Dinucleotide-Results-", Sys.Date(), ".xlsx", sep = "")
    },
    content = function(file) {
      if (!is.null(values$dinucleotideResults)) {
        # Create workbook
        wb <- createWorkbook()
        
        # Add a worksheet
        addWorksheet(wb, "Dinucleotide Results")
        
        # Write data
        writeData(wb, "Dinucleotide Results", values$dinucleotideResults, rowNames = FALSE)
        
        # Style headers
        headerStyle <- createStyle(
          fontColour = "#FFFFFF", fgFill = "#4F81BD",
          halign = "center", valign = "center",
          textDecoration = "bold"
        )
        
        # Apply style to header row
        addStyle(wb, "Dinucleotide Results", headerStyle, rows = 1, cols = 1:ncol(values$dinucleotideResults), gridExpand = TRUE)
        
        # Save workbook
        saveWorkbook(wb, file, overwrite = TRUE)
      } else {
        # Create empty workbook with message
        wb <- createWorkbook()
        addWorksheet(wb, "No Data")
        writeData(wb, "No Data", "No Dinucleotide Shuffling results available. Please generate with Dinucleotide Shuffling first.")
        saveWorkbook(wb, file, overwrite = TRUE)
      }
    }
  )
  
  # Download Dinucleotide results in FASTA format
  output$downloadDinucleotideFASTA <- downloadHandler(
    filename = function() {
      paste("Dinucleotide-Sequences-", Sys.Date(), ".fasta", sep = "")
    },
    content = function(file) {
      if (!is.null(values$dinucleotideResults)) {
        # Create a connection to the file
        con <- file(file, "w")
        
        # Write each sequence in FASTA format
        for (i in 1:nrow(values$dinucleotideResults)) {
          # Header line with sequence ID and seed info
          writeLines(paste0(">", values$dinucleotideResults$OriginalAccession[i], " | Method: Dinucleotide | Seed: ", input$dinucleotideSeed), con)
          
          # Sequence (adding line breaks every 60 characters for readability)
          seq <- values$dinucleotideResults$NegativeSequence[i]
          seq_wrapped <- gsub("(.{60})", "\\1\n", seq)
          writeLines(seq_wrapped, con)
        }
        
        # Close the connection
        close(con)
      } else {
        # Create empty file with message
        writeLines("No Dinucleotide Shuffling results available. Please generate with Dinucleotide Shuffling first.", file)
      }
    }
  )
  
  # Run Random Shuffling method when button is clicked
  observeEvent(input$runRandomBtn, {
    req(values$fastaData)
    
    # Show processing notification 
    notification <- showNotification("Generating negative RNA sequences with Random Shuffling...", type = "message", duration = NULL)
    
    # Get the seed from input
    seed <- input$randomSeed
    
    # Run Random Shuffling with seed (with parallel processing support)
    generatedSeqs <- randomShuffle(values$fastaData$Sequence, 
                                   seed = seed,
                                   use_parallel = input$enableParallel)
    
    # Create result data frame using helper function
    values$randomResults <- create_results_dataframe(
      values$fastaData$Sequence, 
      generatedSeqs, 
      values$fastaData$Accession
    )
    
    # Remove notification
    removeNotification(notification)
    
    # Show success notification
    showNotification("Random shuffled sequences generated successfully", type = "message")
  })
  
  # Display Random Shuffling results
  output$randomOutput <- renderDT({
    req(values$randomResults)
    datatable(values$randomResults, 
              options = list(scrollX = TRUE, pageLength = 5),
              rownames = FALSE)
  })
  
  # Download handlers for Random Shuffling
  output$downloadRandomCSV <- downloadHandler(
    filename = function() {
      paste("Random-Results-", Sys.Date(), ".csv", sep = "")
    },
    content = function(file) {
      if (!is.null(values$randomResults)) {
        write.csv(values$randomResults, file, row.names = FALSE)
      } else {
        # Create empty dataframe with message
        df <- data.frame(Message = "No Random Shuffling results available. Please generate with Random Shuffling first.")
        write.csv(df, file, row.names = FALSE)
      }
    }
  )
  
  output$downloadRandomExcel <- downloadHandler(
    filename = function() {
      paste("Random-Results-", Sys.Date(), ".xlsx", sep = "")
    },
    content = function(file) {
      if (!is.null(values$randomResults)) {
        # Create workbook
        wb <- createWorkbook()
        
        # Add a worksheet
        addWorksheet(wb, "Random Results")
        
        # Write data
        writeData(wb, "Random Results", values$randomResults, rowNames = FALSE)
        
        # Style headers
        headerStyle <- createStyle(
          fontColour = "#FFFFFF", fgFill = "#4F81BD",
          halign = "center", valign = "center",
          textDecoration = "bold"
        )
        
        # Apply style to header row
        addStyle(wb, "Random Results", headerStyle, rows = 1, cols = 1:ncol(values$randomResults), gridExpand = TRUE)
        
        # Save workbook
        saveWorkbook(wb, file, overwrite = TRUE)
      } else {
        # Create empty workbook with message
        wb <- createWorkbook()
        addWorksheet(wb, "No Data")
        writeData(wb, "No Data", "No Random Shuffling results available. Please generate with Random Shuffling first.")
        saveWorkbook(wb, file, overwrite = TRUE)
      }
    }
  )
  
  # Download Random results in FASTA format
  output$downloadRandomFASTA <- downloadHandler(
    filename = function() {
      paste("Random-Sequences-", Sys.Date(), ".fasta", sep = "")
    },
    content = function(file) {
      if (!is.null(values$randomResults)) {
        # Create a connection to the file
        con <- file(file, "w")
        
        # Write each sequence in FASTA format
        for (i in 1:nrow(values$randomResults)) {
          # Header line with sequence ID and seed info
          writeLines(paste0(">", values$randomResults$OriginalAccession[i], " | Method: Random | Seed: ", input$randomSeed), con)
          
          # Sequence (adding line breaks every 60 characters for readability)
          seq <- values$randomResults$NegativeSequence[i]
          seq_wrapped <- gsub("(.{60})", "\\1\n", seq)
          writeLines(seq_wrapped, con)
        }
        
        # Close the connection
        close(con)
      } else {
        # Create empty file with message
        writeLines("No Random Shuffling results available. Please generate with Random Shuffling first.", file)
      }
    }
  )
  
  # Run all generation methods
  observeEvent(input$generateAllBtn, {
    req(values$fastaData)
    
    # Show processing notification
    notification <- showNotification("Generating sequences with all methods...", type = "message", duration = NULL)
    
    # Get parameters
    shifting_size <- input$nernaShiftingSize
    dinucleotide_seed <- input$dinucleotideSeed
    random_seed <- input$randomSeed
    
    # Step 1: Generate NeRNA sequences (with parallel processing)
    tryCatch({
      nerna_seqs <- generateNeRNA(
        values$fastaData$Sequence, 
        shifting_size = shifting_size,
        seed = input$nernaSeed,
        rnafold_results = values$rnafoldResults,
        use_parallel = input$enableParallel
      )
      
      # Create result data frame using helper function
      values$nernaResults <- create_results_dataframe(
        values$fastaData$Sequence, 
        nerna_seqs, 
        values$fastaData$Accession,
        additional_cols = list(ShiftAmount = shifting_size)
      )
    }, error = function(e) {
      showNotification(paste("Error in NeRNA generation:", e$message), type = "error")
    })
    
    # Step 2: Generate Dinucleotide Shuffling sequences (with parallel processing)
    tryCatch({
      dinu_seqs <- dinucleotideShuffle(values$fastaData$Sequence, 
                                       seed = dinucleotide_seed,
                                       use_parallel = input$enableParallel)
      
      # Create result data frame using helper function
      values$dinucleotideResults <- create_results_dataframe(
        values$fastaData$Sequence, 
        dinu_seqs, 
        values$fastaData$Accession
      )
    }, error = function(e) {
      showNotification(paste("Error in Dinucleotide Shuffling:", e$message), type = "error")
    })
    
    # Step 3: Generate Random Shuffling sequences (with parallel processing)
    tryCatch({
      random_seqs <- randomShuffle(values$fastaData$Sequence, 
                                   seed = random_seed,
                                   use_parallel = input$enableParallel)
      
      # Create result data frame using helper function
      values$randomResults <- create_results_dataframe(
        values$fastaData$Sequence, 
        random_seqs, 
        values$fastaData$Accession
      )
    }, error = function(e) {
      showNotification(paste("Error in Random Shuffling:", e$message), type = "error")
    })
    
    # Remove notification
    removeNotification(notification)
    
    # Show success notification
    showNotification("All methods completed successfully", type = "message")
  })
  
  # Download all results as CSV (combined)
  output$downloadAllCSV <- downloadHandler(
    filename = function() {
      paste("All-Negative-Results-", Sys.Date(), ".csv", sep = "")
    },
    content = function(file) {
      # Create a combined dataframe
      combined_data <- data.frame()
      
      # Add NeRNA results if available
      if (!is.null(values$nernaResults)) {
        nerna_data <- values$nernaResults
        nerna_data$Method <- "NeRNA"
        combined_data <- rbind(combined_data, nerna_data)
      }
      
      # Add Dinucleotide results if available
      if (!is.null(values$dinucleotideResults)) {
        dinu_data <- values$dinucleotideResults
        dinu_data$Method <- "Dinucleotide"
        combined_data <- rbind(combined_data, dinu_data)
      }
      
      # Add Random results if available
      if (!is.null(values$randomResults)) {
        random_data <- values$randomResults
        random_data$Method <- "Random"
        combined_data <- rbind(combined_data, random_data)
      }
      
      if (nrow(combined_data) > 0) {
        write.csv(combined_data, file, row.names = FALSE)
      } else {
        # No data available
        df <- data.frame(Message = "No results available. Please generate negative sequences first.")
        write.csv(df, file, row.names = FALSE)
      }
    }
  )
  
  # Download all results as Excel
  output$downloadAllExcel <- downloadHandler(
    filename = function() {
      paste("All-Negative-Results-", Sys.Date(), ".xlsx", sep = "")
    },
    content = function(file) {
      # Create workbook
      wb <- createWorkbook()
      
      # Header style
      headerStyle <- createStyle(
        fontColour = "#FFFFFF", fgFill = "#4F81BD",
        halign = "center", valign = "center",
        textDecoration = "bold"
      )
      
      # Check and add NeRNA results
      if (!is.null(values$nernaResults)) {
        addWorksheet(wb, "NeRNA Results")
        writeData(wb, "NeRNA Results", values$nernaResults, rowNames = FALSE)
        addStyle(wb, "NeRNA Results", headerStyle, rows = 1, cols = 1:ncol(values$nernaResults), gridExpand = TRUE)
      }
      
      # Check and add Dinucleotide results
      if (!is.null(values$dinucleotideResults)) {
        addWorksheet(wb, "Dinucleotide Results")
        writeData(wb, "Dinucleotide Results", values$dinucleotideResults, rowNames = FALSE)
        addStyle(wb, "Dinucleotide Results", headerStyle, rows = 1, cols = 1:ncol(values$dinucleotideResults), gridExpand = TRUE)
      }
      
      # Check and add Random results
      if (!is.null(values$randomResults)) {
        addWorksheet(wb, "Random Results")
        writeData(wb, "Random Results", values$randomResults, rowNames = FALSE)
        addStyle(wb, "Random Results", headerStyle, rows = 1, cols = 1:ncol(values$randomResults), gridExpand = TRUE)
      }
      
      if (is.null(values$nernaResults) && is.null(values$dinucleotideResults) && is.null(values$randomResults)) {
        # No data available
        addWorksheet(wb, "No Data")
        writeData(wb, "No Data", "No results available. Please generate negative sequences first.")
      }
      
      # Save workbook
      saveWorkbook(wb, file, overwrite = TRUE)
    }
  )
  
  
  # Update sequence selection dropdown
  observe({
    if (!is.null(values$fastaData)) {
      updateSelectInput(session, "comparisonSequence", 
                       choices = setNames(values$fastaData$Accession, values$fastaData$Accession))
    }
  })
  
  # Calculate statistics for a specific sequence
  calculateSequenceStats <- function(sequence, structure = NULL) {
    # Length
    length <- nchar(sequence)
    
    # GC Content
    gc_count <- sum(gregexpr("[GC]", sequence)[[1]] > 0)
    gc_content <- round(gc_count / length * 100, 1)
    
    # Paired and unpaired bases are reported only for a valid RNAfold
    # dot-bracket structure. Missing structures must remain missing.
    paired_bases <- NA_real_
    unpaired_bases <- NA_real_

    if (has_valid_structure(structure, length)) {
      structure_chars <- strsplit(structure, "", fixed = TRUE)[[1]]
      paired_bases <- sum(structure_chars %in% c("(", ")"))
      unpaired_bases <- sum(structure_chars == ".")
    }
    
    return(list(
      Length = length,
      PairedBases = paired_bases,
      UnpairedBases = unpaired_bases,
      GCContent = gc_content
    ))
  }

  # RNAfold-derived values must only be displayed when RNAfold returned a
  # valid dot-bracket structure or a finite MFE value. Missing calculations
  # are represented as NULL/NA rather than estimated or substituted values.
  has_valid_structure <- function(structure, expected_length = NULL) {
    if (is.null(structure) || length(structure) != 1 || is.na(structure) || !nzchar(structure)) {
      return(FALSE)
    }

    if (!grepl("^[().]+$", structure)) {
      return(FALSE)
    }

    structure_chars <- strsplit(structure, "", fixed = TRUE)[[1]]
    balance <- cumsum(ifelse(structure_chars == "(", 1L, ifelse(structure_chars == ")", -1L, 0L)))
    if (any(balance < 0L) || tail(balance, 1) != 0L) {
      return(FALSE)
    }

    if (!is.null(expected_length) && nchar(structure) != expected_length) {
      return(FALSE)
    }

    TRUE
  }

  has_valid_mfe <- function(mfe) {
    !is.null(mfe) && length(mfe) == 1 && !is.na(mfe) && is.finite(mfe)
  }
  
  # Helper function to calculate similarity metrics and create results data frame
  create_results_dataframe <- function(original_sequences, generated_sequences, accessions, additional_cols = NULL) {
    # Calculate structural similarity metrics for each sequence
    similarity_metrics <- lapply(1:length(original_sequences), function(i) {
      metrics <- calculate_structural_similarity(original_sequences[i], generated_sequences[i])
      return(metrics)
    })
    
    # Extract metrics for data frame
    seq_identity <- sapply(similarity_metrics, function(x) x$sequence_identity)
    dinuc_similarity <- sapply(similarity_metrics, function(x) x$dinucleotide_similarity)
    gc_similarity <- sapply(similarity_metrics, function(x) x$gc_similarity)
    
    # Create base data frame
    result_df <- data.frame(
      OriginalAccession = accessions,
      OriginalSequence = original_sequences,
      NegativeSequence = generated_sequences,
      SequenceIdentity = round(seq_identity, 2),
      DinucleotideSimilarity = round(dinuc_similarity, 2),
      GCSimilarity = round(gc_similarity, 2),
      stringsAsFactors = FALSE
    )
    
    # Add additional columns if provided
    if (!is.null(additional_cols)) {
      for (col_name in names(additional_cols)) {
        result_df[[col_name]] <- additional_cols[[col_name]]
      }
    }
    
    return(result_df)
  }
  
  # Generate summary statistics
  generateSummaryStats <- function() {
    req(values$fastaData)

    safe_mean <- function(x, digits = 1) {
      x <- suppressWarnings(as.numeric(x))
      x <- x[is.finite(x)]
      if (length(x) == 0) return(NA_real_)
      round(mean(x), digits)
    }

    structure_from_cache <- function(cache, index, expected_length) {
      if (is.null(cache) || length(cache) < index || is.null(cache[[index]])) return(NULL)
      structure <- cache[[index]]$structure
      if (has_valid_structure(structure, expected_length)) structure else NULL
    }

    mfe_from_cache <- function(cache, index) {
      if (is.null(cache) || length(cache) < index || is.null(cache[[index]])) return(NA_real_)
      mfe <- cache[[index]]$mfe
      if (has_valid_mfe(mfe)) as.numeric(mfe) else NA_real_
    }

    summarize_method <- function(sequences, structures = NULL, mfe_values = NULL, cache = NULL) {
      if (is.null(sequences) || length(sequences) == 0) {
        return(c(MFE = NA, Length = NA, Paired = NA, Unpaired = NA, GC = NA))
      }

      stats <- lapply(seq_along(sequences), function(i) {
        structure <- if (!is.null(cache)) {
          structure_from_cache(cache, i, nchar(sequences[i]))
        } else if (!is.null(structures) && length(structures) >= i &&
                   has_valid_structure(structures[i], nchar(sequences[i]))) {
          structures[i]
        } else {
          NULL
        }
        calculateSequenceStats(sequences[i], structure)
      })

      if (!is.null(cache)) {
        mfe_values <- vapply(seq_along(sequences), function(i) {
          mfe_from_cache(cache, i)
        }, numeric(1))
      }

      c(
        MFE = safe_mean(mfe_values, 2),
        Length = safe_mean(vapply(stats, function(x) x$Length, numeric(1))),
        Paired = safe_mean(vapply(stats, function(x) x$PairedBases, numeric(1))),
        Unpaired = safe_mean(vapply(stats, function(x) x$UnpairedBases, numeric(1))),
        GC = safe_mean(vapply(stats, function(x) x$GCContent, numeric(1)))
      )
    }

    original_summary <- summarize_method(
      values$fastaData$Sequence,
      structures = if (!is.null(values$rnafoldResults)) values$rnafoldResults$SecondaryStructure else NULL,
      mfe_values = if (!is.null(values$rnafoldResults)) values$rnafoldResults$MFE else NULL
    )
    nerna_summary <- summarize_method(
      if (!is.null(values$nernaResults)) values$nernaResults$NegativeSequence else NULL,
      cache = values$nernaStructures
    )
    dinucleotide_summary <- summarize_method(
      if (!is.null(values$dinucleotideResults)) values$dinucleotideResults$NegativeSequence else NULL,
      cache = values$dinucStructures
    )
    random_summary <- summarize_method(
      if (!is.null(values$randomResults)) values$randomResults$NegativeSequence else NULL,
      cache = values$randomStructures
    )

    summaries <- rbind(original_summary, nerna_summary, dinucleotide_summary, random_summary)
    data.frame(
      Method = c("Original", "NeRNA", "Dinucleotide", "Random"),
      Avg.MFE = summaries[, "MFE"],
      Avg.Length = summaries[, "Length"],
      Avg.Paired.Bases = summaries[, "Paired"],
      Avg.Unpaired.Bases = summaries[, "Unpaired"],
      Avg.GC.Content.... = summaries[, "GC"],
      stringsAsFactors = FALSE,
      row.names = NULL
    )
  }
  
  # Generate summary statistics table
  output$summaryStatsTable <- renderDT({
    summary_df <- generateSummaryStats()
    
    datatable(summary_df, 
              options = list(
                dom = 't',  # Shows just the table without search, pagination, etc.
                ordering = FALSE,
                columnDefs = list(list(className = 'dt-center', targets = 1:5))
              ),
              rownames = FALSE)
  })
  
  # Compare selected sequence
  observeEvent(input$compareBtn, {
    req(values$fastaData, input$comparisonSequence)
    
    # Get the selected sequence by accession (now it's a string)
    selected_acc <- input$comparisonSequence
    
    # Find the index
    seq_idx <- which(values$fastaData$Accession == selected_acc)
    if (length(seq_idx) == 0) return(NULL)
    seq_idx <- seq_idx[1]
    
    # Extract data for the selected sequence
    original_seq <- values$fastaData$Sequence[seq_idx]
    original_acc <- values$fastaData$Accession[seq_idx]
    
    # Get original structure and MFE if available
    original_structure <- if (!is.null(values$rnafoldResults) && nrow(values$rnafoldResults) > 0) {
      idx <- which(values$rnafoldResults$Accession == selected_acc)
      if (length(idx) > 0) values$rnafoldResults$SecondaryStructure[idx[1]] else NULL
    } else { NULL }
    
    original_mfe <- if (!is.null(values$rnafoldResults) && nrow(values$rnafoldResults) > 0) {
      idx <- which(values$rnafoldResults$Accession == selected_acc)
      if (length(idx) > 0) values$rnafoldResults$MFE[idx[1]] else NA
    } else { NA }
    
    # Get generated sequences by accession
    nerna_seq <- if (!is.null(values$nernaResults) && nrow(values$nernaResults) > 0) {
      idx <- which(values$nernaResults$OriginalAccession == selected_acc)
      if (length(idx) > 0) values$nernaResults$NegativeSequence[idx[1]] else NULL
    } else { NULL }
    
    dinu_seq <- if (!is.null(values$dinucleotideResults) && nrow(values$dinucleotideResults) > 0) {
      idx <- which(values$dinucleotideResults$OriginalAccession == selected_acc)
      if (length(idx) > 0) values$dinucleotideResults$NegativeSequence[idx[1]] else NULL
    } else { NULL }
    
    random_seq <- if (!is.null(values$randomResults) && nrow(values$randomResults) > 0) {
      idx <- which(values$randomResults$OriginalAccession == selected_acc)
      if (length(idx) > 0) values$randomResults$NegativeSequence[idx[1]] else NULL
    } else { NULL }
    
    # Get RNAfold results from cache or calculate if not available
    get_rnafold_results <- function(seq, acc, cache) {
      if (is.null(seq)) return(list(structure = NULL, mfe = NA))
      
      # First, check if result exists in cache
      if (!is.null(cache) && acc %in% names(cache)) {
        message("Using cached structure for: ", acc)
        cached_result <- cache[[acc]]
        return(list(
          structure = if (has_valid_structure(cached_result$structure, nchar(seq))) cached_result$structure else NULL,
          mfe = if (has_valid_mfe(cached_result$mfe)) cached_result$mfe else NA_real_
        ))
      }
      
      # If not in cache, calculate using RNAfold
      tryCatch({
        message("Running RNAfold for sequence: ", substr(seq, 1, 20), "... (not in cache)")
        result <- quick_rnafold(seq)
        valid_structure <- has_valid_structure(result$structure, nchar(seq))
        valid_mfe <- has_valid_mfe(result$mfe)
        message(
          "RNAfold result - Structure: ",
          if (valid_structure) substr(result$structure, 1, 20) else "NOT CALCULATED",
          ", MFE: ",
          if (valid_mfe) result$mfe else "NOT CALCULATED"
        )
        return(list(
          structure = if (valid_structure) result$structure else NULL,
          mfe = if (valid_mfe) result$mfe else NA_real_
        ))
      }, error = function(e) {
        message("RNAfold error for sequence: ", substr(seq, 1, 20), " - Error: ", e$message)
        return(list(
          structure = NULL,
          mfe = NA_real_
        ))
      })
    }
    
    # Get secondary structures and MFE values (from cache if available)
    nerna_estimates <- get_rnafold_results(nerna_seq, selected_acc, values$nernaStructures)
    dinu_estimates <- get_rnafold_results(dinu_seq, selected_acc, values$dinucStructures)
    random_estimates <- get_rnafold_results(random_seq, selected_acc, values$randomStructures)
    
    nerna_structure <- nerna_estimates$structure
    nerna_mfe <- nerna_estimates$mfe
    dinu_structure <- dinu_estimates$structure
    dinu_mfe <- dinu_estimates$mfe
    random_structure <- random_estimates$structure
    random_mfe <- random_estimates$mfe
    
    # Calculate statistics with updated secondary structures
    original_stats <- calculateSequenceStats(original_seq, original_structure)
    nerna_stats <- if (!is.null(nerna_seq)) calculateSequenceStats(nerna_seq, nerna_structure) else NULL
    dinu_stats <- if (!is.null(dinu_seq)) calculateSequenceStats(dinu_seq, dinu_structure) else NULL
    random_stats <- if (!is.null(random_seq)) calculateSequenceStats(random_seq, random_structure) else NULL
    
    # Create data for display
    values$selectedComparison <- list(
      Accession = original_acc,
      Original = list(
        Sequence = original_seq,
        Structure = original_structure,
        MFE = original_mfe,
        Stats = original_stats
      ),
      NeRNA = list(
        Sequence = nerna_seq,
        Structure = nerna_structure,
        MFE = nerna_mfe,
        Stats = nerna_stats
      ),
      Dinucleotide = list(
        Sequence = dinu_seq,
        Structure = dinu_structure,
        MFE = dinu_mfe,
        Stats = dinu_stats
      ),
      Random = list(
        Sequence = random_seq,
        Structure = random_structure,
        MFE = random_mfe,
        Stats = random_stats
      )
    )
  })
  
  # Display selected sequence comparison
  output$selectedSequenceComparison <- renderUI({
    req(values$selectedComparison)
    
    comp <- values$selectedComparison
    
    tagList(
      h4(paste("Selected Sequence:", comp$Accession)),
      
      div(style = "padding: 10px; background-color: #f8f9fa; border-radius: 5px; margin-bottom: 15px;",
        h5("Original Sequence:"),
        p(style = "word-break: break-all; font-family: monospace;", comp$Original$Sequence),
        if (has_valid_structure(comp$Original$Structure, nchar(comp$Original$Sequence))) {
          div(
            h5("RNAfold-predicted Original Structure:"),
            p(style = "word-break: break-all; font-family: monospace;", comp$Original$Structure),
            p(paste("RNAfold MFE:", if (has_valid_mfe(comp$Original$MFE)) comp$Original$MFE else "Not calculated"))
          )
        } else {
          p(style = "color: #856404; font-weight: bold;", "RNAfold structure and MFE: Not calculated")
        }
      ),
      
      if (!is.null(comp$NeRNA$Sequence)) {
        div(style = "padding: 10px; background-color: #d4edda; border-radius: 5px; margin-bottom: 15px;",
          h5("NeRNA Generated Sequence:"),
          p(style = "word-break: break-all; font-family: monospace;", comp$NeRNA$Sequence),
          if (has_valid_structure(comp$NeRNA$Structure, nchar(comp$NeRNA$Sequence))) {
            div(
              h5("RNAfold-predicted Secondary Structure:"),
              p(style = "word-break: break-all; font-family: monospace;", comp$NeRNA$Structure),
              p(paste("RNAfold MFE:", if (has_valid_mfe(comp$NeRNA$MFE)) comp$NeRNA$MFE else "Not calculated"))
            )
          } else {
            p(style = "color: #856404; font-weight: bold;", "RNAfold structure and MFE: Not calculated")
          }
        )
      },
      
      if (!is.null(comp$Dinucleotide$Sequence)) {
        div(style = "padding: 10px; background-color: #cce5ff; border-radius: 5px; margin-bottom: 15px;",
          h5("Dinucleotide Shuffling Generated Sequence:"),
          p(style = "word-break: break-all; font-family: monospace;", comp$Dinucleotide$Sequence),
          if (has_valid_structure(comp$Dinucleotide$Structure, nchar(comp$Dinucleotide$Sequence))) {
            div(
              h5("RNAfold-predicted Secondary Structure:"),
              p(style = "word-break: break-all; font-family: monospace;", comp$Dinucleotide$Structure),
              p(paste("RNAfold MFE:", if (has_valid_mfe(comp$Dinucleotide$MFE)) comp$Dinucleotide$MFE else "Not calculated"))
            )
          } else {
            p(style = "color: #856404; font-weight: bold;", "RNAfold structure and MFE: Not calculated")
          }
        )
      },
      
      if (!is.null(comp$Random$Sequence)) {
        div(style = "padding: 10px; background-color: #fff3cd; border-radius: 5px; margin-bottom: 15px;",
          h5("Random Shuffling Generated Sequence:"),
          p(style = "word-break: break-all; font-family: monospace;", comp$Random$Sequence),
          if (has_valid_structure(comp$Random$Structure, nchar(comp$Random$Sequence))) {
            div(
              h5("RNAfold-predicted Secondary Structure:"),
              p(style = "word-break: break-all; font-family: monospace;", comp$Random$Structure),
              p(paste("RNAfold MFE:", if (has_valid_mfe(comp$Random$MFE)) comp$Random$MFE else "Not calculated"))
            )
          } else {
            p(style = "color: #856404; font-weight: bold;", "RNAfold structure and MFE: Not calculated")
          }
        )
      }
    )
  })
  
  # Switch to Overview Comparison tab when the button is clicked
  observeEvent(input$showOverviewComparisonBtn, {
    updateTabItems(session, "sidebar", "overview_comparison")
    
    # Generate summary data if not already done
    if (is.null(values$comparisonSummary)) {
      generateMethodComparisonSummary()
    }
  })
  
  # Switch to Sequence Specific Comparison (Result Comparison) tab when button is clicked
  observeEvent(input$gotoSequenceComparisonBtn, {
    updateTabItems(session, "sidebar", "comparison")
  })
  
  # Run RNAfold for all sequences (Original + All Negatives)
  observeEvent(input$runAllRNAfoldBtn, {
    req(values$fastaData)
    
    # Show progress indicator (using JavaScript)
    session$sendCustomMessage("showElement", "rnafoldProgress")
    
    tryCatch({
      total_sequences <- 0
      completed_sequences <- 0
      
      # Calculate total number of sequences to process
      total_sequences <- nrow(values$fastaData)
      if (!is.null(values$nernaResults)) total_sequences <- total_sequences + nrow(values$nernaResults)
      if (!is.null(values$dinucleotideResults)) total_sequences <- total_sequences + nrow(values$dinucleotideResults)
      if (!is.null(values$randomResults)) total_sequences <- total_sequences + nrow(values$randomResults)
      
      showNotification(
        sprintf("Calculating secondary structures for %d sequences (Original + Negatives)...", total_sequences),
        type = "message",
        duration = NULL,
        id = "rnafoldAllNotif"
      )
      
      # Process Original sequences (if not already done)
      if (is.null(values$rnafoldResults)) {
        for (i in 1:nrow(values$fastaData)) {
          result <- quick_rnafold(values$fastaData$Sequence[i])
          # Store in rnafoldResults format
          if (i == 1) {
            values$rnafoldResults <- data.frame(
              Accession = values$fastaData$Accession[i],
              Sequence = values$fastaData$Sequence[i],
              SecondaryStructure = result$structure,
              MFE = result$mfe,
              stringsAsFactors = FALSE
            )
          } else {
            values$rnafoldResults <- rbind(
              values$rnafoldResults,
              data.frame(
                Accession = values$fastaData$Accession[i],
                Sequence = values$fastaData$Sequence[i],
                SecondaryStructure = result$structure,
                MFE = result$mfe,
                stringsAsFactors = FALSE
              )
            )
          }
          completed_sequences <- completed_sequences + 1
        }
      } else {
        completed_sequences <- completed_sequences + nrow(values$fastaData)
      }
      
      # Process NeRNA sequences
      if (!is.null(values$nernaResults)) {
        values$nernaStructures <- list()
        for (i in 1:nrow(values$nernaResults)) {
          seq <- values$nernaResults$NegativeSequence[i]
          acc <- values$nernaResults$OriginalAccession[i]
          result <- quick_rnafold(seq)
          values$nernaStructures[[acc]] <- result
          completed_sequences <- completed_sequences + 1
        }
      }
      
      # Process Dinucleotide sequences
      if (!is.null(values$dinucleotideResults)) {
        values$dinucStructures <- list()
        for (i in 1:nrow(values$dinucleotideResults)) {
          seq <- values$dinucleotideResults$NegativeSequence[i]
          acc <- values$dinucleotideResults$OriginalAccession[i]
          result <- quick_rnafold(seq)
          values$dinucStructures[[acc]] <- result
          completed_sequences <- completed_sequences + 1
        }
      }
      
      # Process Random sequences
      if (!is.null(values$randomResults)) {
        values$randomStructures <- list()
        for (i in 1:nrow(values$randomResults)) {
          seq <- values$randomResults$NegativeSequence[i]
          acc <- values$randomResults$OriginalAccession[i]
          result <- quick_rnafold(seq)
          values$randomStructures[[acc]] <- result
          completed_sequences <- completed_sequences + 1
        }
      }
      
      # Hide progress (using JavaScript)
      session$sendCustomMessage("hideElement", "rnafoldProgress")
      
      # Remove processing notification
      removeNotification("rnafoldAllNotif")
      
      count_valid_cache_results <- function(cache, sequences) {
        if (is.null(cache) || length(cache) == 0) return(0L)
        sum(vapply(seq_along(cache), function(i) {
          result <- cache[[i]]
          has_valid_structure(result$structure, nchar(sequences[i])) &&
            has_valid_mfe(result$mfe)
        }, logical(1)))
      }

      original_valid <- if (!is.null(values$rnafoldResults)) {
        sum(vapply(seq_len(nrow(values$rnafoldResults)), function(i) {
          has_valid_structure(
            values$rnafoldResults$SecondaryStructure[i],
            nchar(values$rnafoldResults$Sequence[i])
          ) && has_valid_mfe(values$rnafoldResults$MFE[i])
        }, logical(1)))
      } else 0L

      nerna_valid <- if (!is.null(values$nernaResults)) {
        count_valid_cache_results(values$nernaStructures, values$nernaResults$NegativeSequence)
      } else 0L
      dinuc_valid <- if (!is.null(values$dinucleotideResults)) {
        count_valid_cache_results(values$dinucStructures, values$dinucleotideResults$NegativeSequence)
      } else 0L
      random_valid <- if (!is.null(values$randomResults)) {
        count_valid_cache_results(values$randomStructures, values$randomResults$NegativeSequence)
      } else 0L

      valid_results <- original_valid + nerna_valid + dinuc_valid + random_valid
      failed_results <- total_sequences - valid_results

      showNotification(
        if (failed_results == 0) {
          sprintf("RNAfold successfully calculated %d of %d sequences.", valid_results, total_sequences)
        } else {
          sprintf(
            "RNAfold calculated %d of %d sequences; %d results are marked Not calculated.",
            valid_results, total_sequences, failed_results
          )
        },
        type = if (failed_results == 0) "message" else "warning",
        duration = if (failed_results == 0) 5 else 10
      )
      
      # Update comparison summary with fresh data
      values$comparisonSummary <- generateMethodComparisonSummary()
      
      # Update status output
      output$rnafoldAllStatus <- renderUI({
        tags$div(
          style = "background-color: #d4edda; border-left: 4px solid #28a745; padding: 10px; margin-top: 10px;",
          tags$p(style = "color: #155724; margin: 0;",
            icon("check-circle"),
            tags$strong(sprintf(" Valid RNAfold results: %d of %d sequences", valid_results, total_sequences)),
            tags$br(),
            tags$small(
              "Original: ", original_valid,
              " | NeRNA: ", nerna_valid,
              " | Dinucleotide: ", dinuc_valid,
              " | Random: ", random_valid,
              if (failed_results > 0) paste0(" | Not calculated: ", failed_results) else ""
            )
          )
        )
      })
      
    }, error = function(e) {
      session$sendCustomMessage("hideElement", "rnafoldProgress")
      removeNotification("rnafoldAllNotif")
      showNotification(
        paste("Error calculating structures:", e$message),
        type = "error",
        duration = 10
      )
    })
  })
  
  # Navigate to Upload tab when Get Started button is clicked
  observeEvent(input$getStartedBtn, {
    updateTabItems(session, "sidebar", "upload")
  })
  
  # Navigate to Upload tab from other pages when no data is loaded
  observeEvent(input$goToUploadBtn, {
    updateTabItems(session, "sidebar", "upload")
  })
  
  observeEvent(input$goToUploadBtn2, {
    updateTabItems(session, "sidebar", "upload")
  })
  
  observeEvent(input$goToUploadBtn3, {
    updateTabItems(session, "sidebar", "upload")
  })
  
  observeEvent(input$goToUploadBtn4, {
    updateTabItems(session, "sidebar", "upload")
  })
  
  observeEvent(input$goToUploadBtn5, {
    updateTabItems(session, "sidebar", "upload")
  })
  
  # Helper function to get data for visualizations
  getComparisonData <- function() {
    req(values$selectedComparison)
    comp <- values$selectedComparison
    
    # MFE values returned by RNAfold. Missing values remain NA and are never
    # replaced with fixed examples or estimates.
    mfe_values <- c(
      if (has_valid_mfe(comp$Original$MFE)) comp$Original$MFE else NA_real_,
      if (has_valid_mfe(comp$NeRNA$MFE)) comp$NeRNA$MFE else NA_real_,
      if (has_valid_mfe(comp$Dinucleotide$MFE)) comp$Dinucleotide$MFE else NA_real_,
      if (has_valid_mfe(comp$Random$MFE)) comp$Random$MFE else NA_real_
    )
    
    # Sequence lengths
    length_values <- c(
      comp$Original$Stats$Length,
      if (!is.null(comp$NeRNA$Stats)) comp$NeRNA$Stats$Length else NA_real_,
      if (!is.null(comp$Dinucleotide$Stats)) comp$Dinucleotide$Stats$Length else NA_real_,
      if (!is.null(comp$Random$Stats)) comp$Random$Stats$Length else NA_real_
    )
    
    # Paired and unpaired bases are calculated only from real RNAfold
    # dot-bracket structures. No percentage-based fallback is used.
    paired_values <- c(
      if (has_valid_structure(comp$Original$Structure, nchar(comp$Original$Sequence))) comp$Original$Stats$PairedBases else NA_real_,
      if (!is.null(comp$NeRNA$Stats) && has_valid_structure(comp$NeRNA$Structure, nchar(comp$NeRNA$Sequence))) comp$NeRNA$Stats$PairedBases else NA_real_,
      if (!is.null(comp$Dinucleotide$Stats) && has_valid_structure(comp$Dinucleotide$Structure, nchar(comp$Dinucleotide$Sequence))) comp$Dinucleotide$Stats$PairedBases else NA_real_,
      if (!is.null(comp$Random$Stats) && has_valid_structure(comp$Random$Structure, nchar(comp$Random$Sequence))) comp$Random$Stats$PairedBases else NA_real_
    )
    
    unpaired_values <- c(
      if (!is.na(paired_values[1])) comp$Original$Stats$UnpairedBases else NA_real_,
      if (!is.na(paired_values[2])) comp$NeRNA$Stats$UnpairedBases else NA_real_,
      if (!is.na(paired_values[3])) comp$Dinucleotide$Stats$UnpairedBases else NA_real_,
      if (!is.na(paired_values[4])) comp$Random$Stats$UnpairedBases else NA_real_
    )
    
    # GC content values
    gc_values <- c(
      comp$Original$Stats$GCContent,
      if (!is.null(comp$NeRNA$Stats)) comp$NeRNA$Stats$GCContent else NA_real_,
      if (!is.null(comp$Dinucleotide$Stats)) comp$Dinucleotide$Stats$GCContent else NA_real_,
      if (!is.null(comp$Random$Stats)) comp$Random$Stats$GCContent else NA_real_
    )
    
    # Calculate pairing ratio (percentage of paired bases)
    pairing_ratio <- 100 * paired_values / (paired_values + unpaired_values)
    
    # Method names
    methods <- c("Original", "NeRNA", "Dinucleotide", "Random")
    
    # Colors for consistent visualization with lower alpha
    colors <- c("#17a2b8", "#28a745", "#007bff", "#ffc107")
    
    return(list(
      methods = methods,
      colors = colors,
      mfe = mfe_values,
      length = length_values,
      paired = paired_values,
      unpaired = unpaired_values,
      gc = gc_values,
      pairing_ratio = pairing_ratio
    ))
  }
  
  # Radar plot visualization with reduced alpha values
  output$radarVisualization <- renderPlot({
    data <- getComparisonData()

    complete_methods <- !is.na(data$mfe) & !is.na(data$length) &
      !is.na(data$paired) & !is.na(data$gc)

    if (!any(complete_methods)) {
      plot.new()
      title("RNA Sequence Features")
      text(0.5, 0.55, "Not calculated", cex = 1.2, font = 2)
      text(0.5, 0.45, "Run RNAfold to obtain MFE and pairing values.", cex = 0.9)
      return(invisible(NULL))
    }

    data$methods <- data$methods[complete_methods]
    data$colors <- data$colors[complete_methods]
    data$mfe <- data$mfe[complete_methods]
    data$length <- data$length[complete_methods]
    data$paired <- data$paired[complete_methods]
    data$gc <- data$gc[complete_methods]
    
    # Normalize the data for radar plot
    norm_mfe <- abs(data$mfe) / max(abs(data$mfe))
    norm_length <- data$length / max(data$length)
    norm_paired <- data$paired / max(data$paired)
    norm_gc <- data$gc / 100
    
    # Check if fmsb package is installed, if not install it
    if (!requireNamespace("fmsb", quietly = TRUE)) {
      install.packages("fmsb")
      library(fmsb)
    } else {
      library(fmsb)
    }
    
    # Create data frame for radar plot
    radar_data <- data.frame(
      MFE = norm_mfe,
      Length = norm_length,
      Paired_Bases = norm_paired,
      GC_Content = norm_gc
    )
    
    # Add max and min rows required by fmsb
    radar_data <- rbind(rep(1, 4), rep(0, 4), radar_data)
    
    # Set row names
    rownames(radar_data) <- c("Max", "Min", data$methods)
    
    # Colors for each method with lower alpha values
    radar_colors <- data$colors
    
    # Plot parameters
    par(mar = c(1, 1, 2, 6))
    
    # Create radar plot
    radarchart(
      radar_data,
      axistype = 1,
      pcol = radar_colors,
      pfcol = sapply(radar_colors, function(x) adjustcolor(x, alpha.f = 0.3)),
      plwd = 2,
      cglcol = "grey",
      cglty = 1,
      axislabcol = "grey30",
      caxislabels = seq(0, 1, 0.25),
      cglwd = 0.8,
      vlcex = 0.9,
      title = "RNA Sequence Features"
    )
    
    # Add legend
    legend(
      "topright",
      legend = data$methods,
      col = radar_colors,
      lwd = 2,
      bty = "n",
      cex = 0.8,
      inset = c(-0.1, 0),
      xpd = TRUE
    )
  })
  
  # Spider chart visualization - alternative view
  output$spiderVisualization <- renderPlot({
    data <- getComparisonData()
    
    # Calculate nucleotide counts for all sequences
    calculateBaseCounts <- function(sequences) {
      base_counts <- list(A = 0, G = 0, U = 0, C = 0)

      if (!is.null(sequences)) {
        for (seq in sequences) {
          # Count each nucleotide
          a_count <- sum(strsplit(seq, "")[[1]] == "A")
          g_count <- sum(strsplit(seq, "")[[1]] == "G")
          u_count <- sum(strsplit(seq, "")[[1]] == "U") + sum(strsplit(seq, "")[[1]] == "T")  # Count U and T
          c_count <- sum(strsplit(seq, "")[[1]] == "C")
          
          # Add to totals
          base_counts$A <- base_counts$A + a_count
          base_counts$G <- base_counts$G + g_count
          base_counts$U <- base_counts$U + u_count
          base_counts$C <- base_counts$C + c_count
        }
      }
      
      return(base_counts)
    }
    
    # Get base counts for each method
    original_counts <- calculateBaseCounts(values$fastaData$Sequence)
    nerna_counts <- calculateBaseCounts(if(!is.null(values$nernaResults)) values$nernaResults$NegativeSequence else NULL)
    dinu_counts <- calculateBaseCounts(if(!is.null(values$dinucleotideResults)) values$dinucleotideResults$NegativeSequence else NULL)
    random_counts <- calculateBaseCounts(if(!is.null(values$randomResults)) values$randomResults$NegativeSequence else NULL)
    
    # Create matrix for grouped bar chart
    counts_matrix <- matrix(
      c(
        original_counts$A, nerna_counts$A, dinu_counts$A, random_counts$A,
        original_counts$G, nerna_counts$G, dinu_counts$G, random_counts$G,
        original_counts$U, nerna_counts$U, dinu_counts$U, random_counts$U,
        original_counts$C, nerna_counts$C, dinu_counts$C, random_counts$C
      ),
      nrow = 4,
      byrow = TRUE
    )
    
    # Row and column names
    rownames(counts_matrix) <- c("A", "G", "U", "C")
    colnames(counts_matrix) <- data$methods
    
    # Colors for nucleotides
    nuc_colors <- c("darkgreen", "blue", "red", "orange")
    
    # Plot parameters
    par(mar = c(5, 5, 4, 10))  # Increase right margin for legend
    
    # Create grouped bar chart
    barplot(
      counts_matrix,
      beside = TRUE,
      col = nuc_colors,
      main = "Nucleotide Counts Comparison",
      xlab = "Method",
      ylab = "Number of Nucleotides",
      legend.text = rownames(counts_matrix),
      args.legend = list(x = "right", inset = c(-0.20, 0), xpd = TRUE),
      ylim = c(0, max(counts_matrix) * 1.1),
      las = 1
    )
    
    # Add grid lines
    grid(nx = NA, ny = NULL, lty = 2, col = "lightgray")
    
    # Add explanation
    mtext("Comparison of nucleotide counts across different RNA generation methods", 
          side = 3, line = 0.5, cex = 0.8)
  })
  
  # MFE visualization
  output$mfeVisualization <- renderPlot({
    data <- getComparisonData()

    available <- !is.na(data$mfe)
    plot_values <- ifelse(available, abs(data$mfe), 0)
    y_max <- if (any(available)) max(plot_values[available]) * 1.25 else 1
    if (!is.finite(y_max) || y_max <= 0) y_max <- 1
    
    par(mar = c(6, 4, 4, 2) + 0.1)
    bp <- barplot(
      plot_values,
      names.arg = data$methods,
      col = ifelse(available, data$colors, "gray85"),
      main = "Minimum Free Energy (MFE)",
      ylab = "Absolute MFE Value",
      border = "white",
      ylim = c(0, y_max),
      las = 2
    )
    
    text(
      bp,
      ifelse(available, plot_values + y_max * 0.04, y_max * 0.05),
      labels = ifelse(available, format(data$mfe, trim = TRUE), "Not calculated"),
      cex = 0.85
    )
    
    # Add grid lines
    grid(nx = NA, ny = NULL, lty = 2, col = "gray80")
  })
  
  # Sequence length visualization
  output$lengthVisualization <- renderPlot({
    data <- getComparisonData()
    
    par(mar = c(6, 4, 4, 2) + 0.1)
    bp <- barplot(
      data$length, 
      names.arg = data$methods,
      col = data$colors,
      main = "Sequence Length",
      ylab = "Number of Nucleotides",
      border = "white",
      ylim = c(0, max(data$length) * 1.2),
      las = 2
    )
    
    # Add length values on top of bars
    text(
      bp,
      data$length + max(data$length) * 0.05,
      labels = data$length,
      cex = 0.9
    )
    
    # Add grid lines
    grid(nx = NA, ny = NULL, lty = 2, col = "gray80")
  })
  
  # Base pairing visualization
  output$pairingVisualization <- renderPlot({
    data <- getComparisonData()

    available <- !is.na(data$paired) & !is.na(data$unpaired)
    
    # Create stacked bar data
    stacked_data <- rbind(
      ifelse(available, data$paired, 0),
      ifelse(available, data$unpaired, 0)
    )
    rownames(stacked_data) <- c("Paired", "Unpaired")

    totals <- colSums(stacked_data)
    y_max <- if (any(available)) max(totals[available]) * 1.2 else 1
    if (!is.finite(y_max) || y_max <= 0) y_max <- 1
    
    par(mar = c(6, 4, 4, 8) + 0.1)
    bp <- barplot(
      stacked_data, 
      names.arg = data$methods,
      col = c("#28a745", "#dc3545"),
      main = "Base Pairing Distribution",
      ylab = "Number of Bases",
      border = "white",
      ylim = c(0, y_max),
      las = 2,
      legend.text = rownames(stacked_data),
      args.legend = list(
        x = "topright", 
        bty = "n", 
        cex = 0.8,
        inset = c(-0.15, 0),
        xpd = TRUE
      )
    )
    
    # Add percentage labels for paired bases
    paired_percent <- ifelse(
      available,
      round(data$paired / (data$paired + data$unpaired) * 100, 1),
      NA_real_
    )
    text(
      bp,
      ifelse(available, pmax(data$paired / 2, y_max * 0.05), y_max * 0.05),
      labels = ifelse(available, paste0(paired_percent, "%"), "Not calculated"),
      cex = 0.85
    )
  })
  
  # GC content visualization
  output$gcVisualization <- renderPlot({
    data <- getComparisonData()
    
    # Create GC/AT stacked data
    gc_data <- data$gc
    at_data <- 100 - gc_data
    stacked_data <- rbind(gc_data, at_data)
    rownames(stacked_data) <- c("GC Content", "AT Content")
    
    par(mar = c(6, 4, 4, 8) + 0.1)
    bp <- barplot(
      stacked_data, 
      names.arg = data$methods,
      col = c("#17a2b8", "#fd7e14"),
      main = "GC Content Distribution",
      ylab = "Percentage (%)",
      border = "white",
      ylim = c(0, 120),
      las = 2,
      legend.text = rownames(stacked_data),
      args.legend = list(
        x = "topright", 
        bty = "n", 
        cex = 0.8,
        inset = c(-0.15, 0),
        xpd = TRUE
      )
    )
    
    # Add percentage labels for GC content
    text(
      bp,
      gc_data / 2,
      labels = paste0(gc_data, "%"),
      cex = 0.9
    )
  })
  
  # Generate comprehensive method comparison summary
  generateMethodComparisonSummary <- function() {
    if (is.null(values$fastaData)) return(NULL)

    methods <- c("Original", "NeRNA", "Dinucleotide", "Random")
    metrics <- c("Avg. Sequence Length",
                 "Avg. GC Content (%)",
                 "Avg. MFE",
                 "Avg. Paired Bases",
                 "Avg. Unpaired Bases",
                 "Pairing Ratio (%)")

    summary_df <- data.frame(
      Metric = metrics,
      Original = numeric(length(metrics)),
      NeRNA = numeric(length(metrics)),
      Dinucleotide = numeric(length(metrics)),
      Random = numeric(length(metrics))
    )

    summarize_rnafold_outputs <- function(structures, mfes, sequences) {
      empty_summary <- list(
        avg_mfe = NA_real_,
        avg_paired = NA_real_,
        avg_unpaired = NA_real_,
        avg_pairing_ratio = NA_real_
      )

      if (is.null(structures) || is.null(mfes) ||
          length(structures) != length(sequences) ||
          length(mfes) != length(sequences) || length(sequences) == 0) {
        return(empty_summary)
      }

      valid_structures <- vapply(seq_along(sequences), function(i) {
        has_valid_structure(structures[i], nchar(sequences[i]))
      }, logical(1))

      if (!all(valid_structures)) return(empty_summary)

      paired_counts <- vapply(structures, function(structure) {
        sum(strsplit(structure, "", fixed = TRUE)[[1]] %in% c("(", ")"))
      }, numeric(1))
      total_lengths <- nchar(structures)
      unpaired_counts <- total_lengths - paired_counts

      list(
        avg_mfe = if (all(vapply(mfes, has_valid_mfe, logical(1)))) round(mean(mfes), 2) else NA_real_,
        avg_paired = round(mean(paired_counts), 1),
        avg_unpaired = round(mean(unpaired_counts), 1),
        avg_pairing_ratio = round(mean(paired_counts / total_lengths * 100), 1)
      )
    }

    summarize_structure_cache <- function(cache, sequences) {
      if (is.null(cache) || length(cache) != length(sequences)) {
        return(summarize_rnafold_outputs(NULL, NULL, sequences))
      }

      structures <- vapply(cache, function(result) {
        if (!is.null(result$structure) && length(result$structure) == 1) result$structure else NA_character_
      }, character(1))
      mfes <- vapply(cache, function(result) {
        if (has_valid_mfe(result$mfe)) as.numeric(result$mfe) else NA_real_
      }, numeric(1))

      summarize_rnafold_outputs(structures, mfes, sequences)
    }

    ## --- Original ---
    original_seqs <- values$fastaData$Sequence
    summary_df$Original[1] <- round(mean(nchar(original_seqs)), 1)
    gc_percentages <- sapply(original_seqs, function(seq) {
      gc_count <- sum(gregexpr("[GC]", seq)[[1]] > 0)
      (gc_count / nchar(seq)) * 100
    })
    summary_df$Original[2] <- round(mean(gc_percentages), 1)

    if (!is.null(values$rnafoldResults)) {
      original_structure_summary <- summarize_rnafold_outputs(
        values$rnafoldResults$SecondaryStructure,
        values$rnafoldResults$MFE,
        values$rnafoldResults$Sequence
      )
      summary_df$Original[3] <- original_structure_summary$avg_mfe
      summary_df$Original[4] <- original_structure_summary$avg_paired
      summary_df$Original[5] <- original_structure_summary$avg_unpaired
      summary_df$Original[6] <- original_structure_summary$avg_pairing_ratio
    } else {
      summary_df$Original[3:6] <- NA
    }

    ## helper: ortalama GC ve uzunluk
    avg_gc <- function(seqs) {
      if (is.null(seqs) || length(seqs) == 0) return(NA_real_)
      round(mean(sapply(seqs, function(seq) {
        gc_count <- sum(gregexpr("[GC]", seq)[[1]] > 0)
        (gc_count / nchar(seq)) * 100
      })), 1)
    }

    ## --- NeRNA ---
    if (!is.null(values$nernaResults)) {
      nerna_seqs <- values$nernaResults$NegativeSequence
      summary_df$NeRNA[1] <- round(mean(nchar(nerna_seqs)), 1)
      summary_df$NeRNA[2] <- avg_gc(nerna_seqs)

      nerna_structure_summary <- summarize_structure_cache(values$nernaStructures, nerna_seqs)
      summary_df$NeRNA[3] <- nerna_structure_summary$avg_mfe
      summary_df$NeRNA[4] <- nerna_structure_summary$avg_paired
      summary_df$NeRNA[5] <- nerna_structure_summary$avg_unpaired
      summary_df$NeRNA[6] <- nerna_structure_summary$avg_pairing_ratio
    } else {
      summary_df$NeRNA[1:6] <- NA
    }

    ## --- Dinucleotide ---
    if (!is.null(values$dinucleotideResults)) {
      dinu_seqs <- values$dinucleotideResults$NegativeSequence
      summary_df$Dinucleotide[1] <- round(mean(nchar(dinu_seqs)), 1)
      summary_df$Dinucleotide[2] <- avg_gc(dinu_seqs)

      dinuc_structure_summary <- summarize_structure_cache(values$dinucStructures, dinu_seqs)
      summary_df$Dinucleotide[3] <- dinuc_structure_summary$avg_mfe
      summary_df$Dinucleotide[4] <- dinuc_structure_summary$avg_paired
      summary_df$Dinucleotide[5] <- dinuc_structure_summary$avg_unpaired
      summary_df$Dinucleotide[6] <- dinuc_structure_summary$avg_pairing_ratio
    } else {
      summary_df$Dinucleotide[1:6] <- NA
    }

    ## --- Random ---
    if (!is.null(values$randomResults)) {
      random_seqs <- values$randomResults$NegativeSequence
      summary_df$Random[1] <- round(mean(nchar(random_seqs)), 1)
      summary_df$Random[2] <- avg_gc(random_seqs)

      random_structure_summary <- summarize_structure_cache(values$randomStructures, random_seqs)
      summary_df$Random[3] <- random_structure_summary$avg_mfe
      summary_df$Random[4] <- random_structure_summary$avg_paired
      summary_df$Random[5] <- random_structure_summary$avg_unpaired
      summary_df$Random[6] <- random_structure_summary$avg_pairing_ratio
    } else {
      summary_df$Random[1:6] <- NA
    }

    values$comparisonSummary <- summary_df
    return(summary_df)
  }
  
  # Render comparison summary table
  output$comparisonSummaryTable <- renderDT({
    if (is.null(values$comparisonSummary)) {
      summary_df <- generateMethodComparisonSummary()
    } else {
      summary_df <- values$comparisonSummary
    }
    
    req(summary_df)
    
    # Format table
    datatable(
      summary_df,
      options = list(
        pageLength = 10,
        dom = 'ft',
        ordering = TRUE
      ),
      rownames = FALSE
    ) %>%
    formatStyle(
      columns = c("Original", "NeRNA", "Dinucleotide", "Random"),
      backgroundColor = styleEqual(c(NA), c("lightgray"))
    )
  })
  
  # Generate comparison barchart visualization
  output$methodComparisonPlot <- renderPlot({
    req(values$comparisonSummary)
    
    summary_df <- values$comparisonSummary
    
    # Create data for barplot visualization
    # We'll visualize effectiveness score for each method
    # Higher score = better negative RNA model
    
    # Effectiveness parameters:
    # 1. MFE difference from original (normalized)
    # 2. Pairing ratio difference from original (normalized)
    # 3. GC content preservation (normalized)
    
    # Calculate normalized metrics
    original_mfe <- abs(summary_df$Original[3])
    original_pairing <- summary_df$Original[6]
    original_gc <- summary_df$Original[2]
    
    # Create scoring metrics 
    mfe_scores <- numeric(3)
    pairing_scores <- numeric(3)
    gc_scores <- numeric(3)
    
    # Calculate scores for each method if data exists
    if (!is.na(original_mfe) && !is.na(original_pairing) && !is.na(original_gc)) {
      # NeRNA
      if (!is.na(summary_df$NeRNA[3])) {
        mfe_diff <- abs((abs(summary_df$NeRNA[3]) - original_mfe) / original_mfe)
        mfe_scores[1] <- 1 - min(mfe_diff, 1) # Lower difference = higher score
        
        pairing_diff <- abs(summary_df$NeRNA[6] - original_pairing) / original_pairing
        pairing_scores[1] <- 1 - min(pairing_diff, 1)
        
        gc_diff <- abs(summary_df$NeRNA[2] - original_gc) / original_gc
        gc_scores[1] <- 1 - min(gc_diff, 1)
      }
      
      # Dinucleotide
      if (!is.na(summary_df$Dinucleotide[3])) {
        mfe_diff <- abs((abs(summary_df$Dinucleotide[3]) - original_mfe) / original_mfe)
        mfe_scores[2] <- 1 - min(mfe_diff, 1)
        
        pairing_diff <- abs(summary_df$Dinucleotide[6] - original_pairing) / original_pairing
        pairing_scores[2] <- 1 - min(pairing_diff, 1)
        
        gc_diff <- abs(summary_df$Dinucleotide[2] - original_gc) / original_gc
        gc_scores[2] <- 1 - min(gc_diff, 1)
      }
      
      # Random
      if (!is.na(summary_df$Random[3])) {
        mfe_diff <- abs((abs(summary_df$Random[3]) - original_mfe) / original_mfe)
        mfe_scores[3] <- 1 - min(mfe_diff, 1)
        
        pairing_diff <- abs(summary_df$Random[6] - original_pairing) / original_pairing
        pairing_scores[3] <- 1 - min(pairing_diff, 1)
        
        gc_diff <- abs(summary_df$Random[2] - original_gc) / original_gc
        gc_scores[3] <- 1 - min(gc_diff, 1)
      }
    }
    
    # Scale scores to percentages (0-100)
    mfe_scores <- mfe_scores * 40       # MFE counts for 40% of total
    pairing_scores <- pairing_scores * 30  # Pairing counts for 30% of total
    gc_scores <- gc_scores * 30         # GC content counts for 30% of total
    
    # Create matrix for stacked bar chart
    stacked_data <- rbind(mfe_scores, pairing_scores, gc_scores)
    rownames(stacked_data) <- c("MFE Similarity", "Pairing Ratio", "GC Content")
    
    # Create plot
    methods <- c("NeRNA", "Dinucleotide", "Random")
    
    # Colors for components
    component_colors <- c("#fd7e14", "#17a2b8", "#6f42c1")
    
    # Plot data
    par(mar = c(6, 4, 4, 8) + 0.1)
    
    # Create stacked barplot
    bp <- barplot(
      stacked_data,
      names.arg = methods,
      col = component_colors,
      main = "Method Effectiveness Score",
      ylab = "Score (0-100)",
      ylim = c(0, 100),
      las = 2,
      border = NA,
      legend.text = rownames(stacked_data),
      args.legend = list(
        x = "topright", 
        bty = "n", 
        cex = 0.8,
        inset = c(-0.15, 0),
        xpd = TRUE
      )
    )
    
    # Add total score labels on top of bars
    total_scores <- colSums(stacked_data)
    text(
      bp,
      total_scores + 3,
      labels = paste0("Total: ", round(total_scores, 1)),
      cex = 0.9,
      font = 2
    )
    
    # Add grid lines
    abline(h = seq(0, 100, by = 20), lty = 3, col = "gray80")
  })
  
  # Download handlers for comparison summary
  output$downloadComparisonCSV <- downloadHandler(
    filename = function() {
      paste("nerna-method-comparison-summary-", format(Sys.Date(), "%Y%m%d"), ".csv", sep = "")
    },
    content = function(file) {
      if (is.null(values$comparisonSummary)) {
        values$comparisonSummary <- generateMethodComparisonSummary()
      }
      write.csv(values$comparisonSummary, file, row.names = FALSE)
    }
  )
  
  output$downloadComparisonExcel <- downloadHandler(
    filename = function() {
      paste("nerna-method-comparison-summary-", format(Sys.Date(), "%Y%m%d"), ".xlsx", sep = "")
    },
    content = function(file) {
      if (is.null(values$comparisonSummary)) {
        values$comparisonSummary <- generateMethodComparisonSummary()
      }
      
      if (!requireNamespace("openxlsx", quietly = TRUE)) {
        install.packages("openxlsx")
      }
      library(openxlsx)
      
      wb <- createWorkbook()
      addWorksheet(wb, "Method Comparison")
      writeData(wb, "Method Comparison", values$comparisonSummary)
      
      # Add styling
      headerStyle <- createStyle(
        fontColour = "#FFFFFF", fgFill = "#4472C4",
        halign = "center", valign = "center",
        textDecoration = "bold"
      )
      
      addStyle(wb, "Method Comparison", headerStyle, rows = 1, cols = 1:5, gridExpand = TRUE)
      
      saveWorkbook(wb, file, overwrite = TRUE)
    }
  )
  
  # ============================================================
  # REAL-TIME STRUCTURE VIEWER
  # ============================================================
  
  # Update structure viewer dropdown when RNAfold results are available
  observe({
    if (!is.null(values$rnafoldResults) && nrow(values$rnafoldResults) > 0) {
      choices <- setNames(values$rnafoldResults$Accession, values$rnafoldResults$Accession)
      updateSelectInput(session, "structureViewerSeq", choices = choices)
    }
  })
  
  # Render structure viewer (text)
  output$structureViewer <- renderText({
    req(input$structureViewerSeq)
    req(values$rnafoldResults)
    
    # Find selected sequence
    idx <- which(values$rnafoldResults$Accession == input$structureViewerSeq)
    
    if (length(idx) > 0) {
      seq <- values$rnafoldResults$Sequence[idx[1]]
      struct <- values$rnafoldResults$SecondaryStructure[idx[1]]
      mfe <- values$rnafoldResults$MFE[idx[1]]

      if (has_valid_structure(struct, nchar(seq))) {
        render_structure_text(seq, struct, if (has_valid_mfe(mfe)) mfe else NULL)
      } else {
        "RNAfold structure and MFE: Not calculated"
      }
    } else {
      "No structure available for selected sequence."
    }
  })
  
  # Render 3D structure plot
  output$structure3DPlot <- renderPlotly({
    req(input$structureViewerSeq)
    req(values$rnafoldResults)
    
    # Find selected sequence
    idx <- which(values$rnafoldResults$Accession == input$structureViewerSeq)
    
    if (length(idx) > 0) {
      seq <- values$rnafoldResults$Sequence[idx[1]]
      struct <- values$rnafoldResults$SecondaryStructure[idx[1]]
      mfe <- values$rnafoldResults$MFE[idx[1]]

      if (has_valid_structure(struct, nchar(seq))) {
        create_3d_structure_plot(seq, struct, if (has_valid_mfe(mfe)) mfe else NULL)
      } else {
        NULL
      }
    } else {
      NULL
    }
  })
  
  # ============================================================
  # 4-WAY STRUCTURE COMPARISON VIEWER
  # ============================================================
  
  # Render Original structure (4-way comparison)
  output$comp4OriginalPlot <- renderPlotly({
    req(input$comparisonSequence)
    req(values$rnafoldResults)
    req(input$compareBtn)
    
    # Use isolate to prevent reactive dependency on selectedComparison
    idx <- which(values$rnafoldResults$Accession == input$comparisonSequence)
    
    if (length(idx) > 0) {
      seq <- values$rnafoldResults$Sequence[idx[1]]
      struct <- values$rnafoldResults$SecondaryStructure[idx[1]]
      mfe <- values$rnafoldResults$MFE[idx[1]]

      if (has_valid_structure(struct, nchar(seq))) {
        create_3d_structure_plot(seq, struct, if (has_valid_mfe(mfe)) mfe else NULL)
      } else {
        plotly_empty()
      }
    } else {
      plotly_empty()
    }
  })
  
  output$comp4OriginalText <- renderText({
    req(input$comparisonSequence)
    req(values$rnafoldResults)
    req(input$compareBtn)
    
    idx <- which(values$rnafoldResults$Accession == input$comparisonSequence)
    
    if (length(idx) > 0) {
      seq <- values$rnafoldResults$Sequence[idx[1]]
      struct <- values$rnafoldResults$SecondaryStructure[idx[1]]
      mfe <- values$rnafoldResults$MFE[idx[1]]

      if (has_valid_structure(struct, nchar(seq))) {
        paired <- sum(strsplit(struct, "", fixed = TRUE)[[1]] %in% c("(", ")"))
        total <- nchar(struct)

        paste0("Structure: ", struct, "\n",
               "MFE: ", if (has_valid_mfe(mfe)) paste(mfe, "kcal/mol") else "Not calculated", "\n",
               "Paired: ", paired, "/", total, " (", round(paired/total*100, 1), "%)")
      } else {
        "RNAfold structure and MFE: Not calculated"
      }
    } else {
      "No structure data available"
    }
  })
  
  # Render NeRNA structure (4-way comparison)
  output$comp4NeRNAPlot <- renderPlotly({
    req(input$compareBtn)
    req(values$selectedComparison)
    
    comp <- values$selectedComparison
    if (!is.null(comp$NeRNA$Sequence) &&
        has_valid_structure(comp$NeRNA$Structure, nchar(comp$NeRNA$Sequence))) {
      create_3d_structure_plot(comp$NeRNA$Sequence, comp$NeRNA$Structure, comp$NeRNA$MFE)
    } else {
      NULL
    }
  })
  
  output$comp4NeRNAText <- renderText({
    req(input$compareBtn)
    req(values$selectedComparison)
    
    comp <- values$selectedComparison
    if (is.null(comp$NeRNA$Sequence)) {
      "NeRNA sequence not generated yet"
    } else if (has_valid_structure(comp$NeRNA$Structure, nchar(comp$NeRNA$Sequence))) {
      structure <- comp$NeRNA$Structure
      mfe <- comp$NeRNA$MFE
      paired <- sum(strsplit(structure, "")[[1]] %in% c("(", ")"))
      total <- nchar(structure)
      
      paste0("Structure: ", structure, "\n",
             "MFE: ", if (has_valid_mfe(mfe)) paste(mfe, "kcal/mol") else "Not calculated", "\n",
             "Paired: ", paired, "/", total, " (", round(paired/total*100, 1), "%)")
    } else {
      "RNAfold structure and MFE: Not calculated"
    }
  })
  
  # Render Dinucleotide structure (4-way comparison)
  output$comp4DinucPlot <- renderPlotly({
    req(input$compareBtn)
    req(values$selectedComparison)
    
    comp <- values$selectedComparison
    if (!is.null(comp$Dinucleotide$Sequence) &&
        has_valid_structure(comp$Dinucleotide$Structure, nchar(comp$Dinucleotide$Sequence))) {
      create_3d_structure_plot(comp$Dinucleotide$Sequence, comp$Dinucleotide$Structure, comp$Dinucleotide$MFE)
    } else {
      NULL
    }
  })
  
  output$comp4DinucText <- renderText({
    req(input$compareBtn)
    req(values$selectedComparison)
    
    comp <- values$selectedComparison
    if (is.null(comp$Dinucleotide$Sequence)) {
      "Dinucleotide sequence not generated yet"
    } else if (has_valid_structure(comp$Dinucleotide$Structure, nchar(comp$Dinucleotide$Sequence))) {
      structure <- comp$Dinucleotide$Structure
      mfe <- comp$Dinucleotide$MFE
      paired <- sum(strsplit(structure, "")[[1]] %in% c("(", ")"))
      total <- nchar(structure)
      
      paste0("Structure: ", structure, "\n",
             "MFE: ", if (has_valid_mfe(mfe)) paste(mfe, "kcal/mol") else "Not calculated", "\n",
             "Paired: ", paired, "/", total, " (", round(paired/total*100, 1), "%)")
    } else {
      "RNAfold structure and MFE: Not calculated"
    }
  })
  
  # Render Random structure (4-way comparison)
  output$comp4RandomPlot <- renderPlotly({
    req(input$compareBtn)
    req(values$selectedComparison)
    
    comp <- values$selectedComparison
    if (!is.null(comp$Random$Sequence) &&
        has_valid_structure(comp$Random$Structure, nchar(comp$Random$Sequence))) {
      create_3d_structure_plot(comp$Random$Sequence, comp$Random$Structure, comp$Random$MFE)
    } else {
      NULL
    }
  })
  
  output$comp4RandomText <- renderText({
    req(input$compareBtn)
    req(values$selectedComparison)
    
    comp <- values$selectedComparison
    if (is.null(comp$Random$Sequence)) {
      "Random sequence not generated yet"
    } else if (has_valid_structure(comp$Random$Structure, nchar(comp$Random$Sequence))) {
      structure <- comp$Random$Structure
      mfe <- comp$Random$MFE
      paired <- sum(strsplit(structure, "")[[1]] %in% c("(", ")"))
      total <- nchar(structure)
      
      paste0("Structure: ", structure, "\n",
             "MFE: ", if (has_valid_mfe(mfe)) paste(mfe, "kcal/mol") else "Not calculated", "\n",
             "Paired: ", paired, "/", total, " (", round(paired/total*100, 1), "%)")
    } else {
      "RNAfold structure and MFE: Not calculated"
    }
  })
  
  # Structure comparison summary table
  output$structureComparisonTable <- renderTable({
    req(input$compareBtn)
    req(values$selectedComparison)
    
    comp <- values$selectedComparison
    
    make_structure_row <- function(method, item) {
      sequence_available <- !is.null(item$Sequence)
      structure_available <- sequence_available &&
        has_valid_structure(item$Structure, nchar(item$Sequence))

      if (!structure_available) {
        return(data.frame(
          Method = method,
          MFE = NA_real_,
          Paired_Bases = NA_real_,
          Unpaired_Bases = NA_real_,
          Paired_Percent = NA_real_,
          RNAfold_Status = if (sequence_available) "Not calculated" else "Sequence not generated",
          stringsAsFactors = FALSE
        ))
      }

      paired <- sum(strsplit(item$Structure, "", fixed = TRUE)[[1]] %in% c("(", ")"))
      total <- nchar(item$Structure)
      data.frame(
        Method = method,
        MFE = if (has_valid_mfe(item$MFE)) item$MFE else NA_real_,
        Paired_Bases = paired,
        Unpaired_Bases = total - paired,
        Paired_Percent = round(paired / total * 100, 1),
        RNAfold_Status = if (has_valid_mfe(item$MFE)) "Calculated" else "Structure calculated; MFE unavailable",
        stringsAsFactors = FALSE
      )
    }

    comparison_data <- do.call(rbind, list(
      make_structure_row("Original", comp$Original),
      make_structure_row("NeRNA", comp$NeRNA),
      make_structure_row("Dinucleotide", comp$Dinucleotide),
      make_structure_row("Random", comp$Random)
    ))
    
    return(comparison_data)
  }, striped = TRUE, hover = TRUE, bordered = TRUE, spacing = 's')
  
  # ============================================================
  # COMPREHENSIVE EXCEL EXPORT
  # ============================================================
  
  # Comprehensive Excel Export
  output$downloadComprehensiveExcel <- downloadHandler(
    filename = function() {
      paste0("NeRNA_Comprehensive_Report_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".xlsx")
    },
    content = function(file) {
      tryCatch({
        # Load required library
        if (!requireNamespace("openxlsx", quietly = TRUE)) {
          install.packages("openxlsx", repos = "https://cran.rstudio.com/")
        }
        library(openxlsx)
        
        # Prepare data
        message("=== SERVER: PREPARING EXCEL DATA ===")
        original_data <- values$fastaData
        message("Original data: ", if(is.null(original_data)) "NULL" else paste("class =", class(original_data), "nrow =", nrow(original_data)))
        
        # Get all method results separately and convert to data frames if needed
        nerna_data <- values$nernaResults
        message("NeRNA data (before conversion): ", if(is.null(nerna_data)) "NULL" else paste("class =", paste(class(nerna_data), collapse=",")))
        if (!is.null(nerna_data) && is.list(nerna_data) && !is.data.frame(nerna_data)) {
          message("Converting NeRNA data from list to data frame...")
          nerna_data <- tryCatch(as.data.frame(nerna_data, stringsAsFactors = FALSE), 
                                  error = function(e) { message("Failed: ", e$message); NULL })
        }
        message("NeRNA data (after conversion): ", if(is.null(nerna_data)) "NULL" else paste("class =", paste(class(nerna_data), collapse=","), "nrow =", nrow(nerna_data)))
        
        dinuc_data <- values$dinucleotideResults
        message("Dinuc data (before conversion): ", if(is.null(dinuc_data)) "NULL" else paste("class =", paste(class(dinuc_data), collapse=",")))
        if (!is.null(dinuc_data) && is.list(dinuc_data) && !is.data.frame(dinuc_data)) {
          message("Converting Dinuc data from list to data frame...")
          dinuc_data <- tryCatch(as.data.frame(dinuc_data, stringsAsFactors = FALSE), 
                                  error = function(e) { message("Failed: ", e$message); NULL })
        }
        message("Dinuc data (after conversion): ", if(is.null(dinuc_data)) "NULL" else paste("class =", paste(class(dinuc_data), collapse=","), "nrow =", nrow(dinuc_data)))
        
        random_data <- values$randomResults
        message("Random data (before conversion): ", if(is.null(random_data)) "NULL" else paste("class =", paste(class(random_data), collapse=",")))
        if (!is.null(random_data) && is.list(random_data) && !is.data.frame(random_data)) {
          message("Converting Random data from list to data frame...")
          random_data <- tryCatch(as.data.frame(random_data, stringsAsFactors = FALSE), 
                                  error = function(e) { message("Failed: ", e$message); NULL })
        }
        message("Random data (after conversion): ", if(is.null(random_data)) "NULL" else paste("class =", paste(class(random_data), collapse=","), "nrow =", nrow(random_data)))
        
        # Get RNAfold results if available
        rnafold_results <- values$rnafoldResults
        message("RNAfold results: ", if(is.null(rnafold_results)) "NULL" else paste("nrow =", nrow(rnafold_results)))
        
        # Combine all for similarity metrics
        all_negatives <- list()
        if (!is.null(nerna_data) && is.data.frame(nerna_data) && nrow(nerna_data) > 0) {
          all_negatives[[length(all_negatives) + 1]] <- nerna_data
        }
        if (!is.null(dinuc_data) && is.data.frame(dinuc_data) && nrow(dinuc_data) > 0) {
          all_negatives[[length(all_negatives) + 1]] <- dinuc_data
        }
        if (!is.null(random_data) && is.data.frame(random_data) && nrow(random_data) > 0) {
          all_negatives[[length(all_negatives) + 1]] <- random_data
        }
        
        similarity_metrics <- NULL
        if (length(all_negatives) > 0) {
          # Try to combine - if columns don't match, just use first one
          similarity_metrics <- tryCatch({
            do.call(rbind, all_negatives)
          }, error = function(e) {
            all_negatives[[1]]
          })
        }
        
        # Export with separate method data
        message("Calling export_comprehensive_results...")
        export_comprehensive_results(original_data, nerna_data, dinuc_data, random_data, 
                                     similarity_metrics, rnafold_results, file)
        message("Excel export completed successfully!")
      }, error = function(e) {
        # Write detailed error message to file instead of HTM
        error_msg <- paste(
          "Error generating Excel file:",
          e$message,
          "\n\nStack trace:",
          paste(deparse(e$call), collapse = "\n"),
          "\n\nPlease check the console for more details."
        )
        writeLines(error_msg, file)
        message("Excel export error: ", e$message)
        message("Stack trace: ", deparse(e$call))
      })
    }
  )
  
  # ============================================================
  # SESSION MANAGEMENT
  # ============================================================
  
  # Save session
  output$saveSession <- downloadHandler(
    filename = function() {
      paste0("NeRNA_Session_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".rds")
    },
    content = function(file) {
      # Save session data
      success <- save_session(values, file)
      
      if (success) {
        showNotification("Session saved successfully!", type = "message", duration = 5)
      } else {
        showNotification("Error saving session!", type = "error", duration = 5)
      }
    }
  )
  
  # Load session
  observeEvent(input$loadSessionBtn, {
    req(input$loadSessionFile)
    
    # Load session data
    session_obj <- load_session(input$loadSessionFile$datapath)
    
    if (!is.null(session_obj)) {
      # Restore session data
      values$fastaData <- session_obj$fasta_data
      values$rnafoldResults <- session_obj$rnafold_results
      values$nernaResults <- session_obj$nerna_results
      values$dinucleotideResults <- session_obj$dinucleotide_results
      values$randomResults <- session_obj$random_results
      values$rnaType <- session_obj$rna_type
      values$dataLoaded <- TRUE
      
      showNotification(
        paste("Session loaded successfully! Loaded from:", 
              format(session_obj$timestamp, "%Y-%m-%d %H:%M:%S")),
        type = "message",
        duration = 10
      )
    } else {
      showNotification("Error loading session!", type = "error", duration = 5)
    }
  })
  
  # ============================================================
  # BATCH PROCESSING
  # ============================================================
  
  # Load batch FASTA files
  observeEvent(input$loadBatchBtn, {
    req(input$batchFastaFiles)
    
    # Process batch files
    file_paths <- input$batchFastaFiles$datapath
    
    withProgress(message = 'Processing batch files...', value = 0, {
      combined_data <- process_batch_fasta(file_paths)
      
      if (!is.null(combined_data) && nrow(combined_data) > 0) {
        values$fastaData <- combined_data
        values$dataLoaded <- TRUE
        
        showNotification(
          paste("Loaded", nrow(combined_data), "sequences from", 
                length(file_paths), "files!"),
          type = "message",
          duration = 10
        )
      } else {
        showNotification("Error loading batch files!", type = "error", duration = 5)
      }
    })
  })
  
  # ============================================================
  # BATCH FASTA DOWNLOAD HANDLER
  # ============================================================
  
  # Download all FASTA sequences (Batch Operations)
  output$downloadBatchFASTA <- downloadHandler(
    filename = function() {
      paste0("Batch_All_sequences_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".fasta")
    },
    content = function(file) {
      fasta_content <- character()
      
      # Add original sequences
      if (!is.null(values$fastaData) && nrow(values$fastaData) > 0) {
        for (i in 1:nrow(values$fastaData)) {
          fasta_content <- c(fasta_content, 
                           paste0(">", values$fastaData$Accession[i], "_Original"),
                           values$fastaData$Sequence[i])
        }
      }
      
      # Add NeRNA sequences
      if (!is.null(values$nernaResults) && nrow(values$nernaResults) > 0) {
        for (i in 1:nrow(values$nernaResults)) {
          fasta_content <- c(fasta_content, 
                           paste0(">", values$nernaResults$OriginalAccession[i], "_NeRNA"),
                           values$nernaResults$NegativeSequence[i])
        }
      }
      
      # Add Dinucleotide sequences
      if (!is.null(values$dinucleotideResults) && nrow(values$dinucleotideResults) > 0) {
        for (i in 1:nrow(values$dinucleotideResults)) {
          fasta_content <- c(fasta_content, 
                           paste0(">", values$dinucleotideResults$OriginalAccession[i], "_Dinucleotide"),
                           values$dinucleotideResults$NegativeSequence[i])
        }
      }
      
      # Add Random sequences
      if (!is.null(values$randomResults) && nrow(values$randomResults) > 0) {
        for (i in 1:nrow(values$randomResults)) {
          fasta_content <- c(fasta_content, 
                           paste0(">", values$randomResults$OriginalAccession[i], "_Random"),
                           values$randomResults$NegativeSequence[i])
        }
      }
      
      if (length(fasta_content) > 0) {
        writeLines(fasta_content, file)
      } else {
        writeLines("No sequences available", file)
      }
    }
  )
  
  # ============================================================
  # FASTA Download Handlers for Result Comparison
  # ============================================================
  
  # Download NeRNA FASTA sequences
  output$downloadComparisonNeRNAFASTA <- downloadHandler(
    filename = function() {
      paste0("NeRNA_sequences_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".fasta")
    },
    content = function(file) {
      if (!is.null(values$nernaResults) && nrow(values$nernaResults) > 0) {
        # Create FASTA content
        fasta_content <- character()
        for (i in 1:nrow(values$nernaResults)) {
          fasta_content <- c(fasta_content, 
                           paste0(">", values$nernaResults$OriginalAccession[i], "_NeRNA"),
                           values$nernaResults$NegativeSequence[i])
        }
        writeLines(fasta_content, file)
      } else {
        writeLines("No NeRNA sequences available", file)
      }
    }
  )
  
  # Download Dinucleotide FASTA sequences
  output$downloadDinucFASTA <- downloadHandler(
    filename = function() {
      paste0("Dinucleotide_sequences_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".fasta")
    },
    content = function(file) {
      if (!is.null(values$dinucleotideResults) && nrow(values$dinucleotideResults) > 0) {
        # Create FASTA content
        fasta_content <- character()
        for (i in 1:nrow(values$dinucleotideResults)) {
          fasta_content <- c(fasta_content, 
                           paste0(">", values$dinucleotideResults$OriginalAccession[i], "_Dinucleotide"),
                           values$dinucleotideResults$NegativeSequence[i])
        }
        writeLines(fasta_content, file)
      } else {
        writeLines("No Dinucleotide sequences available", file)
      }
    }
  )
  
  # Download Random FASTA sequences
  output$downloadComparisonRandomFASTA <- downloadHandler(
    filename = function() {
      paste0("Random_sequences_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".fasta")
    },
    content = function(file) {
      if (!is.null(values$randomResults) && nrow(values$randomResults) > 0) {
        # Create FASTA content
        fasta_content <- character()
        for (i in 1:nrow(values$randomResults)) {
          fasta_content <- c(fasta_content, 
                           paste0(">", values$randomResults$OriginalAccession[i], "_Random"),
                           values$randomResults$NegativeSequence[i])
        }
        writeLines(fasta_content, file)
      } else {
        writeLines("No Random sequences available", file)
      }
    }
  )
  
  # Download All FASTA sequences (Original + All Methods)
  output$downloadAllFASTA <- downloadHandler(
    filename = function() {
      paste0("All_sequences_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".fasta")
    },
    content = function(file) {
      fasta_content <- character()
      
      # Add original sequences
      if (!is.null(values$fastaData) && nrow(values$fastaData) > 0) {
        for (i in 1:nrow(values$fastaData)) {
          fasta_content <- c(fasta_content, 
                           paste0(">", values$fastaData$Accession[i], "_Original"),
                           values$fastaData$Sequence[i])
        }
      }
      
      # Add NeRNA sequences
      if (!is.null(values$nernaResults) && nrow(values$nernaResults) > 0) {
        for (i in 1:nrow(values$nernaResults)) {
          fasta_content <- c(fasta_content, 
                           paste0(">", values$nernaResults$OriginalAccession[i], "_NeRNA"),
                           values$nernaResults$NegativeSequence[i])
        }
      }
      
      # Add Dinucleotide sequences
      if (!is.null(values$dinucleotideResults) && nrow(values$dinucleotideResults) > 0) {
        for (i in 1:nrow(values$dinucleotideResults)) {
          fasta_content <- c(fasta_content, 
                           paste0(">", values$dinucleotideResults$OriginalAccession[i], "_Dinucleotide"),
                           values$dinucleotideResults$NegativeSequence[i])
        }
      }
      
      # Add Random sequences
      if (!is.null(values$randomResults) && nrow(values$randomResults) > 0) {
        for (i in 1:nrow(values$randomResults)) {
          fasta_content <- c(fasta_content, 
                           paste0(">", values$randomResults$OriginalAccession[i], "_Random"),
                           values$randomResults$NegativeSequence[i])
        }
      }
      
      if (length(fasta_content) > 0) {
        writeLines(fasta_content, file)
      } else {
        writeLines("No sequences available", file)
      }
    }
  )
  
  # Sequence-Specific Excel Download Handler
  output$downloadSequenceSpecificExcel <- downloadHandler(
    filename = function() {
      selected_seq <- input$comparisonSequence
      if (!is.null(selected_seq) && selected_seq != "") {
        paste0("Sequence_Specific_Report_", selected_seq, "_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".xlsx")
      } else {
        paste0("Sequence_Specific_Report_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".xlsx")
      }
    },
    content = function(file) {
      tryCatch({
        # Check if comparison data is available
        if (is.null(values$selectedComparison)) {
          stop("No comparison data available. Please select a sequence and click 'Compare' first.")
        }
        
        # Get the selected sequence data
        selected_acc <- input$comparisonSequence
        if (is.null(selected_acc) || selected_acc == "") {
          stop("No sequence selected. Please select a sequence first.")
        }
        
        # Find the original sequence data
        orig_idx <- which(values$fastaData$Accession == selected_acc)
        if (length(orig_idx) == 0) {
          stop("Selected sequence not found in original data.")
        }
        
        orig_seq <- values$fastaData$Sequence[orig_idx[1]]
        
        # Get original structure and MFE with safety checks
        rnafold_idx <- which(values$rnafoldResults$Accession == selected_acc)
        if (length(rnafold_idx) > 0) {
          orig_structure <- values$rnafoldResults$SecondaryStructure[rnafold_idx[1]]
          orig_mfe <- values$rnafoldResults$MFE[rnafold_idx[1]]
        } else {
          orig_structure <- "N/A"
          orig_mfe <- NA
        }
        
        # Get comparison data
        comp_data <- values$selectedComparison
        
        # Ensure openxlsx package is loaded
        if (!requireNamespace("openxlsx", quietly = TRUE)) {
          install.packages("openxlsx", repos = "https://cran.rstudio.com/")
        }
        library(openxlsx)
        
        # Create workbook
        wb <- openxlsx::createWorkbook()
        
        # Sheet 1: Sequence Information
        addWorksheet(wb, "Sequence Information")
        seq_info <- data.frame(
          Property = c("Accession", "Original Sequence", "Sequence Length", "Original Structure", "Original MFE"),
          Value = c(selected_acc, orig_seq, nchar(orig_seq), orig_structure, orig_mfe),
          stringsAsFactors = FALSE
        )
        writeData(wb, "Sequence Information", seq_info)
        
        # Sheet 2: Method Comparison
        addWorksheet(wb, "Method Comparison")
        method_comp <- data.frame(
          Method = c("Original", "NeRNA", "Dinucleotide", "Random"),
          Sequence = c(
            orig_seq,
            if (!is.null(comp_data$NeRNA$Sequence)) comp_data$NeRNA$Sequence else "N/A",
            if (!is.null(comp_data$Dinucleotide$Sequence)) comp_data$Dinucleotide$Sequence else "N/A",
            if (!is.null(comp_data$Random$Sequence)) comp_data$Random$Sequence else "N/A"
          ),
          Structure = c(
            orig_structure,
            if (!is.null(comp_data$NeRNA$Structure)) comp_data$NeRNA$Structure else "N/A",
            if (!is.null(comp_data$Dinucleotide$Structure)) comp_data$Dinucleotide$Structure else "N/A",
            if (!is.null(comp_data$Random$Structure)) comp_data$Random$Structure else "N/A"
          ),
          MFE = c(
            orig_mfe,
            if (!is.null(comp_data$NeRNA$MFE)) comp_data$NeRNA$MFE else NA,
            if (!is.null(comp_data$Dinucleotide$MFE)) comp_data$Dinucleotide$MFE else NA,
            if (!is.null(comp_data$Random$MFE)) comp_data$Random$MFE else NA
          ),
          GC_Content = c(
            round(sum(strsplit(orig_seq, "")[[1]] %in% c("G", "C")) / nchar(orig_seq) * 100, 2),
            if (!is.null(comp_data$NeRNA$Sequence)) round(sum(strsplit(comp_data$NeRNA$Sequence, "")[[1]] %in% c("G", "C")) / nchar(comp_data$NeRNA$Sequence) * 100, 2) else NA,
            if (!is.null(comp_data$Dinucleotide$Sequence)) round(sum(strsplit(comp_data$Dinucleotide$Sequence, "")[[1]] %in% c("G", "C")) / nchar(comp_data$Dinucleotide$Sequence) * 100, 2) else NA,
            if (!is.null(comp_data$Random$Sequence)) round(sum(strsplit(comp_data$Random$Sequence, "")[[1]] %in% c("G", "C")) / nchar(comp_data$Random$Sequence) * 100, 2) else NA
          ),
          stringsAsFactors = FALSE
        )
        writeData(wb, "Method Comparison", method_comp)
        
        # Sheet 3: Similarity Metrics
        addWorksheet(wb, "Similarity Metrics")
        
        # Get similarity metrics from the results data frames
        get_similarity_metrics <- function(method_name, selected_acc) {
          tryCatch({
            if (method_name == "NeRNA" && !is.null(values$nernaResults) && is.data.frame(values$nernaResults) && nrow(values$nernaResults) > 0) {
              idx <- which(values$nernaResults$OriginalAccession == selected_acc)
              if (length(idx) > 0) {
                return(list(
                  sequence_identity = if(!is.na(values$nernaResults$SequenceIdentity[idx[1]])) values$nernaResults$SequenceIdentity[idx[1]] else NA,
                  dinucleotide_similarity = if(!is.na(values$nernaResults$DinucleotideSimilarity[idx[1]])) values$nernaResults$DinucleotideSimilarity[idx[1]] else NA,
                  gc_similarity = if(!is.na(values$nernaResults$GCSimilarity[idx[1]])) values$nernaResults$GCSimilarity[idx[1]] else NA
                ))
              }
            } else if (method_name == "Dinucleotide" && !is.null(values$dinucleotideResults) && is.data.frame(values$dinucleotideResults) && nrow(values$dinucleotideResults) > 0) {
              idx <- which(values$dinucleotideResults$OriginalAccession == selected_acc)
              if (length(idx) > 0) {
                return(list(
                  sequence_identity = if(!is.na(values$dinucleotideResults$SequenceIdentity[idx[1]])) values$dinucleotideResults$SequenceIdentity[idx[1]] else NA,
                  dinucleotide_similarity = if(!is.na(values$dinucleotideResults$DinucleotideSimilarity[idx[1]])) values$dinucleotideResults$DinucleotideSimilarity[idx[1]] else NA,
                  gc_similarity = if(!is.na(values$dinucleotideResults$GCSimilarity[idx[1]])) values$dinucleotideResults$GCSimilarity[idx[1]] else NA
                ))
              }
            } else if (method_name == "Random" && !is.null(values$randomResults) && is.data.frame(values$randomResults) && nrow(values$randomResults) > 0) {
              idx <- which(values$randomResults$OriginalAccession == selected_acc)
              if (length(idx) > 0) {
                return(list(
                  sequence_identity = if(!is.na(values$randomResults$SequenceIdentity[idx[1]])) values$randomResults$SequenceIdentity[idx[1]] else NA,
                  dinucleotide_similarity = if(!is.na(values$randomResults$DinucleotideSimilarity[idx[1]])) values$randomResults$DinucleotideSimilarity[idx[1]] else NA,
                  gc_similarity = if(!is.na(values$randomResults$GCSimilarity[idx[1]])) values$randomResults$GCSimilarity[idx[1]] else NA
                ))
              }
            }
          }, error = function(e) {
            # Silent error handling
          })
          
          # Return default values if not found or error
          return(list(
            sequence_identity = NA,
            dinucleotide_similarity = NA,
            gc_similarity = NA
          ))
        }
        
        # Get similarity metrics for each method
        nerna_sim <- get_similarity_metrics("NeRNA", selected_acc)
        dinuc_sim <- get_similarity_metrics("Dinucleotide", selected_acc)
        random_sim <- get_similarity_metrics("Random", selected_acc)
        
        sim_metrics <- data.frame(
          Metric = c("Sequence Identity", "Dinucleotide Similarity", "GC Similarity"),
          NeRNA = c(
            nerna_sim$sequence_identity,
            nerna_sim$dinucleotide_similarity,
            nerna_sim$gc_similarity
          ),
          Dinucleotide = c(
            dinuc_sim$sequence_identity,
            dinuc_sim$dinucleotide_similarity,
            dinuc_sim$gc_similarity
          ),
          Random = c(
            random_sim$sequence_identity,
            random_sim$dinucleotide_similarity,
            random_sim$gc_similarity
          ),
          stringsAsFactors = FALSE
        )
        writeData(wb, "Similarity Metrics", sim_metrics)
        
        # Sheet 4: Structure Analysis
        addWorksheet(wb, "Structure Analysis")
        # Helper function to calculate structure stats
        calc_structure_stats <- function(structure, sequence) {
          if (is.null(sequence) || !has_valid_structure(structure, nchar(sequence))) {
            return(list(paired = NA_real_, unpaired = NA_real_, percentage = NA_real_))
          }
          structure_chars <- strsplit(structure, "", fixed = TRUE)[[1]]
          paired <- sum(structure_chars %in% c("(", ")"))
          unpaired <- sum(structure_chars == ".")
          total <- nchar(structure)
          percentage <- round(paired / total * 100, 2)
          return(list(paired = paired, unpaired = unpaired, percentage = percentage))
        }
        
        orig_stats <- calc_structure_stats(orig_structure, orig_seq)
        nerna_stats <- calc_structure_stats(comp_data$NeRNA$Structure, comp_data$NeRNA$Sequence)
        dinuc_stats <- calc_structure_stats(comp_data$Dinucleotide$Structure, comp_data$Dinucleotide$Sequence)
        random_stats <- calc_structure_stats(comp_data$Random$Structure, comp_data$Random$Sequence)
        
        structure_analysis <- data.frame(
          Method = c("Original", "NeRNA", "Dinucleotide", "Random"),
          Paired_Bases = c(
            orig_stats$paired,
            nerna_stats$paired,
            dinuc_stats$paired,
            random_stats$paired
          ),
          Unpaired_Bases = c(
            orig_stats$unpaired,
            nerna_stats$unpaired,
            dinuc_stats$unpaired,
            random_stats$unpaired
          ),
          Paired_Percentage = c(
            orig_stats$percentage,
            nerna_stats$percentage,
            dinuc_stats$percentage,
            random_stats$percentage
          ),
          stringsAsFactors = FALSE
        )
        writeData(wb, "Structure Analysis", structure_analysis)
        
        # Save workbook
        openxlsx::saveWorkbook(wb, file, overwrite = TRUE)
        
      }, error = function(e) {
        # Write detailed error message to file
        error_msg <- paste(
          "Error generating sequence-specific Excel report:",
          e$message,
          "\n\nStack trace:",
          paste(deparse(e$call), collapse = "\n"),
          "\n\nDebug info:",
          "selected_acc:", selected_acc,
          "orig_seq length:", if(!is.null(orig_seq)) nchar(orig_seq) else "NULL",
          "orig_structure:", if(!is.null(orig_structure)) substr(orig_structure, 1, 50) else "NULL",
          "comp_data null:", is.null(comp_data)
        )
        writeLines(error_msg, file)
      })
    }
  )
  
  # Download handler for sequence batches (ZIP file)
  output$downloadBatchesZip <- downloadHandler(
    filename = function() {
      paste0("NeRNA_batches_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".zip")
    },
    content = function(file) {
      req(values$sequenceBatches)
      
      tryCatch({
        # Get base name from original file
        base_name <- if (!is.null(input$fastaFile)) {
          tools::file_path_sans_ext(basename(input$fastaFile$name))
        } else {
          "sequences"
        }
        
        # Create ZIP file with all batches
        zip_path <- create_batch_zip(values$sequenceBatches, base_name)
        
        if (!is.null(zip_path) && file.exists(zip_path)) {
          file.copy(zip_path, file, overwrite = TRUE)
          showNotification(
            paste("Downloaded", length(values$sequenceBatches), "batch files in ZIP archive."),
            type = "message",
            duration = 5
          )
        } else {
          stop("Failed to create ZIP file")
        }
      }, error = function(e) {
        showNotification(
          paste("Error creating batch files:", e$message),
          type = "error",
          duration = 10
        )
        # Write error to file
        writeLines(paste("Error:", e$message), file)
      })
    }
  )

}

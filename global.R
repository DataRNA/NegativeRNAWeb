library(shiny)
library(shinydashboard)
library(DT)
library(Biostrings)
library(seqinr)

# For radar charts in the comparison tab
if (!requireNamespace("fmsb", quietly = TRUE)) {
  install.packages("fmsb")
}
library(fmsb)

# For parallel processing
if (!requireNamespace("future", quietly = TRUE)) {
  install.packages("future")
}
if (!requireNamespace("future.apply", quietly = TRUE)) {
  install.packages("future.apply")
}
if (!requireNamespace("progressr", quietly = TRUE)) {
  install.packages("progressr")
}
library(future)
library(future.apply)
library(progressr)

# For nice visualization
if (!requireNamespace("RColorBrewer", quietly = TRUE)) {
  install.packages("RColorBrewer")
}
library(RColorBrewer)

# For Excel export
if (!requireNamespace("openxlsx", quietly = TRUE)) {
  install.packages("openxlsx")
}
library(openxlsx)

# For 3D structure visualization
if (!requireNamespace("plotly", quietly = TRUE)) {
  install.packages("plotly")
}
library(plotly)

# ============================================================
# RNAFOLD EXECUTABLE DISCOVERY
# ============================================================

get_rnafold_executable <- function() {
  # First check environment variable
  env_path <- Sys.getenv("RNAFOLD_PATH", unset = "")
  if (env_path != "" && file.exists(env_path)) {
    return(normalizePath(env_path, winslash = "/", mustWork = FALSE))
  }
  
  # On Unix/Linux (including Docker), prioritize system PATH lookup
  if (.Platform$OS.type == "unix") {
    system_rnafold <- Sys.which("RNAfold")
    if (system_rnafold != "" && file.exists(system_rnafold)) {
      return(system_rnafold)
    }
    # Also try /usr/bin/RNAfold explicitly (common location)
    if (file.exists("/usr/bin/RNAfold")) {
      return("/usr/bin/RNAfold")
    }
  }
  
  # Then check local utils folder (for Windows or local builds)
  candidates <- unique(c(
    if (.Platform$OS.type == "windows") file.path("utils", "RNAfold.exe") else file.path("utils", "RNAfold"),
    file.path("utils", "RNAfold.exe"),
    file.path("utils", "RNAfold")
  ))
  
  for (candidate in candidates) {
    if (!is.null(candidate) && candidate != "" && file.exists(candidate)) {
      return(normalizePath(candidate, winslash = "/", mustWork = FALSE))
    }
  }
  
  # Final fallback: try system PATH again (for runtime lookup)
  system_rnafold <- Sys.which("RNAfold")
  if (system_rnafold != "") {
    return(system_rnafold)
  }
  
  # Last resort: return command name and hope it's in PATH
  return(if (.Platform$OS.type == "windows") "RNAfold.exe" else "RNAfold")
}

# ============================================================
# PARALLEL PROCESSING SETUP AND UTILITIES
# ============================================================

# Initialize parallel processing
setup_parallel_processing <- function(n_cores = NULL) {
  tryCatch({
    # Detect available cores
    available_cores <- future::availableCores()
    
    # If n_cores not specified, use 75% of available cores
    if (is.null(n_cores)) {
      n_cores <- max(1, floor(available_cores * 0.75))
    } else {
      n_cores <- min(n_cores, available_cores)
    }
    
    message(sprintf("Setting up parallel processing with %d cores (out of %d available)", 
                    n_cores, available_cores))
    
    # Setup multicore or multisession based on OS
    if (.Platform$OS.type == "unix") {
      # Unix/Linux/Mac: use multicore (faster, shared memory)
      plan(multicore, workers = n_cores)
    } else {
      # Windows: use multisession (separate R sessions)
      plan(multisession, workers = n_cores)
    }
    
    return(n_cores)
  }, error = function(e) {
    message("Error setting up parallel processing: ", e$message)
    message("Falling back to sequential processing")
    plan(sequential)
    return(1)
  })
}

# Get current parallel plan info
get_parallel_info <- function() {
  current_plan <- class(plan())[1]
  n_workers <- nbrOfWorkers()
  
  return(list(
    plan = current_plan,
    workers = n_workers,
    is_parallel = n_workers > 1
  ))
}

# Disable parallel processing (return to sequential)
disable_parallel_processing <- function() {
  plan(sequential)
  message("Parallel processing disabled. Using sequential processing.")
}

# Function to read FASTA file and convert to data frame
readFastaToTable <- function(fastaFile, rnaType = NULL) {
  if (is.null(fastaFile)) {
    return(NULL)
  }
  
  # Check if it's RNA or DNA
  firstLines <- readLines(fastaFile, n = 10)
  sequences <- character(0)
  for (i in 1:length(firstLines)) {
    if (!startsWith(firstLines[i], ">")) {
      sequences <- c(sequences, firstLines[i])
      if (length(sequences) >= 2) break
    }
  }
  
  # Check if the sequence contains 'T' (DNA) instead of 'U' (RNA)
  isDNA <- any(grepl("T", sequences, ignore.case = FALSE))
  
  # Try to read the file
  fasta_data <- try({
    readRNAStringSet(fastaFile)
  }, silent = TRUE)
  
  if (inherits(fasta_data, "try-error")) {
    # If RNA reading fails, try with DNA as fallback
    fasta_data <- readDNAStringSet(fastaFile)
  }
  
  # Create data frame with Accession and Sequence
  df <- data.frame(
    Accession = names(fasta_data),
    Sequence = as.character(fasta_data),
    stringsAsFactors = FALSE
  )
  
  # Clean accession IDs (remove everything after first space)
  df$Accession <- sub(" .*", "", df$Accession)
  
  # Add a flag to indicate if the data is DNA
  attr(df, "isDNA") <- isDNA
  
  # Clean non-standard bases for tRNA (keep only A, U, G, C, T)
  if (!is.null(rnaType) && tolower(rnaType) == "trna") {
    original_seq_count <- nrow(df)
    cleaned_sequences <- character(nrow(df))
    non_standard_bases_found <- FALSE
    
    for (i in 1:nrow(df)) {
      original_seq <- df$Sequence[i]
      # Remove any character that is not A, U, G, C, T (case insensitive)
      cleaned_seq <- gsub("[^AUGCTaugct]", "", original_seq)
      cleaned_sequences[i] <- cleaned_seq
      
      if (nchar(cleaned_seq) != nchar(original_seq)) {
        non_standard_bases_found <- TRUE
      }
    }
    
    # Update sequences with cleaned versions
    df$Sequence <- cleaned_sequences
    
    # Remove sequences that became empty after cleaning
    empty_seqs <- nchar(df$Sequence) == 0
    if (any(empty_seqs)) {
      message(sprintf("Removed %d sequences that became empty after cleaning non-standard bases.", 
                      sum(empty_seqs)))
      df <- df[!empty_seqs, ]
    }
    
    if (non_standard_bases_found) {
      message(sprintf("tRNA sequences cleaned: non-standard bases (not A, U, G, C, T) removed from %d sequences.", 
                      original_seq_count))
    }
  }
  
  # Filter sequences by length (max 1000 nucleotides for system stability)
  max_length <- 1000
  original_count <- nrow(df)
  sequence_lengths <- nchar(df$Sequence)
  
  # Store filtered sequences info before filtering
  filtered_sequences <- df[sequence_lengths > max_length, ]
  
  df <- df[sequence_lengths <= max_length, ]
  filtered_count <- original_count - nrow(df)
  
  if (filtered_count > 0) {
    message(sprintf("Filtered %d sequences longer than %d nucleotides for system stability.", 
                    filtered_count, max_length))
    message(sprintf("Remaining sequences: %d", nrow(df)))
  }
  
  # Store original count and filtered sequences
  attr(df, "original_count") <- original_count
  attr(df, "filtered_count") <- filtered_count
  attr(df, "filtered_sequences") <- filtered_sequences
  
  # Check if sequences exceed maximum allowed (250)
  max_sequences <- 250
  if (nrow(df) > max_sequences) {
    attr(df, "exceeds_limit") <- TRUE
    attr(df, "total_sequences") <- nrow(df)
    attr(df, "max_sequences") <- max_sequences
  } else {
    attr(df, "exceeds_limit") <- FALSE
  }
  
  return(df)
}

# Function to split sequences into batches of 250
split_sequences_to_batches <- function(fasta_data, batch_size = 250) {
  total_sequences <- nrow(fasta_data)
  num_batches <- ceiling(total_sequences / batch_size)
  
  # Get attributes from original data
  filtered_sequences <- attr(fasta_data, "filtered_sequences")
  filtered_count <- attr(fasta_data, "filtered_count")
  original_count <- attr(fasta_data, "original_count")
  
  batches <- list()
  
  for (i in 1:num_batches) {
    start_idx <- (i - 1) * batch_size + 1
    end_idx <- min(i * batch_size, total_sequences)
    
    batch <- fasta_data[start_idx:end_idx, ]
    
    # Copy attributes to each batch (for README generation)
    attr(batch, "filtered_sequences") <- filtered_sequences
    attr(batch, "filtered_count") <- filtered_count
    attr(batch, "original_count") <- original_count
    
    batches[[i]] <- batch
  }
  
  return(batches)
}

# Function to write batch FASTA files
write_batch_fasta_files <- function(batches, output_dir, base_name = "sequences") {
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  
  file_paths <- character(length(batches))
  
  for (i in 1:length(batches)) {
    batch <- batches[[i]]
    file_name <- sprintf("%s_batch_%d_of_%d.fasta", base_name, i, length(batches))
    file_path <- file.path(output_dir, file_name)
    
    # Write FASTA file
    fasta_content <- character(nrow(batch) * 2)
    for (j in 1:nrow(batch)) {
      fasta_content[(j - 1) * 2 + 1] <- paste0(">", batch$Accession[j])
      fasta_content[(j - 1) * 2 + 2] <- batch$Sequence[j]
    }
    
    writeLines(fasta_content, file_path)
    file_paths[i] <- file_path
  }
  
  return(file_paths)
}

# Function to create a single ZIP file with all batches
create_batch_zip <- function(batches, base_name = "sequences") {
  # Create unique temp directory for this batch
  temp_dir <- tempdir()
  batch_id <- format(Sys.time(), "%Y%m%d_%H%M%S")
  batch_dir <- file.path(temp_dir, paste0("batches_", batch_id))
  
  # Clean up if directory already exists
  if (dir.exists(batch_dir)) {
    unlink(batch_dir, recursive = TRUE)
  }
  dir.create(batch_dir, recursive = TRUE)
  
  message("Creating batch files in: ", batch_dir)
  
  # Write all batch files
  file_paths <- write_batch_fasta_files(batches, batch_dir, base_name)
  
  # Verify files were created
  message("Created ", length(file_paths), " batch files")
  for (fp in file_paths) {
    if (file.exists(fp)) {
      message("  - ", basename(fp), " (", file.size(fp), " bytes)")
    } else {
      message("  - ERROR: ", basename(fp), " was not created!")
    }
  }
  
  # Create README file
  readme_path <- file.path(batch_dir, "README.txt")
  readme_content <- c(
    "NeRNA Batch Processing",
    "======================",
    "",
    sprintf("Your FASTA file contained %d sequences, which exceeds the maximum allowed (%d).", 
            sum(sapply(batches, nrow)), 250),
    sprintf("The sequences have been split into %d batches of up to 250 sequences each.", 
            length(batches)),
    "",
    "Instructions:",
    "1. Extract all FASTA files from this ZIP archive",
    "2. Upload each batch file separately to NeRNA",
    "3. Process each batch individually",
    "4. Combine results from all batches for your final analysis",
    "",
    "Batch Information:",
    "-------------------"
  )
  
  for (i in 1:length(batches)) {
    batch_info <- sprintf("Batch %d: %d sequences (%s to %s)", 
                          i, 
                          nrow(batches[[i]]),
                          batches[[i]]$Accession[1],
                          batches[[i]]$Accession[nrow(batches[[i]])])
    readme_content <- c(readme_content, batch_info)
  }
  
  readme_content <- c(readme_content, 
                      "",
                      "Generated by NeRNA",
                      sprintf("Date: %s", Sys.time()))
  
  # Add filtered sequences information if available
  filtered_sequences <- attr(batches[[1]], "filtered_sequences")
  filtered_count <- attr(batches[[1]], "filtered_count")
  
  if (!is.null(filtered_sequences) && !is.null(filtered_count) && filtered_count > 0) {
    readme_content <- c(
      readme_content,
      "",
      "",
      "================================================================================",
      "FILTERED SEQUENCES (>1000 nucleotides)",
      "================================================================================",
      "",
      sprintf("The following %d sequences were filtered out because they exceed", filtered_count),
      "the maximum allowed length of 1000 nucleotides for system stability:",
      "",
      "Seq#  | Accession                          | Length (nt) | Status",
      "------|---------------------------------------|-------------|--------"
    )
    
    for (i in 1:nrow(filtered_sequences)) {
      seq_info <- sprintf("%-5d | %-37s | %11d | Filtered",
                         i,
                         filtered_sequences$Accession[i],
                         nchar(filtered_sequences$Sequence[i]))
      readme_content <- c(readme_content, seq_info)
    }
    
    readme_content <- c(
      readme_content,
      "",
      "Note: These sequences were automatically excluded from processing.",
      "If you need to analyze these sequences, please:",
      "  1. Split them into smaller segments (≤1000 nt each), or",
      "  2. Use external tools for long sequence analysis",
      ""
    )
  }
  
  writeLines(readme_content, readme_path)
  message("Created README.txt (", file.size(readme_path), " bytes)")
  
  # Create a separate FASTA file for filtered sequences if any exist
  if (!is.null(filtered_sequences) && !is.null(filtered_count) && filtered_count > 0) {
    filtered_fasta_path <- file.path(batch_dir, "FILTERED_sequences.fasta")
    filtered_fasta_content <- character(nrow(filtered_sequences) * 2)
    
    for (i in 1:nrow(filtered_sequences)) {
      filtered_fasta_content[(i - 1) * 2 + 1] <- paste0(">", 
                                                        filtered_sequences$Accession[i],
                                                        " | Length: ",
                                                        nchar(filtered_sequences$Sequence[i]),
                                                        " nt (>1000) | FILTERED")
      filtered_fasta_content[(i - 1) * 2 + 2] <- filtered_sequences$Sequence[i]
    }
    
    writeLines(filtered_fasta_content, filtered_fasta_path)
    message("Created FILTERED_sequences.fasta with ", filtered_count, " sequences (", 
            file.size(filtered_fasta_path), " bytes)")
  }
  
  # Create ZIP file
  zip_name <- paste0(base_name, "_batches.zip")
  zip_path <- file.path(temp_dir, zip_name)
  
  # Remove old zip if exists
  if (file.exists(zip_path)) {
    unlink(zip_path)
  }
  
  # Try different methods to create ZIP
  zip_success <- FALSE
  
  # Method 1: Use utils::zip
  tryCatch({
    message("Attempting to create ZIP using utils::zip...")
    all_files <- c(file_paths, readme_path)
    result <- utils::zip(
      zipfile = zip_path,
      files = all_files,
      flags = "-r9Xj"  # -j: junk paths (store files without directory structure)
    )
    
    if (result == 0 && file.exists(zip_path) && file.size(zip_path) > 0) {
      message("ZIP created successfully with utils::zip")
      zip_success <- TRUE
    }
  }, error = function(e) {
    message("utils::zip failed: ", e$message)
  })
  
  # Method 2: Use system2 with PowerShell (Windows)
  if (!zip_success && .Platform$OS.type == "windows") {
    tryCatch({
      message("Attempting to create ZIP using PowerShell...")
      
      # PowerShell command to create ZIP
      ps_cmd <- sprintf(
        'Compress-Archive -Path "%s\\*" -DestinationPath "%s" -Force',
        normalizePath(batch_dir, winslash = "/"),
        normalizePath(zip_path, winslash = "/", mustWork = FALSE)
      )
      
      result <- system2("powershell", 
                       args = c("-Command", ps_cmd),
                       stdout = TRUE,
                       stderr = TRUE)
      
      if (file.exists(zip_path) && file.size(zip_path) > 0) {
        message("ZIP created successfully with PowerShell")
        zip_success <- TRUE
      }
    }, error = function(e) {
      message("PowerShell compression failed: ", e$message)
    })
  }
  
  # Method 3: Use system zip command (Unix/Linux/Mac)
  if (!zip_success && .Platform$OS.type == "unix") {
    tryCatch({
      message("Attempting to create ZIP using system zip command...")
      current_dir <- getwd()
      setwd(batch_dir)
      
      result <- system2("zip",
                       args = c("-r", zip_path, "."),
                       stdout = TRUE,
                       stderr = TRUE)
      
      setwd(current_dir)
      
      if (file.exists(zip_path) && file.size(zip_path) > 0) {
        message("ZIP created successfully with system zip")
        zip_success <- TRUE
      }
    }, error = function(e) {
      message("System zip failed: ", e$message)
      setwd(current_dir)
    })
  }
  
  # Check final result
  if (zip_success && file.exists(zip_path)) {
    message("Final ZIP file: ", zip_path, " (", file.size(zip_path), " bytes)")
    return(zip_path)
  } else {
    message("ERROR: Failed to create ZIP file")
    return(NULL)
  }
}

# Function to run RNAfold with parallel processing support
runRNAfold <- function(fastaFile, rnaType = "mirna", use_parallel = TRUE) {
  if (is.null(fastaFile)) {
    return(NULL)
  }
  
  # Read FASTA file to count sequences
  fasta_data <- try(readLines(fastaFile), silent = TRUE)
  if (inherits(fasta_data, "try-error")) {
    return(NULL)
  }
  
  # Count sequences (lines starting with >)
  n_sequences <- sum(grepl("^>", fasta_data))
  message(sprintf("Processing %d sequences with RNAfold", n_sequences))
  
  # Decide whether to use parallel processing
  # Only use parallel if we have many sequences and it's enabled
  use_parallel <- use_parallel && n_sequences >= 10
  
  if (use_parallel) {
    message("Using PARALLEL processing for RNAfold")
    return(runRNAfold_parallel(fastaFile, rnaType))
  } else {
    message("Using SEQUENTIAL processing for RNAfold")
    return(runRNAfold_sequential(fastaFile, rnaType))
  }
}

# Sequential RNAfold (original implementation)
runRNAfold_sequential <- function(fastaFile, rnaType = "mirna") {
  if (is.null(fastaFile)) {
    return(NULL)
  }
  
  # Create temp directory for input and output
  tempDir <- tempdir()
  tempFile <- file.path(tempDir, "input.fasta")
  
  # Copy input file to temp location
  file.copy(fastaFile, tempFile, overwrite = TRUE)
  
  # Path to RNAfold executable (cross-platform)
  rnafoldPath <- get_rnafold_executable()
  
  # Set arguments based on RNA type
  args <- c("-i", tempFile, "--noPS")
  
  # Add --circ argument for circRNA
  if (rnaType == "circrna") {
    args <- c(args, "--circ")
    message("Using circular RNA mode with --circ argument")
  }
  
  # Run RNAfold with the appropriate arguments
  message("Running RNAfold with command: ", rnafoldPath, " ", paste(args, collapse = " "))
  
  # Change working directory to tempDir to capture output there
  oldwd <- getwd()
  setwd(tempDir)
  
  cmd_result <- system2(rnafoldPath, 
                     args = args, 
                     wait = TRUE,
                     stdout = TRUE,
                     stderr = TRUE)
  
  # Log command output
  message("RNAfold command output: ", paste(cmd_result, collapse = "\n"))
  
  # Return to original working directory
  setwd(oldwd)
  
  # New parsing algorithm for RNAfold output
  if (length(cmd_result) > 0) {
    # Initialize variables
    accessions <- c()
    sequences <- c()
    structures <- c()
    mfe_values <- c()
    
    i <- 1
    while (i <= length(cmd_result)) {
      line <- cmd_result[i]
      
      # Check if this is a header line
      if (startsWith(line, ">")) {
        accession <- sub("^>", "", line)
        
        # Check if there are at least two more lines (sequence and structure)
        if (i + 1 <= length(cmd_result) && i + 2 <= length(cmd_result)) {
          # Extract sequence
          sequence <- cmd_result[i + 1]
          
          # Extract structure and MFE
          structure_line <- cmd_result[i + 2]
          # Split by space followed by parenthesis
          structure_parts <- strsplit(structure_line, " \\(")[[1]]
          
          if (length(structure_parts) >= 2) {
            structure <- structure_parts[1]
            # Extract MFE value (remove closing parenthesis)
            mfe <- as.numeric(gsub("\\)", "", structure_parts[2]))
            
            # Add to our data
            accessions <- c(accessions, accession)
            sequences <- c(sequences, sequence)
            structures <- c(structures, structure)
            mfe_values <- c(mfe_values, mfe)
          }
        }
        
        # Skip to next header (3 lines per entry)
        i <- i + 3
      } else {
        # Skip unrecognized line
        i <- i + 1
      }
    }
    
    if (length(accessions) > 0) {
      # Create dataframe
      result_df <- data.frame(
        Accession = accessions,
        Sequence = sequences,
        SecondaryStructure = structures,
        MFE = mfe_values,
        stringsAsFactors = FALSE
      )
      return(result_df)
    } else {
      message("Failed to parse RNAfold output")
      return(NULL)
    }
  } else {
    return(NULL)
  }
}

# Parallel RNAfold implementation
runRNAfold_parallel <- function(fastaFile, rnaType = "mirna") {
  tryCatch({
    # Read and parse FASTA file
    fasta_lines <- readLines(fastaFile)
    
    # Parse sequences
    sequences_list <- list()
    current_acc <- NULL
    current_seq <- ""
    
    for (line in fasta_lines) {
      if (startsWith(line, ">")) {
        # Save previous sequence if exists
        if (!is.null(current_acc) && nchar(current_seq) > 0) {
          sequences_list[[length(sequences_list) + 1]] <- list(
            accession = current_acc,
            sequence = current_seq
          )
        }
        # Start new sequence
        current_acc <- sub("^>", "", line)
        current_seq <- ""
      } else {
        current_seq <- paste0(current_seq, trimws(line))
      }
    }
    
    # Add last sequence
    if (!is.null(current_acc) && nchar(current_seq) > 0) {
      sequences_list[[length(sequences_list) + 1]] <- list(
        accession = current_acc,
        sequence = current_seq
      )
    }
    
    n_sequences <- length(sequences_list)
    message(sprintf("Parsed %d sequences for parallel processing", n_sequences))
    
    # Setup parallel processing
    setup_parallel_processing()
    
    # Process sequences in parallel with progress
    message("Starting parallel RNAfold processing...")
    
    results_list <- future_lapply(sequences_list, function(seq_data) {
      # Run quick_rnafold for each sequence
      result <- quick_rnafold(seq_data$sequence)
      
      return(list(
        Accession = seq_data$accession,
        Sequence = seq_data$sequence,
        SecondaryStructure = result$structure,
        MFE = result$mfe
      ))
    }, future.seed = TRUE)
    
    # Convert to data frame
    result_df <- data.frame(
      Accession = sapply(results_list, function(x) x$Accession),
      Sequence = sapply(results_list, function(x) x$Sequence),
      SecondaryStructure = sapply(results_list, function(x) x$SecondaryStructure),
      MFE = sapply(results_list, function(x) x$MFE),
      stringsAsFactors = FALSE
    )
    
    message(sprintf("Parallel RNAfold completed: %d sequences processed", nrow(result_df)))
    
    return(result_df)
    
  }, error = function(e) {
    message("Error in parallel RNAfold: ", e$message)
    message("Falling back to sequential processing...")
    return(runRNAfold_sequential(fastaFile, rnaType))
  })
}

# NeRNA Method - Correct Implementation with Structural Transformation and Octal Shifting
generateNeRNA <- function(sequences, shifting_size = 1, seed = NULL, rnafold_results = NULL, use_parallel = TRUE) {
  # Set seed if provided
  if (!is.null(seed)) {
    set.seed(seed)
  }
  
  # Decide whether to use parallel processing
  n_sequences <- length(sequences)
  use_parallel <- use_parallel && n_sequences >= 10
  
  if (use_parallel) {
    message(sprintf("Using PARALLEL processing for NeRNA (%d sequences)", n_sequences))
    return(generateNeRNA_parallel(sequences, shifting_size, seed, rnafold_results))
  } else {
    message(sprintf("Using SEQUENTIAL processing for NeRNA (%d sequences)", n_sequences))
    return(generateNeRNA_sequential(sequences, shifting_size, seed, rnafold_results))
  }
}

# Sequential NeRNA (original implementation)
generateNeRNA_sequential <- function(sequences, shifting_size = 1, seed = NULL, rnafold_results = NULL) {
  # Set seed if provided
  if (!is.null(seed)) {
    set.seed(seed)
  }
  
  # Generate NeRNA for a single sequence
  generate_single_nerna <- function(seq) {
    # Step 1: Get secondary structure (use RNAfold if available, otherwise use simple pattern)
    structure <- get_secondary_structure(seq)
    
    # DEBUG: Print all steps
    cat("=== DEBUG NeRNA ===\n")
    cat("Original:", seq, "\n")
    cat("Structure:", structure, "\n")
    
    # Use the correct rna_octal_rotate_once function
    result <- rna_octal_rotate_once(seq, structure, shifting_size)
    
    cat("Encoded bits:", result$encoded_bits, "\n")
    cat("Bit string length:", nchar(result$encoded_bits), "\n")
    cat("Rotated bits:", result$rotated_bits, "\n")
    cat("Result:", result$decoded_sequence, "\n")
    cat("Result length:", nchar(result$decoded_sequence), "vs Original length:", nchar(seq), "\n")
    cat("==================\n")
    
    return(result$decoded_sequence)
  }
  
  # Get secondary structure - use RNAfold results if available
  get_secondary_structure <- function(seq) {
    # Try to find the structure in RNAfold results
    if (!is.null(rnafold_results) && nrow(rnafold_results) > 0) {
      # Find matching sequence in RNAfold results
      seq_idx <- which(rnafold_results$Sequence == seq)
      if (length(seq_idx) > 0) {
        return(rnafold_results$SecondaryStructure[seq_idx[1]])
      }
    }
    
    # If no RNAfold results available, return error
    stop("RNAfold results are required for NeRNA algorithm. Please run secondary structure analysis first.")
  }
  
  # Apply structural transformation
  apply_structural_transformation <- function(seq, structure) {
    chars <- strsplit(seq, "")[[1]]
    struct_chars <- strsplit(structure, "")[[1]]
    
    # Ensure same length
    if (length(chars) != length(struct_chars)) {
      return(seq)  # Return original if lengths don't match
    }
    
    result <- character(length(chars))
    
    for (i in 1:length(chars)) {
      if (struct_chars[i] == ".") {
        # Unpaired base - convert to double letter (e.g., A -> Aa)
        result[i] <- paste0(chars[i], tolower(chars[i]))
      } else {
        # Paired base - keep as is
        result[i] <- chars[i]
      }
    }
    
    return(paste0(result, collapse = ""))
  }
  
  # Convert nucleotide sequence to octal representation
  convert_to_octal <- function(seq) {
    # Parse the sequence to handle double letter bases
    bases <- parse_sequence_with_double_letters(seq)
    
    octal_sequence <- character()
    
    for (base in bases) {
      octal_code <- get_octal_code(base)
      if (!is.null(octal_code)) {
        octal_sequence <- c(octal_sequence, octal_code)
      }
    }
    
    return(paste0(octal_sequence, collapse = ""))
  }
  
  # Parse sequence to correctly handle double letter bases
  parse_sequence_with_double_letters <- function(seq) {
    # Remove spaces
    seq <- gsub(" ", "", seq)
    
    bases <- character()
    i <- 1
    
    while (i <= nchar(seq)) {
      char <- substr(seq, i, i)
      
      # Check if next character is lowercase (double letter base)
      if (i < nchar(seq) && substr(seq, i + 1, i + 1) == tolower(char)) {
        # This is a double letter base (unpaired)
        bases <- c(bases, paste0(char, tolower(char)))
        i <- i + 2  # Skip both characters
      } else {
        # Regular base (paired)
        bases <- c(bases, char)
        i <- i + 1
      }
    }
    
    return(bases)
  }
  
  # Get octal code for a nucleotide base
  get_octal_code <- function(base) {
    octal_map <- list(
      "A" = "000",    # Paired A
      "G" = "001",    # Paired G
      "C" = "010",    # Paired C
      "U" = "100",    # Paired U
      "Aa" = "011",   # Unpaired A
      "Gg" = "110",   # Unpaired G
      "Cc" = "101",   # Unpaired C
      "Uu" = "111"    # Unpaired U
    )
    
    return(octal_map[[base]])
  }
  
  # 0/1 bit dizisini sola döndür (k bit)
  rotate_bits_left <- function(bitstr, k = 1) {
    s <- gsub("\\s+", "", bitstr)
    stopifnot(grepl("^[01]+$", s))
    n <- nchar(s)
    if (n == 0) return(s)
    k <- k %% n
    if (k == 0) return(s)
    paste0(substr(s, k + 1, n), substr(s, 1, k))
  }
  
  # Sekans + dot-bracket -> 3-bit kod -> TEK SEFER rotate(k) -> bazlara çöz
  rna_octal_rotate_once <- function(sequence, structure, k) {
    # --- sözlükler ---
    octal_to_base <- c("000"="A","001"="G","010"="C","100"="U",
                       "011"="a","110"="g","101"="c","111"="u")
    paired_map    <- c(A="000", G="001", C="010", U="100")
    unpaired_map  <- c(A="011", G="110", C="101", U="111")
    
    # --- yardımcılar ---
    clean_seq <- function(seq) {
      s <- toupper(gsub("\\s+", "", seq))
      chartr("T", "U", s)
    }
    clean_struct <- function(struct) {
      s <- gsub("\\s+", "", struct)
      if (grepl("[^().]", s)) stop("Structure yalnızca '.', '(' ve ')' içermeli (MFE olmamalı).")
      s
    }
    encode_bits <- function(sequence, structure) {
      seq <- clean_seq(sequence)
      str <- clean_struct(structure)
      if (nchar(seq) != nchar(str)) {
        stop(sprintf("Uzunluklar farklı: sekans=%d, yapı=%d", nchar(seq), nchar(str)))
      }
      bases   <- strsplit(seq, "", fixed = TRUE)[[1]]
      structc <- strsplit(str, "",  fixed = TRUE)[[1]]
      if (!all(bases %in% c("A","C","G","U"))) {
        bad <- unique(bases[!bases %in% c("A","C","G","U")])
        stop(sprintf("Geçersiz baz(lar): %s", paste(bad, collapse=",")))
      }
      paired <- structc %in% c("(", ")")
      codes <- character(length(bases))
      codes[paired]  <- paired_map[bases[paired]]
      codes[!paired] <- unpaired_map[bases[!paired]]
      paste0(codes, collapse = "")
    }
    decode_bits_to_bases <- function(bitstr) {
      n_bits <- nchar(bitstr)
      if (n_bits %% 3 != 0) stop("Bit uzunluğu 3'ün katı olmalı.")
      i <- seq(1, n_bits, by = 3)
      trip <- substring(bitstr, i, i + 2)
      if (!all(trip %in% names(octal_to_base))) {
        bad <- unique(trip[!trip %in% names(octal_to_base)])
        stop(sprintf("Sözlükte olmayan üçlü(ler): %s", paste(bad, collapse=",")))
      }
      # Convert all to uppercase
      result <- paste0(unname(octal_to_base[trip]), collapse = "")
      toupper(result)
    }
    
    # --- k doğrulama: 1..L (L = sekans uzunluğu, baz sayısı) ---
    L <- nchar(clean_seq(sequence))
    if (length(k) != 1 || is.na(as.integer(k))) stop("k tek bir tam sayı olmalı.")
    k <- as.integer(k)
    if (k < 1L || k > L) stop(sprintf("k 1 ile sekans uzunluğu (%d) arasında olmalı.", L))
    
    # --- akış ---
    bits0   <- encode_bits(sequence, structure)   # uzunluk = 3*L bit
    bitsRot <- rotate_bits_left(bits0, k)         # k bit sola döndür
    decSeq  <- decode_bits_to_bases(bitsRot)      # üçlüleri bazlara çöz
    
    list(
      k = k,
      encoded_bits     = bits0,
      rotated_bits     = bitsRot,
      decoded_sequence = decSeq
    )
  }
  
  # Apply bit-level circular shifting using the correct method
  apply_bit_circular_shift <- function(octal_sequence, shift) {
    # Convert octal sequence to bit string
    bit_string <- convert_octal_to_bits(octal_sequence)
    
    # Apply circular shift using the correct method
    shifted_bits <- rotate_bits_left(bit_string, shift)
    
    # Convert back to octal
    shifted_octal <- convert_bits_to_octal(shifted_bits)
    
    return(shifted_octal)
  }
  
  # Convert octal sequence to bit string
  convert_octal_to_bits <- function(octal_sequence) {
    # Split into 3-bit groups and convert to binary
    n <- nchar(octal_sequence)
    if (n == 0) return("")
    
    bit_string <- ""
    
    for (i in seq(1, n, by = 3)) {
      octal_group <- substr(octal_sequence, i, i + 2)
      if (nchar(octal_group) == 3) {
        bit_group <- convert_octal_group_to_bits(octal_group)
        bit_string <- paste0(bit_string, bit_group)
      }
    }
    
    return(bit_string)
  }
  
  # Convert 3-bit octal group to binary
  convert_octal_group_to_bits <- function(octal_group) {
    # Convert octal to binary (3 bits)
    decimal <- strtoi(octal_group, base = 8)
    
    # Convert decimal to 3-bit binary
    if (decimal == 0) return("000")
    if (decimal == 1) return("001")
    if (decimal == 2) return("010")
    if (decimal == 3) return("011")
    if (decimal == 4) return("100")
    if (decimal == 5) return("101")
    if (decimal == 6) return("110")
    if (decimal == 7) return("111")
    
    return("000")  # Fallback
  }
  
  # Convert bit string back to octal
  convert_bits_to_octal <- function(bit_string) {
    n <- nchar(bit_string)
    if (n == 0) return("")
    
    # Pad with zeros to make length multiple of 3
    remainder <- n %% 3
    if (remainder != 0) {
      padding <- 3 - remainder
      bit_string <- paste0(bit_string, paste(rep("0", padding), collapse = ""))
      n <- nchar(bit_string)
    }
    
    octal_sequence <- ""
    
    for (i in seq(1, n, by = 3)) {
      bit_group <- substr(bit_string, i, i + 2)
      if (nchar(bit_group) == 3) {
        octal_group <- convert_bit_group_to_octal(bit_group)
        octal_sequence <- paste0(octal_sequence, octal_group)
      }
    }
    
    return(octal_sequence)
  }
  
  # Convert 3-bit group to octal
  convert_bit_group_to_octal <- function(bit_group) {
    # Convert binary to octal
    decimal <- strtoi(bit_group, base = 2)
    
    # Convert decimal to octal
    if (decimal == 0) return("000")
    if (decimal == 1) return("001")
    if (decimal == 2) return("010")
    if (decimal == 3) return("011")
    if (decimal == 4) return("100")
    if (decimal == 5) return("101")
    if (decimal == 6) return("110")
    if (decimal == 7) return("111")
    
    return("000")  # Fallback
  }
  
  # Convert octal sequence back to nucleotide sequence
  convert_from_octal <- function(octal_sequence) {
    # Split into 3-digit groups
    n <- nchar(octal_sequence)
    if (n == 0) return("")
    
    bases <- character()
    
    for (i in seq(1, n, by = 3)) {
      octal_group <- substr(octal_sequence, i, i + 2)
      if (nchar(octal_group) == 3) {
        base <- get_base_from_octal(octal_group)
        if (!is.null(base)) {
          # Convert double letter bases back to single letter
          # Aa -> A, Gg -> G, Cc -> C, Uu -> U
          if (base == "Aa") clean_base <- "A"
          else if (base == "Gg") clean_base <- "G"
          else if (base == "Cc") clean_base <- "C"
          else if (base == "Uu") clean_base <- "U"
          else clean_base <- base  # Already single letter
          
          bases <- c(bases, clean_base)
        }
      }
    }
    
    # Join bases without spaces
    return(paste0(bases, collapse = ""))
  }
  
  # Get nucleotide base from octal code
  get_base_from_octal <- function(octal_code) {
    octal_to_base <- list(
      "000" = "A",    # Paired A
      "001" = "G",    # Paired G
      "010" = "C",    # Paired C
      "100" = "U",    # Paired U
      "011" = "Aa",   # Unpaired A
      "110" = "Gg",   # Unpaired G
      "101" = "Cc",   # Unpaired C
      "111" = "Uu"    # Unpaired U
    )
    
    return(octal_to_base[[octal_code]])
  }
  
  # Apply the NeRNA algorithm to all sequences
  nerna_seqs <- lapply(sequences, generate_single_nerna)
  
  return(unlist(nerna_seqs))
}

# Parallel NeRNA implementation
generateNeRNA_parallel <- function(sequences, shifting_size = 1, seed = NULL, rnafold_results = NULL) {
  tryCatch({
    # Setup parallel processing
    setup_parallel_processing()
    
    # Create a helper function that has access to rnafold_results
    process_sequence_nerna <- function(seq_idx) {
      seq <- sequences[seq_idx]
      
      # Get secondary structure from RNAfold results
      if (!is.null(rnafold_results) && nrow(rnafold_results) > 0) {
        # Find matching sequence in RNAfold results
        seq_match_idx <- which(rnafold_results$Sequence == seq)
        
        if (length(seq_match_idx) > 0) {
          structure <- rnafold_results$SecondaryStructure[seq_match_idx[1]]
        } else {
          stop("RNAfold results are required for NeRNA algorithm. Sequence not found in results.")
        }
      } else {
        stop("RNAfold results are required for NeRNA algorithm.")
      }
      
      # Apply NeRNA algorithm
      result <- rna_octal_rotate_once(seq, structure, shifting_size)
      
      return(result$decoded_sequence)
    }
    
    # Process sequences in parallel
    message("Starting parallel NeRNA processing...")
    
    nerna_seqs <- future_sapply(1:length(sequences), process_sequence_nerna, 
                                 future.seed = if(!is.null(seed)) seed else TRUE)
    
    message(sprintf("Parallel NeRNA completed: %d sequences processed", length(nerna_seqs)))
    
    return(nerna_seqs)
    
  }, error = function(e) {
    message("Error in parallel NeRNA: ", e$message)
    message("Falling back to sequential processing...")
    return(generateNeRNA_sequential(sequences, shifting_size, seed, rnafold_results))
  })
}

# Test function for debugging NeRNA
test_nerna_debug <- function() {
  # Test with simple sequence
  test_seq <- "AUGC"
  test_structure <- "(..)"
  
  cat("=== NeRNA Debug Test (New System) ===\n")
  cat("Original sequence:", test_seq, "\n")
  cat("Structure:", test_structure, "\n")
  
  # Apply structural transformation
  transformed <- apply_structural_transformation(test_seq, test_structure)
  cat("Transformed:", transformed, "\n")
  
  # Parse the transformed sequence
  parsed_bases <- parse_sequence_with_double_letters(transformed)
  cat("Parsed bases:", paste(parsed_bases, collapse = " "), "\n")
  
  # Test each base conversion
  for (base in parsed_bases) {
    octal_code <- get_octal_code(base)
    cat("Base:", base, "-> Octal:", octal_code, "\n")
  }
  
  # Convert to octal
  octal_seq <- convert_to_octal(transformed)
  cat("Octal:", octal_seq, "\n")
  
  # Test octal to bits conversion
  bit_string <- convert_octal_to_bits(octal_seq)
  cat("Bit string:", bit_string, "\n")
  cat("Bit string length:", nchar(bit_string), "\n")
  
  # Test each octal group conversion
  cat("Testing octal group conversions:\n")
  for (i in seq(1, nchar(octal_seq), by = 3)) {
    octal_group <- substr(octal_seq, i, i + 2)
    if (nchar(octal_group) == 3) {
      bit_group <- convert_octal_group_to_bits(octal_group)
      cat("Octal group:", octal_group, "-> Bit group:", bit_group, "\n")
    }
  }
  
  # Test bit shift manually
  n_bits <- nchar(bit_string)
  if (n_bits > 0) {
    # Manual shift: move first bit to end
    shifted_bits <- paste0(substr(bit_string, 2, n_bits), substr(bit_string, 1, 1))
    cat("Manual shifted bits:", shifted_bits, "\n")
    cat("Shifted bits length:", nchar(shifted_bits), "\n")
    
    # Convert back to octal manually
    manual_octal <- convert_bits_to_octal(shifted_bits)
    cat("Manual shifted octal:", manual_octal, "\n")
  }
  
  # Apply bit shift using function
  shifted_octal <- apply_bit_circular_shift(octal_seq, 1)
  cat("Function shifted octal:", shifted_octal, "\n")
  
  # Convert back to nucleotides
  result <- convert_from_octal(shifted_octal)
  cat("Result:", result, "\n")
  
  return(result)
}

# Calculate structural similarity between original and NeRNA sequences
calculate_structural_similarity <- function(original_seq, nerna_seq) {
  # Calculate various similarity metrics
  metrics <- list()
  
  # 1. Sequence identity
  metrics$sequence_identity <- calculate_sequence_identity(original_seq, nerna_seq)
  
  # 2. Dinucleotide composition similarity
  metrics$dinucleotide_similarity <- calculate_dinucleotide_similarity(original_seq, nerna_seq)
  
  # 3. GC content similarity
  metrics$gc_similarity <- calculate_gc_similarity(original_seq, nerna_seq)
  
  # 4. Structural motif preservation
  metrics$motif_preservation <- calculate_motif_preservation(original_seq, nerna_seq)
  
  # 5. Overall structural similarity score
  metrics$overall_similarity <- calculate_overall_similarity(metrics)
  
  return(metrics)
}

# Calculate sequence identity percentage
calculate_sequence_identity <- function(seq1, seq2) {
  chars1 <- strsplit(seq1, "")[[1]]
  chars2 <- strsplit(seq2, "")[[1]]
  
  if (length(chars1) != length(chars2)) {
    return(0)
  }
  
  identical_positions <- sum(chars1 == chars2)
  return(identical_positions / length(chars1) * 100)
}

# Calculate dinucleotide composition similarity
calculate_dinucleotide_similarity <- function(seq1, seq2) {
  # Get dinucleotide frequencies
  freq1 <- get_dinucleotide_frequencies(seq1)
  freq2 <- get_dinucleotide_frequencies(seq2)
  
  # Calculate cosine similarity
  common_dinucs <- intersect(names(freq1), names(freq2))
  
  if (length(common_dinucs) == 0) return(0)
  
  dot_product <- sum(freq1[common_dinucs] * freq2[common_dinucs])
  norm1 <- sqrt(sum(freq1^2))
  norm2 <- sqrt(sum(freq2^2))
  
  if (norm1 == 0 || norm2 == 0) return(0)
  
  return((dot_product / (norm1 * norm2)) * 100)
}

# Get dinucleotide frequencies
get_dinucleotide_frequencies <- function(seq) {
    chars <- strsplit(seq, "")[[1]]
    n <- length(chars)
    
  if (n < 2) return(numeric(0))
  
  dinucs <- paste0(chars[-n], chars[-1])
  freq_table <- table(dinucs)
  return(as.numeric(freq_table) / length(dinucs))
}

# Calculate GC content similarity
calculate_gc_similarity <- function(seq1, seq2) {
  gc1 <- sum(gregexpr("[GC]", seq1)[[1]] > 0) / nchar(seq1) * 100
  gc2 <- sum(gregexpr("[GC]", seq2)[[1]] > 0) / nchar(seq2) * 100
  
  # Return similarity as percentage (100 - difference)
  return(max(0, 100 - abs(gc1 - gc2)))
}

# Calculate structural motif preservation (simplified for NeRNA)
calculate_motif_preservation <- function(seq1, seq2) {
  # Simplified motif preservation calculation
  # For NeRNA, we'll use a simple length-based similarity
  len1 <- nchar(seq1)
  len2 <- nchar(seq2)
  
  if (len1 == 0 || len2 == 0) return(0)
  
  # Calculate length similarity
  length_similarity <- 100 - abs(len1 - len2) / max(len1, len2) * 100
  
  # For NeRNA, we assume some structural preservation due to the algorithm
  # This is a simplified metric
  return(max(0, length_similarity * 0.7))
}

# Calculate overall similarity score
calculate_overall_similarity <- function(metrics) {
  # Weighted average of all metrics
  weights <- c(0.3, 0.25, 0.2, 0.25)  # Weights for each metric
  scores <- c(
    metrics$sequence_identity,
    metrics$dinucleotide_similarity,
    metrics$gc_similarity,
    metrics$motif_preservation
  )
  
  return(sum(weights * scores))
}

# Dinükleotit-korumalı karıştırma (Altschul & Erickson) — Euler izi ile
dinucleotideShuffle <- function(sequences, seed = NULL, use_parallel = TRUE) {
  if (!is.null(seed)) set.seed(seed)
  
  # Decide whether to use parallel processing
  n_sequences <- length(sequences)
  use_parallel <- use_parallel && n_sequences >= 10
  
  if (use_parallel) {
    message(sprintf("Using PARALLEL processing for Dinucleotide Shuffling (%d sequences)", n_sequences))
    return(dinucleotideShuffle_parallel(sequences, seed))
  } else {
    message(sprintf("Using SEQUENTIAL processing for Dinucleotide Shuffling (%d sequences)", n_sequences))
    return(dinucleotideShuffle_sequential(sequences, seed))
  }
}

# Sequential Dinucleotide Shuffling (original implementation)
dinucleotideShuffle_sequential <- function(sequences, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)

  # Yardımcı: tek bir sekansı karıştır
  shuffle_one <- function(seq) {
    # Temizlik ve doğrulamalar
    s <- toupper(gsub("\\s+", "", seq))
    s <- chartr("T", "U", s)
    if (!grepl("^[ACGU]+$", s)) {
      stop("Dizide yalnızca A,C,G,U (veya T->U dönüştürülebilir) olmalı.")
    }
    n <- nchar(s)
    if (n < 2) return(s)

    chars <- strsplit(s, "", fixed = TRUE)[[1]]
    
    # DEBUG: Print original composition
    cat("=== DEBUG Dinucleotide Shuffling (Hierholzer) ===\n")
    cat("Original:", seq, "\n")
    cat("Original GC count:", sum(chars %in% c("G", "C")), "\n")
    cat("Original GC ratio:", sum(chars %in% c("G", "C")) / n, "\n")
    
    # Dinükleotit sayımları
    dinucs <- paste0(chars[-n], chars[-1])
    all_dinucs <- as.vector(outer(c("A","C","G","U"), c("A","C","G","U"), paste0))
    counts <- integer(length(all_dinucs)); names(counts) <- all_dinucs
    tab <- table(dinucs); counts[names(tab)] <- as.integer(tab)
    
    cat("Dinucleotide counts:", paste(names(counts)[counts > 0], "=", counts[counts > 0], collapse = ", "), "\n")

    # Çoklu graf adjacency (çıktı kenar hedefleri listesi)
    nodes <- c("A","C","G","U")
    adj <- setNames(vector("list", length(nodes)), nodes)
    for (d in names(counts)) {
      cnt <- counts[[d]]
      if (cnt > 0) {
        from <- substr(d,1,1); to <- substr(d,2,2)
        adj[[from]] <- c(adj[[from]], rep(to, cnt))
      }
    }

    # Euler izi için başlangıç: özgün dizinin ilk harfi uygundur
    start <- chars[1]

    # Hierholzer algoritması
    stack <- c(start)
    circuit <- character(0)

    while (length(stack) > 0) {
      v <- stack[length(stack)]
      if (length(adj[[v]]) > 0) {
        # Rastgele bir çıkış kenarı seç (çoklu kenarlar arasından)
        idx <- sample.int(length(adj[[v]]), 1)
        w <- adj[[v]][idx]
        # Kenarı tüket
        adj[[v]] <- adj[[v]][-idx]
        stack <- c(stack, w)
      } else {
        # Geri sarma: devreye düğümü ekle
        circuit <- c(circuit, v)
        stack <- stack[-length(stack)]
      }
    }

    # Hierholzer devresi ters döner; yolu elde etmek için çevir
    path <- rev(circuit)

    # Yol uzunluğu n düğüm olmalı (n-1 kenar) → sekans uzunluğu n
    if (length(path) != n) {
      stop("Euler izi uzunluğu beklenenle eşleşmiyor (içsel tutarsızlık).")
    }

    shuffled_seq <- paste0(path, collapse = "")
    cat("Shuffled:", shuffled_seq, "\n")
    cat("Shuffled GC count:", sum(path %in% c("G", "C")), "\n")
    cat("Shuffled GC ratio:", sum(path %in% c("G", "C")) / n, "\n")
    cat("=====================================\n")
    
    return(shuffled_seq)
  }

  vapply(sequences, shuffle_one, FUN.VALUE = character(1))
}

# Parallel Dinucleotide Shuffling implementation
dinucleotideShuffle_parallel <- function(sequences, seed = NULL) {
  tryCatch({
    # Setup parallel processing
    setup_parallel_processing()
    
    # Define the shuffle function (same as sequential)
    shuffle_one <- function(seq) {
      # Temizlik ve doğrulamalar
      s <- toupper(gsub("\\s+", "", seq))
      s <- chartr("T", "U", s)
      if (!grepl("^[ACGU]+$", s)) {
        stop("Dizide yalnızca A,C,G,U (veya T->U dönüştürülebilir) olmalı.")
      }
      n <- nchar(s)
      if (n < 2) return(s)

      chars <- strsplit(s, "", fixed = TRUE)[[1]]
      
      # Dinükleotit sayımları
      dinucs <- paste0(chars[-n], chars[-1])
      all_dinucs <- as.vector(outer(c("A","C","G","U"), c("A","C","G","U"), paste0))
      counts <- integer(length(all_dinucs)); names(counts) <- all_dinucs
      tab <- table(dinucs); counts[names(tab)] <- as.integer(tab)

      # Çoklu graf adjacency
      nodes <- c("A","C","G","U")
      adj <- setNames(vector("list", length(nodes)), nodes)
      for (d in names(counts)) {
        cnt <- counts[[d]]
        if (cnt > 0) {
          from <- substr(d,1,1); to <- substr(d,2,2)
          adj[[from]] <- c(adj[[from]], rep(to, cnt))
        }
      }

      # Euler izi için başlangıç
      start <- chars[1]

      # Hierholzer algoritması
      stack <- c(start)
      circuit <- character(0)

      while (length(stack) > 0) {
        v <- stack[length(stack)]
        if (length(adj[[v]]) > 0) {
          idx <- sample.int(length(adj[[v]]), 1)
          w <- adj[[v]][idx]
          adj[[v]] <- adj[[v]][-idx]
          stack <- c(stack, w)
        } else {
          circuit <- c(circuit, v)
          stack <- stack[-length(stack)]
        }
      }

      path <- rev(circuit)

      if (length(path) != n) {
        stop("Euler izi uzunluğu beklenenle eşleşmiyor (içsel tutarsızlık).")
      }

      return(paste0(path, collapse = ""))
    }
    
    # Process sequences in parallel
    message("Starting parallel Dinucleotide Shuffling...")
    
    shuffled_seqs <- future_sapply(sequences, shuffle_one, 
                                    future.seed = if(!is.null(seed)) seed else TRUE,
                                    USE.NAMES = FALSE)
    
    message(sprintf("Parallel Dinucleotide Shuffling completed: %d sequences processed", length(shuffled_seqs)))
    
    return(shuffled_seqs)
    
  }, error = function(e) {
    message("Error in parallel Dinucleotide Shuffling: ", e$message)
    message("Falling back to sequential processing...")
    return(dinucleotideShuffle_sequential(sequences, seed))
  })
}

# Random Shuffling Method with seed parameter
randomShuffle <- function(sequences, seed = NULL, use_parallel = TRUE) {
  # Set seed if provided
  if (!is.null(seed)) {
    set.seed(seed)
  }
  
  # Decide whether to use parallel processing
  n_sequences <- length(sequences)
  use_parallel <- use_parallel && n_sequences >= 10
  
  if (use_parallel) {
    message(sprintf("Using PARALLEL processing for Random Shuffling (%d sequences)", n_sequences))
    return(randomShuffle_parallel(sequences, seed))
  } else {
    message(sprintf("Using SEQUENTIAL processing for Random Shuffling (%d sequences)", n_sequences))
    return(randomShuffle_sequential(sequences, seed))
  }
}

# Sequential Random Shuffling (original implementation)
randomShuffle_sequential <- function(sequences, seed = NULL) {
  # Set seed if provided
  if (!is.null(seed)) {
    set.seed(seed)
  }
  
  shuffled <- lapply(sequences, function(seq) {
    chars <- strsplit(seq, "")[[1]]
    paste0(sample(chars), collapse = "")
  })
  
  # Make sure to return a character vector, not a list
  return(unlist(shuffled))
}

# Parallel Random Shuffling implementation
randomShuffle_parallel <- function(sequences, seed = NULL) {
  tryCatch({
    # Setup parallel processing
    setup_parallel_processing()
    
    # Shuffle function
    shuffle_one <- function(seq) {
      chars <- strsplit(seq, "")[[1]]
      paste0(sample(chars), collapse = "")
    }
    
    # Process sequences in parallel
    message("Starting parallel Random Shuffling...")
    
    shuffled_seqs <- future_sapply(sequences, shuffle_one, 
                                    future.seed = if(!is.null(seed)) seed else TRUE,
                                    USE.NAMES = FALSE)
    
    message(sprintf("Parallel Random Shuffling completed: %d sequences processed", length(shuffled_seqs)))
    
    return(shuffled_seqs)
    
  }, error = function(e) {
    message("Error in parallel Random Shuffling: ", e$message)
    message("Falling back to sequential processing...")
    return(randomShuffle_sequential(sequences, seed))
  })
}

# Calculate GC Content
calculateGC <- function(sequence) {
  bases <- strsplit(sequence, "")[[1]]
  gc_count <- sum(bases %in% c("G", "C"))
  gc_percentage <- (gc_count / length(bases)) * 100
  return(round(gc_percentage, 1))
}

# Calculate paired bases from structure notation
calculatePaired <- function(structure) {
  if (is.null(structure)) return(list(paired = 0, unpaired = 0))
  
  paired_count <- sum(grepl("[\\(\\)]", strsplit(structure, "")[[1]]))
  total_length <- nchar(structure)
  unpaired_count <- total_length - paired_count
  
  return(list(paired = paired_count, unpaired = unpaired_count))
}

# ============================================================
# COMPREHENSIVE EXCEL EXPORT
# ============================================================

# Comprehensive Excel export with detailed analysis
export_comprehensive_results <- function(original_data, nerna_data = NULL, dinuc_data = NULL, random_data = NULL, similarity_metrics, rnafold_results = NULL, filename) {
  tryCatch({
    message("=== EXCEL EXPORT DEBUG START ===")
    message("Original data: ", if(is.null(original_data)) "NULL" else paste("nrow =", nrow(original_data)))
    message("NeRNA data: ", if(is.null(nerna_data)) "NULL" else paste("nrow =", nrow(nerna_data)))
    message("Dinuc data: ", if(is.null(dinuc_data)) "NULL" else paste("nrow =", nrow(dinuc_data)))
    message("Random data: ", if(is.null(random_data)) "NULL" else paste("nrow =", nrow(random_data)))
    message("Similarity metrics: ", if(is.null(similarity_metrics)) "NULL" else paste("nrow =", nrow(similarity_metrics)))
    message("RNAfold results: ", if(is.null(rnafold_results)) "NULL" else paste("nrow =", nrow(rnafold_results)))
    
    # Convert lists to data frames if needed
    message("=== CONVERTING DATA STRUCTURES ===")
    if (!is.null(nerna_data) && is.list(nerna_data) && !is.data.frame(nerna_data)) {
      message("Converting nerna_data from list to data frame...")
      nerna_data <- tryCatch(as.data.frame(nerna_data, stringsAsFactors = FALSE), 
                              error = function(e) { message("Failed: ", e$message); NULL })
    }
    if (!is.null(dinuc_data) && is.list(dinuc_data) && !is.data.frame(dinuc_data)) {
      message("Converting dinuc_data from list to data frame...")
      dinuc_data <- tryCatch(as.data.frame(dinuc_data, stringsAsFactors = FALSE), 
                              error = function(e) { message("Failed: ", e$message); NULL })
    }
    if (!is.null(random_data) && is.list(random_data) && !is.data.frame(random_data)) {
      message("Converting random_data from list to data frame...")
      random_data <- tryCatch(as.data.frame(random_data, stringsAsFactors = FALSE), 
                              error = function(e) { message("Failed: ", e$message); NULL })
    }
    message("Conversions complete!")
    
    wb <- createWorkbook()
    
    # ============================================================
    # SHEET 1: EXECUTIVE SUMMARY
    # ============================================================
    addWorksheet(wb, "Executive Summary")
    
    # Calculate statistics safely
    total_original <- if (!is.null(original_data)) nrow(original_data) else 0
    total_nerna <- if (!is.null(nerna_data)) nrow(nerna_data) else 0
    total_dinuc <- if (!is.null(dinuc_data)) nrow(dinuc_data) else 0
    total_random <- if (!is.null(random_data)) nrow(random_data) else 0
    
    avg_orig_length <- if (!is.null(original_data) && nrow(original_data) > 0) {
      round(mean(nchar(as.character(original_data$Sequence))), 2)
    } else NA
    
    avg_orig_gc <- if (!is.null(original_data) && nrow(original_data) > 0) {
      round(mean(sapply(as.character(original_data$Sequence), function(s) {
        chars <- strsplit(s, "")[[1]]
        sum(chars %in% c("G", "C")) / length(chars) * 100
      })), 2)
    } else NA
    
    avg_seq_identity <- if (!is.null(similarity_metrics) && "SequenceIdentity" %in% names(similarity_metrics)) {
      round(mean(as.numeric(similarity_metrics$SequenceIdentity), na.rm = TRUE), 2)
    } else NA
    
    avg_dinuc_sim <- if (!is.null(similarity_metrics) && "DinucleotideSimilarity" %in% names(similarity_metrics)) {
      round(mean(as.numeric(similarity_metrics$DinucleotideSimilarity), na.rm = TRUE), 2)
    } else NA
    
    avg_gc_sim <- if (!is.null(similarity_metrics) && "GCSimilarity" %in% names(similarity_metrics)) {
      round(mean(as.numeric(similarity_metrics$GCSimilarity), na.rm = TRUE), 2)
    } else NA
    
    avg_motif_pres <- if (!is.null(similarity_metrics) && "MotifPreservation" %in% names(similarity_metrics)) {
      round(mean(as.numeric(similarity_metrics$MotifPreservation), na.rm = TRUE), 2)
    } else NA
    
    avg_overall_sim <- if (!is.null(similarity_metrics) && "OverallSimilarity" %in% names(similarity_metrics)) {
      round(mean(as.numeric(similarity_metrics$OverallSimilarity), na.rm = TRUE), 2)
    } else NA
    
    summary_stats <- data.frame(
      Metric = c("Report Generated", 
                 "Total Original Sequences",
                 "Total NeRNA Sequences",
                 "Total Dinucleotide Sequences",
                 "Total Random Sequences",
                 "Avg Original Length (nt)",
                 "Min Original Length (nt)", 
                 "Max Original Length (nt)",
                 "Avg Original GC%"),
      Value = c(format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
                total_original,
                total_nerna,
                total_dinuc,
                total_random,
                avg_orig_length,
                if (!is.null(original_data) && nrow(original_data) > 0) min(nchar(as.character(original_data$Sequence))) else NA,
                if (!is.null(original_data) && nrow(original_data) > 0) max(nchar(as.character(original_data$Sequence))) else NA,
                avg_orig_gc),
      stringsAsFactors = FALSE
    )
    writeData(wb, "Executive Summary", summary_stats)
    
    # Add header style
    headerStyle <- createStyle(fontSize = 12, fontColour = "#FFFFFF", 
                              halign = "center", fgFill = "#4F81BD",
                              textDecoration = "bold", border = "TopBottomLeftRight")
    addStyle(wb, "Executive Summary", headerStyle, rows = 1, cols = 1:2, gridExpand = TRUE)
    
    # ============================================================
    # SHEET 2: ORIGINAL SEQUENCES WITH DETAILED STATS
    # ============================================================
    if (!is.null(original_data) && is.data.frame(original_data) && nrow(original_data) > 0) {
      addWorksheet(wb, "Original Sequences")
      
      # Add detailed statistics for each sequence
      orig_detailed <- data.frame(
        Accession = original_data$Accession,
        Sequence = as.character(original_data$Sequence),
        Length = nchar(as.character(original_data$Sequence)),
        GC_Content = sapply(as.character(original_data$Sequence), function(s) {
          chars <- strsplit(s, "")[[1]]
          round(sum(chars %in% c("G", "C")) / length(chars) * 100, 2)
        }),
        stringsAsFactors = FALSE
      )
      
      # Add secondary structure if available
      if (!is.null(rnafold_results) && is.data.frame(rnafold_results) && nrow(rnafold_results) > 0) {
        orig_detailed$Structure <- sapply(original_data$Accession, function(acc) {
          idx <- which(rnafold_results$Accession == acc)
          if (length(idx) > 0) {
            as.character(rnafold_results$Structure[idx[1]])
          } else {
            "N/A"
          }
        })
        orig_detailed$MFE <- sapply(original_data$Accession, function(acc) {
          idx <- which(rnafold_results$Accession == acc)
          if (length(idx) > 0) {
            rnafold_results$MFE[idx[1]]
          } else {
            NA
          }
        })
      } else {
        orig_detailed$Structure <- "N/A"
        orig_detailed$MFE <- NA
      }
      
      writeData(wb, "Original Sequences", orig_detailed)
      addStyle(wb, "Original Sequences", headerStyle, rows = 1, cols = 1:ncol(orig_detailed), gridExpand = TRUE)
    }
    
    # ============================================================
    # SHEET 3: NeRNA SEQUENCES
    # ============================================================
    if (!is.null(nerna_data) && is.data.frame(nerna_data) && nrow(nerna_data) > 0) {
      addWorksheet(wb, "NeRNA Sequences")
      
      seq_col_name <- if ("NegativeSequence" %in% names(nerna_data)) {
        "NegativeSequence"
      } else if ("Sequence" %in% names(nerna_data)) {
        "Sequence"
      } else {
        NA_character_
      }
      
      if (!is.null(seq_col_name) && !is.na(seq_col_name) && seq_col_name %in% names(nerna_data)) {
        nerna_detailed <- data.frame(
          OriginalAccession = if ("OriginalAccession" %in% names(nerna_data)) nerna_data$OriginalAccession else NA,
          NegativeSequence = as.character(nerna_data[[seq_col_name]]),
          Length = nchar(as.character(nerna_data[[seq_col_name]])),
          GC_Content = sapply(as.character(nerna_data[[seq_col_name]]), function(s) {
            chars <- strsplit(s, "")[[1]]
            round(sum(chars %in% c("G", "C")) / length(chars) * 100, 2)
          }),
          ShiftAmount = if ("ShiftAmount" %in% names(nerna_data)) nerna_data$ShiftAmount else NA,
          stringsAsFactors = FALSE
        )
        
        # Add secondary structure by calling RNAfold
        message("Predicting secondary structures for NeRNA sequences...")
        nerna_detailed$Structure <- sapply(as.character(nerna_data[[seq_col_name]]), function(seq) {
          result <- tryCatch({
            quick_rnafold(seq)
          }, error = function(e) {
            list(structure = "N/A", mfe = NA)
          })
          result$structure
        })
        
        nerna_detailed$MFE <- sapply(as.character(nerna_data[[seq_col_name]]), function(seq) {
          result <- tryCatch({
            quick_rnafold(seq)
          }, error = function(e) {
            list(structure = "N/A", mfe = NA)
          })
          result$mfe
        })
        
        writeData(wb, "NeRNA Sequences", nerna_detailed)
        if (nrow(nerna_detailed) > 0) {
          addStyle(wb, "NeRNA Sequences", headerStyle, rows = 1, cols = 1:ncol(nerna_detailed), gridExpand = TRUE)
        }
      }
    }
    
    # ============================================================
    # SHEET 4: DINUCLEOTIDE SEQUENCES
    # ============================================================
    if (!is.null(dinuc_data) && nrow(dinuc_data) > 0) {
      addWorksheet(wb, "Dinucleotide Sequences")
      
      seq_col_name <- if ("NegativeSequence" %in% names(dinuc_data)) {
        "NegativeSequence"
      } else if ("Sequence" %in% names(dinuc_data)) {
        "Sequence"
      } else {
        NA_character_
      }
      
      if (!is.null(seq_col_name) && !is.na(seq_col_name) && seq_col_name %in% names(dinuc_data)) {
        dinuc_detailed <- data.frame(
          OriginalAccession = if ("OriginalAccession" %in% names(dinuc_data)) dinuc_data$OriginalAccession else NA,
          NegativeSequence = as.character(dinuc_data[[seq_col_name]]),
          Length = nchar(as.character(dinuc_data[[seq_col_name]])),
          GC_Content = sapply(as.character(dinuc_data[[seq_col_name]]), function(s) {
            chars <- strsplit(s, "")[[1]]
            round(sum(chars %in% c("G", "C")) / length(chars) * 100, 2)
          }),
          stringsAsFactors = FALSE
        )
        
        # Add secondary structure by calling RNAfold
        message("Predicting secondary structures for Dinucleotide sequences...")
        dinuc_detailed$Structure <- sapply(as.character(dinuc_data[[seq_col_name]]), function(seq) {
          result <- tryCatch({
            quick_rnafold(seq)
          }, error = function(e) {
            list(structure = "N/A", mfe = NA)
          })
          result$structure
        })
        
        dinuc_detailed$MFE <- sapply(as.character(dinuc_data[[seq_col_name]]), function(seq) {
          result <- tryCatch({
            quick_rnafold(seq)
          }, error = function(e) {
            list(structure = "N/A", mfe = NA)
          })
          result$mfe
        })
        
        writeData(wb, "Dinucleotide Sequences", dinuc_detailed)
        if (nrow(dinuc_detailed) > 0) {
          addStyle(wb, "Dinucleotide Sequences", headerStyle, rows = 1, cols = 1:ncol(dinuc_detailed), gridExpand = TRUE)
        }
      }
    }
    
    # ============================================================
    # SHEET 5: RANDOM SEQUENCES
    # ============================================================
    if (!is.null(random_data) && is.data.frame(random_data) && nrow(random_data) > 0) {
      addWorksheet(wb, "Random Sequences")
      
      seq_col_name <- if ("NegativeSequence" %in% names(random_data)) {
        "NegativeSequence"
      } else if ("Sequence" %in% names(random_data)) {
        "Sequence"
      } else {
        NA_character_
      }
      
      if (!is.null(seq_col_name) && !is.na(seq_col_name) && seq_col_name %in% names(random_data)) {
        random_detailed <- data.frame(
          OriginalAccession = if ("OriginalAccession" %in% names(random_data)) random_data$OriginalAccession else NA,
          NegativeSequence = as.character(random_data[[seq_col_name]]),
          Length = nchar(as.character(random_data[[seq_col_name]])),
          GC_Content = sapply(as.character(random_data[[seq_col_name]]), function(s) {
            chars <- strsplit(s, "")[[1]]
            round(sum(chars %in% c("G", "C")) / length(chars) * 100, 2)
          }),
          stringsAsFactors = FALSE
        )
        
        # Add secondary structure by calling RNAfold
        message("Predicting secondary structures for Random sequences...")
        random_detailed$Structure <- sapply(as.character(random_data[[seq_col_name]]), function(seq) {
          result <- tryCatch({
            quick_rnafold(seq)
          }, error = function(e) {
            list(structure = "N/A", mfe = NA)
          })
          result$structure
        })
        
        random_detailed$MFE <- sapply(as.character(random_data[[seq_col_name]]), function(seq) {
          result <- tryCatch({
            quick_rnafold(seq)
          }, error = function(e) {
            list(structure = "N/A", mfe = NA)
          })
          result$mfe
        })
        
        writeData(wb, "Random Sequences", random_detailed)
        if (nrow(random_detailed) > 0) {
          addStyle(wb, "Random Sequences", headerStyle, rows = 1, cols = 1:ncol(random_detailed), gridExpand = TRUE)
        }
      }
    }
    
    # ============================================================
    # SHEET 7: SIMILARITY METRICS (ALL METRICS)
    # ============================================================
    if (!is.null(similarity_metrics) && is.data.frame(similarity_metrics) && nrow(similarity_metrics) > 0) {
      addWorksheet(wb, "Similarity Metrics")
      writeData(wb, "Similarity Metrics", similarity_metrics)
      addStyle(wb, "Similarity Metrics", headerStyle, rows = 1, cols = 1:ncol(similarity_metrics), gridExpand = TRUE)
    }
    
    # ============================================================
    # SHEET 5: SECONDARY STRUCTURES (IF AVAILABLE)
    # ============================================================
    if (!is.null(rnafold_results) && nrow(rnafold_results) > 0) {
      addWorksheet(wb, "Secondary Structures")
      writeData(wb, "Secondary Structures", rnafold_results)
      addStyle(wb, "Secondary Structures", headerStyle, rows = 1, cols = 1:ncol(rnafold_results), gridExpand = TRUE)
    }
    
    # ============================================================
    # SHEET 6: DINUCLEOTIDE COMPOSITION
    # ============================================================
    if (!is.null(original_data) && is.data.frame(original_data) && nrow(original_data) > 0) {
      addWorksheet(wb, "Dinucleotide Composition")
      
      # Calculate dinucleotide counts for each sequence
      dinuc_comp_list <- lapply(1:nrow(original_data), function(i) {
        seq <- as.character(original_data$Sequence[i])
        chars <- strsplit(seq, "")[[1]]
        n <- length(chars)
        
        if (n < 2) {
          # Return zeros for sequences too short
          return(data.frame(
            Accession = original_data$Accession[i],
            AA = 0, AC = 0, AG = 0, AU = 0,
            CA = 0, CC = 0, CG = 0, CU = 0,
            GA = 0, GC = 0, GG = 0, GU = 0,
            UA = 0, UC = 0, UG = 0, UU = 0,
            stringsAsFactors = FALSE
          ))
        }
        
        # Get dinucleotides
        dinucs <- paste0(chars[-n], chars[-1])
        dinuc_table <- table(factor(dinucs, levels = c("AA", "AC", "AG", "AU",
                                                        "CA", "CC", "CG", "CU",
                                                        "GA", "GC", "GG", "GU",
                                                        "UA", "UC", "UG", "UU")))
        
        # Create data frame
        result <- data.frame(
          Accession = original_data$Accession[i],
          stringsAsFactors = FALSE
        )
        
        for (dn in names(dinuc_table)) {
          result[[dn]] <- as.numeric(dinuc_table[dn])
        }
        
        return(result)
      })
      
      dinuc_comp_table <- do.call(rbind, dinuc_comp_list)
      writeData(wb, "Dinucleotide Composition", dinuc_comp_table)
      addStyle(wb, "Dinucleotide Composition", headerStyle, rows = 1, cols = 1:ncol(dinuc_comp_table), gridExpand = TRUE)
    }
    
    
    # ============================================================
    # SHEET 9: METHOD COMPARISON
    # ============================================================
    addWorksheet(wb, "Method Comparison")
    
    method_comparison <- data.frame(
      Method = character(),
      Count = numeric(),
      Avg_Length = numeric(),
      Avg_GC = numeric(),
      stringsAsFactors = FALSE
    )
    
    # Add Original stats
    if (!is.null(original_data) && is.data.frame(original_data) && nrow(original_data) > 0) {
      seq_col_name <- if ("Sequence" %in% names(original_data)) {
        "Sequence"
      } else {
        NA_character_
      }
      
      if (!is.null(seq_col_name) && !is.na(seq_col_name) && seq_col_name %in% names(original_data)) {
        avg_gc <- mean(sapply(as.character(original_data[[seq_col_name]]), function(s) {
          chars <- strsplit(s, "")[[1]]
          sum(chars %in% c("G", "C")) / length(chars) * 100
        }))
        
        method_comparison <- rbind(method_comparison, data.frame(
          Method = "Original",
          Count = nrow(original_data),
          Avg_Length = round(mean(nchar(as.character(original_data[[seq_col_name]]))), 2),
          Avg_GC = round(avg_gc, 2),
          stringsAsFactors = FALSE
        ))
      }
    }
    
    # Add NeRNA stats
    message("=== NERNA DEBUG ===")
    if (!is.null(nerna_data)) {
      message("NeRNA data class: ", paste(class(nerna_data), collapse = ", "))
      
      # If nerna_data is a list, try to convert to data frame
      if (is.list(nerna_data) && !is.data.frame(nerna_data)) {
        message("NeRNA data is a list, attempting to convert to data frame...")
        tryCatch({
          nerna_data <- as.data.frame(nerna_data, stringsAsFactors = FALSE)
          message("Conversion successful!")
        }, error = function(e) {
          message("Conversion failed: ", e$message)
        })
      }
    }
    
    if (!is.null(nerna_data) && is.data.frame(nerna_data) && nrow(nerna_data) > 0) {
      seq_col_name <- if ("NegativeSequence" %in% names(nerna_data)) {
        "NegativeSequence"
      } else if ("Sequence" %in% names(nerna_data)) {
        "Sequence"
      } else {
        NA_character_
      }
      
      if (!is.null(seq_col_name) && seq_col_name %in% names(nerna_data)) {
        avg_gc <- mean(sapply(as.character(nerna_data[[seq_col_name]]), function(s) {
          chars <- strsplit(s, "")[[1]]
          sum(chars %in% c("G", "C")) / length(chars) * 100
        }))
        
        method_comparison <- rbind(method_comparison, data.frame(
          Method = "NeRNA",
          Count = nrow(nerna_data),
          Avg_Length = round(mean(nchar(as.character(nerna_data[[seq_col_name]]))), 2),
          Avg_GC = round(avg_gc, 2),
          stringsAsFactors = FALSE
        ))
      }
    }
    
    # Add Dinucleotide stats
    message("=== DINUCLEOTIDE DEBUG ===")
    message("Dinuc data is null: ", is.null(dinuc_data))
    if (!is.null(dinuc_data)) {
      message("Dinuc data class: ", paste(class(dinuc_data), collapse = ", "))
      message("Dinuc data length: ", length(dinuc_data))
      
      # If dinuc_data is a list, try to convert to data frame
      if (is.list(dinuc_data) && !is.data.frame(dinuc_data)) {
        message("Dinuc data is a list, attempting to convert to data frame...")
        tryCatch({
          dinuc_data <- as.data.frame(dinuc_data, stringsAsFactors = FALSE)
          message("Conversion successful!")
        }, error = function(e) {
          message("Conversion failed: ", e$message)
        })
      }
      
      message("Dinuc data nrow: ", if(is.null(nrow(dinuc_data))) "NULL" else nrow(dinuc_data))
      message("Dinuc data names: ", if(is.null(names(dinuc_data))) "NULL" else paste(names(dinuc_data), collapse = ", "))
    }
    
    if (!is.null(dinuc_data) && is.data.frame(dinuc_data) && nrow(dinuc_data) > 0) {
      message("Dinuc data is valid data frame with rows: ", nrow(dinuc_data))
      
      seq_col_name <- if ("NegativeSequence" %in% names(dinuc_data)) {
        "NegativeSequence"
      } else if ("Sequence" %in% names(dinuc_data)) {
        "Sequence"
      } else {
        NA_character_
      }
      
      message("Dinuc seq_col_name: ", seq_col_name)
      message("seq_col_name is null: ", is.null(seq_col_name))
      message("seq_col_name is na: ", is.na(seq_col_name))
      message("seq_col_name in names: ", seq_col_name %in% names(dinuc_data))
      
      if (!is.null(seq_col_name) && !is.na(seq_col_name) && seq_col_name %in% names(dinuc_data)) {
        message("Dinuc sequences to process: ", length(dinuc_data[[seq_col_name]]))
        
        gc_values <- sapply(as.character(dinuc_data[[seq_col_name]]), function(s) {
          chars <- strsplit(s, "")[[1]]
          gc_count <- sum(chars %in% c("G", "C"))
          total_count <- length(chars)
          if (total_count > 0) {
            (gc_count / total_count) * 100
          } else {
            0
          }
        })
        
        message("Dinuc gc_values: ", paste(gc_values, collapse = ", "))
        avg_gc <- mean(gc_values)
        message("Dinuc avg_gc calculated: ", avg_gc)
        
        message("Dinuc avg_gc: ", avg_gc)
        
        length_values <- nchar(as.character(dinuc_data[[seq_col_name]]))
        message("Dinuc length_values: ", paste(length_values, collapse = ", "))
        avg_length <- round(mean(length_values), 2)
        message("Dinuc avg_length: ", avg_length)
        
        new_row <- data.frame(
          Method = "Dinucleotide",
          Count = nrow(dinuc_data),
          Avg_Length = avg_length,
          Avg_GC = round(avg_gc, 2),
          stringsAsFactors = FALSE
        )
        
        message("Dinuc new_row: ", paste(names(new_row), collapse = ", "))
        message("Dinuc new_row values: ", paste(unlist(new_row), collapse = ", "))
        
        message("method_comparison before rbind: ", nrow(method_comparison))
        message("method_comparison names: ", paste(names(method_comparison), collapse = ", "))
        
        method_comparison <- rbind(method_comparison, new_row)
        
        message("method_comparison after rbind: ", nrow(method_comparison))
      }
    }
    
    # Add Random stats
    message("=== RANDOM DEBUG ===")
    if (!is.null(random_data)) {
      message("Random data class: ", paste(class(random_data), collapse = ", "))
      
      # If random_data is a list, try to convert to data frame
      if (is.list(random_data) && !is.data.frame(random_data)) {
        message("Random data is a list, attempting to convert to data frame...")
        tryCatch({
          random_data <- as.data.frame(random_data, stringsAsFactors = FALSE)
          message("Conversion successful!")
        }, error = function(e) {
          message("Conversion failed: ", e$message)
        })
      }
    }
    
    if (!is.null(random_data) && is.data.frame(random_data) && nrow(random_data) > 0) {
      seq_col_name <- if ("NegativeSequence" %in% names(random_data)) {
        "NegativeSequence"
      } else if ("Sequence" %in% names(random_data)) {
        "Sequence"
      } else {
        NA_character_
      }
      
      if (!is.null(seq_col_name) && !is.na(seq_col_name) && seq_col_name %in% names(random_data)) {
        avg_gc <- mean(sapply(as.character(random_data[[seq_col_name]]), function(s) {
          chars <- strsplit(s, "")[[1]]
          sum(chars %in% c("G", "C")) / length(chars) * 100
        }))
        
        method_comparison <- rbind(method_comparison, data.frame(
          Method = "Random",
          Count = nrow(random_data),
          Avg_Length = round(mean(nchar(as.character(random_data[[seq_col_name]]))), 2),
          Avg_GC = round(avg_gc, 2),
          stringsAsFactors = FALSE
        ))
      }
    }
    
    if (nrow(method_comparison) > 0) {
      writeData(wb, "Method Comparison", method_comparison)
      addStyle(wb, "Method Comparison", headerStyle, rows = 1, cols = 1:ncol(method_comparison), gridExpand = TRUE)
    }
    
    # ============================================================
    # SHEET 10: ANALYSIS METADATA
    # ============================================================
    addWorksheet(wb, "Analysis Metadata")
    
    metadata <- data.frame(
      Property = c("Analysis Date", "Tool Name", "Tool Version",
                   "Max Sequence Length", "Total Sequences Analyzed",
                   "RNAfold Used", "Export Format"),
      Value = c(
        format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
        "NeRNA - Negative RNA Generator",
        "1.0.0",
        "1000 nucleotides",
        total_original,
        if (!is.null(rnafold_results)) "Yes" else "No",
        "Excel (.xlsx)"
      ),
      stringsAsFactors = FALSE
    )
    
    writeData(wb, "Analysis Metadata", metadata)
    addStyle(wb, "Analysis Metadata", headerStyle, rows = 1, cols = 1:2, gridExpand = TRUE)
    
    # Save workbook
    saveWorkbook(wb, filename, overwrite = TRUE)
    
    return(TRUE)
  }, error = function(e) {
    message("Error in export_comprehensive_results: ", e$message)
    message("Stack trace: ", deparse(e$call))
    stop("Error exporting to Excel: ", e$message)
  })
}

# ============================================================
# METHOD COMPARISON SUMMARY
# ============================================================
# NOTE: The generateMethodComparisonSummary function is now defined
# inside server.R where it has access to the reactive values object

# ============================================================
# REAL-TIME STRUCTURE VIEWER
# ============================================================

# Create arc diagram data for structure visualization
create_arc_diagram_data <- function(sequence, structure) {
  tryCatch({
    # Parse structure to find base pairs
    chars <- strsplit(structure, "")[[1]]
    n <- length(chars)
    
    # Find matching pairs
    stack <- list()
    pairs <- list()
    
    for (i in 1:n) {
      if (chars[i] == "(") {
        stack[[length(stack) + 1]] <- i
      } else if (chars[i] == ")") {
        if (length(stack) > 0) {
          left <- stack[[length(stack)]]
          stack <- stack[-length(stack)]
          pairs[[length(pairs) + 1]] <- c(left, i)
        }
      }
    }
    
    # Create data frame for visualization
    if (length(pairs) > 0) {
      pairs_df <- data.frame(
        start = sapply(pairs, function(p) p[1]),
        end = sapply(pairs, function(p) p[2]),
        stringsAsFactors = FALSE
      )
      
      return(list(
        sequence = sequence,
        structure = structure,
        pairs = pairs_df,
        length = n
      ))
    } else {
      return(NULL)
    }
  }, error = function(e) {
    message("Error creating arc diagram: ", e$message)
    return(NULL)
  })
}

# Render structure as text-based visualization
render_structure_text <- function(sequence, structure, mfe = NULL) {
  tryCatch({
    # Create a simple text-based visualization
    seq_chars <- strsplit(sequence, "")[[1]]
    struct_chars <- strsplit(structure, "")[[1]]
    
    # Create position numbers
    positions <- sprintf("%3d", 1:length(seq_chars))
    
    # Format output
    output <- paste0(
      "Position: ", paste(positions, collapse = " "), "\n",
      "Sequence: ", paste(sprintf("%3s", seq_chars), collapse = " "), "\n",
      "Structure:", paste(sprintf("%3s", struct_chars), collapse = " ")
    )
    
    if (!is.null(mfe)) {
      output <- paste0(output, "\n", "MFE: ", mfe, " kcal/mol")
    }
    
    # Add statistics
    paired <- sum(struct_chars %in% c("(", ")"))
    unpaired <- sum(struct_chars == ".")
    
    output <- paste0(output, "\n\n",
                    "Statistics:\n",
                    "  Paired bases: ", paired, " (", round(paired/length(struct_chars)*100, 1), "%)\n",
                    "  Unpaired bases: ", unpaired, " (", round(unpaired/length(struct_chars)*100, 1), "%)")
    
    return(output)
  }, error = function(e) {
    return("Error rendering structure")
  })
}

# Create 3D arc diagram for RNA secondary structure
create_3d_structure_plot <- function(sequence, structure, mfe = NULL) {
  tryCatch({
    library(plotly)
    
    seq_chars <- strsplit(sequence, "")[[1]]
    struct_chars <- strsplit(structure, "")[[1]]
    n <- length(seq_chars)
    
    # Find base pairs
    stack <- list()
    pairs <- list()
    
    for (i in 1:n) {
      if (struct_chars[i] == "(") {
        stack[[length(stack) + 1]] <- i
      } else if (struct_chars[i] == ")") {
        if (length(stack) > 0) {
          left <- stack[[length(stack)]]
          stack <- stack[-length(stack)]
          pairs[[length(pairs) + 1]] <- c(left, i)
        }
      }
    }
    
    # Create base positions along x-axis
    x_pos <- 1:n
    y_base <- rep(0, n)
    
    # Color bases by type
    base_colors <- sapply(seq_chars, function(b) {
      if (b == "A") return("red")
      if (b == "U") return("blue")
      if (b == "G") return("green")
      if (b == "C") return("orange")
      return("gray")
    })
    
    # Create plot
    p <- plot_ly()
    
    # Add base positions (sequence)
    p <- add_trace(p,
                   x = x_pos,
                   y = y_base,
                   type = 'scatter',
                   mode = 'markers+text',
                   marker = list(size = 12, color = base_colors, 
                                line = list(color = 'black', width = 1)),
                   text = seq_chars,
                   textposition = 'bottom center',
                   textfont = list(size = 10, color = 'black'),
                   name = 'Bases',
                   hoverinfo = 'text',
                   hovertext = paste0("Position: ", 1:n, "<br>Base: ", seq_chars, 
                                     "<br>Structure: ", struct_chars))
    
    # Add arcs for base pairs
    if (length(pairs) > 0) {
      for (pair in pairs) {
        i <- pair[1]
        j <- pair[2]
        
        # Create arc
        arc_x <- seq(i, j, length.out = 50)
        arc_height <- (j - i) / 4
        arc_y <- arc_height * sin(pi * (arc_x - i) / (j - i))
        
        p <- add_trace(p,
                      x = arc_x,
                      y = arc_y,
                      type = 'scatter',
                      mode = 'lines',
                      line = list(color = 'rgba(100, 100, 200, 0.5)', width = 2),
                      showlegend = FALSE,
                      hoverinfo = 'text',
                      hovertext = paste0("Base pair: ", i, "-", j, 
                                        " (", seq_chars[i], "-", seq_chars[j], ")"))
      }
    }
    
    # Layout
    title_text <- paste0("RNA Secondary Structure")
    if (!is.null(mfe)) {
      title_text <- paste0(title_text, " (MFE: ", mfe, " kcal/mol)")
    }
    
    p <- layout(p,
                title = title_text,
                xaxis = list(title = "Position", showgrid = FALSE),
                yaxis = list(title = "Structure", showgrid = FALSE, zeroline = TRUE),
                hovermode = 'closest',
                showlegend = TRUE,
                legend = list(x = 0.1, y = 0.9))
    
    return(p)
  }, error = function(e) {
    message("Error creating 3D structure plot: ", e$message)
    return(NULL)
  })
}

# Quick RNAfold for single sequence (for comparison viewer)
quick_rnafold <- function(sequence) {
  tryCatch({
    # Run RNAfold directly with sequence
    rnafold_path <- get_rnafold_executable()
    
    # Create input for RNAfold (FASTA format)
    input_text <- paste0(">temp\n", sequence)
    
    output <- system2(rnafold_path, 
                     input = input_text,
                     stdout = TRUE, 
                     stderr = TRUE)
    
    # Debug: print raw output
    message("RNAfold raw output length: ", length(output))
    if (length(output) > 0) {
      for (i in 1:min(5, length(output))) {
        message("  Line ", i, ": ", output[i])
      }
    }
    
    # Parse output - RNAfold returns:
    # >temp
    # sequence
    # structure (MFE)
    if (length(output) >= 3) {
      structure_line <- output[3]
      message("Structure line: '", structure_line, "'")
      
      # Extract structure and MFE
      # The line should be like: "(((...))) (-12.34)"
      # Split by space to separate structure and MFE
      parts <- strsplit(trimws(structure_line), "\\s+")[[1]]
      message("Parts: ", paste(parts, collapse = " | "))
      
      if (length(parts) >= 2) {
        structure <- parts[1]
        
        # Extract MFE from the last part which should be like "(-12.34)"
        mfe_part <- parts[length(parts)]
        message("MFE part: '", mfe_part, "'")
        
        mfe_match <- regmatches(mfe_part, regexpr("-?[0-9.]+", mfe_part))
        message("MFE match: ", if(length(mfe_match) > 0) mfe_match[1] else "NONE")
        
        if (length(mfe_match) > 0) {
          mfe <- as.numeric(mfe_match[1])
        } else {
          mfe <- NA
        }
        
        return(list(structure = structure, mfe = mfe))
      }
    }
    
    # Fallback
    return(list(structure = paste(rep(".", nchar(sequence)), collapse = ""), mfe = NA))
  }, error = function(e) {
    message("Error running quick RNAfold: ", e$message)
    return(list(structure = paste(rep(".", nchar(sequence)), collapse = ""), mfe = NA))
  })
}

# ============================================================
# SESSION MANAGEMENT
# ============================================================

# Save session data
save_session <- function(session_data, filename) {
  tryCatch({
    # Create session object
    session_obj <- list(
      timestamp = Sys.time(),
      version = "1.0.0",
      fasta_data = session_data$fastaData,
      rnafold_results = session_data$rnafoldResults,
      nerna_results = session_data$nernaResults,
      dinucleotide_results = session_data$dinucleotideResults,
      random_results = session_data$randomResults,
      rna_type = session_data$rnaType
    )
    
    # Save as RDS file
    saveRDS(session_obj, filename)
    
    return(TRUE)
  }, error = function(e) {
    message("Error saving session: ", e$message)
    return(FALSE)
  })
}

# Load session data
load_session <- function(filename) {
  tryCatch({
    # Load RDS file
    session_obj <- readRDS(filename)
    
    # Validate session object
    if (!is.list(session_obj)) {
      stop("Invalid session file format")
    }
    
    return(session_obj)
  }, error = function(e) {
    message("Error loading session: ", e$message)
    return(NULL)
  })
}

# ============================================================
# BATCH PROCESSING
# ============================================================

# Process multiple FASTA files
process_batch_fasta <- function(fasta_files, max_length = 1000) {
  tryCatch({
    all_data <- list()
    
    for (i in 1:length(fasta_files)) {
      file_path <- fasta_files[i]
      
      # Read FASTA file
      data <- readFastaToTable(file_path)
      
      if (!is.null(data) && nrow(data) > 0) {
        # Add source file information
        data$SourceFile <- basename(file_path)
        data$FileIndex <- i
        
        all_data[[i]] <- data
      }
    }
    
    # Combine all data
    if (length(all_data) > 0) {
      combined_data <- do.call(rbind, all_data)
      return(combined_data)
    } else {
      return(NULL)
    }
  }, error = function(e) {
    message("Error in batch processing: ", e$message)
    return(NULL)
  })
}

# Batch negative generation
batch_generate_negatives <- function(sequences, method = "nerna", shifting_size = 23, seed = NULL, rnafold_results = NULL) {
  tryCatch({
    if (method == "nerna") {
      return(generateNeRNA(sequences, shifting_size, seed, rnafold_results))
    } else if (method == "dinucleotide") {
      return(dinucleotideShuffle(sequences, seed))
    } else if (method == "random") {
      return(randomShuffle(sequences, seed))
    } else {
      stop("Unknown method")
    }
  }, error = function(e) {
    message("Error in batch generation: ", e$message)
    return(NULL)
  })
} 
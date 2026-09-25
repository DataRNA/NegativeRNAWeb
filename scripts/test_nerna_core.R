expressions <- parse(file = "global.R")
required_functions <- c(
  "rna_octal_rotate_once",
  "generateNeRNA_sequential",
  "generateNeRNA_parallel"
)

for (expression in expressions) {
  if (
    is.call(expression) &&
      identical(as.character(expression[[1]]), "<-") &&
      as.character(expression[[2]]) %in% required_functions
  ) {
    eval(expression, envir = .GlobalEnv)
  }
}

# Lightweight stand-ins allow the parallel worker path to be tested without
# starting background processes or loading the complete Shiny application.
setup_parallel_processing <- function() 1L
future_sapply <- function(X, FUN, future.seed = TRUE) sapply(X, FUN)

sequences <- c("AUGC", "GGAAUUCC")
rnafold_results <- data.frame(
  Sequence = sequences,
  SecondaryStructure = c("(..)", "((....))"),
  stringsAsFactors = FALSE
)

sequential <- generateNeRNA_sequential(
  sequences,
  shifting_size = 1L,
  rnafold_results = rnafold_results
)
parallel <- generateNeRNA_parallel(
  sequences,
  shifting_size = 1L,
  rnafold_results = rnafold_results
)

stopifnot(
  identical(unname(sequential), unname(parallel)),
  identical(nchar(sequential), nchar(sequences)),
  all(grepl("^[ACGU]+$", sequential)),
  inherits(try(rna_octal_rotate_once("AUGC", "((..", 1L), silent = TRUE), "try-error")
)

cat("NeRNA sequential and parallel paths agree.\n")

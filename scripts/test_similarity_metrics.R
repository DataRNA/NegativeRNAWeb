expressions <- parse(file = "global.R")
required_functions <- c(
  "calculate_structural_similarity",
  "calculate_sequence_identity",
  "calculate_dinucleotide_similarity",
  "get_dinucleotide_frequencies",
  "calculate_gc_similarity"
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

frequencies <- get_dinucleotide_frequencies("AACGUU")
stopifnot(
  length(frequencies) == 16L,
  !is.null(names(frequencies)),
  length(unique(names(frequencies))) == 16L
)

metrics <- calculate_structural_similarity("AACGUU", "AACGUU")
stopifnot(
  identical(
    sort(names(metrics)),
    sort(c("sequence_identity", "dinucleotide_similarity", "gc_similarity"))
  ),
  abs(metrics$dinucleotide_similarity - 100) < 1e-9
)

cat("Similarity functions OK; fields:", paste(names(metrics), collapse = ", "), "\n")

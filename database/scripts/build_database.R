source_path <- "database/source/variants.csv"
output_path <- "Scores_FECH_UROD.csv"

thresholds <- c(
  PP2 = 0.908,
  SIFT = 0.05,
  Consurf = 7,
  MA = 3.5,
  PROVEAN = -2.5,
  FATHMM = -2.0,
  REVEL = 0.7,
  CADD = 25,
  AlphaMissense_score = 0.5,
  BayesDel = 0.07
)

# TRUE means a value at or above the threshold supports a damaging prediction.
# FALSE means a value at or below the threshold supports one.
damaging_when_high <- c(
  PP2 = TRUE,
  SIFT = FALSE,
  Consurf = TRUE,
  MA = TRUE,
  PROVEAN = FALSE,
  FATHMM = FALSE,
  REVEL = TRUE,
  CADD = TRUE,
  AlphaMissense_score = TRUE,
  BayesDel = TRUE
)

if (!file.exists(source_path)) {
  stop("Database source not found: ", source_path)
}

variants <- read.csv(
  source_path,
  fileEncoding = "UTF-8",
  na.strings = c("NA", ""),
  stringsAsFactors = FALSE,
  check.names = FALSE
)

missing_columns <- setdiff(names(thresholds), names(variants))
if (length(missing_columns) > 0) {
  stop("Missing predictor columns: ", paste(missing_columns, collapse = ", "))
}

for (predictor in names(thresholds)) {
  original <- variants[[predictor]]
  numeric_values <- suppressWarnings(as.numeric(original))
  invalid <- is.na(numeric_values) & !is.na(original) &
    nzchar(trimws(as.character(original)))

  if (any(invalid)) {
    stop(
      "Non-numeric values found in ", predictor, " at source rows: ",
      paste(which(invalid) + 1L, collapse = ", ")
    )
  }

  variants[[predictor]] <- numeric_values
}

predictor_pass <- vapply(
  names(thresholds),
  function(predictor) {
    values <- variants[[predictor]]
    if (damaging_when_high[[predictor]]) {
      values >= thresholds[[predictor]]
    } else {
      values <= thresholds[[predictor]]
    }
  },
  logical(nrow(variants))
)

valid_predictors <- rowSums(!is.na(predictor_pass))
damaging_predictors <- rowSums(predictor_pass, na.rm = TRUE)

variants[["Consensus_Score"]] <- ifelse(
  valid_predictors == 0,
  NA_integer_,
  as.integer(round(damaging_predictors / valid_predictors * 100))
)

if (any(
  variants[["Consensus_Score"]] < 0 |
    variants[["Consensus_Score"]] > 100,
  na.rm = TRUE
)) {
  stop("Generated consensus scores fall outside 0–100.")
}

source_lines <- readLines(source_path, encoding = "UTF-8", warn = FALSE)
if (length(source_lines) != nrow(variants) + 1L) {
  stop("Source CSV contains embedded line breaks or an unexpected row count.")
}

score_text <- ifelse(
  is.na(variants[["Consensus_Score"]]),
  "NA",
  as.character(variants[["Consensus_Score"]])
)

output_lines <- c(
  paste0(source_lines[[1]], ",Consensus_Score"),
  paste0(source_lines[-1], ",", score_text)
)
writeLines(output_lines, output_path, useBytes = TRUE)

message("Generated ", output_path, " from ", source_path)
message("Rows: ", nrow(variants))
message("Rows scored: ", sum(!is.na(variants[["Consensus_Score"]])))
message("Rows without predictor data: ", sum(valid_predictors == 0))

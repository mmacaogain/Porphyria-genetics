library(shiny)
library(DT)

porphyria_data <- read.csv(
  "Scores_FECH_UROD.csv",
  fileEncoding = "UTF-8",
  na.strings = c("NA", ""),
  stringsAsFactors = FALSE,
  check.names = FALSE
)

display_columns <- c(
  "Gene",
  "GRCh37",
  "HGVS_nucleotide",
  "HGVS_protein",
  "ClinVarID",
  "Variant_source",
  "PP2",
  "Consurf",
  "PROVEAN",
  "FATHMM",
  "CADD",
  "AlphaMissense_score",
  "BayesDel",
  "REVEL",
  "SIFT",
  "MA",
  "Consensus_Score"
)

porphyria_display <- porphyria_data[display_columns]
names(porphyria_display)[
  names(porphyria_display) == "AlphaMissense_score"
] <- "AM"
names(porphyria_display)[
  names(porphyria_display) == "Consensus_Score"
] <- "PAS"

extended_predictors <- c("PP2", "Consurf", "PROVEAN", "FATHMM", "CADD")
extended_column_indices <- match(
  extended_predictors,
  names(porphyria_display)
) - 1L

predictor_thresholds <- c(
  PP2 = 0.908,
  Consurf = 7,
  PROVEAN = -2.5,
  FATHMM = -2.0,
  CADD = 25,
  AM = 0.5,
  BayesDel = 0.07,
  REVEL = 0.7,
  SIFT = 0.05,
  MA = 3.5
)

damaging_when_high <- c(
  PP2 = TRUE,
  Consurf = TRUE,
  PROVEAN = FALSE,
  FATHMM = FALSE,
  CADD = TRUE,
  AM = TRUE,
  BayesDel = TRUE,
  REVEL = TRUE,
  SIFT = FALSE,
  MA = TRUE
)

threshold_background <- function(threshold, high_is_damaging) {
  crossed_colour <- "#f6c1c1"
  neutral_colour <- "#ffffff"

  if (high_is_damaging) {
    # styleInterval() assigns equality to the lower interval. Move the cut
    # imperceptibly lower so values equal to the threshold are highlighted.
    cutoff <- threshold -
      (.Machine$double.eps * max(1, abs(threshold)) * 8)
    DT::styleInterval(cutoff, c(neutral_colour, crossed_colour))
  } else {
    DT::styleInterval(threshold, c(crossed_colour, neutral_colour))
  }
}

server <- function(input, output, session) {
  output$mytable <- DT::renderDT({
    table <- DT::datatable(
      porphyria_display,
      rownames = FALSE,
      options = list(
        pageLength = 25,
        scrollX = TRUE,
        columnDefs = list(
          list(
            targets = extended_column_indices,
            visible = FALSE
          )
        )
      )
    )

    for (predictor in names(predictor_thresholds)) {
      table <- DT::formatStyle(
        table,
        predictor,
        backgroundColor = threshold_background(
          predictor_thresholds[[predictor]],
          damaging_when_high[[predictor]]
        )
      )
    }

    DT::formatStyle(
      table,
      "PAS",
      backgroundColor = DT::styleInterval(
        seq(10, 90, by = 10),
        c(
          "#fff5f0", "#fee0d2", "#fcbba1", "#fc9272", "#fb6a4a",
          "#ef3b2c", "#cb181d", "#a50f15", "#800026", "#67000d"
        )
      ),
      color = DT::styleInterval(60, c("#3b0a0a", "#ffffff")),
      fontWeight = DT::styleInterval(70, c("normal", "bold"))
    )
  }, server = TRUE)

  extended_visible <- reactiveVal(FALSE)
  table_proxy <- DT::dataTableProxy("mytable")

  observeEvent(input$toggle_extended_predictors, {
    show_extended <- !extended_visible()
    extended_visible(show_extended)

    if (show_extended) {
      DT::showCols(table_proxy, extended_column_indices)
    } else {
      DT::hideCols(table_proxy, extended_column_indices)
    }

    updateActionButton(
      session,
      "toggle_extended_predictors",
      label = if (show_extended) {
        "Hide extended predictors"
      } else {
        "Extended predictor list"
      }
    )
  })
}

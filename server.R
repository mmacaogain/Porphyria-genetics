library(shiny)
library(DT)

porphyria_data <- read.csv(
  "Scores_FECH_UROD.csv",
  fileEncoding = "UTF-8",
  na.strings = c("NA", ""),
  stringsAsFactors = FALSE,
  check.names = FALSE
)

server <- function(input, output, session) {
  output$mytable <- DT::renderDT({
    DT::datatable(
      porphyria_data,
      rownames = FALSE,
      filter = "top",
      options = list(
        pageLength = 25,
        scrollX = TRUE
      )
    ) |>
      DT::formatStyle(
        "Consensus_Score",
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
}

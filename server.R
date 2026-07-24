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
          c(25, 50, 75),
          c("#f1f5f9", "#dbeafe", "#93c5fd", "#2563eb")
        ),
        color = DT::styleInterval(75, c("#1f2937", "#ffffff")),
        fontWeight = DT::styleInterval(75, c("normal", "bold"))
      )
  }, server = TRUE)
}

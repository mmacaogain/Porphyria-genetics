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
    )
  }, server = TRUE)
}

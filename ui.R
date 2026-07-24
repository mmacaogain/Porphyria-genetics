library(shiny)
library(DT)

ui <- fluidPage(
  titlePanel("PorphyriaDB"),
  tags$p(
    "Browse and search missense variants associated with porphyria. ",
    "Use the filters beneath each column heading to narrow the results."
  ),
  tags$div(
    role = "region",
    `aria-label` = "Porphyria missense variant database",
    DT::DTOutput("mytable")
  )
)

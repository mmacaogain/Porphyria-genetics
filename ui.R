library(shiny)
library(DT)

ui <- fluidPage(
  titlePanel("PorphyriaDB"),
  tags$p("Browse and search missense variants associated with porphyria."),
  tags$div(
    style = paste(
      "display: flex;",
      "align-items: center;",
      "gap: 0.75rem;",
      "margin-bottom: 1rem;"
    ),
    actionButton(
      "toggle_extended_predictors",
      "Extended predictor list",
      `aria-controls` = "mytable"
    ),
    tags$span("PAS incorporates all available predictors.")
  ),
  tags$div(
    role = "region",
    `aria-label` = "Porphyria missense variant database",
    DT::DTOutput("mytable")
  )
)

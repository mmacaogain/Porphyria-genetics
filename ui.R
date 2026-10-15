library(shiny)
library(DT)

ui <- fluidPage(
  titlePanel("PorphyriaDB"),
  tags$p("Browse and search missense variants associated with porphyria."),
  tags$div(
    style = paste(
      "display: flex;",
      "align-items: center;",
      "justify-content: space-between;",
      "gap: 1rem;",
      "flex-wrap: wrap;",
      "margin-bottom: 1rem;"
    ),
    tags$div(
      style = paste(
        "display: flex;",
        "align-items: center;",
        "gap: 0.75rem;",
        "flex-wrap: wrap;"
      ),
      actionButton(
        "toggle_extended_predictors",
        "Extended predictor list",
        `aria-controls` = "mytable"
      ),
      tags$span("PAS incorporates all available predictors.")
    ),
    downloadButton(
      "download_dataset",
      "Download dataset (CSV)",
      icon = icon("download"),
      title = "Download the complete PorphyriaDB variant predictor dataset"
    )
  ),
  tags$div(
    role = "region",
    `aria-label` = "Porphyria missense variant database",
    DT::DTOutput("mytable")
  )
)

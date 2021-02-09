library(DT)
library(shiny)

PorData<-read.csv("Scores.csv", row.names = 1)
server <- function(input, output) {
  output$mytable = DT::renderDataTable({
    PorData
  })
}

#shinyApp(ui, server)


#needed to update various packages to deploy to shiny.io
#install.packages('xtable')
#install.packages('DT')
#install.packages("base64enc")
#install.packages('shiny')
#install.packages('htmltools')

library(DT)
library(shiny)
require(devtools) #added devtools to try and ensure DT:: command was followed


PorData<-read.csv("Scores.csv", row.names = 1)
server <- function(input, output) {
  output$mytable = DT::renderDataTable({
    PorData
  })
}

#shinyApp(ui, server)

#search box is still not working...
#DataTables warning: table id=DataTables_Table_0 - Error in tolower(x): invalid multibyte string 100 
#
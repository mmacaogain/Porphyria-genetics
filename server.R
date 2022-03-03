#needed to update various packages to deploy to shiny.io
#install.packages('xtable')
#install.packages('DT')
#install.packages("base64enc")
#install.packages('shiny')
#install.packages('htmltools')

library(DT)
library(tidyverse)
library(shiny)
require(devtools) #added devtools to try and ensure DT:: command was followed

#PorData<-readr::read_csv("Scores.csv", encoding="UTF-8") 
PorData<-read.csv("Scores.csv",encoding="UTF-8")
server <- function(input, output) {
  output$mytable = DT::renderDataTable({
    PorData
  })
}

#shinyApp(ui, server)

#search box is still not working...
#DataTables warning: table id=DataTables_Table_0 - Error in tolower(x): invalid multibyte string 100 
#Solved: https://github.com/rstudio/DT/issues/497
#Solved2: There was actually just some wierd apple symbol in one of the Chen et al reference cells. Original warning "DataTables warning: table id=DataTables_Table_0 - Error in tolower(x): invalid input 'Chen et al. 2019 PMID:<f0>30740734' in 'utf8towcs'"
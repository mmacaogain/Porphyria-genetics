library(DT)
#library(shiny)

PorData<-read.csv("scores.csv", row.names = 1)
ui<-basicPage(
  h2("PorphyriaDB - Biochemical Genetics Laboratory, St. James's Hospital "),
  DT::dataTableOutput("mytable")
)

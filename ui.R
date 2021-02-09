library(DT)
library(shiny)

PorData<-read.csv("Scores.csv", row.names = 1)
basicPage(
  h2("PorphyriaDB - Biochemical Genetics Laboratory, St. James's Hospital "),
  DT::dataTableOutput("mytable")
)


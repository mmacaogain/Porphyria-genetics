library(DT)
#library(shiny)
#PorData<-read.csv("Scores.csv",encoding="UTF-8")
#PorData<-readr::read_csv("scores.csv")
PorData<-read.csv("Scores.csv", row.names = 1)
ui<-basicPage(
  h2("PorphyriaDB - Biochemical Genetics Laboratory, St. James's Hospital "),
  DT::dataTableOutput("mytable")
)

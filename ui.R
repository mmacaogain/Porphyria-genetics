library(DT)
#library(shiny)
#PorData<-read.csv("Scores.csv",encoding="UTF-8")
#PorData<-readr::read_csv("scores.csv")
PorData<-read.csv("Scores.csv", row.names = 1)
ui<-basicPage(
  h2("PorphyriaDB: a genetic database of missense variants associated with acute intermittent porphyrias"),
  DT::dataTableOutput("mytable")
)

# PorphyriaDB

PorphyriaDB is a searchable database of missense variants associated with
porphyria. The Shiny application displays the curated variants in
`Scores_FECH_UROD.csv`.

The database is curated by members of the Biochemical Genetics Laboratory at
St. James's Hospital, Dublin, for the National Centre for the Investigation and
Diagnosis of Porphyria.

## Run locally

Install the two runtime dependencies:

```r
install.packages(c("shiny", "DT"))
```

Then run the application from the repository root:

```r
shiny::runApp()
```

The deployed application is available at
<https://porphyria.shinyapps.io/PorphyriaDB_2/>.

# PorphyriaDB database workflow

This folder separates database development from the Shiny application.

## Structure

- `source/variants.csv` is the canonical curated input. Predictor values and
  variant annotations are edited here.
- `scripts/build_database.R` validates the source, calculates predictor
  agreement, and generates the app database at `Scores_FECH_UROD.csv`.
- `analysis/` is reserved for reproducible exploratory or publication
  analyses. Analysis outputs are not required to run the Shiny app.

The root `Scores_FECH_UROD.csv` is generated. Do not edit it directly.

## Making a database change

1. Create a database-specific branch.
2. Edit `database/source/variants.csv`.
3. From the repository root, run:

   ```powershell
   Rscript database/scripts/build_database.R
   ```

4. Review and commit both the source and generated CSV changes.
5. Run the Shiny app and confirm the table before deployment.

App presentation changes should be made on separate branches and should
normally touch `server.R`, `ui.R`, or future content files—not the database
source.

## Predictor agreement

`Consensus_Score` is the percentage of available predictors whose values cross
their configured damaging threshold. Predictors with missing values are
excluded from the denominator. If a variant has no predictor values, its
consensus score remains missing.

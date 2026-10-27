# data-raw

Raw inputs and the scripts that turn them into the package's data. Not part of the built package.

| Path | What |
|---|---|
| `Maine.Ref/maine_referendum_2008_2025.csv` | Yes / No votes on Maine's 72 statewide questions, Nov 2008–2025, by reporting area (462). Source of `maine`. |
| `Maine.Ref/maine_ballot_measures_november_2008_2025.csv` | Year, number, subject, summary and statewide result of each question. Source of `maine_questions`. |
| `Maine.Ref/referendum_row_key.csv` | Member municipalities and county of each reporting area. |
| `Maine.Ref/maine_2008_2019_referendums.csv`, `maine_ballot_questions.txt` | The earlier 2008–2019 file (423 towns) used by version 0.1.0. |
| `Mass.Ref/`, `Multnomah.Ref/`, `Portland.Maine.Ref/` | Massachusetts town-level referendums (2000–2020), Multnomah County, Oregon (1995–2021), and Portland, Maine, from version 0.1.0. Read them with `bbmm_data_from_wide()`. |
| `maine.R` | Builds `data/maine.rda`, `data/maine_questions.rda`, `inst/extdata/maine_towns.gpkg`. |
| `tutorial_results.R` | Runs the analyses shown in the tutorial and saves them to `inst/extdata/maine_tutorial.rds` (about 1.5 hours on 10 threads). |

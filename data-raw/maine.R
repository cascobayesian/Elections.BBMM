# data-raw/maine.R -- builds data/maine.rda, data/maine_questions.rda and inst/extdata/maine_towns.gpkg.
#   Rscript data-raw/maine.R            (from the package root; needs pkgload and sf)
# Votes: Maine Secretary of State tabulations by municipality, November 2008-2025, aggregated to the
# 462 reporting areas that report statewide referendums (Maine.Ref/maine_referendum_2008_2025.csv;
# Maine.Ref/referendum_row_key.csv lists each area's member municipalities and county).
# Questions: Maine.Ref/maine_ballot_measures_november_2008_2025.csv.
# Geometry: Maine Office of GIS town and township polygons (LAND parts only) dissolved to reporting areas,
# simplified here to ~200 m. Set MAINE_UNITS_GPKG to the unsimplified file to rebuild it.
pkgload::load_all(quiet = TRUE)
raw <- read.csv("data-raw/Maine.Ref/maine_referendum_2008_2025.csv", check.names = FALSE)
key <- read.csv("data-raw/Maine.Ref/referendum_row_key.csv")
qs <- read.csv("data-raw/Maine.Ref/maine_ballot_measures_november_2008_2025.csv", check.names = FALSE)

maine <- bbmm_data_from_wide(raw, id = "TOWN", yes_suffix = "_YES", no_suffix = "_NO")
stopifnot(nrow(maine$X) == 462, ncol(maine$X) == 72, identical(key$TOWN, maine$units))
maine_questions <- data.frame(question_id = maine$questions, year = qs$year, number = qs$`question number`,
                              label = sprintf("%d Q%d: %s", qs$year, qs$`question number`, qs$subject),
                              summary = qs$`ten word summary`, statewide_yes = qs$`overall resulting percentage`)
stopifnot(identical(maine_questions$question_id, sprintf("REFERENDUM_%d_%d", qs$year, qs$`question number`)))
maine$labels <- maine_questions$label
maine$unit_info <- data.frame(unit_id = key$TOWN, name = key$reporting_area, county = key$counties)
usethis_save <- function(obj, name) {
  e <- new.env(); assign(name, obj, e)
  save(list = name, envir = e, file = file.path("data", paste0(name, ".rda")), compress = "xz")
}
usethis_save(maine, "maine"); usethis_save(maine_questions, "maine_questions")

src <- Sys.getenv("MAINE_UNITS_GPKG", "../../output/table2/maine/units.gpkg")
if (file.exists(src)) {
  library(sf)
  u <- st_read(src, layer = "units", quiet = TRUE)
  stopifnot(setequal(u$unit_id, maine$units))
  u <- u[match(maine$units, u$unit_id), ]
  g <- st_simplify(st_make_valid(u), preserveTopology = TRUE, dTolerance = 200)
  g <- st_cast(st_make_valid(g), "MULTIPOLYGON")
  g <- g[, "unit_id"]
  out <- "inst/extdata/maine_towns.gpkg"
  unlink(out); st_write(g, out, layer = "maine_towns", quiet = TRUE)
  cat(sprintf("geometry: %d units, %.0f KB\n", nrow(g), file.size(out) / 1024))
}
cat(sprintf("maine: %d x %d; %s\n", nrow(maine$X), ncol(maine$X), paste(file.size(list.files("data", full.names = TRUE)), collapse = " ")))

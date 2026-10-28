# bbmm_data: the count panel every fitting function takes.

#' Build a BBMM data set from Yes and total counts
#'
#' Every fitting function takes a `bbmm_data` object: two \eqn{N \times Q}
#' integer matrices of Yes votes (`X`) and two-option totals, Yes + No (`M`),
#' for \eqn{N} units and \eqn{Q} questions. A cell with no votes -- a question
#' not on the ballot in that unit, a unit that did not exist yet -- is `NA` and
#' contributes nothing to the likelihood.
#'
#' @param yes Matrix (or data frame) of Yes counts, units in rows.
#' @param total Matrix of the same shape: Yes + No counts.
#' @param units Unit identifiers; default the row names of `yes`, else `1:N`.
#' @param questions Question identifiers; default the column names of `yes`.
#' @param labels Optional short labels for the questions, same order (used in
#'   tables and figures).
#' @param unit_info Optional data frame with one row per unit (same order) of
#'   extra columns kept for reference, e.g. a display name or county.
#' @return An object of class `bbmm_data`: a list with `X`, `M` (integer
#'   matrices with `NA` for empty cells), `units`, `questions`, `labels` and
#'   `unit_info`.
#' @seealso [bbmm_data_from_wide()] for the common one-row-per-unit layout.
#' @examples
#' yes   <- matrix(c(30, 10, 25, 5, 40, 12), 3, 2)
#' total <- matrix(c(50, 40, 45, 20, 60, 30), 3, 2)
#' d <- bbmm_data(yes, total, units = c("A", "B", "C"), questions = c("Q1", "Q2"))
#' d
#' @export
bbmm_data <- function(yes, total, units = NULL, questions = NULL, labels = NULL, unit_info = NULL) {
  X <- as.matrix(yes); M <- as.matrix(total)
  if (!identical(dim(X), dim(M))) stop("`yes` and `total` must have the same dimensions.", call. = FALSE)
  if (nrow(X) < 2 || ncol(X) < 1) stop("Need at least two units and one question.", call. = FALSE)
  units <- as.character(units %||% rownames(X) %||% seq_len(nrow(X)))
  questions <- as.character(questions %||% colnames(X) %||% paste0("Q", seq_len(ncol(X))))
  if (length(units) != nrow(X) || anyDuplicated(units)) stop("`units` must be unique, one per row.", call. = FALSE)
  if (length(questions) != ncol(X) || anyDuplicated(questions))
    stop("`questions` must be unique, one per column.", call. = FALSE)
  labels <- as.character(labels %||% questions)
  if (length(labels) != ncol(X)) stop("`labels` must have one entry per question.", call. = FALSE)
  miss <- is.na(X) | is.na(M) | M == 0
  ok <- !miss
  if (any(X[ok] < 0 | X[ok] > M[ok] | X[ok] != round(X[ok])))
    stop("Counts must be whole numbers with 0 <= yes <= total.", call. = FALSE)
  if (any(M[ok] > .Machine$integer.max)) stop("Counts exceed the integer range.", call. = FALSE)
  X[miss] <- NA; M[miss] <- NA
  X <- matrix(as.integer(X), nrow(X), dimnames = list(units, questions))
  M <- matrix(as.integer(M), nrow(M), dimnames = list(units, questions))
  empty_u <- rowSums(!miss) == 0; empty_q <- colSums(!miss) == 0
  if (any(empty_u)) warning(sum(empty_u), " unit(s) have no votes on any question.", call. = FALSE)
  if (any(empty_q)) warning(sum(empty_q), " question(s) have no votes in any unit.", call. = FALSE)
  if (!is.null(unit_info)) {
    unit_info <- as.data.frame(unit_info)
    if (nrow(unit_info) != nrow(X)) stop("`unit_info` must have one row per unit.", call. = FALSE)
  }
  structure(list(X = X, M = M, units = units, questions = questions, labels = labels, unit_info = unit_info),
            class = "bbmm_data")
}

#' Build a BBMM data set from a wide table of Yes and No columns
#'
#' Reads the layout most election offices publish: one row per unit, and for
#' each question a Yes column and a No column whose names share a stem, e.g.
#' `Q_1_Y` and `Q_1_N`, or `REFERENDUM_2008_1_YES` and `REFERENDUM_2008_1_NO`.
#'
#' @param df A data frame.
#' @param id Name of the column holding unit identifiers.
#' @param yes_suffix,no_suffix Suffixes that mark the Yes and No columns.
#' @param labels Optional question labels, in the order the Yes columns appear.
#' @return A [bbmm_data] object; question identifiers are the column stems.
#' @examples
#' df <- data.frame(town = c("A", "B", "C"), Q1_Y = c(30, 10, 25), Q1_N = c(20, 30, 20),
#'                  Q2_Y = c(5, 40, 12), Q2_N = c(15, 20, 18))
#' bbmm_data_from_wide(df, id = "town")
#' @export
bbmm_data_from_wide <- function(df, id, yes_suffix = "_Y", no_suffix = "_N", labels = NULL) {
  df <- as.data.frame(df)
  if (!id %in% names(df)) stop("Column `", id, "` not found.", call. = FALSE)
  ends <- function(s, suf) substr(s, nchar(s) - nchar(suf) + 1L, nchar(s)) == suf
  ycols <- names(df)[ends(names(df), yes_suffix)]
  stems <- substr(ycols, 1L, nchar(ycols) - nchar(yes_suffix))
  ncols <- paste0(stems, no_suffix)
  if (!length(ycols)) stop("No columns end in `", yes_suffix, "`.", call. = FALSE)
  if (!all(ncols %in% names(df)))
    stop("No column found for: ", paste(ncols[!ncols %in% names(df)], collapse = ", "), call. = FALSE)
  yes <- as.matrix(df[ycols]); no <- as.matrix(df[ncols])
  info <- df[setdiff(names(df), c(ycols, ncols))]
  bbmm_data(yes, yes + no, units = df[[id]], questions = stems, labels = labels,
            unit_info = if (ncol(info) > 1) info else NULL)
}

#' @export
print.bbmm_data <- function(x, ...) {
  obs <- !is.na(x$M)
  cat(sprintf("BBMM data: %d units x %d questions, %s votes; %.1f%% of cells empty\n",
              nrow(x$X), ncol(x$X), format(sum(as.numeric(x$M), na.rm = TRUE), big.mark = ","),
              100 * mean(!obs)))
  cat(sprintf("Median votes per cell %s; pooled Yes share %.3f\n",
              format(stats::median(x$M[obs]), big.mark = ","), sum(as.numeric(x$X), na.rm = TRUE) /
                sum(as.numeric(x$M), na.rm = TRUE)))
  invisible(x)
}

# Internal: check the argument is a bbmm_data object.
check_data <- function(data) {
  if (!inherits(data, "bbmm_data")) stop("`data` must be a bbmm_data object (see bbmm_data()).", call. = FALSE)
  data
}

`%||%` <- function(a, b) if (is.null(a)) b else a

#' Maine statewide referendums by town, 2008--2025
#'
#' Yes and No votes on Maine's 72 statewide ballot questions of November
#' 2008--2025 (citizen initiatives, people's vetoes, bonds and constitutional
#' amendments), in the 462 reporting areas that report them. Most reporting
#' areas are single towns; a few combine municipalities that vote together.
#'
#' @format A [bbmm_data] object: `X` and `M` are 462 x 72 integer matrices of
#'   Yes votes and Yes + No votes; `units` the reporting-area keys (as in
#'   [maine_towns()]); `labels` short question labels; `unit_info` a data frame
#'   with `unit_id`, `name` and `county`.
#' @source Maine Secretary of State, Bureau of Corporations, Elections and
#'   Commissions, tabulations of referendum results by municipality,
#'   aggregated to reporting areas.
#' @seealso [maine_questions], [maine_towns()]
#' @examples
#' maine
#' head(maine_questions)
"maine"

#' The Maine ballot questions, 2008--2025
#'
#' One row per question in [maine], in the same order.
#'
#' @format A data frame with 72 rows:
#' \describe{
#'   \item{question_id}{Column name in `maine$X`.}
#'   \item{year, number}{Election year (November) and question number.}
#'   \item{label}{Short label, e.g. "2009 Q1: Same-sex marriage repeal".}
#'   \item{summary}{One-line description.}
#'   \item{statewide_yes}{Statewide Yes percentage.}
#' }
#' @source Maine Secretary of State.
"maine_questions"

#' Town boundaries for the Maine data
#'
#' Reads the reporting-area polygons that match `maine$units`, simplified for
#' plotting (about 200 m tolerance). Requires the \pkg{sf} package.
#'
#' @return An `sf` data frame with columns `unit_id` and `geometry`
#'   (NAD83 / UTM zone 19N).
#' @source Maine GeoLibrary municipal boundaries, dissolved to referendum
#'   reporting areas.
#' @examples
#' if (requireNamespace("sf", quietly = TRUE)) {
#'   towns <- maine_towns()
#'   nrow(towns)
#' }
#' @export
maine_towns <- function() {
  if (!requireNamespace("sf", quietly = TRUE)) stop("maine_towns() needs the sf package.", call. = FALSE)
  sf::st_read(system.file("extdata", "maine_towns.gpkg", package = "Elections.BBMM"), quiet = TRUE)
}

#' Subset a BBMM data set
#'
#' @param x A [bbmm_data] object.
#' @param i,j Units (rows) and questions (columns) to keep: indices, names or
#'   logical vectors.
#' @param ... Unused.
#' @return A [bbmm_data] object.
#' @examples
#' maine[1:50, 1:10]
#' @export
`[.bbmm_data` <- function(x, i, j, ...) {
  if (missing(i)) i <- seq_len(nrow(x$X))
  if (missing(j)) j <- seq_len(ncol(x$X))
  i <- if (is.character(i)) match(i, x$units) else seq_len(nrow(x$X))[i]
  j <- if (is.character(j)) match(j, x$questions) else seq_len(ncol(x$X))[j]
  if (anyNA(i) || anyNA(j)) stop("Unknown units or questions.", call. = FALSE)
  bbmm_data(x$X[i, j, drop = FALSE], x$M[i, j, drop = FALSE], units = x$units[i], questions = x$questions[j],
            labels = x$labels[j], unit_info = if (!is.null(x$unit_info)) x$unit_info[i, , drop = FALSE])
}

#' @export
dim.bbmm_data <- function(x) dim(x$X)

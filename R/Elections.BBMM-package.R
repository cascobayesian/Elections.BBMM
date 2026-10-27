#' Elections.BBMM: voting blocs from referendum returns
#'
#' Fits the beta-binomial mixture model (BBMM) to a panel of referendum (or
#' other two-option) returns: \eqn{N} units -- precincts, towns, counties --
#' by \eqn{Q} questions. Units in the same bloc share a profile of Yes
#' probabilities across all questions; the beta-binomial adds over-dispersion
#' within a bloc.
#'
#' Three engines, all with a likelihood power \eqn{\zeta} (tempering, or the
#' coarsened posterior of Miller and Dunson 2019) that keeps very large vote
#' counts from splitting the electorate into ever more blocs:
#' \describe{
#'   \item{[bbmm_em()], [bbmm_em_select()]}{Tempered EM at a fixed
#'     dispersion: fast point estimates, multiple random starts, BIC and
#'     minimum-share rules for the number of blocs.}
#'   \item{[bbmm_mcmc()]}{Fixed-K Markov chain Monte Carlo at power \eqn{\zeta}
#'     (split-merge, collapsed allocation updates, random-walk updates of the
#'     profiles and per-bloc dispersions).}
#'   \item{[bbmm_ss()], [bbmm_evidence()]}{Stepping-stone sequential Monte Carlo:
#'     the marginal likelihood \eqn{p(x \mid K)} for each \eqn{K}, and with a
#'     prior on \eqn{K} the posterior on the number of blocs, plus a weighted
#'     particle sample from the posterior at each \eqn{K}.}
#' }
#' Posterior summaries ([bloc_summary()], [distinctive_questions()],
#' [question_fit()]) and figures ([plot_bloc_map()], [plot_bloc_profiles()],
#' [plot_js_pca()], [plot_evidence()], [plot_coclustering()], [plot_trace()])
#' work on the output of any of the three engines.
#'
#' The package ships Maine's 72 statewide referendums of 2008--2025 by town
#' ([maine]) with town boundaries ([maine_towns()]).
#'
#' @keywords internal
#' @useDynLib Elections.BBMM, .registration = TRUE
#' @importFrom Rcpp sourceCpp
#' @importFrom ggplot2 .data
#' @importFrom stats aggregate median quantile rbeta rgamma runif rnorm rexp sd var lm
#' @importFrom utils head tail
"_PACKAGE"

# Figures. All return ggplot objects, so they can be modified, combined or saved with ggplot2::ggsave().

#' Colours for blocs
#'
#' A qualitative palette (Okabe--Ito first, then Paul Tol's muted set) for up
#' to 20 blocs, avoiding a red/blue pairing that would read as partisan.
#'
#' @param n Number of blocs.
#' @return Character vector of colours.
#' @examples
#' bloc_colours(5)
#' @export
bloc_colours <- function(n) {
  pal <- c("#E69F00", "#56B4E9", "#009E73", "#CC79A7", "#0072B2", "#D55E00", "#F0E442", "#999999",
           "#332288", "#88CCEE", "#44AA99", "#117733", "#999933", "#DDCC77", "#CC6677", "#882255",
           "#AA4499", "#661100", "#6699CC", "#000000")
  if (n > length(pal)) pal <- grDevices::hcl.colors(n, "Dark 3")
  pal[seq_len(n)]
}

bloc_scale <- function(summary, aesthetic = "fill", name = "Bloc") {
  b <- summary$blocs
  labs <- sprintf("%d  (%d units, mean Yes %.0f%%)", b$bloc, b$n_units, 100 * b$mean_yes)
  ggplot2::scale_discrete_manual(aesthetic, values = stats::setNames(bloc_colours(nrow(b)), b$bloc),
                                 labels = labs, name = name, drop = FALSE)
}

#' Map the blocs
#'
#' Draws each unit coloured by its bloc, or shaded by how certainly it belongs.
#'
#' @param summary A [bloc_summary] object.
#' @param geometry An `sf` data frame of unit polygons (or points).
#' @param id Column of `geometry` holding the unit identifiers used in the data.
#' @param fill `"bloc"` or `"certainty"`.
#' @param title Plot title.
#' @return A ggplot object.
#' @examples
#' \donttest{
#' if (requireNamespace("sf", quietly = TRUE)) {
#'   fit <- bbmm_em(maine, K = 4, zeta = 0.3, n_starts = 2, seed = 1)
#'   plot_bloc_map(bloc_summary(fit, maine), maine_towns(), id = "unit_id")
#' }
#' }
#' @export
plot_bloc_map <- function(summary, geometry, id = "unit_id", fill = c("bloc", "certainty"), title = NULL) {
  if (!requireNamespace("sf", quietly = TRUE)) stop("plot_bloc_map() needs the sf package.", call. = FALSE)
  fill <- match.arg(fill)
  a <- summary$assignment
  g <- geometry[match(a$unit, as.character(geometry[[id]])), ]
  if (anyNA(match(a$unit, as.character(geometry[[id]]))))
    warning(sum(is.na(match(a$unit, as.character(geometry[[id]])))), " unit(s) have no geometry.", call. = FALSE)
  g <- g[!is.na(g[[id]]), ]
  a <- a[match(as.character(g[[id]]), a$unit), ]
  g$bloc <- factor(a$bloc, levels = summary$blocs$bloc); g$certainty <- a$certainty
  p <- ggplot2::ggplot(g) + ggplot2::theme_void(base_size = 11) +
    ggplot2::labs(title = title %||% if (fill == "bloc") "Voting blocs" else "Certainty of bloc membership",
                  subtitle = summary$method)
  if (fill == "bloc") {
    p + ggplot2::geom_sf(ggplot2::aes(fill = .data$bloc), colour = "white", linewidth = 0.05) + bloc_scale(summary)
  } else {
    p + ggplot2::geom_sf(ggplot2::aes(fill = .data$certainty), colour = "white", linewidth = 0.05) +
      ggplot2::scale_fill_viridis_c(name = "Certainty", limits = c(0, 1), option = "magma", direction = -1)
  }
}

#' Heatmap of bloc Yes shares by question
#'
#' One row per question, one column per bloc; cells show the bloc's
#' ballot-weighted Yes share (or its gap from all units) with the value printed.
#'
#' @param summary A [bloc_summary] object.
#' @param value `"yes_share"`, `"gap"` (difference from all units) or `"p_mean"` (model).
#' @param questions Optional subset of question identifiers or labels to show.
#' @param digits Decimal places of the printed percentages (`NA` to omit them).
#' @return A ggplot object.
#' @examples
#' small <- maine[1:120, 1:20]
#' s <- bloc_summary(bbmm_em(small, K = 3, zeta = 0.3, n_starts = 3, seed = 1), small)
#' plot_bloc_profiles(s)
#' @export
plot_bloc_profiles <- function(summary, value = c("yes_share", "gap", "p_mean"), questions = NULL, digits = 0) {
  value <- match.arg(value)
  q <- summary$questions
  if (!is.null(questions)) q <- q[q$question %in% questions | q$label %in% questions, ]
  q$label <- factor(q$label, levels = rev(unique(q$label)))
  q$bloc <- factor(q$bloc)
  q$v <- q[[value]]
  p <- ggplot2::ggplot(q, ggplot2::aes(.data$bloc, .data$label, fill = .data$v)) +
    ggplot2::geom_tile(colour = "white", linewidth = 0.3) +
    ggplot2::theme_minimal(base_size = 10) +
    ggplot2::theme(panel.grid = ggplot2::element_blank()) +
    ggplot2::labs(x = "Bloc", y = NULL, subtitle = summary$method,
                  title = switch(value, yes_share = "Yes share by bloc", gap = "Yes share by bloc, minus all units",
                                 p_mean = "Model Yes probability by bloc"))
  if (!is.na(digits)) p <- p + ggplot2::geom_text(ggplot2::aes(label = round(100 * .data$v, digits)), size = 2.6,
                                                  colour = ifelse(abs(q$v - if (value == "gap") 0 else 0.5) > 0.25,
                                                                  "white", "grey15"))
  if (value == "gap") {
    p + ggplot2::scale_fill_gradient2(name = "Gap", low = "#5E3C99", mid = "grey97", high = "#E66101",
                                      labels = function(x) paste0(round(100 * x), " pts"))
  } else {
    p + ggplot2::scale_fill_gradient2(name = "Yes", low = "#5E3C99", mid = "grey97", high = "#E66101",
                                      midpoint = 0.5, limits = c(0, 1), labels = function(x) paste0(100 * x, "%"))
  }
}

#' Units in Jensen--Shannon principal coordinates, coloured by bloc
#'
#' Treats each unit's Yes share on each question as a Bernoulli distribution
#' (with a Jeffreys adjustment, (Yes + 0.5)/(N + 1)), sums the Jensen--Shannon
#' divergence between two units over questions, takes the square root as a
#' distance, and plots the leading principal coordinates (classical
#' multidimensional scaling) in pairs. Units missing a question are dropped
#' (after dropping questions missing in more than `max_missing` of units).
#' If the blocs are slices of a single axis, they will line up along PC1;
#' blocs separated on PC2 or PC3 differ in ways one left--right score misses.
#'
#' @param data A [bbmm_data] object.
#' @param summary A [bloc_summary] object for `data`.
#' @param n_pc Principal coordinates to compute (pairs among the first three are plotted).
#' @param max_missing Drop questions missing in more than this share of units first.
#' @param orient Optional numeric vector, one value per unit (e.g. a party's
#'   vote share): PC1 is flipped so that it rises with it.
#' @return A ggplot object (three facets: PC1 v PC2, PC1 v PC3, PC2 v PC3) with
#'   the coordinates, percent of the total each axis carries, and the share of
#'   each axis explained by bloc (R squared) in the facet titles.
#' @examples
#' small <- maine[1:120, 1:20]
#' s <- bloc_summary(bbmm_em(small, K = 3, zeta = 0.3, n_starts = 3, seed = 1), small)
#' plot_js_pca(small, s)
#' @export
plot_js_pca <- function(data, summary, n_pc = 3L, max_missing = 0.02, orient = NULL) {
  check_data(data)
  js <- js_pcoa(data, n_pc = max(3L, n_pc), max_missing = max_missing)
  co <- js$coords
  bloc <- summary$assignment$bloc[match(rownames(co), summary$assignment$unit)]
  if (!is.null(orient)) {
    o <- orient[match(rownames(co), data$units)]
    if (stats::cor(co[, 1], o, use = "complete.obs") < 0) co[, 1] <- -co[, 1]
  }
  r2 <- vapply(1:3, function(j) summary(lm(co[, j] ~ factor(bloc)))$r.squared, 1)
  ax <- sprintf("PC%d (%.1f%%, bloc R\u00b2 %.2f)", 1:3, js$pct[1:3], r2)
  pairs <- list(c(1, 2), c(1, 3), c(2, 3))
  long <- do.call(rbind, lapply(pairs, function(pr) data.frame(
    panel = sprintf("x: %s\ny: %s", ax[pr[1]], ax[pr[2]]), x = co[, pr[1]], y = co[, pr[2]],
    bloc = factor(bloc, levels = summary$blocs$bloc))))
  long$panel <- factor(long$panel, levels = unique(long$panel))
  ctr <- stats::aggregate(cbind(x, y) ~ panel + bloc, long, stats::median)
  ggplot2::ggplot(long, ggplot2::aes(.data$x, .data$y, colour = .data$bloc)) +
    ggplot2::geom_point(size = 0.9, alpha = 0.6) +
    ggplot2::geom_label(data = ctr, ggplot2::aes(label = .data$bloc), size = 3.2, fontface = "bold",
                        show.legend = FALSE, fill = grDevices::adjustcolor("white", 0.85)) +
    ggplot2::facet_wrap(~panel, scales = "free", nrow = 1) +
    bloc_scale(summary, "colour") +
    ggplot2::guides(colour = ggplot2::guide_legend(override.aes = list(size = 3, alpha = 1))) +
    ggplot2::theme_minimal(base_size = 10) + ggplot2::theme(legend.position = "bottom", aspect.ratio = 1) +
    ggplot2::labs(x = NULL, y = NULL, title = "Units in Jensen-Shannon principal coordinates",
                  subtitle = sprintf("%d units complete on %d questions; number = bloc median", nrow(co), js$n_questions))
}

#' Jensen--Shannon principal coordinates of the units
#'
#' The computation behind [plot_js_pca()].
#'
#' @inheritParams plot_js_pca
#' @return A list: `coords` (units x `n_pc`), `values` (eigenvalues), `pct`
#'   (percent of the trace), `n_questions`, `dropped_questions`, `dropped_units`.
#' @examples
#' js <- js_pcoa(maine[1:100, 1:20])
#' round(js$pct, 1)
#' @export
js_pcoa <- function(data, n_pc = 3L, max_missing = 0.02) {
  check_data(data)
  miss <- is.na(data$M)
  qk <- colMeans(miss) <= max_missing
  rk <- rowSums(miss[, qk, drop = FALSE]) == 0
  if (sum(rk) < 3) stop("Fewer than three units are complete on the retained questions.", call. = FALSE)
  p <- (data$X[rk, qk, drop = FALSE] + 0.5) / (data$M[rk, qk, drop = FALSE] + 1)
  ent <- function(z) -(z * log(z) + (1 - z) * log(1 - z))
  N <- nrow(p)
  D2 <- matrix(0, N, N)
  for (q in seq_len(ncol(p))) {
    pq <- p[, q]
    D2 <- D2 + ent(outer(pq, pq, "+") / 2) - outer(ent(pq), ent(pq), "+") / 2
  }
  D2[D2 < 0] <- 0; diag(D2) <- 0
  B <- -0.5 * D2; B <- sweep(B, 1, rowMeans(B)); B <- sweep(B, 2, colMeans(B))
  k <- min(n_pc, N - 1)
  e <- if (N > 2000 && requireNamespace("RSpectra", quietly = TRUE)) RSpectra::eigs_sym(B, k, which = "LA") else {
    ee <- eigen(B, symmetric = TRUE); list(values = ee$values[1:k], vectors = ee$vectors[, 1:k, drop = FALSE]) }
  coords <- sweep(e$vectors, 2, sqrt(pmax(e$values, 0)), "*")
  dimnames(coords) <- list(data$units[rk], paste0("PC", seq_len(k)))
  list(coords = coords, values = e$values, pct = 100 * e$values / sum(diag(B)), n_questions = sum(qk),
       dropped_questions = data$questions[!qk], dropped_units = data$units[!rk])
}

#' Plot the posterior on the number of blocs
#'
#' Three panels from a [bbmm_evidence()] run: the log evidence by \eqn{K}
#' (relative to the best, with replicate means), the posterior \eqn{p(K \mid
#' x)} over mixture components, and \eqn{p(T \mid x)} over occupied blocs.
#'
#' @param evidence A `bbmm_evidence` object.
#' @return A ggplot object.
#' @examples
#' tiny <- maine[1:40, 1:6]
#' ev <- bbmm_evidence(tiny, K = 1:3, zeta = 0.3, seed = 1, progress = FALSE,
#'                     control = list(n_particles = 50, cess_target = 0.5))
#' plot_evidence(ev)
#' @export
plot_evidence <- function(evidence) {
  stopifnot(inherits(evidence, "bbmm_evidence"))
  pk <- evidence$posterior_K; runs <- evidence$runs
  top <- max(pk$logZ)
  d <- rbind(data.frame(panel = "log p(x | K), relative to best (nats)", n = runs$K, v = runs$logZ - top, kind = "pt"),
             data.frame(panel = "log p(x | K), relative to best (nats)", n = pk$K, v = pk$logZ - top, kind = "line"),
             data.frame(panel = "p(K | x): mixture components", n = pk$K, v = pk$posterior, kind = "bar"))
  if (!is.null(evidence$posterior_T))
    d <- rbind(d, data.frame(panel = "p(T | x): occupied blocs", n = evidence$posterior_T$T,
                             v = evidence$posterior_T$posterior, kind = "bar"))
  d$panel <- factor(d$panel, levels = unique(d$panel))
  ggplot2::ggplot(d, ggplot2::aes(.data$n, .data$v)) +
    ggplot2::geom_col(data = d[d$kind == "bar", ], fill = "grey35", width = 0.7) +
    ggplot2::geom_line(data = d[d$kind == "line", ], colour = "grey55") +
    ggplot2::geom_point(data = d[d$kind == "pt", ], size = 1.8) +
    ggplot2::facet_wrap(~panel, scales = "free_y", nrow = 1) +
    ggplot2::scale_x_continuous(breaks = function(l) seq(ceiling(l[1]), floor(l[2]))) +
    ggplot2::theme_minimal(base_size = 10) + ggplot2::theme(panel.grid.minor = ggplot2::element_blank()) +
    ggplot2::labs(x = "K or T", y = NULL, title = sprintf("Number of blocs, stepping-stone SMC (zeta = %g)", evidence$zeta))
}

#' Plot EM model-selection criteria
#'
#' BIC and the smallest bloc's share of units against \eqn{K}, from
#' [bbmm_em_select()]. Dashed lines mark the minimum-share thresholds.
#'
#' @param selection A `bbmm_em_select` object.
#' @param min_share Thresholds to draw.
#' @return A ggplot object.
#' @examples
#' small <- maine[1:120, 1:20]
#' plot_em_selection(bbmm_em_select(small, K = 1:4, zeta = 0.3, n_starts = 2, seed = 1))
#' @export
plot_em_selection <- function(selection, min_share = c(0.03, 0.05)) {
  stopifnot(inherits(selection, "bbmm_em_select"))
  t <- selection$table
  d <- rbind(data.frame(panel = "BIC (lower is better)", K = t$K, v = t$BIC),
             data.frame(panel = "Smallest bloc, share of units", K = t$K, v = t$min_share))
  hl <- data.frame(panel = "Smallest bloc, share of units", v = min_share)
  ggplot2::ggplot(d, ggplot2::aes(.data$K, .data$v)) + ggplot2::geom_line(colour = "grey55") +
    ggplot2::geom_point(size = 1.8) +
    ggplot2::geom_hline(data = hl, ggplot2::aes(yintercept = .data$v), linetype = 2, colour = "grey40") +
    ggplot2::facet_wrap(~panel, scales = "free_y") +
    ggplot2::scale_x_continuous(breaks = t$K) +
    ggplot2::theme_minimal(base_size = 10) + ggplot2::theme(panel.grid.minor = ggplot2::element_blank()) +
    ggplot2::labs(x = "K", y = NULL, title = sprintf("Tempered EM over K (zeta = %g)", selection$zeta))
}

#' Plot the posterior similarity matrix
#'
#' The probability that two units share a bloc, with units ordered by their
#' representative bloc. Sharp blocks on the diagonal mean a well-determined
#' partition; grey off-diagonal bands show units that move between blocs.
#'
#' @param summary A [bloc_summary] from a sampler (with `psm`).
#' @param max_units Draw at most this many units (a random subset, kept in bloc order).
#' @return A ggplot object.
#' @examples
#' tiny <- maine[1:60, 1:10]
#' fit <- bbmm_mcmc(tiny, K = 3, zeta = 0.3, n_iter = 100, n_burn = 50, seed = 1, progress = FALSE)
#' plot_coclustering(bloc_summary(fit, tiny))
#' @export
plot_coclustering <- function(summary, max_units = 1500L) {
  if (is.null(summary$psm)) stop("This summary has no similarity matrix (EM fits do not).", call. = FALSE)
  a <- summary$assignment
  o <- order(a$bloc, -a$certainty)
  if (length(o) > max_units) o <- o[sort(sample.int(length(o), max_units))]   # keeps the bloc order
  P <- summary$psm[o, o]
  n <- length(o)
  d <- data.frame(i = rep(seq_len(n), n), j = rep(seq_len(n), each = n), v = as.vector(P))
  edges <- cumsum(table(a$bloc[o]))
  ggplot2::ggplot(d, ggplot2::aes(.data$i, .data$j, fill = .data$v)) + ggplot2::geom_raster() +
    ggplot2::geom_vline(xintercept = edges + 0.5, colour = "white", linewidth = 0.2) +
    ggplot2::geom_hline(yintercept = edges + 0.5, colour = "white", linewidth = 0.2) +
    ggplot2::scale_fill_gradient(low = "white", high = "black", limits = c(0, 1), name = "P(same bloc)") +
    ggplot2::scale_y_reverse() + ggplot2::coord_equal() + ggplot2::theme_void(base_size = 10) +
    ggplot2::labs(title = "Posterior similarity of units", subtitle = "units ordered by representative bloc")
}

#' Trace plots for bbmm_mcmc
#'
#' The complete-data log-likelihood and the number of occupied blocs at
#' each saved iteration, one line per chain.
#'
#' @param fit A `bbmm_mcmc` object.
#' @return A ggplot object.
#' @examples
#' tiny <- maine[1:60, 1:10]
#' fit <- bbmm_mcmc(tiny, K = 3, zeta = 0.3, n_iter = 100, n_burn = 50, seed = 1, progress = FALSE)
#' plot_trace(fit)
#' @export
plot_trace <- function(fit) {
  if (!all(c("ll", "T", "chain", "iteration") %in% names(fit))) stop("`fit` must be a bbmm_mcmc object.", call. = FALSE)
  d <- rbind(data.frame(panel = "log-likelihood", it = fit$iteration, v = fit$ll, chain = factor(fit$chain)),
             data.frame(panel = "occupied blocs", it = fit$iteration, v = fit$T, chain = factor(fit$chain)))
  ggplot2::ggplot(d, ggplot2::aes(.data$it, .data$v, colour = .data$chain)) + ggplot2::geom_line(linewidth = 0.3) +
    ggplot2::facet_wrap(~panel, scales = "free_y", ncol = 1) + ggplot2::theme_minimal(base_size = 10) +
    ggplot2::scale_colour_manual(values = bloc_colours(fit$n_chains)) +
    ggplot2::labs(x = "iteration", y = NULL, colour = "Chain",
                  title = sprintf("Fixed-K MCMC traces (K = %d, zeta = %g)", fit$K, fit$zeta))
}

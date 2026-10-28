# data-raw/tutorial_results.R -- the runs behind the tutorial (vignettes/Elections.BBMM.Rmd, TUTORIAL.md),
# saved slimmed to inst/extdata/maine_tutorial.rds so the tutorial can show real results without
# rerunning them.  Rscript data-raw/tutorial_results.R [threads]   (from the package root)
# Run time on an Apple M4 Pro (10 threads): EM about 1 minute; stepping-stone K = 2..10 about
# 1.5 hours; MCMC a few minutes. Per-K SS fits are cached in data-raw/tutorial_cache/.
pkgload::load_all(quiet = TRUE)
thr <- as.integer(commandArgs(TRUE)[1]); if (is.na(thr)) thr <- 8L
cache <- "data-raw/tutorial_cache"; dir.create(cache, showWarnings = FALSE)
cached <- function(name, expr) { f <- file.path(cache, paste0(name, ".rds"))
  if (file.exists(f)) return(readRDS(f)); v <- expr; saveRDS(v, f); v }

sel <- cached("em_select", bbmm_em_select(maine, K = 1:12, zeta = 0.3, n_starts = 10, seed = 2026, cores = thr))
K_em <- sel$selected[["min_share_5%"]]
fits <- lapply(2:10, function(K) cached(sprintf("ss_K%02d", K),
  bbmm_ss(maine, K = K, zeta = 0.3, seed = 2026 + K, n_threads = thr, progress = TRUE)))
reps <- lapply(c(6L, 10L), function(K) cached(sprintf("ss_K%02d_rep2", K),
  bbmm_ss(maine, K = K, zeta = 0.3, seed = 4026 + K, n_threads = thr, progress = TRUE)))
ev <- evidence_from_fits(c(fits, reps))
T_mode <- ev$posterior_T$T[which.max(ev$posterior_T$posterior)]
mc <- cached("mcmc", bbmm_mcmc(maine, K = T_mode, zeta = 0.3, n_iter = 2000, n_burn = 1000, thin = 2,
                               n_chains = 4, seed = 2026, n_threads = min(4L, thr), progress = TRUE))

s_em <- bloc_summary(sel$fits[[as.character(K_em)]], maine)
s_ss <- bloc_summary(ev, maine)
s_mc <- bloc_summary(mc, maine)
slim_sel <- sel; slim_sel$fits <- slim_sel$fits[as.character(K_em)]
slim_ev <- ev; slim_ev$fits <- NULL
slim_mc <- mc[c("ll", "T", "chain", "iteration", "acceptance", "K", "zeta", "n_iter", "n_burn", "thin", "n_chains",
                "seconds", "seed", "init")]
class(slim_mc) <- "bbmm_mcmc_trace"
round_psm <- function(s) { s$psm <- round(s$psm, 3); s }
res <- list(sel = slim_sel, ev = slim_ev, mc = slim_mc, s_em = s_em, s_ss = round_psm(s_ss), s_mc = round_psm(s_mc),
            ss_seconds = vapply(fits, `[[`, 1, "seconds"), ss_steps = vapply(fits, `[[`, 1L, "n_steps"),
            mc_diag = mcmc_diagnostics(mc),
            mc_ari_ss = vapply(tapply(seq_along(mc$chain), mc$chain, max), function(i)
              adjusted_rand(mc$z[i, ], s_ss$assignment$bloc), 1), threads = thr, machine = "Apple M4 Pro, 14 cores", date = Sys.Date())
saveRDS(res, "inst/extdata/maine_tutorial.rds", compress = "xz")
cat(sprintf("saved: %.0f KB\n", file.size("inst/extdata/maine_tutorial.rds") / 1024))

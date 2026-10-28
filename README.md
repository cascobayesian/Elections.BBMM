# Elections.BBMM: identifying voting blocs from referendum data

The Elections.BBMM package fits a beta-binomial mixture model (BBMM) to longitudinal referendum data. The model identifies **voting blocs**: long-term patterns of vote similarity across a set of political districts (precincts, municipalities, counties). Districts in the same bloc share a profile of Yes probabilities across every question on the ballot, and the beta-binomial allows for district-to-district variation within a bloc. The results can be used to infer:

* the distribution on the number of voting blocs;
* which bloc each district belongs to, and how certain that assignment is;
* each bloc's expected Yes share, with uncertainty, on every question, and how tightly its districts cluster around it;
* which questions set each bloc apart, and which questions cut across the blocs (questions that behave unusually from the model's point of view);
* an information-theoretic (Jensen–Shannon) view of how districts differ and where the blocs sit in that space; and
* the co-clustering structure of the districts: which ones always vote together and which move between blocs.

This kind of analysis might appeal to political scientists, sociologists or historians. Why might you want to perform it?

* Identify political voting blocs and their presence in a set of communities;
* Use voting blocs to identify long-term political and cultural structures and dynamics;
* Identify questions that sharply distinguish voting blocs; and
* Identify questions that run counter to the inferred blocs.

## What is in the package

Three ways to fit the model, all with a **likelihood power ζ** (tempering, or the coarsened posterior of Miller and Dunson 2019). With hundreds of votes per district, the ordinary likelihood (ζ = 1) treats every small, persistent difference as decisive and keeps splitting the electorate into more blocs; ζ < 1 asks for structure that holds up to that kind of misspecification. Smaller ζ gives fewer, broader blocs.

| Function | Method | Speed | Use it for |
|---|---|---|---|
| `bbmm_em()`, `bbmm_em_select()` | Tempered EM at a fixed dispersion ν, many random starts; BIC and minimum-bloc-share rules for K | Seconds to minutes | A first look; choosing K quickly; large data |
| `bbmm_mcmc()` | Fixed-K MCMC on the tempered posterior: split–merge, collapsed allocation updates, random-walk updates of the bloc profiles and per-bloc dispersions | Minutes | Uncertainty about the partition and the bloc profiles at a chosen K |
| `bbmm_ss()`, `bbmm_evidence()` | Stepping-stone sequential Monte Carlo: the marginal likelihood p(x \| K) for each K, the posterior on the number of blocs, and a posterior sample at every K | Minutes to hours per K | The number of blocs, with a full Bayesian account |

Then, for the output of any of the three:

| Function | Gives |
|---|---|
| `bloc_summary()` | Unit assignments and certainty, bloc sizes, bloc-by-question Yes shares with intervals; for the samplers, the posterior similarity matrix and a representative partition |
| `distinctive_questions()`, `question_fit()` | The questions that set each bloc apart; how much of each question's variation the blocs explain |
| `plot_bloc_map()` | Map of the blocs (or of assignment certainty) |
| `plot_bloc_profiles()` | Heatmap of Yes shares by bloc and question |
| `plot_js_pca()`, `js_pcoa()` | Districts in Jensen–Shannon principal coordinates, coloured by bloc |
| `plot_evidence()`, `plot_em_selection()` | Posterior on the number of blocs; EM selection criteria |
| `plot_coclustering()`, `plot_trace()`, `mcmc_diagnostics()` | Posterior similarity; MCMC traces and R-hat |

The package includes test data: Maine's 72 statewide referendums of 2008–2025 by town (`maine`, `maine_questions`) with town boundaries (`maine_towns()`). The raw files behind it, and the earlier Massachusetts, Multnomah County and Portland data sets, are in [`data-raw/`](data-raw).

## Installation

```r
# install.packages("remotes")
remotes::install_github("cascobayesian/Elections.BBMM", build_vignettes = TRUE)
```

The package compiles C++ code, so you need a compiler: Xcode command-line tools on macOS (`xcode-select --install`), Rtools on Windows, `build-essential` on Linux. There are no other system requirements; maps need the **sf** package.

## Quick start

```r
library(Elections.BBMM)
maine                                    # 462 towns x 72 questions

# Fast: tempered EM over K = 1..12, choose K by the 3% rule
sel <- bbmm_em_select(maine, K = 1:12, zeta = 0.3, n_starts = 10, cores = 4)
fit <- sel$fits[[as.character(sel$selected[["min_share_3%"]])]]
s <- bloc_summary(fit, maine)
plot_bloc_map(s, maine_towns())
plot_bloc_profiles(s)
plot_js_pca(maine, s)

# Slow, fully Bayesian: evidence for each K by stepping-stone SMC (about an hour on 8 threads)
ev <- bbmm_evidence(maine, K = 2:10, zeta = 0.3, n_threads = 8)
plot_evidence(ev)
s_ss <- bloc_summary(ev, maine)
```

## Tutorial

**[TUTORIAL.md](TUTORIAL.md)** walks through installation, the data format, fitting Maine with EM and with the stepping-stone sampler (with timings), and reading the posterior summaries. The same material is the package vignette: `vignette("Elections.BBMM")`.

## Model

For unit *i* in bloc *k*, the Yes count on question *q* is BetaBinomial(*M*<sub>iq</sub>, ν<sub>k</sub> *p*<sub>kq</sub>, ν<sub>k</sub>(1 − *p*<sub>kq</sub>)): *p*<sub>kq</sub> is the bloc's expected Yes share and ν<sub>k</sub> its concentration. Questions not on a unit's ballot contribute nothing. Priors in the Bayesian samplers: *p*<sub>kq</sub> ~ Beta(10, 10); ν<sub>k</sub> ~ Gamma(4, rate 0.04) truncated to ν ≥ 50, so no bloc can become a catch-all for districts that fit nowhere else; bloc weights ~ Dirichlet(1, …, 1), integrated out; K ~ Geometric(0.1) on 1–50 for the posterior on K. The target is prior × likelihood<sup>ζ</sup>. EM uses one fixed ν (default 100) and the same Beta prior on *p*.

The stepping-stone code has been checked against exact evidences computed by enumeration and quadrature on small data sets (errors below 0.07 nats), and its kernel against exact tempered posteriors over partitions.

## Acknowledgements

Nathan Kellerman contributed to the EM algorithm as part of an undergraduate honors thesis, *A Rapid Algorithm for Inferring Latent Mixture Structure in Replicate Social Science Data* (Bowdoin College, 2026). Aaron Kaufman (New York University Abu Dhabi) helped develop the tools for exploring the model results.

## Citation

Kaufman, A. and O'Brien, J. D. The structure of conflict in local politics. Working paper.

## License

GPL (≥ 3).

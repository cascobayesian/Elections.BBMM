# Elections.BBMM 1.0.0

A rewrite. The birth-death sampler of version 0.1.0 is replaced by three engines with a common
likelihood power `zeta` (tempering / coarsened posterior):

* `bbmm_em()`, `bbmm_em_select()`: tempered EM at a fixed dispersion, with multiple starts and
  BIC and minimum-share rules for the number of blocs.
* `bbmm_mcmc()`: fixed-K MCMC on the tempered posterior with per-bloc dispersions and
  Jain-Neal split-merge moves.
* `bbmm_ss()`, `bbmm_evidence()`, `evidence_from_fits()`: stepping-stone sequential Monte Carlo
  for the marginal likelihood at each K and the posterior on the number of blocs.

Posterior summaries (`bloc_summary()`, `distinctive_questions()`, `question_fit()`) and figures
(`plot_bloc_map()`, `plot_bloc_profiles()`, `plot_js_pca()`, `plot_evidence()`,
`plot_coclustering()`, `plot_trace()`, `plot_em_selection()`) work on all three.

Test data: Maine's statewide referendums of 2008-2025 by town (`maine`, `maine_questions`,
`maine_towns()`).

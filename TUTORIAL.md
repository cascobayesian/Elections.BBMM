Tutorial: voting blocs in Maine with Elections.BBMM
================

- [1. Installing the package and getting the test
  data](#1-installing-the-package-and-getting-the-test-data)
  - [Your own data](#your-own-data)
- [2. The model in brief](#2-the-model-in-brief)
- [3. Fast: tempered EM](#3-fast-tempered-em)
- [4. Slow: the stepping-stone
  sampler](#4-slow-the-stepping-stone-sampler)
  - [Fixed-K MCMC](#fixed-k-mcmc)
- [5. Posterior summaries and how to read
  them](#5-posterior-summaries-and-how-to-read-them)
  - [Where the blocs are](#where-the-blocs-are)
  - [What the blocs vote for](#what-the-blocs-vote-for)
  - [Questions the blocs do not
    explain](#questions-the-blocs-do-not-explain)
  - [The shape of the blocs: Jensen–Shannon principal
    coordinates](#the-shape-of-the-blocs-jensenshannon-principal-coordinates)
  - [How settled is the partition?](#how-settled-is-the-partition)
  - [Comparing the methods](#comparing-the-methods)
- [6. Where to go next](#6-where-to-go-next)
  - [References](#references)

This tutorial fits the beta-binomial mixture model (BBMM) to the test
data that ship with the package – Maine’s 72 statewide referendums of
2008–2025 in its 462 towns – in two ways: a fast tempered EM fit, and
the slow but fully Bayesian stepping-stone sampler. It then walks
through the posterior summaries and how to read them.

Code chunks marked **long** were run once on an Apple M4 Pro laptop (12
threads); their results are loaded from the package so that this
document builds quickly. Running them yourself gives the same numbers
(the seeds are set), on the time scales quoted.

# 1. Installing the package and getting the test data

Install from GitHub. The package compiles C++ code, so you need a
compiler: the Xcode command-line tools on macOS
(`xcode-select --install`), Rtools on Windows, or `build-essential` on
Linux. Maps also need the **sf** package.

``` r
install.packages(c("remotes", "sf"))
remotes::install_github("cascobayesian/Elections.BBMM", build_vignettes = TRUE)
library(Elections.BBMM)
```

The test data are part of the package: `maine` (the votes),
`maine_questions` (what each question was) and `maine_towns()` (town
boundaries).

``` r
maine
#> BBMM data: 462 units x 72 questions, 37,279,962 votes; 0.1% of cells empty
#> Median votes per cell 504; pooled Yes share 0.555
head(maine_questions[, c("question_id", "label", "statewide_yes")])
#>         question_id                                         label statewide_yes
#> 1 REFERENDUM_2008_1         2008 Q1: Dirigo Health funding repeal         64.65
#> 2 REFERENDUM_2008_2                 2008 Q2: Oxford County casino         46.16
#> 3 REFERENDUM_2008_3            2008 Q3: Water infrastructure bond         50.27
#> 4 REFERENDUM_2009_1             2009 Q1: Same-sex marriage repeal         52.90
#> 5 REFERENDUM_2009_2               2009 Q2: Vehicle tax reductions         25.89
#> 6 REFERENDUM_2009_3 2009 Q3: School district consolidation repeal         41.68
maine$X[1:3, 1:4]   # Yes votes
#>         REFERENDUM_2008_1 REFERENDUM_2008_2 REFERENDUM_2008_3 REFERENDUM_2009_1
#> ABBOT                 296               190               117               242
#> ACTON                 930               699               537               564
#> ADDISON               403               161               177               324
maine$M[1:3, 1:4]   # Yes + No votes
#>         REFERENDUM_2008_1 REFERENDUM_2008_2 REFERENDUM_2008_3 REFERENDUM_2009_1
#> ABBOT                 389               392               384               339
#> ACTON                1335              1340              1299              1020
#> ADDISON               610               613               599               479
```

A `bbmm_data` object is two matrices with one row per unit (here, a town
or a group of towns that report together) and one column per question:
Yes votes `X` and two-option totals `M`. Blank ballots are left out. A
cell with no votes – a question that was not on that unit’s ballot – is
`NA` and simply contributes nothing.

## Your own data

Election offices usually publish one row per unit with a Yes and a No
column per question. The raw Maine file in the repository has that
layout, and `bbmm_data_from_wide()` reads it:

``` r
url <- paste0("https://raw.githubusercontent.com/cascobayesian/Elections.BBMM/main/",
              "data-raw/Maine.Ref/maine_referendum_2008_2025.csv")
raw <- read.csv(url, check.names = FALSE)
names(raw)[1:5]
#> "TOWN" "REFERENDUM_2008_1_YES" "REFERENDUM_2008_1_NO" "REFERENDUM_2008_2_YES" ...
my_data <- bbmm_data_from_wide(raw, id = "TOWN", yes_suffix = "_YES", no_suffix = "_NO")
```

With two matrices of Yes counts and totals, use `bbmm_data(yes, total)`
instead. Units must be the same across all questions, so if district
lines moved over the years, put every election on one set of units
first.

# 2. The model in brief

Each unit belongs to one bloc. A bloc *k* has a Yes probability $p_{kq}$
for every question *q*; a unit in bloc *k* casts its $M_{iq}$ votes on
question *q* as a beta-binomial draw with mean $p_{kq}$ and
concentration $\nu_k$, so units in a bloc scatter around the bloc’s
profile, more tightly when $\nu_k$ is large. The blocs are found from
all questions at once: two towns are in the same bloc when they vote
alike across the board, not on any one question.

**Tempering.** A town casts hundreds or thousands of votes on each of 72
questions, so the ordinary likelihood treats every small, persistent
difference between towns as decisive, and the number of blocs keeps
growing with the data. All three fitting methods therefore raise the
likelihood to a power $\zeta \le 1$ (the coarsened posterior of Miller
and Dunson, 2019): roughly, they look for structure that would survive
if the data were a fraction $\zeta$ of their size. Smaller $\zeta$ gives
fewer, broader blocs. This tutorial uses the untempered model,
$\zeta = 1$, throughout, so you can see what the full likelihood does
with data this large; the same analysis at $\zeta = 0.3$ finds five
blocs instead of sixteen. Treat $\zeta$ as a choice about the grain of
the analysis, and check that the conclusions you care about hold at
neighbouring values.

# 3. Fast: tempered EM

`bbmm_em_select()` fits the model by expectation-maximisation for each
*K* in a range, from several random starts each, holding the
concentration fixed at $\nu = 100$ for every bloc.

``` r
# long: about 1 minute on 12 cores, 4 minutes on 4
sel <- bbmm_em_select(maine, K = 1:16, zeta = 1, n_starts = 10, seed = 2026, cores = 12)
```

``` r
sel
#> BBMM tempered EM over K = 1-16 (zeta = 1, nu = 100; 462 units x 72 questions)
#>   K  loglik occupied min_share    BIC converged n_at_best
#>   1 -180637        1 1.0000000 362023      TRUE        10
#>   2 -163445        2 0.3225108 328400      TRUE        10
#>   3 -158441        3 0.1601732 319152      TRUE        10
#>   4 -156536        4 0.1385281 316102      TRUE         9
#>   5 -155502        5 0.1212121 314795      TRUE         4
#>   6 -154571        6 0.0844156 313693      TRUE         2
#>   7 -153845        7 0.0454545 313001      TRUE         1
#>   8 -153406        8 0.0411255 312881      TRUE         1
#>   9 -153086        9 0.0411255 313003      TRUE         1
#>  10 -152890       10 0.0238095 313370      TRUE         1
#>  11 -152547       11 0.0324675 313444      TRUE         1
#>  12 -152421       12 0.0108225 313953      TRUE         1
#>  13 -152029       13 0.0259740 313928      TRUE         1
#>  14 -151778       14 0.0238095 314186      TRUE         1
#>  15 -151611       15 0.0173160 314614      TRUE         1
#>  16 -151629       15 0.0259740 314648      TRUE         1
#> 
#> Selected number of blocs:
#>          BIC min_share_3% min_share_5% 
#>            8           11            6
```

Each row is the best of the 10 starts at that *K*: its tempered
log-likelihood, the number of components that hold at least one town
(`occupied`), the smallest bloc’s share of towns, BIC, and how many
starts reached that same best fit (`n_at_best`). When few starts agree,
the partition at that *K* is one of several nearly equally good ones,
and conclusions about individual small blocs should be held loosely.

Two rules for choosing the number of blocs are reported. BIC balances
fit against the number of parameters; with tempered likelihoods it is
only a heuristic. The minimum-share rule takes the largest *K* at which
every bloc still holds at least 3% (or 5%) of the units – it stops
before the model starts carving off blocs of a handful of towns.

``` r
plot_em_selection(sel)
```

<img src="TUTORIAL_files/figure-gfm/em-plot-1.png" alt="" width="100%" />

On Maine the rules disagree: BIC picks 8 blocs, the 5% rule 6 and the 3%
rule 11. The log-likelihood keeps rising all the way to *K* = 16 – the
untempered likelihood always finds something to gain from another bloc –
and from *K* = 5 on fewer than half of the starts reach the best fit;
from *K* = 7 on only one does, so the extra blocs are not well
determined. We carry the 5% fit forward as a coarse summary; the
stepping-stone sampler in the next section addresses the number of blocs
directly.

Pick a fit and summarise it:

``` r
K_em <- sel$selected[["min_share_5%"]]
fit_em <- sel$fits[[as.character(K_em)]]
s_em <- bloc_summary(fit_em, maine)
s_em
#> Bloc summary: tempered EM, K = 6 (6 occupied), zeta = 1, fixed nu = 100 
#>  bloc units share    ballots mean_yes  nu certainty
#>     1    39 0.084  8,270,147    0.636 100     0.999
#>     2    62 0.134  8,303,343    0.583 100     0.998
#>     3    97 0.210 10,538,699    0.546 100     0.995
#>     4    40 0.087  1,589,005    0.540 100     0.997
#>     5   131 0.284  6,840,626    0.507 100     0.988
#>     6    93 0.201  1,738,142    0.461 100     0.992
```

``` r
plot_bloc_map(s_em, maine_towns())
```

<img src="TUTORIAL_files/figure-gfm/em-map-1.png" alt="" width="100%" />

To fit a single *K*, use `bbmm_em(maine, K = 6, zeta = 1)`. EM also
returns each town’s membership probabilities (`fit_em$gamma`); at this
$\zeta$ almost all are close to 0 or 1, which says more about the
fixed-$\nu$ model than about how certain the partition is.

# 4. Slow: the stepping-stone sampler

EM gives one best partition at each *K* and no measure of uncertainty,
and choosing *K* rests on rules of thumb. The Bayesian version of the
model puts priors on everything – including a separate concentration
$\nu_k$ for each bloc – and the stepping-stone sampler computes the
**marginal likelihood** $p(x \mid K)$ for each *K*. With a prior on *K*
(geometric, favouring fewer blocs) that gives the posterior probability
of each number of blocs, and every run also leaves a weighted sample of
partitions from the posterior at its *K*.

The sampler moves a population of 200 “particles” from the prior to the
tempered posterior through a few hundred intermediate distributions, so
it is far slower than EM.

> **How long will this take?** On Maine (462 towns x 72 questions) at
> $\zeta = 1$, one *K* took from about 4 minutes (*K* = 2) to about 17
> minutes (*K* = 16) on 12 threads of an otherwise idle machine, and *K*
> = 2 to 20 took about four hours in all, plus an hour for two
> replicates. Time grows with units x questions, with *K*, and with
> $\zeta$ (a larger $\zeta$ means a longer path from prior to
> posterior). A state’s precincts (thousands of units) can take hours to
> days per *K*. Start with a range of *K* around the EM answer, extend
> it until the evidence turns down, use as many threads as you can spare
> (`n_threads`), and for big data run each *K* as a separate job and
> combine them with `evidence_from_fits()`.

``` r
# long: about 4 hours on 12 threads (plus two extra replicates, 1 hour)
ev <- bbmm_evidence(maine, K = 2:20, zeta = 1, n_threads = 12, seed = 2026)

# or, one K at a time (e.g. as separate jobs), then combine:
fits <- lapply(2:20, function(K) bbmm_ss(maine, K = K, zeta = 1, seed = 2026 + K, n_threads = 12))
# a second replicate at K = 16 and K = 18, to see the Monte Carlo error
reps <- lapply(c(16, 18), function(K) bbmm_ss(maine, K = K, zeta = 1, seed = 4026 + K, n_threads = 12))
ev <- evidence_from_fits(c(fits, reps))
ev$runs        # one row per run
```

``` r
ev
#> BBMM stepping-stone evidence, zeta = 1: 21 run(s)
#>   K      logZ    sd posterior
#>   2 -23950454    NA  0.00e+00
#>   3 -23947026    NA  0.00e+00
#>   4 -23945584    NA  0.00e+00
#>   5 -23944609    NA  0.00e+00
#>   6 -23944131    NA  0.00e+00
#>   7 -23943566    NA  0.00e+00
#>   8 -23943358    NA  0.00e+00
#>   9 -23943116    NA 6.78e-278
#>  10 -23942951    NA 4.68e-206
#>  11 -23942790    NA 1.72e-136
#>  12 -23942738    NA 5.99e-114
#>  13 -23942632    NA  5.47e-68
#>  14 -23942609    NA  5.35e-58
#>  15 -23942583    NA  9.04e-47
#>  16 -23942510 14.75  7.10e-15
#>  17 -23942500    NA  1.80e-10
#>  18 -23942477 36.39  1.00e+00
#>  19 -23942515    NA  4.37e-17
#>  20 -23942577    NA  2.78e-44
#> Posterior on occupied blocs T: 16: 0.97, 17: 0.03
```

``` r
plot_evidence(ev)
```

<img src="TUTORIAL_files/figure-gfm/ss-plot-1.png" alt="" width="100%" />

Reading this output:

- `logZ` is the log marginal likelihood at each *K*, in nats (with the
  mean over replicates where there are several, and their `sd`).
  Differences matter, not levels: a gap of 5 nats is a factor of about
  150 in posterior odds.
- The evidence climbs steeply to about *K* = 8, keeps climbing more
  slowly to *K* = 16, and is flat from *K* = 16 to 18 before turning
  down at 19 and 20. Without tempering, the data support a lot of blocs:
  each extra bloc buys hundreds of nats up to about ten, and tens of
  nats after that.
- Near the top, the differences are the size of the Monte Carlo error.
  The two runs at *K* = 16 differ by 21 nats and the two at *K* = 18
  by 51. The stepping-stone estimate is noisy and biased low when
  particles lag behind, so one lucky run can take nearly all of
  $p(K \mid x)$, as *K* = 18 does here. Do not read $p(K \mid x)$
  literally when the gaps are this small; run replicates (`n_rep` in
  `bbmm_evidence()`) at the competing values.
- $p(K \mid x)$ is over mixture *components*, and components can stay
  empty. The robust quantity is $p(T \mid x)$, the posterior on the
  number of *occupied* blocs, averaged over *K*. The runs at *K* = 16 to
  20 all end with 15 to 17 occupied blocs, and the runs that dominate
  the evidence leave two of their eighteen components empty, so
  $p(T = 16 \mid x)$ = 0.97: at $\zeta = 1$ the answer to “how many
  blocs?” is sixteen. (At $\zeta = 0.3$ it is five.)

Settings live in `bbmm_control()`; the defaults are the ones used here
(200 particles, conditional ESS 0.9, floor $\nu_k \ge 50$).
`bbmm_rejuvenate()` enlarges a finished run’s particle sample before
summarising it.

## Fixed-K MCMC

When the number of blocs is settled, `bbmm_mcmc()` samples the same
posterior at that *K* with ordinary Markov chains, each started from one
of the best EM fits. It is the quickest way to put uncertainty on one
partition.

``` r
# long: about 3 minutes for 4 chains x 3,000 iterations on 4 threads
mc <- bbmm_mcmc(maine, K = 16, zeta = 1, n_iter = 2000, n_burn = 1000, thin = 2,
                n_chains = 4, seed = 2026, n_threads = 4)
mcmc_diagnostics(mc)
plot_trace(mc)
```

<img src="TUTORIAL_files/figure-gfm/mcmc-plot-1.png" alt="" width="100%" />

Within a chain the traces look like stationary noise, but the chains do
not agree: the split-$\hat R$ of the log-likelihood is 2.6, far above
the 1.01 one would want.

``` r
round(res$mc_diag$mean_ll - max(res$mc_diag$mean_ll))   # mean log-likelihood per chain, vs the best
#>    1    2    3    4 
#>  -49 -150    0 -114
round(res$mc_diag$ari, 2)                               # agreement between final partitions
#>         chain_1 chain_2 chain_3 chain_4
#> chain_1    1.00    0.56    0.57    0.53
#> chain_2    0.56    1.00    0.52    0.61
#> chain_3    0.57    0.52    1.00    0.60
#> chain_4    0.53    0.61    0.60    1.00
round(res$mc_ari_ss, 2)                                 # each chain's final partition vs the SS one
#>    1    2    3    4 
#> 0.65 0.67 0.63 0.61
```

Each chain settles in a different partition, with mean log-likelihoods
up to 150 apart, and any two chains agree only moderately (adjusted Rand
0.52–0.61). They agree about as well with the stepping-stone partition
(0.61–0.67). With sixteen blocs, many of them neighbours that differ on
a handful of questions, the untempered posterior has many modes that
differ in where dozens of towns go. The allocation updates move one town
at a time, and on data this large the split–merge moves are almost never
accepted (`mc$acceptance` shows 0 for `split` and `merge`), so each
chain explores the neighbourhood of the EM optimum it started from. Use
`bbmm_mcmc()` to describe uncertainty around a mode, run several chains,
and always check `mcmc_diagnostics()`; finding the mode is the job of EM
with many starts or of the stepping-stone sampler, whose particles move
gradually from the prior. A smaller $\zeta$ also helps: at $\zeta = 0.3$
three of four chains agree.

# 5. Posterior summaries and how to read them

`bloc_summary()` turns any fit into the same tables. For the samplers it
first computes the **posterior similarity matrix** – for each pair of
towns, the posterior probability that they are in the same bloc – and
picks the sampled partition that agrees with it best (least-squares
clustering). Blocs are numbered by their average Yes share, highest
first.

``` r
s_ss <- bloc_summary(ev, maine)   # pools the runs at every K, weighted by p(K | x)
```

``` r
s_ss
#> Bloc summary: stepping-stone SMC pooled over K = 2-20 by p(K | x), zeta = 1 
#>  bloc units share   ballots mean_yes  nu certainty
#>     1    23 0.050 5,711,160    0.648  98     0.814
#>     2    33 0.071 4,069,547    0.605 127     0.906
#>     3    10 0.022 3,709,343    0.586 192     0.860
#>     4    36 0.078 3,034,031    0.575 146     0.904
#>     5    15 0.032   444,928    0.570 166     1.000
#>     6    44 0.095 6,481,813    0.553 193     0.979
#>     7    31 0.067 2,268,336    0.537 218     0.900
#>     8    19 0.041   851,197    0.527 223     0.884
#>     9    47 0.102 3,103,754    0.523 207     0.961
#>    10    12 0.026 1,075,121    0.519 199     0.999
#>    11    36 0.078 2,048,313    0.511 189     0.858
#>    12    25 0.054   484,801    0.510 121     0.969
#>    13    26 0.056 1,695,563    0.500 224     0.977
#>    14    35 0.076 1,272,453    0.476 286     0.929
#>    15    22 0.048   299,626    0.476 165     0.929
#>    16    48 0.104   729,976    0.448 122     0.939
#> Posterior on the number of occupied blocs: 16: 0.97, 17: 0.03
```

To see what a bloc is, start with its largest towns:

``` r
a <- s_ss$assignment
votes <- rowSums(maine$M, na.rm = TRUE)
for (b in s_ss$blocs$bloc) {
  members <- a$unit[a$bloc == b]
  top <- head(members[order(-votes[members])], 8)
  cat(b, ": ", paste(maine$unit_info$name[match(top, maine$units)], collapse = ", "), "\n", sep = "")
}
#> 1: Portland, South Portland, Brunswick, York, Cape Elizabeth, Yarmouth, Kittery, Orono
#> 2: Falmouth, Kennebunk, Waterville, Freeport, Cumberland, Bath, Belfast, Eliot
#> 3: Scarborough, Saco, Biddeford, Gorham, Westbrook, Wells, Topsham, Old Orchard Beach
#> 4: Bangor, Augusta, Alton + Edinburg + Old Town + Penobscot Nation Voting District, Farmington, Boothbay, Bowdoinham, Readfield, Bethel
#> 5: Madawaska, Fort Kent, Van Buren, Limestone, Eagle Lake + Winterville Plantation, Frenchville, Saint Agatha, Stockholm
#> 6: Lewiston, Auburn, Windham, Sanford, Standish, Gray, Buxton, Winslow
#> 7: Brewer, Hampden, Ellsworth, Waldoboro, Bucksport, Woolwich, Gouldsboro + Steuben, Union
#> 8: Presque Isle, Caribou, Houlton, Fort Fairfield, Mapleton, Washburn, Easton, Woodland
#> 9: Lisbon, Skowhegan, Waterboro, Fairfield, Lebanon, Durham, Jay, Weld + Wilton
#> 10: New Gloucester, Poland, Norway, Paris, Oxford, Mechanic Falls, Harrison, Otisfield
#> 11: Hermon, Millinocket, Glenburn, Orrington, Dover-Foxcroft, Winterport, Cherryfield + Franklin + Milbridge, Holden
#> 12: Calais, Cooper + East Machias, Baileyville, Addison, Jonesport, Machiasport, Perry, Princeton
#> 13: Turner, Sidney, Sabattus, Greene, West Gardiner, Litchfield, Bowdoin, Pittston
#> 14: Lincoln, Carmel, Levant, Eddington, Albion, Saint Albans, Corinna, Palmyra
#> 15: Mars Hill, Patten, Littleton, Island Falls, Linneus, Blaine, Monticello, Danforth
#> 16: Mattawamkeag + Medway + Sherman + Stacyville, Corinth, Wales, Cary Plantation + Hodgdon, Charleston, Greenbush, Atkinson + Sebec, Jackman
```

Read together with the map and the distinctive questions below, the
sixteen blocs fall into a few families (bloc numbers in brackets; the
descriptions are for this tutorial, not established names):

- **Greater Portland and the coast** \[1–3\]: Portland and the college
  and coastal towns (1, the most liberal bloc); the affluent coastal
  suburbs with Waterville, Bath and Belfast (2); and the southern
  suburban cities – Scarborough, Saco, Biddeford, Westbrook (3). All
  three are most distinctive on the gun-safety questions.
- **Service centres** \[4, 6, 7\]: Bangor, Augusta, Old Town and
  Farmington (4); Lewiston–Auburn, Windham and Sanford (6), set apart by
  the Oxford County casino votes; Brewer, Ellsworth and the midcoast
  (7).
- **Aroostook and Washington counties** \[5, 8, 12, 15\]: the St. John
  Valley (5, Madawaska, Fort Kent, Van Buren – the most certain bloc,
  with its own votes on same-sex marriage and the flag); central
  Aroostook (8, Presque Isle, Caribou, Houlton); Washington County (12,
  Calais and the Machias area); southern Aroostook (15, Mars Hill,
  Patten, Island Falls). The last two voted heavily to repeal
  school-district consolidation in 2009.
- **Small-town southern and central Maine** \[9–11, 13\]: including the
  Oxford Hills (10, Norway, Paris, Oxford, Poland – the casino’s host
  area, sharply against the racetrack and Lewiston casino questions) and
  the Bangor suburbs (11).
- **The rural interior** \[14, 16\]: small towns of Penobscot,
  Piscataquis and Somerset counties and beyond, most opposed to the
  gun-safety measures and most supportive of voter ID.

The columns of `s_ss$blocs`:

- `n_units`, `share`, `ballots`: how many towns, what fraction of towns,
  and how many two-option votes they cast in all. Blocs of small towns
  can hold many towns but few voters.
- `mean_yes`: the bloc’s Yes share averaged over questions – a rough
  left–right ordering, not a description.
- `nu`: the bloc’s concentration (posterior median). Low values mean a
  loose bloc whose towns differ a good deal; high values a tight one.
- `certainty`: the average, over the bloc’s towns, of their posterior
  probability of sharing a bloc with the bloc’s other towns. Near 1 is a
  solid bloc; around 0.5 a bloc whose edges are uncertain.

## Where the blocs are

``` r
plot_bloc_map(s_ss, maine_towns())
```

<img src="TUTORIAL_files/figure-gfm/map-1.png" alt="" width="100%" />

Blocs are defined only by how towns vote, not by where they are, so
geographic coherence is a finding, not an assumption. Town-level
certainty shows where the boundaries are firm:

``` r
plot_bloc_map(s_ss, maine_towns(), fill = "certainty")
```

<img src="TUTORIAL_files/figure-gfm/map-certainty-1.png" alt="" width="100%" />

`s_ss$assignment` has each town’s bloc and certainty; `s_ss$membership`
gives every town’s average similarity to each bloc, so for an uncertain
town you can see which bloc it is torn between:

``` r
a <- s_ss$assignment
uncertain <- head(a[order(a$certainty), ], 5)
cbind(uncertain, round(s_ss$membership[uncertain$unit, ], 2))
#>                          unit bloc certainty bloc_1 bloc_2 bloc_3 bloc_4 bloc_5
#> 329 PLEASANT_RIDGE_PLANTATION    8 0.2253447   0.00   0.00   0.00   0.00      0
#> 179                   HAMMOND   14 0.2617102   0.00   0.00   0.00   0.00      0
#> 96                    CLINTON   11 0.2726409   0.00   0.00   0.00   0.00      0
#> 286                NEW_SWEDEN    7 0.2813770   0.00   0.00   0.00   0.00      0
#> 297            NORTH_YARMOUTH    3 0.3189948   0.02   0.34   0.32   0.38      0
#>     bloc_6 bloc_7 bloc_8 bloc_9 bloc_10 bloc_11 bloc_12 bloc_13 bloc_14 bloc_15
#> 329      0   0.00   0.23   0.55       0    0.21    0.01    0.01    0.00    0.01
#> 179      0   0.00   0.00   0.00       0    0.00    0.00    0.00    0.26    0.71
#> 96       0   0.00   0.03   0.39       0    0.27    0.00    0.33    0.00    0.00
#> 286      0   0.28   0.67   0.00       0    0.01    0.00    0.00    0.00    0.00
#> 297      0   0.00   0.00   0.00       0    0.00    0.00    0.00    0.00    0.00
#>     bloc_16
#> 329    0.00
#> 179    0.01
#> 96     0.00
#> 286    0.00
#> 297    0.00
```

## What the blocs vote for

The heatmap shows each bloc’s Yes share on every question (from the
data, ballot-weighted). Read across a row to see how one question
divided the blocs, and down a column for a bloc’s profile.

``` r
plot_bloc_profiles(s_ss)
```

<img src="TUTORIAL_files/figure-gfm/profiles-1.png" alt="" width="100%" />

Showing each bloc’s gap from the state as a whole makes the distinctive
questions stand out:

``` r
plot_bloc_profiles(s_ss, value = "gap")
```

<img src="TUTORIAL_files/figure-gfm/gaps-1.png" alt="" width="100%" />

`distinctive_questions()` lists them:

``` r
distinctive_questions(s_ss, n = 3)
#>    bloc          question                                             label
#> 1     1 REFERENDUM_2016_3                2016 Q3: Firearm background checks
#> 2     1 REFERENDUM_2025_2           2025 Q2: Extreme risk protection orders
#> 3     1 REFERENDUM_2025_1 2025 Q1: Voter identification and absentee voting
#> 4     2 REFERENDUM_2016_3                2016 Q3: Firearm background checks
#> 5     2 REFERENDUM_2025_2           2025 Q2: Extreme risk protection orders
#> 6     2 REFERENDUM_2024_5                   2024 Q5: State flag replacement
#> 7     3 REFERENDUM_2016_3                2016 Q3: Firearm background checks
#> 8     3 REFERENDUM_2025_2           2025 Q2: Extreme risk protection orders
#> 9     3 REFERENDUM_2009_3     2009 Q3: School district consolidation repeal
#> 10    4 REFERENDUM_2010_1                     2010 Q1: Oxford County casino
#> 11    4 REFERENDUM_2011_2                  2011 Q2: Racetrack slot machines
#> 12    4 REFERENDUM_2015_3                      2015 Q3: Transportation bond
#> 13    5 REFERENDUM_2009_1                 2009 Q1: Same-sex marriage repeal
#> 14    5 REFERENDUM_2024_5                   2024 Q5: State flag replacement
#> 15    5 REFERENDUM_2012_1                        2012 Q1: Same-sex marriage
#> 16    6 REFERENDUM_2008_2                     2008 Q2: Oxford County casino
#> 17    6 REFERENDUM_2010_1                     2010 Q1: Oxford County casino
#> 18    6 REFERENDUM_2011_3                          2011 Q3: Lewiston casino
#> 19    7 REFERENDUM_2008_2                     2008 Q2: Oxford County casino
#> 20    7 REFERENDUM_2016_3                2016 Q3: Firearm background checks
#> 21    7 REFERENDUM_2010_1                     2010 Q1: Oxford County casino
#> 22    8 REFERENDUM_2009_1                 2009 Q1: Same-sex marriage repeal
#> 23    8 REFERENDUM_2024_5                   2024 Q5: State flag replacement
#> 24    8 REFERENDUM_2012_1                        2012 Q1: Same-sex marriage
#> 25    9 REFERENDUM_2016_3                2016 Q3: Firearm background checks
#> 26    9 REFERENDUM_2025_2           2025 Q2: Extreme risk protection orders
#> 27    9 REFERENDUM_2025_1 2025 Q1: Voter identification and absentee voting
#> 28   10 REFERENDUM_2011_2                  2011 Q2: Racetrack slot machines
#> 29   10 REFERENDUM_2011_3                          2011 Q3: Lewiston casino
#> 30   10 REFERENDUM_2008_2                     2008 Q2: Oxford County casino
#> 31   11 REFERENDUM_2016_3                2016 Q3: Firearm background checks
#> 32   11 REFERENDUM_2025_2           2025 Q2: Extreme risk protection orders
#> 33   11 REFERENDUM_2025_1 2025 Q1: Voter identification and absentee voting
#> 34   12 REFERENDUM_2009_3     2009 Q3: School district consolidation repeal
#> 35   12 REFERENDUM_2010_1                     2010 Q1: Oxford County casino
#> 36   12 REFERENDUM_2008_2                     2008 Q2: Oxford County casino
#> 37   13 REFERENDUM_2025_2           2025 Q2: Extreme risk protection orders
#> 38   13 REFERENDUM_2016_3                2016 Q3: Firearm background checks
#> 39   13 REFERENDUM_2025_1 2025 Q1: Voter identification and absentee voting
#> 40   14 REFERENDUM_2025_2           2025 Q2: Extreme risk protection orders
#> 41   14 REFERENDUM_2016_3                2016 Q3: Firearm background checks
#> 42   14 REFERENDUM_2025_1 2025 Q1: Voter identification and absentee voting
#> 43   15 REFERENDUM_2009_3     2009 Q3: School district consolidation repeal
#> 44   15 REFERENDUM_2025_1 2025 Q1: Voter identification and absentee voting
#> 45   15 REFERENDUM_2025_2           2025 Q2: Extreme risk protection orders
#> 46   16 REFERENDUM_2025_2           2025 Q2: Extreme risk protection orders
#> 47   16 REFERENDUM_2025_1 2025 Q1: Voter identification and absentee voting
#> 48   16 REFERENDUM_2016_3                2016 Q3: Firearm background checks
#>    yes_share all_share         gap
#> 1  0.7430548 0.4803974  0.26265733
#> 2  0.8423558 0.6277386  0.21461723
#> 3  0.1525831 0.3594147 -0.20683168
#> 4  0.6150749 0.4803974  0.13467748
#> 5  0.7518238 0.6277386  0.12408520
#> 6  0.5592854 0.4412227  0.11806265
#> 7  0.6128318 0.4803974  0.13243435
#> 8  0.7458267 0.6277386  0.11808809
#> 9  0.3089472 0.4167205 -0.10777324
#> 10 0.4555931 0.5044110 -0.04881785
#> 11 0.4031508 0.4503344 -0.04718365
#> 12 0.7726032 0.7273352  0.04526801
#> 13 0.7181237 0.5289603  0.18916335
#> 14 0.2544742 0.4412227 -0.18674849
#> 15 0.3438608 0.5260023 -0.18214155
#> 16 0.5509809 0.4618037  0.08917727
#> 17 0.5714054 0.5044110  0.06699445
#> 18 0.4293910 0.3654093  0.06398165
#> 19 0.3754697 0.4618037 -0.08633395
#> 20 0.3971874 0.4803974 -0.08321000
#> 21 0.4258294 0.5044110 -0.07858159
#> 22 0.7121462 0.5289603  0.18318587
#> 23 0.2649731 0.4412227 -0.17624963
#> 24 0.3541294 0.5260023 -0.17187294
#> 25 0.3241545 0.4803974 -0.15624297
#> 26 0.4853734 0.6277386 -0.14236516
#> 27 0.4765685 0.3594147  0.11715381
#> 28 0.2683579 0.4503344 -0.18197651
#> 29 0.1876873 0.3654093 -0.17772202
#> 30 0.6088978 0.4618037  0.14709412
#> 31 0.3197266 0.4803974 -0.16067083
#> 32 0.4719458 0.6277386 -0.15579283
#> 33 0.5032830 0.3594147  0.14386831
#> 34 0.6531267 0.4167205  0.23640622
#> 35 0.2777632 0.5044110 -0.22664778
#> 36 0.2353894 0.4618037 -0.22641421
#> 37 0.4420792 0.6277386 -0.18565940
#> 38 0.2964694 0.4803974 -0.18392798
#> 39 0.5412753 0.3594147  0.18186060
#> 40 0.3816254 0.6277386 -0.24611316
#> 41 0.2374016 0.4803974 -0.24299579
#> 42 0.5872312 0.3594147  0.22781642
#> 43 0.7201990 0.4167205  0.30347855
#> 44 0.6253019 0.3594147  0.26588720
#> 45 0.3640483 0.6277386 -0.26369027
#> 46 0.3061540 0.6277386 -0.32158459
#> 47 0.6508041 0.3594147  0.29138941
#> 48 0.2008245 0.4803974 -0.27957292
```

The bloc-level parameters from the model, with 90% intervals, are in
`s_ss$questions` (`p_mean`, `p_lo`, `p_hi`) next to the empirical shares
(`yes_share`). The two should agree closely; if they do not, the bloc is
not well described by a single profile.

``` r
head(s_ss$questions[, c("bloc", "label", "yes_share", "p_mean", "p_lo", "p_hi")])
#>   bloc                                         label yes_share    p_mean
#> 1    1         2008 Q1: Dirigo Health funding repeal 0.5565439 0.5318419
#> 2    1                 2008 Q2: Oxford County casino 0.4057344 0.3630177
#> 3    1            2008 Q3: Water infrastructure bond 0.6355820 0.6236447
#> 4    1             2009 Q1: Same-sex marriage repeal 0.3228870 0.3216191
#> 5    1               2009 Q2: Vehicle tax reductions 0.2783155 0.2698416
#> 6    1 2009 Q3: School district consolidation repeal 0.3413734 0.4034721
#>        p_lo      p_hi
#> 1 0.5086269 0.5542335
#> 2 0.3445964 0.3810912
#> 3 0.6028763 0.6419015
#> 4 0.3029876 0.3413916
#> 5 0.2493335 0.2883374
#> 6 0.3846210 0.4245340
```

**Naming blocs.** The model returns numbers, not names. Give a bloc a
name only after looking at its distinctive questions, its map, and
something the votes do not contain – demographics, partisanship, local
knowledge. Two kinds of town can vote alike for different reasons and
still share a bloc.

## Questions the blocs do not explain

`question_fit()` reports, for each question, how much of the
town-to-town variation in Yes share the blocs account for
(ballot-weighted $R^2$) and how far apart the blocs are. Questions near
the top divided towns along lines the blocs do not follow – often local
interests such as a casino site, a transmission line or a school
consolidation – and are worth reading on their own.

``` r
qf <- question_fit(s_ss, maine)
head(qf, 6)
#>            question                                         label all_share
#> 1 REFERENDUM_2009_2               2009 Q2: Vehicle tax reductions 0.2588771
#> 2 REFERENDUM_2021_1   2021 Q1: Electric transmission restrictions 0.5917978
#> 3 REFERENDUM_2018_1                  2018 Q1: Universal home care 0.3700421
#> 4 REFERENDUM_2023_2 2023 Q2: Foreign government campaign spending 0.8632294
#> 5 REFERENDUM_2021_3         2021 Q3: Constitutional right to food 0.6081520
#> 6 REFERENDUM_2023_1        2023 Q1: Voter approval of public debt 0.6548835
#>   spread_data spread_model        r2
#> 1  0.06504863   0.07226177 0.3146891
#> 2  0.23277121   0.23201510 0.3560318
#> 3  0.17408974   0.17682783 0.3683838
#> 4  0.06351992   0.05923398 0.3805039
#> 5  0.17174391   0.17966253 0.4039053
#> 6  0.10859217   0.10041875 0.4367419
tail(qf, 3)
#>             question                                             label
#> 70 REFERENDUM_2025_1 2025 Q1: Voter identification and absentee voting
#> 71 REFERENDUM_2018_3         2018 Q3: Transportation and culverts bond
#> 72 REFERENDUM_2016_6                      2016 Q6: Transportation bond
#>    all_share spread_data spread_model        r2
#> 70 0.3594147   0.4982211    0.4799102 0.9163967
#> 71 0.6763771   0.2750924    0.2864465 0.9175285
#> 72 0.6111483   0.3182071    0.3297193 0.9236774
```

## The shape of the blocs: Jensen–Shannon principal coordinates

`plot_js_pca()` places the towns in a space built only from their votes:
the distance between two towns sums, over questions, the Jensen–Shannon
divergence between their Yes shares; principal coordinates then find the
main axes. If the blocs were slices of a single left–right dimension,
they would line up along PC1 and PC1 would carry nearly everything.
Blocs that separate on PC2 or PC3 differ in ways one ideological score
misses.

``` r
plot_js_pca(maine, s_ss)
```

<img src="TUTORIAL_files/figure-gfm/js-pca-1.png" alt="" width="100%" />

Each facet title gives the share of total variation that axis carries
and the share of that axis explained by bloc membership (bloc $R^2$).

## How settled is the partition?

The posterior similarity matrix, with towns ordered by bloc, shows it
directly: dark squares on the diagonal are blocs whose towns always go
together; grey bands off the diagonal are towns, or whole blocs, that
the posterior sometimes merges.

``` r
plot_coclustering(s_ss)
```

<img src="TUTORIAL_files/figure-gfm/psm-1.png" alt="" width="100%" />

## Comparing the methods

The partitions from EM, the stepping-stone sampler and MCMC can be
compared with the adjusted Rand index (1 = identical, 0 = no more
agreement than chance):

``` r
s_mc <- res$s_mc
c(EM_vs_SS = adjusted_rand(s_em$assignment$bloc, s_ss$assignment$bloc),
  MCMC_vs_SS = adjusted_rand(s_mc$assignment$bloc, s_ss$assignment$bloc),
  EM_vs_MCMC = adjusted_rand(s_em$assignment$bloc, s_mc$assignment$bloc))
#>   EM_vs_SS MCMC_vs_SS EM_vs_MCMC 
#>  0.3471003  0.6578153  0.3367368
table(SS = s_ss$assignment$bloc, MCMC = s_mc$assignment$bloc)
#>     MCMC
#> SS    1  2  3  4  5  6  7  8  9 10 11 12 13 14 15 16
#>   1  14  9  0  0  0  0  0  0  0  0  0  0  0  0  0  0
#>   2   0 14 17  2  0  0  0  0  0  0  0  0  0  0  0  0
#>   3   0  0  0 10  0  0  0  0  0  0  0  0  0  0  0  0
#>   4   0  0  6  0 30  0  0  0  0  0  0  0  0  0  0  0
#>   5   0  0  0  0  0 15  0  0  0  0  0  0  0  0  0  0
#>   6   0  0  0  0  0  0 34  0 10  0  0  0  0  0  0  0
#>   7   0  0  0  0  3  0  0 27  0  1  0  0  0  0  0  0
#>   8   0  0  0  0  0  0  0  0  0 19  0  0  0  0  0  0
#>   9   0  0  0  0  0  0  0  0 34  0  0  1 12  0  0  0
#>   10  0  0  0  0  0  0  0  0  0  0 12  0  0  0  0  0
#>   11  0  0  0  0  0  0  0  3  1 10  0 19  3  0  0  0
#>   12  0  0  0  0  0  0  0  0  0  0  0 19  0  0  6  0
#>   13  0  0  0  0  0  0  0  0  0  0  0  0 26  0  0  0
#>   14  0  0  0  0  0  0  0  0  0  0  0  2  0 33  0  0
#>   15  0  0  0  0  0  0  0  0  0  1  0  0  0  1 20  0
#>   16  0  0  0  0  0  0  0  0  0  0  0  0  0  1  6 41
```

The MCMC representative partition agrees only moderately with the
stepping-stone one – not surprising given the mixing problems above. EM
and the stepping-stone sampler agree even less by the adjusted Rand
index, but for a different reason, which the cross-tabulation shows:

``` r
table(SS = s_ss$assignment$bloc, EM = s_em$assignment$bloc)
#>     EM
#> SS    1  2  3  4  5  6
#>   1  23  0  0  0  0  0
#>   2  16 17  0  0  0  0
#>   3   0 10  0  0  0  0
#>   4   0 33  3  0  0  0
#>   5   0  0  0 15  0  0
#>   6   0  0 44  0  0  0
#>   7   0  2 28  1  0  0
#>   8   0  0  0 16  3  0
#>   9   0  0 11  0 36  0
#>   10  0  0  8  0  4  0
#>   11  0  0  3  2 31  0
#>   12  0  0  0  6 14  5
#>   13  0  0  0  0 26  0
#>   14  0  0  0  0 15 20
#>   15  0  0  0  0  2 20
#>   16  0  0  0  0  0 48
```

Most of the sixteen stepping-stone blocs sit largely inside one of EM’s
six – nine of them have 90% or more of their towns in a single EM bloc –
and the rest straddle two neighbouring EM blocs (SS 14, for example,
splits between EM 5 and 6). EM bloc 1 is SS 1 plus part of SS 2; EM bloc
6 is SS 16 plus most of SS 14 and 15. Largely, the two are views of the
same structure at different grains. The adjusted Rand index penalises
the difference in grain as heavily as genuine disagreement, so read the
cross-tabulation, not only the index.

# 6. Where to go next

- Fit your own data: build a `bbmm_data` object (section 1), run
  `bbmm_em_select()` over a broad range of *K*, then the stepping-stone
  sampler over a narrower range around the EM answer.
- Try neighbouring values of `zeta` and see which blocs persist.
- Every summary is a plain data frame, and every figure a ggplot object,
  so both can be extended with your own covariates and styling.

## References

Del Moral, P., Doucet, A. and Jasra, A. (2006). Sequential Monte Carlo
samplers. *Journal of the Royal Statistical Society B* 68, 411–436.

Dahl, D. B. (2006). Model-based clustering for expression data via a
Dirichlet process mixture model. In *Bayesian Inference for Gene
Expression and Proteomics*, 201–218.

Jain, S. and Neal, R. M. (2004). A split-merge Markov chain Monte Carlo
procedure for the Dirichlet process mixture model. *Journal of
Computational and Graphical Statistics* 13, 158–182.

Miller, J. W. and Dunson, D. B. (2019). Robust Bayesian inference via
coarsening. *Journal of the American Statistical Association* 114,
1113–1125.

O’Brien, J. (2026). A Bayesian mixture model captures temporal and
spatial structure of voting blocs within longitudinal referendum data.
*Political Science Research and Methods*, 1–22.
<doi:10.1017/psrm.2025.10050>.

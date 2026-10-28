## Submission

First submission of Elections.BBMM.

## Test environments

* local macOS (aarch64), R 4.6.1
* GitHub Actions: macOS, Windows and Ubuntu (R release, devel and oldrel-1)

## R CMD check results

0 errors | 0 warnings | 1 note

* New submission.

## Notes for CRAN

* The C++ code uses `std::thread`; examples and tests use at most 2 threads, and every
  user-facing function defaults to `n_threads = 1` / `cores = 1`.
* The samplers use their own seeded random streams (`std::mt19937_64`); the seed is taken from
  R's random-number generator when the user does not supply one, so `set.seed()` makes runs
  reproducible. Functions do not change the user's `.Random.seed` when a `seed` is supplied.
* The long-running analyses in the vignette are precomputed (`inst/extdata/maine_tutorial.rds`);
  the vignette builds in seconds.

# data-raw/render_tutorial.R -- renders the vignette to TUTORIAL.md (+ TUTORIAL_files/) at the repository
# root, for reading on GitHub. Needs the package installed and pandoc on the PATH.
#   Rscript data-raw/render_tutorial.R      (from the package root)
Sys.setenv(BBMM_TUTORIAL_DPI = "120")
file.copy("vignettes/Elections.BBMM.Rmd", "TUTORIAL.Rmd", overwrite = TRUE)
on.exit(unlink("TUTORIAL.Rmd"))
unlink("TUTORIAL_files", recursive = TRUE)
rmarkdown::render("TUTORIAL.Rmd", output_format = rmarkdown::github_document(html_preview = FALSE, toc = TRUE,
                  toc_depth = 2, fig_width = 8, dev = "png"), output_file = "TUTORIAL.md", quiet = TRUE,
                  envir = new.env())
unlink("TUTORIAL.Rmd")

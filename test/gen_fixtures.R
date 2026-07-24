# Regenerate golden fixtures. Run: Rscript test/gen_fixtures.R
suppressMessages({library(ggdist); library(posterior)})
set.seed(20260724)
dir.create("test/fixtures", showWarnings = FALSE, recursive = TRUE)

emit_samples <- function(name, draws) {
  write.csv(data.frame(draw = draws),
            file.path("test/fixtures", paste0(name, "_draws.csv")), row.names = FALSE)
}
emit_pi <- function(name, draws, width, interval_fn, interval_name) {
  pi <- interval_fn(draws, .width = width)   # median_qi / median_hdci / mode_hdi
  out <- data.frame(value = pi$y, lower = pi$ymin, upper = pi$ymax,
                    width = width, interval = interval_name)
  write.csv(out, file.path("test/fixtures", paste0(name, "_", interval_name, "_", width, ".csv")),
            row.names = FALSE)
}

cases <- list(
  normal   = rnorm(4000, 0, 1),
  beta     = rbeta(4000, 2, 8),
  bimodal  = c(rnorm(2500, -3, 0.4), rnorm(2500, 3, 0.4)),
  studentt = rt(4000, df = 3)
)

for (nm in names(cases)) {
  d <- cases[[nm]]
  emit_samples(nm, d)
  for (w in c(0.66, 0.95)) {
    emit_pi(nm, d, w, median_qi,   "qi")
    emit_pi(nm, d, w, median_hdci, "hdci")
  }
}

writeLines(capture.output(sessionInfo()), "test/fixtures/sessionInfo.txt")
cat("fixtures written\n")

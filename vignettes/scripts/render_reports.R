#!/usr/bin/env Rscript
# Render publication Quarto reports (reduced mode by default).

args <- commandArgs(trailingOnly = TRUE)
mode <- if ("--mode" %in% args) {
  args[which(args == "--mode")[1] + 1]
} else {
  "reduced"
}

repo <- Sys.getenv("IMMUNE_HEALTHY_BODY_MAP_ROOT", unset = "")
if (!nzchar(repo)) {
  repo <- normalizePath(file.path("..", ".."), mustWork = FALSE)
}
Sys.setenv(IMMUNE_HEALTHY_BODY_MAP_ROOT = repo)
Sys.setenv(
  IMMUNE_HEALTHY_BODY_MAP_DATA = Sys.getenv(
    "IMMUNE_HEALTHY_BODY_MAP_DATA",
    unset = file.path(repo, "vignettes", "data")
  )
)

reports_dir <- file.path(repo, "vignettes", "reports")
reports <- c(
  "Fig1_umap.qmd",
  "Fig2.qmd",
  "Fig4.qmd",
  "Fig5.qmd",
  "Fig2_supplementary_plots.qmd",
  "Fig4_supplementary_plots.qmd",
  "Fig5_supplementary_plots.qmd"
)

message("Mode: ", mode)
message("Repo: ", repo)
message("Data: ", Sys.getenv("IMMUNE_HEALTHY_BODY_MAP_DATA"))

if (!requireNamespace("quarto", quietly = TRUE)) {
  stop(
    "Install the quarto R package and Quarto CLI to render reports.",
    call. = FALSE
  )
}

results <- lapply(reports, function(f) {
  path <- file.path(reports_dir, f)
  message("\n=== Rendering ", f, " ===")
  ok <- FALSE
  err <- NULL
  tryCatch({
    quarto::quarto_render(path, as_job = FALSE)
    ok <<- TRUE
  }, error = function(e) {
    err <<- conditionMessage(e)
    message("FAILED: ", err)
  })
  data.frame(report = f, ok = ok, error = if (is.null(err)) "" else err)
})

tab <- do.call(rbind, results)
print(tab)
out <- file.path(repo, "vignettes", "manifests", "render_results.csv")
utils::write.csv(tab, out, row.names = FALSE)
message("Wrote ", out)

if (mode == "reduced") {
  message(
    "\nReduced mode: failures due to missing Zenodo/large files are expected\n",
    "until caches cover every chunk. Supplementary plot qmds should succeed\n",
    "if source_tables Excel files are present.\n"
  )
}

quit(save = "no", status = if (all(tab$ok)) 0L else 1L)

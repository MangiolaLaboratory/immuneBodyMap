#!/usr/bin/env Rscript
# Download published Zenodo inputs into IMMUNE_HEALTHY_BODY_MAP_DATA.
#
# Pseudobulk (available now):
#   DOI 10.5281/zenodo.15798373
#   https://zenodo.org/records/15798373
#   file: pseudobulk_se.h5ad
#   md5:  88c71c0fd1d6ce2fe15eccdd7b36110f
#
# sccomp L3/L0 estimates:
#   DOI 10.5281/zenodo.21389126
#   https://doi.org/10.5281/zenodo.21389126
# Download the model record and set SCCOMP_ESTIMATES_DIR to its RDS directory.

args <- commandArgs(trailingOnly = TRUE)
download_pseudobulk <- "--pseudobulk" %in% args || "--all" %in% args

data_root <- Sys.getenv(
  "IMMUNE_HEALTHY_BODY_MAP_DATA",
  unset = file.path(
    Sys.getenv(
      "IMMUNE_HEALTHY_BODY_MAP_ROOT",
      unset = normalizePath(file.path("..", ".."), mustWork = FALSE)
    ),
    "vignettes",
    "data"
  )
)

dest_pb <- file.path(data_root, "zenodo_release", "pseudobulk")
dir.create(dest_pb, recursive = TRUE, showWarnings = FALSE)
out_file <- file.path(dest_pb, "pseudobulk_se.h5ad")

pseudobulk_url <- paste0(
  "https://zenodo.org/records/15798373/files/pseudobulk_se.h5ad?download=1"
)
expected_md5 <- "88c71c0fd1d6ce2fe15eccdd7b36110f"

message(
  "download_publication_release.R\n",
  "Destination root: ", data_root, "\n"
)

if (!download_pseudobulk) {
  message(
    "Pseudobulk Zenodo record: https://doi.org/10.5281/zenodo.15798373\n",
    "  Expected path: ", out_file, "\n",
    "  Expected md5:  ", expected_md5, "\n",
    "\n",
    "Re-run with --pseudobulk to download (~10 GB), or set:\n",
    "  export PSEUDOBULK_H5AD=/path/to/pseudobulk_se.h5ad\n",
    "\n",
    "sccomp L3/L0 estimates: https://doi.org/10.5281/zenodo.21389126\n",
    "  Download the model record and set SCCOMP_ESTIMATES_DIR.\n",
    "See vignettes/manifests/README.md.\n"
  )
  quit(save = "no", status = 0L)
}

if (file.exists(out_file)) {
  message("Already present: ", out_file)
  message("Verify with: md5sum ", out_file, "  # expect ", expected_md5)
  quit(save = "no", status = 0L)
}

message("Downloading pseudobulk_se.h5ad (~10 GB) from Zenodo…")
utils::download.file(pseudobulk_url, destfile = out_file, mode = "wb", quiet = FALSE)
message("Saved: ", out_file)
message("Verify: md5sum should equal ", expected_md5)
message("Then: export PSEUDOBULK_H5AD=", out_file)

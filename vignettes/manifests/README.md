# External inputs for reproduction

Small figure panels, SI workbooks, and code live on GitHub under `vignettes/`.
Large analysis objects are **not** in git. Use this folder to see what else is required.

These resources support *A multi-tissue immune map across 4,240 single-cell
donors resolves asynchronous immune ageing* (target journal: *Nature Aging*).

| Resource | How to obtain | Env override |
|----------|---------------|--------------|
| `pseudobulk_se.h5ad` | [Zenodo 10.5281/zenodo.15798373](https://doi.org/10.5281/zenodo.15798373) (md5 `88c71c0fd1d6ce2fe15eccdd7b36110f`, ~10 GB) | `PSEUDOBULK_H5AD` |
| sccomp L3 + L0 estimate RDS | [Zenodo 10.5281/zenodo.21389127](https://doi.org/10.5281/zenodo.21389127) | `SCCOMP_ESTIMATES_DIR` |
| cellNexus metadata 1.0.12 | Public URL via `get_publication_metadata()` in `R/paths.R` | optional `CELLNEXUS_METADATA` |
| `NEW_CELL_TYPE/_targets/` | **Not deposited** (exemption). Rebuild with `targets_pipeline/` from pseudobulk, or use shipped panel RDS for reduced mode | `NEW_CELL_TYPE_ROOT` |

The plasma-niche input workflow is documented in
[`../reports/plasma_niche_pseudobulk_qc.qmd`](../reports/plasma_niche_pseudobulk_qc.qmd).
It retrieves the selected cellNexus cells and writes six local, gitignored
pseudobulk objects plus QC summaries. These derived objects support plasma niche
transcriptomics and pathway analyses and can be regenerated from the report.

See also `external_large_files.csv` (deposit decisions) and `dataset_manifest.csv`.
Figure → data map: `figure_data_map.csv`.

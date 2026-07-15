# External inputs for reproduction

Small figure panels, SI workbooks, and code live on GitHub under `vignettes/`.
Large analysis objects are **not** in git. Use this folder to see what else is required.

| Resource | How to obtain | Env override |
|----------|---------------|--------------|
| `pseudobulk_se.h5ad` | [Zenodo 10.5281/zenodo.15798373](https://doi.org/10.5281/zenodo.15798373) (md5 `88c71c0fd1d6ce2fe15eccdd7b36110f`, ~10 GB) | `PSEUDOBULK_H5AD` |
| sccomp L3 + L0 estimate RDS | Companion deposit still pending (`<ZENODO_DOI_PENDING>`) | `SCCOMP_ESTIMATES_DIR` |
| cellNexus metadata 1.0.12 | Public URL via `get_publication_metadata()` in `R/paths.R` | optional `CELLNEXUS_METADATA` |
| `NEW_CELL_TYPE/_targets/` | **Not deposited** (exemption). Rebuild with `targets_pipeline/` from pseudobulk, or use shipped panel RDS for reduced mode | `NEW_CELL_TYPE_ROOT` |

See also `external_large_files.csv` (deposit decisions) and `dataset_manifest.csv`.
Figure → data map: `figure_data_map.csv`.

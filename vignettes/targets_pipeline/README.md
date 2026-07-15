# NEW_CELL_TYPE targets pipeline (publication copy)

Portable copy of `HPC_posterior/TAR_SCRIPTS` for the Nature Aging package.
The full `_targets/` object stores are **not** deposited (exemption); this code
reproduces them from the Zenodo primary input `pseudobulk_se.h5ad`.

## Layout

```text
targets_pipeline/
├── batch_run.R              # loop over cell types → V1_<ct>/_targets
├── dynamic_tar_script.R     # writes per–cell-type _targets.R + tar graph
├── functions/
│   ├── edit_covariates.R
│   ├── remove_unwanted_effect.R
│   ├── get_adjusted_matrix.R
│   ├── fit_to_age_monotonic_changes.R
│   ├── fit_to_sex_differ_at_decades.R
│   └── get_brms_tar.R       # read one gene fit from an existing store
└── README.md
```

## Inputs / outputs / intermediates

| Stage | Target name | Role | Size / notes |
|-------|-------------|------|--------------|
| Input | — | `PSEUDOBULK_H5AD` (`pseudobulk_se.h5ad`) | [Zenodo 10.5281/zenodo.15798373](https://doi.org/10.5281/zenodo.15798373); md5 `88c71c0f…` |
| Aux | — | `disease_data_grouped_further.csv` | shipped in `vignettes/data/source_tables/` |
| Aux | — | `edit_covariates.R` | this tree or `rebuttal_CellPress/` |
| Load | `pseudobulk_sample` | Filtered / scaled SE for one cell type | large HDF5-backed SE |
| Side | `reference_sample_*.rds` | Written under `_targets/` store | small |
| Branch | `feature_df` → `estimates_chunk` | Per-gene `brms` ZINB fits | dominant cost; qs objects |
| Branch | `hypothesis_age_monotonic_and_adjust_tissue` | Age contrasts + tissue-adjusted estimates | Fig1 tissue UMAP input |
| Branch | `hypothesis_sex_differ_at_decades` | Sex×decade hypotheses | Fig5-related |
| Branch | `adjust_age` | Age-adjusted expression estimates | Fig1 age UMAP input |

Downstream (outside this pipeline, already in vignettes):

- Fig1 `compute-adjusted-counts` → `adjust_*_list.rds` via `tar_read_raw` of `adjust_age` / tissue adjust columns
- `scripts/rebuild_fig4_gene_predictions.R` → `pred_fixed.rds` / `pred_total.rds` from `estimates_chunk` brms objects

## Configure

```bash
export IMMUNE_HEALTHY_BODY_MAP_ROOT=/path/to/immuneHealthyBodyMap
export PSEUDOBULK_H5AD=/path/to/pseudobulk_se.h5ad
export NEW_CELL_TYPE_ROOT=/path/to/NEW_CELL_TYPE          # created if missing
# optional
export TARGETS_USE_CREW=true                              # false = no Slurm crew
export TARGETS_SLURM_ACCOUNT=saigencir003
export TARGETS_MAKE_NAMES=adjust_age                      # or other target names
export DISEASE_GROUPING_CSV=$IMMUNE_HEALTHY_BODY_MAP_ROOT/vignettes/data/source_tables/disease_data_grouped_further.csv
```

## Run (HPC)

```bash
Rscript $IMMUNE_HEALTHY_BODY_MAP_ROOT/vignettes/targets_pipeline/batch_run.R
```

Each cell type gets `NEW_CELL_TYPE_ROOT/V1_<make.names(ct)>/_targets/`.

## Rebuild Fig4 predictions from an existing store

```bash
export NEW_CELL_TYPE_ROOT=/path/to/NEW_CELL_TYPE
export NEW_PRED_DIR=$IMMUNE_HEALTHY_BODY_MAP_ROOT/vignettes/data/intermediate/fig4
Rscript vignettes/scripts/rebuild_fig4_gene_predictions.R
```

Writes `pred_fixed.rds` / `pred_total.rds` from shipped `new_gene_hyp_df.rds` +
brms objects in the store. Details: `scripts/rebuild_fig4_gene_predictions.md`.

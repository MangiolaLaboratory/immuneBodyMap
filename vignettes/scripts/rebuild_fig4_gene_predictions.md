# Rebuild Fig4 gene-expression ageing predictions

## Core script

```bash
export IMMUNE_HEALTHY_BODY_MAP_ROOT=/path/to/immuneHealthyBodyMap
export NEW_CELL_TYPE_ROOT=/path/to/NEW_CELL_TYPE   # targets stores (exemption)
export NEW_PRED_DIR=$IMMUNE_HEALTHY_BODY_MAP_ROOT/vignettes/data/intermediate/fig4
# optional: export FIG4_PRED_WORKERS=8

Rscript vignettes/scripts/rebuild_fig4_gene_predictions.R
```

Produces:

- `pred_fixed.rds`
- `pred_total.rds`

from `new_gene_hyp_df.rds` (GitHub) + per-gene brms fits in `NEW_CELL_TYPE`.

The fitting pipeline itself is vendored at
`vignettes/targets_pipeline/` (from `HPC_posterior/TAR_SCRIPTS`).

This is the portable form of `Age_Clock/new_pred/pred_fixed_total.R`.

## Tissue extensions (stomach / missing)

These remain analysis-archive scripts (still require HPC inputs):

| Output | Historical script under `AGE_CLOCK_ARCHIVE` |
|--------|-----------------------------------------------|
| `pred_total_stomach.rds` | `TAR_stomach/stomach_pseudobulk_data.R` |
| `pred_total_missing.rds` | `TAR_missing_tissue/missing_tissue_data.R` |

Point `Fig4.qmd` at them with `NEW_PRED_DIR` after they are rebuilt, or keep
the small shipped copies under `vignettes/data/intermediate/fig4/` if already
present for reduced mode.

## Relation to primary inputs

Gene predictions do **not** come from sccomp. They come from the **targets**
gene-expression fits (same store family as Fig1 `adjust_*_list`). With a
targets-cache exemption, deposit only:

1. `pseudobulk_se.h5ad`
2. sccomp L3 / L0 estimates  
3. (exempt) `NEW_CELL_TYPE/_targets/` + this rebuild script

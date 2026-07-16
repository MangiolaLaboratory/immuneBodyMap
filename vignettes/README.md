# Figure reproducibility package

Code and small artefacts supporting the manuscript:

> *A multi-tissue immune map across 4,240 single-cell donors resolves
> asynchronous immune ageing*

**Target journal:** *Nature Aging*

The `Age_Clock` label remains in a few filenames and environment variables as a
legacy analysis-tree name; it is not the manuscript title.

**Reduced mode** (default): render Quarto reports using the shipped
panel RDS and SI workbooks (~60 MB under `data/` + Extended Data plot objects).

**Full recompute** (HPC): start from public checkpoints + optional targets rebuild.

- Pseudobulk: [10.5281/zenodo.15798373](https://doi.org/10.5281/zenodo.15798373)
- sccomp L3/L0 model record: **`<ZENODO_DOI_PENDING>`**

---

## Contents

```text
vignettes/
├── README.md
├── reports/                  # Quarto figure reports (canonical)
├── R/                        # paths.R + SI helpers
├── config/paths.example.env
├── data/
│   ├── source_tables/        # SI Excel + small auxiliaries
│   ├── figure_ready/         # panel ggplot RDS
│   └── intermediate/         # compact caches
├── supplementary_plot_objects/
├── targets_pipeline/         # NEW_CELL_TYPE gene-fit pipeline code
├── scripts/                  # render / rebuild helpers
└── manifests/                # external inputs + figure→data map
```

---

## Quick start (reduced reproduction)

```bash
cd /path/to/immuneHealthyBodyMap
export IMMUNE_HEALTHY_BODY_MAP_ROOT=$PWD
export IMMUNE_HEALTHY_BODY_MAP_DATA=$PWD/vignettes/data

Rscript vignettes/scripts/render_reports.R --mode reduced
```

Copy `config/paths.example.env` to a gitignored `paths.local.env` / `.Renviron`
only if you need overrides.

Reduced mode uses:

1. `data/figure_ready/**`
2. `data/intermediate/**`
3. `data/source_tables/*.xlsx`
4. `supplementary_plot_objects/*.rds`

---

## External inputs (not in git)

| Input | Role |
|-------|------|
| [Zenodo `pseudobulk_se.h5ad`](https://doi.org/10.5281/zenodo.15798373) | Full recompute / targets (md5 `88c71c0f…`) |
| Zenodo sccomp L3 + L0 RDS | Full recompute of composition panels; companion record pending (`<ZENODO_DOI_PENDING>`) |
| cellNexus metadata (cloud) | Auto via `get_publication_metadata()` |
| `NEW_CELL_TYPE/_targets/` | **Exemption** — code in `targets_pipeline/`; reduced mode does not need the store |

Details: `manifests/README.md` and `manifests/external_large_files.csv`.

Download pseudobulk into `data/zenodo_release/pseudobulk/` or set `PSEUDOBULK_H5AD`.
After sccomp estimates are deposited, set `SCCOMP_ESTIMATES_DIR` similarly.

---

## Full HPC recomputation

1. Set env vars from `config/paths.example.env`
2. Optional gene pipeline: `Rscript vignettes/targets_pipeline/batch_run.R`  
   (see `targets_pipeline/README.md`)
3. Gene predictions from an existing store:  
   `Rscript vignettes/scripts/rebuild_fig4_gene_predictions.R`
4. Re-render reports after deleting selected caches under `data/intermediate/`

Requires a recent R stack (sccomp, cellNexus, brms/cmdstanr, targets, tidybulk, …)
and substantial HPC resources for sccomp fits and per-gene brms targets.

---

## Path policy

Manuscript code resolves paths via `R/paths.R` and environment variables.
There are no hard-coded `/hpcfs` paths in `reports/` or `targets_pipeline/`.

Missing large files raise an error that names the file, override variable and,
where applicable, the pending Zenodo record.

---

## Reports

| Report | Reduced readiness |
|--------|-------------------|
| `Fig1_umap.qmd` | High (shipped UMAP RDS) |
| `Fig2.qmd` | Partial (panels + SI; volcano/heatmap rebuild need L3) |
| `Fig4.qmd` | Partial (panels + intermediates; gene preds optional rebuild) |
| `Fig5.qmd` | High for shipped panels |
| `Fig*_supplementary_plots.qmd` | High (SI xlsx + Extended Data RDS) |

Manuscript numbering: Fig1, Fig2, Fig4, Fig5 (some cache files still use
legacy `Fig3_` / `Fig4_` prefixes).

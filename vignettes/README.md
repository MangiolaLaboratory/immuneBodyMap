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
- sccomp L3/L0 model record: [10.5281/zenodo.21389126](https://doi.org/10.5281/zenodo.21389126)

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
│   ├── intermediate/         # compact caches
│   └── processed/            # generated local data (gitignored)
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
| [Zenodo `pseudobulk_se.h5ad`](https://doi.org/10.5281/zenodo.15798373) | Full recompute / targets (md5 `88c71c0fd1d6ce2fe15eccdd7b36110f`) |
| [Zenodo sccomp L3 + L0 RDS](https://doi.org/10.5281/zenodo.21389126) | Full recompute of composition panels |
| cellNexus metadata (cloud) | Auto via `get_publication_metadata()` |
| cellNexus single-cell matrices (cloud/cache) | Plasma-niche pseudobulk generation via `get_single_cell_experiment()` |
| `NEW_CELL_TYPE/_targets/` | **Exemption** — code in `targets_pipeline/`; reduced mode does not need the store |

Details: `manifests/README.md` and `manifests/external_large_files.csv`.

Download pseudobulk into `data/zenodo_release/pseudobulk/` or set `PSEUDOBULK_H5AD`.
Download the sccomp estimates from their Zenodo record and set
`SCCOMP_ESTIMATES_DIR` to the directory containing the L3 and L0 RDS files.

---

## Full HPC recomputation

1. Set env vars from `config/paths.example.env`
2. Plasma-niche inputs:

   `quarto render vignettes/reports/plasma_niche_pseudobulk_qc.qmd`
3. Optional gene pipeline: `Rscript vignettes/targets_pipeline/batch_run.R`

   (see `targets_pipeline/README.md`)
4. Gene predictions from an existing store:

   `Rscript vignettes/scripts/rebuild_fig4_gene_predictions.R`
5. Re-render reports after deleting selected caches under `data/intermediate/`

Requires a recent R stack (sccomp, cellNexus, Seurat, harmony, edgeR, scuttle,
brms/cmdstanr, targets, tidybulk, …) and substantial HPC resources for sccomp
fits, cellNexus single-cell retrieval and per-gene brms targets.

### Plasma niche pseudobulk inputs

[`reports/plasma_niche_pseudobulk_qc.qmd`](reports/plasma_niche_pseudobulk_qc.qmd)
merges the former respiratory/multi-tissue and large-intestine preparation
reports into one parameterized workflow. By default it creates six separate
quality-controlled pseudobulk `SummarizedExperiment` objects:

- respiratory system: plasma, stromal and epithelial;
- large intestine: plasma, stromal and epithelial.

Before pseudobulk aggregation, the report performs Seurat normalization and PCA,
Harmony integration by sample, and single-cell UMAP visualisation by compartment,
tissue and sample. It saves the integrated Seurat object. Interactive Plotly
panels show sample retention, effective library size and sparsity. Set
`CELLNEXUS_CACHE_DIR` for the downloaded cells and optionally set
`PLASMA_NICHE_OUTPUT_DIR`; otherwise outputs go to the gitignored
`data/processed/plasma_niche/` directory. The pseudobulk outputs are inputs for
**plasma niche transcriptomics and pathway analyses**.

---

## Path policy

Manuscript code resolves data paths via `R/paths.R` and environment variables.
Environment-specific R library paths, when needed on HPC, should be configured
locally and are not input-data dependencies.

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
| `plasma_niche_pseudobulk_qc.qmd` | Full-data workflow; requires cellNexus single-cell retrieval |

Manuscript numbering: Fig1, Fig2, Fig4, Fig5 (some cache files still use
legacy `Fig3_` / `Fig4_` prefixes).

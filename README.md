# A multi-tissue immune map across 4,240 single-cell donors resolves asynchronous immune ageing

This repository contains analysis code, figure reports and reproducibility
materials for the manuscript:

> *A multi-tissue immune map across 4,240 single-cell donors resolves
> asynchronous immune ageing*

**Target journal:** *Nature Aging*

The manuscript-facing package focuses on multi-tissue immune composition,
gene-expression ageing trajectories, immune-age clocks and sex-associated
differences in immune ageing. Historical `Age_Clock` names are retained only
where changing them would break paths to archived analyses or cached objects.

## Reproducibility package

Canonical Quarto reports and their supporting documentation live under
[`vignettes/`](vignettes/README.md):

- [`Fig1_umap.qmd`](vignettes/reports/Fig1_umap.qmd) — multi-tissue immune-cell expression UMAPs
- [`Fig2.qmd`](vignettes/reports/Fig2.qmd) — asynchronous ageing of immune-cell composition
- [`Fig4.qmd`](vignettes/reports/Fig4.qmd) — immune-ageing trajectories and clocks
- [`Fig5.qmd`](vignettes/reports/Fig5.qmd) — sex differences in immune ageing
- `Fig*_supplementary_plots.qmd` — supplementary analyses built from source-data workbooks

For a lightweight reproduction using the shipped figure-ready objects:

```bash
export IMMUNE_HEALTHY_BODY_MAP_ROOT=/path/to/immuneHealthyBodyMap
export IMMUNE_HEALTHY_BODY_MAP_DATA=$IMMUNE_HEALTHY_BODY_MAP_ROOT/vignettes/data
Rscript vignettes/scripts/render_reports.R --mode reduced
```

See [`vignettes/README.md`](vignettes/README.md) for the directory layout,
configuration and full HPC recomputation workflow.

## Data and fitted models

Large inputs are not stored in Git:

- Harmonised pseudobulk data: [Zenodo 10.5281/zenodo.15798373](https://doi.org/10.5281/zenodo.15798373)
- Fitted `sccomp` L3 and L0 models: companion Zenodo record pending
- Public cellNexus metadata 1.0.12: resolved by `get_publication_metadata()`
- Per-cell-type gene-model target stores: computational-data exemption requested;
  pipeline code is provided under `vignettes/targets_pipeline/`

The authoritative inventory and deposit decisions are documented in
[`vignettes/manifests/`](vignettes/manifests/README.md).

## Manuscript title

Use the following title verbatim in repository, archive and supplementary
metadata:

> A multi-tissue immune map across 4,240 single-cell donors resolves asynchronous immune ageing

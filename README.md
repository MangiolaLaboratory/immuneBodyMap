This is the code repository for the study

## A body reference map of immune cell composition and communication tracks inflammation and plasticity loss through ageing

### Abstract
The presence of immune cells in non-lymphoid tissues throughout the body is vital in the fight against infections and cancer. Yet, a detailed map of immune cell distribution and interactions throughout the body is lacking. To address this gap, we harmonised and annotated 29 million cells across 12,981 single-cell RNA sequencing samples covering 45 anatomical sites to create a comprehensive compositional and communication map of the healthy immune system. This resource represents a 30-fold increase in cell count and a 10-fold increase in tissue coverage compared to the existing state-of-the-art. We used this resource to model compositional changes in the immune system with age and diversity across demographic groups and to investigate the impact of different technologies using a novel multilevel Bayesian model. Our analysis revealed patterns of progressive tissue-specific inflammation, loss of plasticity with age, and significant differences between ethnicities. This study represents a comprehensive healthy reference for precision medicine, immunology, infectious disease and cancer.

## Code for figures

- [Figure 3](https://github.com/stemangiola/immuneHealthyBodyMap/blob/master/HCA_pipeline/figure_age.R)
- [Figure 4](https://github.com/stemangiola/immuneHealthyBodyMap/blob/master/HCA_pipeline/figure_communication.R)
- [Figure 5](https://github.com/stemangiola/immuneHealthyBodyMap/blob/master/HCA_pipeline/figure_demography.R)
- [Figure 6](https://github.com/stemangiola/immuneHealthyBodyMap/blob/master/HCA_pipeline/figure_tissue.R)
- [Figure 7](https://github.com/stemangiola/immuneHealthyBodyMap/blob/master/HCA_pipeline/figure_assay.R)

## Data

The data used for this study is available at `CuratedAtlasQueryR` version v0.99.2

https://github.com/stemangiola/CuratedAtlasQueryR/releases/tag/v0.99.2


---

## Nature Aging figure reproducibility (Age Clock)

Publication-facing Quarto reports, small reproducibility data, and Zenodo
manifests live under [`vignettes/`](vignettes/README.md).

Configure data roots with `IMMUNE_HEALTHY_BODY_MAP_DATA` (see
`vignettes/config/paths.example.env`). Large intermediate files are **not**
stored in this GitHub repository; see `vignettes/manifests/`.

Pseudobulk: https://doi.org/10.5281/zenodo.15798373
sccomp estimates and other companions: `<ZENODO_DOI_PENDING>`.

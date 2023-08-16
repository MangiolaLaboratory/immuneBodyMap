---
title: "HCA sccomp pipeline"
author: "Stefano Mangiola"
date: "2022-12-13"
output: html_document
---

```{r setup, include=FALSE}
knitr::opts_chunk$set(echo = TRUE)
```

```{r}

Rscript ~/Documents/immuneHealthyBodyMap/HCA_pipeline/create_pipeline.R /vast/projects/cellxgene_curated/metadata_annotated_0.2.3.rds ~/Documents/immuneHealthyBodyMap/sccomp_on_HCA_0.2.3.5
```

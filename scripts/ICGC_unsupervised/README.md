# ICGC_unsupervised

## 2021/04/20

### Added

- Added clustering results followed by data integration by Seurat.
  - Later proven to be negative in survival analysis:
    - P-value: 0.6589839
    - HR: 1.02
## 2021/04/19

### Fixed

- Fixed bugs in the ordering of samples, resulting in the samples in UMAP being almost completely separated by the derived cohort.

## 2021/04/18

### Added

- Add logging functions

## 2021/04/04

### Added

- PCA and UMAP of ICGC data
- Clustering result of ICGC data, resulting in three clusters
- Gene set enrichment anlaysis on cluter-specific marker genes
  - cluster 0: Nectin/Necl trans heterodimerization, angiogenesis, Adenocarcinoma of pancreas, tumor vasculature
  - cluster 1: IL-18 signaling pathway, Inflammation, Immunosuppression
  - cluster 2: blood vessel development, humoral immune response, T cell activation, Immunoregulatory interactions between a Lymphoid and a non-Lymphoid cell

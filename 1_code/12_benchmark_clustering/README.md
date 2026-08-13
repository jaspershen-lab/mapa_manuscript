# Revised clustering benchmark

Run from the `mapa_manuscript` project root:

```sh
Rscript 1_code/12_benchmark_clustering/01_curated_pathway_dataset/01_run_benchmark.R
Rscript 1_code/12_benchmark_clustering/02_informative_GO_term_dataset/01_run_benchmark.R
```

The scans cover:

- PAVER: `minClusterSize = c(1, 2, 5, 10, 15, 20)`,
  `maxCoreScatter = seq(0.01, 1, 0.01)`, coupled
  `minGap = (1 - maxCoreScatter) * 3/4`, and
  `hclust_method = c("ward.D2", "average", "complete")`.
- aPEAR: `clustMethod = c("markov", "hier")`,
  `minClusterSize = c(2, 3, 5, 10)`, and
  `simMethod = c("jaccard", "cosine", "cor")`.
- enrichplot: `layout = c("kk", "fr")`, `min_edge = c(0.1, 0.2, 0.3, 0.4)`,
  K-means clustering, and `nCluster = 2:min(100, n_pathways - 1)`.

MAPA assignments and ARIs are loaded from the existing analyses and are not
recomputed. All ARIs use the complete pathway set. Pathways filtered or left
unassigned by a method are treated as singleton clusters, while native coverage
statistics are retained in the tool-specific result tables.

Implementation notes:

- PAVER clustering is evaluated from the exact clustering core used by
  `generate_themes()` (`cosine_dissimilarity` -> `hclust` ->
  `dynamicTreeCut::cutreeDynamic`). This avoids repeating PAVER's downstream
  MDS/plot construction for every grid point and allows label-0 pathways to be
  handled consistently. The public `generate_themes()` default-call status is
  recorded separately because it fails on the 44-pathway dataset when label 0
  is present.
- aPEAR similarity matrices are computed once per `simMethod`, followed by its
  package clustering implementation. Plot-only cutoffs are not involved.
- enrichplot network coordinates are computed once per `layout`/`min_edge`
  pair. The RNG state immediately after layout generation is reused so every
  K-means result exactly matches a standalone `set.seed(123); emapplot(...)`
  call.

# GO informative-frontier benchmark builder

Run the full main-analysis scan from any working directory:

```bash
python3 mapa_manuscript/1_code/informative_go_set/build_informative_frontier_benchmark.py
```

Inputs default to the two local snapshots in
`mapa/demo_data/GO_term_human_annotation/`. Results default to
`mapa_manuscript/3_data_analysis/informative_go_set/`.

Sensitivity analysis excluding IEA is supported explicitly:

```bash
python3 mapa_manuscript/1_code/informative_go_set/build_informative_frontier_benchmark.py \
  --exclude-iea \
  --output-dir mapa_manuscript/3_data_analysis/informative_go_set_no_iea
```

The implementation is deterministic. It constructs the complete current BP
graph before filtering, propagates annotations over `is_a` and `part_of`, and
uses a fixed-point solver for anchor/member assignment. Every k-specific
directory contains the four required files (`mapa_input.csv`,
`ground_truth.csv`, `module_gene_sets.tsv`, `manifest.json`), `qc_report.json`,
and additional audit tables.

The recommended k is selected only from dataset size and lexical leakage after
removing the nine specified GO-generic words. No clustering or similarity
output from any evaluated method is read by the builder.

# Alternative Qwen model test

Run all scripts from the `mapa_manuscript` project root.

## Design

- Embedding model: `Qwen/Qwen3-Embedding-8B`.
- Clustering: MAPA 3.0.0 defaults, Louvain with cosine-similarity cutoff `0.55`.
- Original embedding comparator: `text-embedding-3-small` matrix in the existing control-data analysis.
- LLM: `Qwen/Qwen3.6-35B-A3B`.
- The Qwen LLM receives the exact prompts saved by the existing GPT run. Its literature retrieval and pathway inputs are therefore fixed to the GPT run; only the generation model changes.
- Qwen/GPT and LLM/expert semantic comparisons use `text-embedding-3-small`, matching the existing annotation-comparison workflow.
- The real/random comparison is paired by functional-module ID.
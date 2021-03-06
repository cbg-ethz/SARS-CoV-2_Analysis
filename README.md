# SARS-CoV-2_Analysis

[![Mega-Linter](https://github.com/cbg-ethz/SARS-CoV-2_Analysis/workflows/Mega-Linter/badge.svg?branch=master)](https://github.com/cbg-ethz/SARS-CoV-2_Analysis/actions?query=workflow%3AMega-Linter+branch%3Amaster)

A Snakemake workflow for large-scale SARS-CoV-2 analyses.

## Usage

First, clone the repository:

```bash
git clone https://github.com/cbg-ethz/SARS-CoV-2_Analysis
```

If needed, search for new samples using `QuerySRAAccessions.ipynb`.

The whole analysis can then be reproduced by executing the following line:

```bash
./run_pipeline.sh
```

This will download the reads, compute coverages, call variants, and plot summary figures.

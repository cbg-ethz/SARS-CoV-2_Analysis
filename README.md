# SARS-CoV-2_Analysis

A Snakemake workflow for large-scale SARS-CoV-2 analyses.

## Usage

First, clone the repository:

```bash
$ git clone --recurse-submodules https://github.com/cbg-ethz/SARS-CoV-2_Analysis
```

If needed, search for new samples using `QuerySRAAccessions.ipynb`.

The whole analysis can then be reproduced by executing the following line:

```bash
$ ./run_pipeline.sh
```

This will download the reads, run V-pipe, and plot the summary figures.

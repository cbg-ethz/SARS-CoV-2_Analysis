# SARS-CoV-2_Analysis

## Usage

First, clone the repository:

```bash
$ git clone --recurse-submodules https://github.com/cbg-ethz/SARS-CoV-2_Analysis
```

The whole analysis can then be reproduced by executing the following line (for paired-end reads only):

```bash
$ ./run_pipeline.sh
```

This will download the reads, run V-pipe, and plot the summary figures.

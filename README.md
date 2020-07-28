# SARS-CoV-2_Analysis

## Usage

First, clone the repository:

```bash
$ git clone --recurse-submodules https://github.com/cbg-ethz/SARS-CoV-2_Analysis
```

The whole analysis can then be reproduced by executing the following lines:

```bash
$ ./run_data_retrieval.sh
$ python prepare_vpipe_input.py
$ ./run_vpipe.sh
$ python gather_vcf_files.py
$ ./run_analysis.sh
```

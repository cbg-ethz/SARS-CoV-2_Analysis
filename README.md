# SARS-CoV-2_Analysis

## Usage

First, clone the repository:

```bash
$ git clone --recurse-submodules https://github.com/cbg-ethz/SARS-CoV-2_Analysis
```

The whole analysis can then be reproduced by executing the following lines (for paired-end reads only):

```bash
$ ./run_data_retrieval.sh
$ python prepare_vpipe_input.py

$ cp ./vpipe.config V-pipe/
$ ./run_vpipe.sh

$ python gather_vcf_files.py
$ ./run_analysis.sh
```

## Notes

The accession numbers of all samples used in our analyses can be found in `./assets`.

# SARS-CoV-2_Analysis

## Usage

The whole analysis can be reproduced by executing the following lines:

```bash
$ ./run_data_retrieval.sh
$ python prepare_vpipe_input.py
$ ./run_vpipe.sh
$ python gather_vcf_files.py
$ ./run_analysis.sh
```

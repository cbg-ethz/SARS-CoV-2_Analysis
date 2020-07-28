#!/usr/bin/env bash

bsub \
  -N \
  -R 'rusage[mem=4000]' \
  -q normal.24h \
  -oo snake.out -eo snake.err \
snakemake \
  --profile lsf \
  -s data_retrieval.smk \
  --use-conda \
  -pr \
  "$@"

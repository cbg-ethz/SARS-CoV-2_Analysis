#!/usr/bin/env bash

bsub \
  -N \
  -R 'rusage[mem=4000]' \
  -q normal.24h \
  -oo snake.analysis.out -eo snake.analysis.err \
snakemake \
  --profile lsf \
  -s analysis.smk \
  --use-conda \
  -pr \
  "$@"

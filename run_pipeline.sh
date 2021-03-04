#!/usr/bin/env bash

bsub \
  -N \
  -R 'rusage[mem=15000]' \
  -q normal.120h \
  -oo snake.out -eo snake.err \
snakemake \
  --profile lsf \
  -pr \
  --use-conda \
  --cores 200 \
  --local-cores 1 \
  --latency-wait 30 \
  --keep-going \
  --show-failed-logs \
  "$@"

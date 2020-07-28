#!/usr/bin/env bash

cluster='bsub -J COVID-vpipe-SRA-{rule} -M {params.mem} -n {threads} -W {params.time} -R rusage[mem={params.mem},scratch={params.scratch}] -e {log.errfile} -o {log.outfile}'

bsub \
  -N \
  -R 'rusage[mem=4000]' \
  -q normal.120h \
  -oo vpipe.out -eo vpipe.err \
snakemake -s vpipe.snake \
  -pr \
  --use-conda \
  --cluster "$cluster" \
  -j 100 \
  $@

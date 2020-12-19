bsub \
  -N \
  -R 'rusage[mem=2000]' \
  -q normal.120h \
  -oo snake.out -eo snake.err \
snakemake \
  --profile lsf \
  -pr \
  --use-conda \
  --cores 100 \
  --local-cores 1 \
  --latency-wait 30 \
  "$@"

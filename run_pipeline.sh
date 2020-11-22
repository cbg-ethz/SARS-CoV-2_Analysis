snakemake \
    -pr \
    -j 1 \
    --use-conda \
    --conda-frontend mamba \
    "$@"

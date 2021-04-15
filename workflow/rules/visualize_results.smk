rule plot_heatmaps:
    input:
        fname='results/variant_calling/vcf_data.csv',
        fname_covariates='results/data_retrieval/results/final_dataframe.csv.gz',
        fname_genes='resources/references/genes.csv',
    output:
        outdir=directory('results/analysis/plots/heatmaps/'),
    script:
        '../scripts/covid_heatmaps.R'


rule plot_histograms:
    input:
        fname='results/variant_calling/vcf_data.csv',
    output:
        fname_entropy_samples='results/analysis/results/log_entropy_samples.csv',
        fname_entropy_positions='results/analysis/results/log_entropy_positions.csv',
        outdir=directory('results/analysis/plots/histograms/'),
    script:
        '../scripts/covid_histograms.R'


rule plot_mutation_histograms:
    input:
        fname_vcf='results/variant_calling/vcf_data.csv',
    output:
        fname='results/analysis/results/mut_samples.csv',
        outdir=directory('results/analysis/plots/mutation_histograms/'),
    script:
        '../scripts/covid_histogram_mut.R'


rule compute_regression:
    input:
        fname_covariates='results/data_retrieval/results/final_dataframe.csv.gz',
        fname_entropy_samples='results/analysis/results/log_entropy_samples.csv',
    output:
        fname_lm_summary='results/analysis/regression/lm_summary.txt',
    script:
        '../scripts/covid_regression.R'


rule snv_coverage_plot:
    input:
        fname_coverage=(
            'results/data_retrieval/coverage/coverage.{accession}.csv.gz'.format(
                accession=config['input']['snv_coverage_plot']['accession']
            )
        ),
        fname_vcf='results/variant_calling/calls/{accession}.vcf'.format(
            accession=config['input']['snv_coverage_plot']['accession']
        ),
        fname_genes='resources/references/genes.csv',
    output:
        fname='results/analysis/plots/snv_coverage_plot.pdf',
    params:
        sample_accession=config['input']['snv_coverage_plot']['accession'],
        snv_highlights=[(23403, 'D614G')],
    script:
        '../scripts/snv_coverage_plot.py'


rule compute_top_positions:
    input:
        fname_entropy_positions='results/analysis/results/log_entropy_positions.csv',
        fname_genes='resources/references/genes.csv',
    output:
        fname_top_deletions='results/analysis/results/top_deletions.txt',
        fname_top_bases='results/analysis/results/top_bases.txt',
        fname_av_top_positions_pdf='results/analysis/plots/av_top_positions.pdf',
        fname_av_top_positions_png='results/analysis/plots/av_top_positions.png',
    script:
        '../scripts/covid_gene_av.R'


rule generate_extended_top_positions_table:
    input:
        fname_top_bases='results/analysis/results/top_bases.txt',
        fname_vcf='results/variant_calling/vcf_data.csv',
    output:
        fname='results/analysis/results/top_bases_extended.txt',
    script:
        '../scripts/covid_top_positions_extended_table.R'


rule generate_extended_top_samples_table:
    input:
        fname_entropy_samples='results/analysis/results/log_entropy_samples.csv',
        fname_vcf='results/variant_calling/vcf_data.csv',
    output:
        fname='results/analysis/results/top_samples_extended.txt',
    script:
        '../scripts/covid_top_samples_extended.R'

localrules: all


rule all:
    input:
        'plots/heatmaps/',
        'plots/histograms/',
        'plots/mutation_histograms/',
        'plots/snv_coverage_plot.pdf',
        'results/top_bases_extended.txt',
        'results/top_samples_extended.txt'


rule combine_se_pe_vcfs:
    input:
        vcf_dir_SE = '../../results/pipeline_run__vpipe_wrapper_SE/vcf_files',
        vcf_dir_PE = '../../results/pipeline_run__vpipe_wrapper_PE/vcf_files'
    output:
        dname = directory('data/vcf_files')
    run:
        import os
        import shutil

        se_files = list(os.scandir(input.vcf_dir_SE))
        pe_files = list(os.scandir(input.vcf_dir_PE))
        all_files = se_files + pe_files

        os.makedirs(output.dname)
        for entry in all_files:
            if not entry.name.endswith('.vcf'):
                continue
            shutil.copy(entry, f'{output.dname}/{entry.name}')


rule gather_vcf_files:
    input:
        vcf_dir = 'data/vcf_files'
    output:
        fname = 'data/vcf_data.csv'
    script:
        '../scripts/get_base_freqs.py'


rule plot_heatmaps:
    input:
        fname = 'data/vcf_data.csv',
        fname_genes = srcdir('references/genes.csv'),
        fname_covariates = srcdir(config['input']['covariates'])
    output:
        outdir = directory('plots/heatmaps/')
    script:
        '../scripts/covid_heatmaps.R'


rule plot_histograms:
    input:
        fname = 'data/vcf_data.csv'
    output:
        fname_entropy_samples = 'results/log_entropy_samples.csv',
        fname_entropy_positions = 'results/log_entropy_positions.csv',
        outdir = directory('plots/histograms/')
    script:
        '../scripts/covid_histograms.R'


rule plot_mutation_histograms:
    input:
        fname_vcf = 'data/vcf_data.csv'
    output:
        fname = 'results/mut_samples.csv',
        outdir = directory('plots/mutation_histograms/')
    script:
        '../scripts/covid_histogram_mut.R'


rule compute_regression:
    input:
        fname_covariates = srcdir(config['input']['covariates']),
        fname_entropy_samples = 'results/log_entropy_samples.csv'
    output:
        fname_lm_summary = 'regression/lm_summary.txt'
    script:
        '../scripts/covid_regression.R'


rule snv_coverage_plot:
    input:
        fname_coverage = srcdir(config['input']['snv_coverage_plot']['coverage']),
        fname_vcf = srcdir(config['input']['snv_coverage_plot']['vcf']),
        fname_genes = srcdir('../../resources/references/genes.csv')
    output:
        fname = 'plots/snv_coverage_plot.pdf'
    params:
        sample_accession = config['input']['snv_coverage_plot']['accession'],
        snv_highlights = [(23403, 'D614G')]
    script:
        '../scripts/snv_coverage_plot.py'


rule compute_top_positions:
    input:
        fname_entropy_positions = 'results/log_entropy_positions.csv',
        fname_genes = srcdir('references/genes.csv')
    output:
        fname_top_deletions = 'results/top_deletions.txt',
        fname_top_bases = 'results/top_bases.txt',
        fname_av_top_positions_pdf = 'plots/av_top_positions.pdf',
        fname_av_top_positions_png = 'plots/av_top_positions.png'
    script:
        '../scripts/covid_gene_av.R'


rule generate_extended_top_positions_table:
    input:
        fname_top_bases = 'results/top_bases.txt',
        fname_vcf = 'data/vcf_data.csv'
    output:
        fname = 'results/top_bases_extended.txt'
    script:
        '../scripts/covid_top_positions_extended_table.R'


rule generate_extended_top_samples_table:
    input:
        fname_entropy_samples = 'results/log_entropy_samples.csv',
        fname_vcf = 'data/vcf_data.csv'
    output:
        fname = 'results/top_samples_extended.txt'
    script:
        '../scripts/covid_top_samples_extended.R'

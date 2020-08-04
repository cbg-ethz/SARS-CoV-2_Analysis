configfile: 'analysis.config.yaml'
workdir: config['workdir']

localrules: all


rule all:
    input:
        'plots/heatmaps/',
        'plots/histograms/',
        'plots/snv_coverage_plot.pdf'


rule gather_vcf_files:
    input:
        vcf_dir = srcdir(config['input']['vcf_directory'])
    output:
        fname = 'data/vcf_data.csv'
    scripts:
        'scripts/get_base_freqs.py'


rule plot_heatmaps:
    input:
        fname = 'data/vcf_data.csv',
        fname_genes = srcdir('references/genes.csv')
    output:
        outdir = directory('plots/heatmaps/')
    script:
        'scripts/covid_heatmaps.R'


rule plot_histograms:
    input:
        fname = 'data/vcf_data.csv'
    output:
        entropy_samples_fname = 'results/log_entropy_samples.csv',
        entropy_positions_fname = 'results/log_entropy_positions.csv',
        outdir = directory('plots/histograms/')
    script:
        'scripts/covid_histograms.R'


rule snv_coverage_plot:
    input:
        fname_coverage = srcdir(config['input']['snv_coverage_plot']['coverage']),
        fname_vcf = srcdir(config['input']['snv_coverage_plot']['vcf']),
        fname_genes = srcdir('references/genes.csv')
    output:
        fname = 'plots/snv_coverage_plot.pdf'
    params:
        sample_accession = config['input']['snv_coverage_plot']['accession']
    script:
        'scripts/snv_coverage_plot.py'

configfile: 'analysis.config.yaml'
workdir: config['workdir']

localrules: all


rule all:
    input:
        'plots/heatmaps/',
        'plots/histograms/'


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

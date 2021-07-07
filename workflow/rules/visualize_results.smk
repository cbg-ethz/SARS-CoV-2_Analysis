rule snv_coverage_plot:
    input:
        fname_coverage=(
            "results/data_retrieval/coverage/coverage.{accession}.csv.gz".format(
                accession=config["input"]["snv_coverage_plot"]["accession"]
            )
        ),
        fname_vcf="results/variant_calling/calls/{accession}.vcf".format(
            accession=config["input"]["snv_coverage_plot"]["accession"]
        ),
        fname_genes="resources/references/genes.csv",
    output:
        fname="results/analysis/plots/snv_coverage_plot.pdf",
    params:
        sample_accession=config["input"]["snv_coverage_plot"]["accession"],
        snv_highlights=[(23403, "D614G")],
    script:
        "../scripts/snv_coverage_plot.py"


rule covid_summarise_entropy:
    input:
        fname_entropy="results/analysis/entropy.csv",
    output:
        fname_persample="results/analysis/entropy_persample.csv",
        fname_perpos="results/analysis/entropy_perposition.csv",
        fname_persamplemuts="results/analysis/entropy_persamplemutations.csv",
        fname_filtered="results/analysis/entropy_filtered.csv",
    script:
        "../scripts/covid_summarise_entropy.R"


rule covid_gene_av:
    input:
        fname_perpos="results/analysis/entropy_perposition.csv",
        fname_genes="resources/references/genes.csv",
    output:
        fname_top_regions="results/analysis/results/top_deletions.txt",
        fname_top_bases="results/analysis/results/top_bases.txt",
        fname_plot_pdf="results/analysis/plots/av_top_positions.pdf",
        fname_plot_png="results/analysis/plots/av_top_positions.png",
    script:
        "../scripts/covid_gene_av.R"


rule covid_compute_entropy:
    input:
        fname="results/variant_calling/vcf_data.csv",
        fname_samples="results/data_retrieval/results/selected_samples.csv",
    output:
        fname="results/analysis/entropy.csv",
    script:
        "../scripts/covid_compute_entropy.R"


rule covid_heatmaps_padded:
    input:
        fname="results/variant_calling/vcf_data.csv",
        fname_genes="resources/references/genes.csv",
    output:
        fname="results/analysis/plots/heatmaps.png",
    script:
        "../scripts/covid_heatmaps_padded.R"


rule plot_histograms:
    input:
        fname_persample="results/analysis/entropy_persample.csv",
        fname_perpos="results/analysis/entropy_perposition.csv",
        fname_persamplemuts="results/analysis/entropy_persamplemutations.csv",
    output:
        outdir=directory("results/analysis/plots/histograms/"),
    resources:
        mem_mb=30000,
    script:
        "../scripts/covid_histograms.R"


rule covid_regression_public:
    input:
        fname="results/variant_calling/vcf_data.csv",
        fname_entropy="results/analysis/entropy.csv",
    output:
        fname="results/analysis/covid.csv",
        fname_hist="results/analysis/plots/time_histogram.png",
        fname_regression="results/analysis/plots/time_regression.png",
    script:
        "../scripts/covid_regression_public.R"


rule covid_top_bases_extended:
    input:
        fname="results/variant_calling/vcf_data.csv",
        fname_samples="results/data_retrieval/results/selected_samples.csv",
        fname_entropy="results/analysis/entropy.csv",
    output:
        fname="results/analysis/top_bases_extended.csv",
    script:
        "../scripts/covid_top_bases_extended.R"


rule covid_top_samples:
    input:
        fname_persample="results/analysis/entropy_persample.csv",
        fname_persamplemuts="results/analysis/entropy_persamplemutations.csv",
    output:
        fname="results/analysis/top_samples.csv",
    script:
        "../scripts/covid_top_samples.R"

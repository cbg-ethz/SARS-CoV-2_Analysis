import itertools

import numpy as np
import pandas as pd

import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib.patches import Patch

from dna_features_viewer import GraphicFeature, GraphicRecord

import vcf


def flip(items, ncol):
    # https://stackoverflow.com/a/10101532/1474740
    return list(itertools.chain(*[items[i::ncol] for i in range(ncol)]))


def convert_vcf(fname):
    """Convert VCF to JSON."""
    tmp = []
    with open(fname) as fd:
        vcf_reader = vcf.Reader(fd)

        for record in vcf_reader:
            tmp.append(
                {
                    'position': record.POS,
                    'reference': record.REF,
                    'variant': [v.sequence for v in record.ALT],
                    'frequency': round(np.mean(
                        [v for k, v in record.INFO.items() if k.startswith('Freq')]
                    ), 3),
                    'posterior': round(1 - 10**(-record.QUAL / 10), 3)
                }
            )

    return pd.DataFrame(tmp)


def main(fname_coverage, fname_snv, accession, snv_highlights, fname_genes, fname_output):
    # read data
    df_coverage = pd.read_csv(fname_coverage, sep='\t')
    df_vcf = convert_vcf(fname_snv)

    df_genes = pd.read_csv(fname_genes)

    # plot
    fig, (ax_genes, ax_data) = plt.subplots(
        nrows=2, ncols=1,
        gridspec_kw={'height_ratios': [1, 5]},
        sharex=True,
        figsize=(10, 5))

    # genes
    features = df_genes.apply(
        lambda x: GraphicFeature(
            start=x.start, end=x.end, color=x.color),
        axis=1).tolist()

    record = GraphicRecord(sequence_length=df_coverage.shape[0], features=features)
    record.plot(ax=ax_genes, with_ruler=False)

    legend_elements = df_genes.apply(
        lambda x: Patch(
            facecolor=x.color, edgecolor='black', label=x.gene_name),
        axis=1).tolist()

    ncol = df_genes.shape[0] // 2 + 1
    ax_genes.legend(
        handles=flip(legend_elements, ncol), loc='upper center',
        ncol=ncol, fontsize=10, frameon=False)

    # rest
    ax_vcf = ax_data
    ax_coverage = ax_vcf.twinx()

    ax_coverage.fill_between(df_coverage.index, df_coverage[f'{accession}-19700101'], alpha=0.3, color='#88bae3')
    ax_coverage.set_ylabel('Coverage', color='#88bae3')
    ax_coverage.tick_params(axis='y', labelcolor='#88bae3')
    ax_coverage.set_ylim(bottom=0)

    ax_vcf.stem(
        df_vcf['position'], df_vcf['frequency'],
        basefmt=' ', linefmt='grey', markerfmt='ok',
        use_line_collection=True)
    ax_vcf.set_xlabel('Genomic position [bp]')
    ax_vcf.set_ylabel('SNV frequency', color='k')
    ax_vcf.tick_params(axis='y', labelcolor='k')
    ax_vcf.set_ylim(bottom=0)

    for pos, name in snv_highlights:
        freq = df_vcf.loc[df_vcf['position'] == pos, 'frequency'].iloc[0]
        ax_vcf.annotate(
            name,
            xy=(pos, freq), xytext=(pos - 3000, freq),
            arrowprops=dict(arrowstyle='->')
        )

    plt.tight_layout()
    plt.savefig(fname_output)


if __name__ == '__main__':
    main(
        snakemake.input.fname_coverage,
        snakemake.input.fname_vcf,
        snakemake.params.sample_accession,
        snakemake.params.snv_highlights,
        snakemake.input.fname_genes,
        snakemake.output.fname)

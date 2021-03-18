import pandas as pd


# setup
localrules: all


accession_list = pd.read_csv(config['sample_accession_file'])['accession'].tolist()


# job grouping requires execution dependent resources
job_resources = {
    'download_fastq': {
        'mem_mb': 1_000,
        'threads': 6
    },
    'vpipe_trim': {
        'mem_mb': 10_000,
        'threads': 1
    },
    'bwa_mem': {
        'mem_mb': 4_000,
        'threads': 8
    }
}
if config['job_grouping_mode']:
    # to be used with '--group-components data_processing=10'
    job_resources = {
        'download_fastq': {
            'mem_mb': 200,
            'threads': 1
        },
        'vpipe_trim': {
            'mem_mb': 200,
            'threads': 1
        },
        'bwa_mem': {
            'mem_mb': 200,
            'threads': 1
        }
    }


# workflow
rule download_fastq:
    output:
        fname_marker = touch('results/data_retrieval/data/{accession}.marker'),
        temp_dir = temp(directory('results/data_retrieval/data/tmp.{accession}'))
    params:
        restart_times = 10
    log:
        outfile = 'logs/download.{accession}.out.log',
        errfile = 'logs/download.{accession}.err.log'
    resources:
        mem_mb = job_resources['download_fastq']['mem_mb']
    threads: job_resources['download_fastq']['threads']
    group: 'data_processing'
    priority: -1
    run:
        import os
        import time
        import subprocess
        from pathlib import Path

        outdir = Path(os.path.dirname(output.fname_marker))

        # delete output files if they already exist
        # because fasterq-dump crashes otherwise
        for path in outdir.glob(f'{wildcards.accession}*.fastq'):
            path.unlink(missing_ok=True)

        # commence download
        counter = 0
        while True:
            if counter > 0:
                # turns out that fasterq-dump crashes for some accessions
                # when using more than 1 thread
                threads = 1

            try:
                shell(
                    'fasterq-dump'
                    ' --threads {threads}'
                    ' -p'  # --progress
                    ' --outdir {outdir}'
                    ' --temp {output.temp_dir}'
                    ' {wildcards.accession}'
                    ' > >(tee {log.outfile}) 2>&1'
                )
            except subprocess.CalledProcessError:
                print('Download process crashed, hopefully this is just a fluke...')
                time.sleep(100)

            # make sure the files were actually created
            available_files = list(outdir.glob(f"{wildcards.accession}*.fastq"))
            if len(available_files) in (1, 2, 3):
                # downloaded SE, PE, varying read number per spot
                break

            # no files were downloaded, retry...
            shell('echo "Download failed, restarting" >> {log.errfile}')
            counter += 1

            if counter > params.restart_times:
                raise RuntimeError(f'Download failed {counter} times')


rule vpipe_trim:
    input:
        fname_marker = 'results/data_retrieval/data/{accession}.marker'
    output:
        touch('results/data_retrieval/trimmed/{accession}.marker')
    params:
        extra = '-ns_max_n 4 -min_qual_mean 30 -trim_qual_left 30 -trim_qual_right 30 -trim_qual_window 10',
        trim_percent_cutoff = .8,
    log:
        outfile = 'logs/trimming.{accession}.out.log',
        errfile = 'logs/trimming.{accession}.err.log'
    resources:
        mem_mb = job_resources['vpipe_trim']['mem_mb']
    threads: job_resources['vpipe_trim']['threads']
    group: 'data_processing'
    conda:
        '../envs/preprocessing.yaml'
    priority: 1
    shell:
        """
        # compute length cutoff
        fname_marker="{input.fname_marker}"
        fastq_fname="${{fname_marker%.marker}}.fastq"
        if [ ! -f "$fastq_fname" ]; then
            # is paired-end
            fastq_fname="${{fname_marker%.marker}}_1.fastq"
        fi

        read_len=$(bioawk -c fastx '{{{{ bases += length($seq); count++ }}}} END {{{{print int(bases/count)}}}}' "$fastq_fname")

        len_cutoff=$(echo "($read_len * {params.trim_percent_cutoff}) / 1" | bc)

        # detect SE/PE read type
        filecount=$(ls results/data_retrieval/data/{wildcards.accession}*.fastq | wc -l)
        case $filecount in
            1)
                # SE reads
                echo "Read type: SE" >> {log.outfile}
                input_spec="-fastq results/data_retrieval/data/{wildcards.accession}.fastq"
                ;;
            2)
                # PE reads
                echo "Read type: PE" >> {log.outfile}
                input_spec="-fastq results/data_retrieval/data/{wildcards.accession}_1.fastq -fastq2 results/data_retrieval/data/{wildcards.accession}_2.fastq"
                ;;
            3)
                # some runs have a variable number of reads per spot.
                # why? how? we might never know.
                # for now, let's pretend everything is a-okay.
                echo "Read type: 'variable number of reads per spot'" >> {log.outfile}
                input_spec="-fastq results/data_retrieval/data/{wildcards.accession}_1.fastq -fastq2 results/data_retrieval/data/{wildcards.accession}_2.fastq"

                rm results/data_retrieval/data/{wildcards.accession}.fastq # there is nothing to see here, walk along
                ;;
            *)
                # oh no
                exit 1
                ;;
        esac

        # do trimming
        prinseq-lite.pl \
            $input_spec \
            {params.extra} \
            -out_format 3 \
            -out_good results/data_retrieval/trimmed/{wildcards.accession} \
            -out_bad null \
            -min_len "$len_cutoff" \
            -log {log.outfile} 2> >(tee {log.errfile} >&2)

        # remove singletons
        rm -f results/data_retrieval/trimmed/{wildcards.accession}*_singletons.fastq

        # if no reads survive, create empty fastq file
        # (such that mapping does not fail)
        trimmedfilecount=$(shopt -s nullglob; files=(results/data_retrieval/trimmed/{wildcards.accession}*.fastq); echo ${{#files[@]}})
        if [ "$trimmedfilecount" -eq "0" ]; then
            echo "No non-singletons survived trimming, creating empty FastQ file" >> {log.outfile}
            touch results/data_retrieval/trimmed/{wildcards.accession}.fastq
        fi
        """


rule bwa_index:
    input:
        config['reference']
    output:
        'results/data_retrieval/references/reference.amb',
        'results/data_retrieval/references/reference.ann',
        'results/data_retrieval/references/reference.bwt',
        'results/data_retrieval/references/reference.pac',
        'results/data_retrieval/references/reference.sa'
    log:
        'logs/bwa_index.log'
    params:
        prefix = 'results/data_retrieval/references/reference'
    wrapper:
        '0.68.0/bio/bwa/index'


rule bwa_mem:
    input:
        fname_marker = 'results/data_retrieval/trimmed/{accession}.marker',
        index = 'results/data_retrieval/references/reference.amb',
        fname_ref = config['reference']
    output:
        fname_cram = 'results/data_retrieval/alignment/{accession}.cram'
    log:
        outfile = 'logs/alignment.{accession}.out.log',
        errfile = 'logs/alignment.{accession}.err.log'
    params:
        index = 'results/data_retrieval/references/reference',
        sort = 'samtools',
        sort_order = 'coordinate',
    resources:
        mem_mb = job_resources['bwa_mem']['mem_mb']
    threads: job_resources['bwa_mem']['threads']
    conda:
        '../envs/alignment.yaml'
    group: 'data_processing'
    priority: 2
    shell:
        """
        # remove (potential) left over files from previous runs
        rm -vf "{output.fname_cram}.tmp."*".bam"

        # run alignment
        (bwa mem \
            -t {threads} \
            {params.index} \
            results/data_retrieval/trimmed/{wildcards.accession}*.fastq \
            | samtools sort \
                --threads {threads} \
                --reference {input.fname_ref} \
                --output-fmt CRAM,embed_ref,use_bzip2,use_lzma,level=9,seqs_per_slice=1000000 \
                -o {output.fname_cram}) \
        > {log.outfile} 2> {log.errfile}

        # delete used fastq files
        if [ {config[data_saver_mode]} == True ]; then
            rm results/data_retrieval/data/{wildcards.accession}*.fastq #results/data_retrieval/data/{wildcards.accession}.marker
            rm results/data_retrieval/trimmed/{wildcards.accession}*.fastq #results/data_retrieval/trimmed/{wildcards.accession}.marker
        fi
        """


rule samtools_index:
    input:
        'results/data_retrieval/alignment/{accession}.cram'
    output:
        'results/data_retrieval/alignment/{accession}.cram.crai'
    group: 'data_processing'
    priority: 3
    wrapper:
        '0.68.0/bio/samtools/index'


rule compute_coverage:
    input:
        fname = 'results/data_retrieval/alignment/{accession}.cram',
        index = 'results/data_retrieval/alignment/{accession}.cram.crai'
    output:
        fname = 'results/data_retrieval/coverage/coverage.{accession}.csv.gz'
    group: 'data_processing'
    priority: 3
    run:
        import pysam

        import numpy as np
        import pandas as pd

        bam = pysam.AlignmentFile(input.fname, 'rc')
        assert bam.has_index()
        assert len(bam.references) == 1
        ref = bam.references[0]

        coverage = np.sum(bam.count_coverage(ref), axis=0)

        pd.DataFrame({
            wildcards.accession: coverage
        }).to_csv(output.fname, index=False)


rule aggregate_results:
    input:
        fname_list = expand(
            'results/data_retrieval/coverage/coverage.{accession}.csv.gz',
            accession=accession_list)
    output:
        fname = 'results/data_retrieval/results/coverage_quantiles.csv.gz',
    benchmark:
        'benchmarks/aggregate_results.benchmark.txt'
    params:
        quantile_list = [0.25, 0.5, 0.75]
    resources:
        mem_mb = 15_000
    run:
        from pathlib import Path
        import pandas as pd
        from tqdm import tqdm

        tmp = []
        for fname in tqdm(input.fname_list):
            accession = Path(fname).stem.split('.')[1]
            df_tmp = pd.read_csv(fname, squeeze=True)

            tmp.append({
                'accession': accession,
                **{f'q{q}': df_tmp.quantile(q=q) for q in params.quantile_list}
            })

        pd.DataFrame(tmp).sort_values('accession').to_csv(output.fname, index=False)


rule plot_coverage_per_locus:
    input:
        unpack(lambda wildcards: {f'fname_coverage_{accession}': f'results/data_retrieval/coverage/coverage.{accession}.csv.gz' for accession in accession_list}),
        fname_selection = 'results/data_retrieval/results/selected_samples.csv'
    output:
        fname = report('results/data_retrieval/plots/coverage_per_locus.pdf', caption='report/empty_caption.rst')
    benchmark:
        'benchmarks/plot_coverage_per_locus.benchmark.txt'
    params:
        target_sample_num = 50_000
    resources:
        mem_mb = 30_000
    run:
        import itertools

        import pandas as pd

        import seaborn as sns
        import matplotlib.pyplot as plt
        from matplotlib.patches import Patch
        from matplotlib.gridspec import GridSpec

        from dna_features_viewer import GraphicFeature, GraphicRecord

        # fix big boom on macOS
        import matplotlib
        matplotlib.use('Agg')

        # helper functions
        def flip(items, ncol):
            # https://stackoverflow.com/a/10101532/1474740
            return list(itertools.chain(*[items[i::ncol] for i in range(ncol)]))

        # read data
        selected_samples = pd.read_csv(input.fname_selection, squeeze=True)
        sample_freq = max(1, selected_samples.shape[0] // params.target_sample_num)

        cov_list = []
        for accession in selected_samples[::sample_freq]:
            cov_list.append(pd.read_csv(input[f'fname_coverage_{accession}']))
        df_cov = pd.concat(cov_list, axis=1)

        df_genes = pd.read_csv('resources/references/genes.csv')
        df_primers = pd.read_csv(
            'resources/references/nCoV-2019.bed', sep='\t', header=None,
            names=['chrom', 'chromStart', 'chromEnd', 'name', 'foo', 'strand'])

        # compute data statistics
        df_stats = pd.DataFrame({
            'lower quartile': df_cov.quantile(q=.25, axis=1),
            'median': df_cov.quantile(q=.5, axis=1),
            'upper quartile': df_cov.quantile(q=.75, axis=1)
        })

        # prepare primers
        df_primers_sub = df_primers[~df_primers['name'].str.contains('alt')].copy()
        df_primers_sub['primer_id'] = df_primers_sub['name'].str.split('_').str[:2].str.join('_')
        df_pos = pd.DataFrame({
            'start': df_primers_sub.groupby('primer_id')['chromStart'].min(),
            'end': df_primers_sub.groupby('primer_id')['chromEnd'].max()
        })

        # plot data
        fig, (ax_genes, ax_coverage) = plt.subplots(
            nrows=2, ncols=1,
            gridspec_kw={'height_ratios': [1, 5]},
            sharex=True,
            figsize=(10, 6))

        # genes
        features = df_genes.apply(
            lambda x: GraphicFeature(
                start=x.start, end=x.end, color=x.color),
            axis=1).tolist()

        record = GraphicRecord(sequence_length=df_cov.shape[0], features=features)
        record.plot(ax=ax_genes, with_ruler=False)

        legend_elements = df_genes.apply(
            lambda x: Patch(
                facecolor=x.color, edgecolor='black', label=x.gene_name),
            axis=1).tolist()

        ncol = df_genes.shape[0] // 2 + 1
        ax_genes.legend(
            handles=flip(legend_elements, ncol), loc='upper center',
            ncol=ncol, fontsize=10, frameon=False)

        # coverage
        ax_coverage.plot(
            df_stats.index,
            df_stats['median'],
            label='median')
        ax_coverage.plot(
            df_stats.index,
            df_stats['lower quartile'],
            alpha=.5,
            label='lower quartile')
        ax_coverage.plot(
            df_stats.index,
            df_stats['upper quartile'],
            alpha=.5,
            label='upper quartile')

        ax_coverage.set_xlabel('Genomic position [bp]')
        ax_coverage.set_ylabel('Per base read count')

        ax_coverage.legend(loc='best')

        # primers
        for row in df_pos.sort_values('start').itertuples():
            ax_coverage.axvline(row.start, 0, .04, color='black', lw=1, ls='solid')
            ax_coverage.axvline(row.end, 0, .04, color='black', lw=1, ls='dashed')

        fig.tight_layout()
        fig.savefig(output.fname)


rule plot_coverage_per_sample:
    input:
        fname = 'results/data_retrieval/results/coverage_quantiles.csv.gz',
    output:
        fname = report('results/data_retrieval/plots/coverage_per_sample.pdf', caption='report/empty_caption.rst')
    benchmark:
        'benchmarks/plot_coverage_per_sample.benchmark.txt'
    params:
        target_sample_num = 1000
    resources:
        mem_mb = 5_000
    run:
        from pathlib import Path

        import pandas as pd

        import seaborn as sns
        import matplotlib.pyplot as plt

        # fix big boom on macOS
        import matplotlib
        matplotlib.use('Agg')

        # read data
        df = pd.read_csv(
            input.fname, index_col='accession',
            low_memory=False
        ).sort_values('q0.5')

        # subsample data for performant plot
        sample_freq = max(1, df.shape[0] // params.target_sample_num)
        df = df.iloc[::sample_freq]

        # plot data
        plt.figure(figsize=(10, 6))

        plt.plot(
            df['q0.5'],
            label='median')
        plt.plot(
            df['q0.25'],
            alpha=.5,
            label='lower quartile')
        plt.plot(
            df['q0.75'],
            alpha=.5,
            label='upper quartile')

        plt.axhline(
            config['thresholds']['median_minium'],
            color='black', ls='dashed', alpha=0.2)
        plt.axhline(
            config['thresholds']['quartile_range'][0],
            color='black', ls='dashed', alpha=0.2)
        plt.axhline(
            config['thresholds']['quartile_range'][1],
            color='black', ls='dashed', alpha=0.2)

        plt.xlabel('SRA Accession')
        plt.ylabel('Per base read count')
        plt.yscale('log')
        plt.tick_params(
            axis='x',
            which='both',
            labelbottom=False)

        plt.legend(loc='best')

        plt.tight_layout()
        plt.savefig(output.fname)


rule retrieve_sra_metadata:
    output:
        fname = 'results/data_retrieval/results/sra_metadata.csv.gz'
    benchmark:
        'benchmarks/retrieve_sra_metadata.benchmark.txt'
    run:
        import json

        import pandas as pd
        from tqdm import tqdm

        from pysradb.sraweb import SRAweb

        # chunks are necessary because `SRAweb` crashes otherwise
        chunk_size = 200

        def chunker(seq, size):
            return list(seq[pos:pos + size] for pos in range(0, len(seq), size))

        # query SRA
        db = SRAweb()

        df_list = []
        for sub_list in tqdm(chunker(accession_list, chunk_size)):
            while True:  # save us from network issues
                try:
                    tmp = db.sra_metadata(sub_list, detailed=True)
                    df_list.append(tmp)
                    break
                except (KeyError, json.JSONDecodeError) as e:
                    print(f'Woopsie ({e}) starting with', sub_list[0])
                    continue

        df_meta = pd.concat(df_list)
        assert set(df_meta['run_accession'].tolist()) == set(accession_list)

        # engineer additional features
        # assert 'location' not in df_meta.columns
        # df_meta['location'] = (
        #     df_meta['geographic location (country and/or sea)'].combine_first(
        #         df_meta['geo_loc_name'].str
        #                                .split(':')
        #                                .str[0]
        #     ).replace({
        #         'Not Applicable': pd.NA,
        #         'not applicable': pd.NA,
        #         'missing': pd.NA,
        #         'not collected': pd.NA
        #     })
        # )

        # save data
        df_meta.to_csv(output.fname, index=False)


rule compute_additional_properties:
    input:
        fname_list = expand(
            'results/data_retrieval/alignment/{accession}.cram',
            accession=accession_list),
        index_list = expand(
            'results/data_retrieval/alignment/{accession}.cram.crai',
            accession=accession_list)
    output:
        fname = 'results/data_retrieval/results/extra_properties.csv.gz'
    benchmark:
        'benchmarks/compute_additional_properties.benchmark.txt'
    resources:
        mem_mb = 20_000
    run:
        import os
        import glob

        import numpy as np
        import pandas as pd

        import pysam
        from tqdm import tqdm

        tmp = []
        for fname in tqdm(input.fname_list):
            accession = os.path.splitext(os.path.basename(fname))[0]
            cram = pysam.AlignmentFile(fname, 'rc')

            read_len = np.mean(
                [read.infer_read_length()
                 for read in cram.fetch(until_eof=True)
                 if read.infer_read_length() is not None]
            ).astype(int)

            tmp.append({
                'accession': accession,
                'avg_read_length': read_len
            })

            cram.close()

        pd.DataFrame(tmp).to_csv(output.fname, index=False)


rule assemble_final_dataframe:
    input:
        fname_quantiles = 'results/data_retrieval/results/coverage_quantiles.csv.gz',
        fname_extra = 'results/data_retrieval/results/extra_properties.csv.gz',
        fname_meta = 'results/data_retrieval/results/sra_metadata.csv.gz'
    output:
        fname = report('results/data_retrieval/results/final_dataframe.csv.gz')
    resources:
        mem_mb = 40_000
    run:
        import pandas as pd

        # read data
        df_quantiles = pd.read_csv(input.fname_quantiles, index_col='accession', low_memory=False)

        df_extra = pd.read_csv(input.fname_extra, index_col='accession')
        df_meta = pd.read_csv(input.fname_meta, index_col='run_accession')

        assert sorted(df_quantiles.index) == sorted(df_extra.index) == sorted(df_meta.index)

        # merge data
        df_merge = pd.concat([
            df_quantiles.rename(columns={
                'q0.25': 'per_base_read_count_lower_quartile',
                'q0.5': 'per_base_read_count_median',
                'q0.75': 'per_base_read_count_upper_quartile'
            }),
            df_extra,
            df_meta
        ], axis=1)

        df_merge.to_csv(output.fname)


rule select_samples:
    input:
        fname = 'results/data_retrieval/results/final_dataframe.csv.gz'
    output:
        fname = report('results/data_retrieval/results/selected_samples.csv', caption='report/empty_caption.rst')
    resources:
        mem_mb = 20_000
    run:
        import pandas as pd

        # read data
        df = pd.read_csv(input.fname, index_col=0)

        # apply thresholds
        df_sub = df[
            (df['per_base_read_count_median'] >= config['thresholds']['median_minium']) &
            (df['per_base_read_count_lower_quartile'] >= config['thresholds']['quartile_range'][0]) &
            (df['per_base_read_count_upper_quartile'] <= config['thresholds']['quartile_range'][1])
        ]

        selection = df_sub.index.tolist()

        # save results
        with open(output.fname, 'w') as fd:
            fd.write('accession\n')
            fd.write('\n'.join(selection))

import pandas as pd


localrules: all


accession_list = pd.read_csv(srcdir('../../'  + config['sample_accession_file']))['accession'].tolist()


def len_cutoff(wildcards, trim_percent_cutoff=.8):
    """Compute minimal length based on average read length."""
    import glob

    import numpy as np
    from snakemake import shell

    available_files = list(sorted(
        glob.glob(f'data/{wildcards.accession}*.fastq')))
    if len(available_files) not in (1, 2, 3):
        raise RuntimeError(
            'Unexpected number of FastQ files for ' +
            f'{wildcards.accession}: {len(available_files)}')

    fname = available_files[0]
    read_len = int(shell.check_output(
        f"bioawk -c fastx '{{{{ bases += length($seq); count++ }}}} END{{{{print int(bases/count)}}}}' {fname}"
    ).rstrip())

    len_cutoff = int(trim_percent_cutoff * read_len)
    return len_cutoff


rule all:
    input:
        'plots/coverage_per_locus.pdf',
        'plots/coverage_per_sample.pdf',
        'results/selected_samples.csv',
        'processed_fastq_data'


rule download_fastq:
    output:
        fname_marker = touch('data/{accession}.marker')
    params:
        restart_times = 10
    log:
        outfile = 'logs/download.{accession}.out.log',
        errfile = 'logs/download.{accession}.err.log'
    resources:
        mem_mb = 5_000
    threads: 6
    group: 'data_processing'
    priority: -1
    run:
        import os
        import time
        import subprocess
        from pathlib import Path

        outdir = Path(os.path.dirname(output.fname_marker))
        tmpdir = Path(os.path.join(outdir, f'tmp.{wildcards.accession}'))

        counter = 0
        while True:
            try:
                shell(
                    'fasterq-dump --threads {threads} --outdir {outdir} --temp {tmpdir} {wildcards.accession} > >(tee {log.outfile}) 2>&1'
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
        fname_marker = 'data/{accession}.marker'
    output:
        touch('trimmed/{accession}.marker')
    params:
        extra = '-ns_max_n 4 -min_qual_mean 30 -trim_qual_left 30 -trim_qual_right 30 -trim_qual_window 10',
        len_cutoff = len_cutoff
    log:
        outfile = 'logs/trimming.{accession}.out.log',
        errfile = 'logs/trimming.{accession}.err.log'
    resources:
        mem_mb = 10_000
    group: 'data_processing'
    conda:
        '../envs/preprocessing.yaml'
    priority: 1
    shell:
        """
        echo "The length cutoff is: {params.len_cutoff}" >> {log.outfile}

        # detect SE/PE read type
        filecount=$(ls data/{wildcards.accession}*.fastq | wc -l)
        case $filecount in
            1)
                # SE reads
                echo "Read type: SE" >> {log.outfile}
                input_spec="-fastq data/{wildcards.accession}.fastq"
                ;;
            2)
                # PE reads
                echo "Read type: PE" >> {log.outfile}
                input_spec="-fastq data/{wildcards.accession}_1.fastq -fastq2 data/{wildcards.accession}_2.fastq"
                ;;
            3)
                # some runs have a variable number of reads per spot.
                # why? how? we might never know.
                # for now, let's pretend everything is a-okay.
                echo "Read type: 'variable number of reads per spot'" >> {log.outfile}
                input_spec="-fastq data/{wildcards.accession}_1.fastq -fastq2 data/{wildcards.accession}_2.fastq"

                rm data/{wildcards.accession}.fastq # there is nothing to see here, walk along
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
            -out_good trimmed/{wildcards.accession} \
            -out_bad null \
            -min_len {params.len_cutoff} \
            -log {log.outfile} 2> >(tee {log.errfile} >&2)

        # remove singletons
        rm -f trimmed/{wildcards.accession}*_singletons.fastq

        # if no reads survive, create empty fastq file
        # (such that mapping does not fail)
        trimmedfilecount=$(shopt -s nullglob; files=(trimmed/{wildcards.accession}*.fastq); echo ${{#files[@]}})
        if [ "$trimmedfilecount" -eq "0" ]; then
            echo "No non-singletons survived trimming, creating empty FastQ file" >> {log.outfile}
            touch trimmed/{wildcards.accession}.fastq
        fi
        """


rule bwa_index:
    input:
        srcdir('../../' + config['reference'])
    output:
        'references/reference.amb',
        'references/reference.ann',
        'references/reference.bwt',
        'references/reference.pac',
        'references/reference.sa'
    log:
        'logs/bwa_index.log'
    params:
        prefix = 'references/reference'
    wrapper:
        '0.68.0/bio/bwa/index'


rule bwa_mem:
    input:
        fname_marker = 'trimmed/{accession}.marker',
        index = 'references/reference.amb',
        fname_ref = srcdir('../../' + config['reference'])
    output:
        fname_cram = 'alignment/{accession}.cram'
    log:
        outfile = 'logs/alignment.{accession}.out.log',
        errfile = 'logs/alignment.{accession}.err.log'
    params:
        index = 'references/reference',
        sort = 'samtools',
        sort_order = 'coordinate',
    resources:
        mem_mb = 16_000
    threads: 8
    conda:
        '../envs/alignment.yaml'
    group: 'data_processing'
    priority: 2
    shell:
        """
        (bwa mem \
            -t {threads} \
            {params.index} \
            trimmed/{wildcards.accession}*.fastq \
            | samtools sort \
                --threads {threads} \
                --reference {input.fname_ref} \
                --output-fmt CRAM,embed_ref,use_bzip2,use_lzma,level=9,seqs_per_slice=1000000 \
                -o {output.fname_cram}) \
        > {log.outfile} 2> {log.errfile}

        # delete used fastq files
        if [ {config[data_saver_mode]} == True ]; then
            rm data/{wildcards.accession}*.fastq data/{wildcards.accession}.marker
            rm trimmed/{wildcards.accession}*.fastq trimmed/{wildcards.accession}.marker
        fi
        """


rule samtools_index:
    input:
        'alignment/{accession}.cram'
    output:
        'alignment/{accession}.cram.crai'
    group: 'data_processing'
    priority: 3
    wrapper:
        '0.68.0/bio/samtools/index'


rule compute_coverage:
    input:
        fname = 'alignment/{accession}.cram',
        index = 'alignment/{accession}.cram.crai'
    output:
        fname = 'coverage/coverage.{accession}.csv.gz'
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
            'coverage/coverage.{accession}.csv.gz',
            accession=accession_list)
    output:
        fname = 'results/coverage.csv.gz',
        fname_stats = report('results/statistics.csv.gz', caption='report/empty_caption.rst'),
        fname_lowquar = 'results/coverage_lowerquartile.csv.gz',
        fname_median = 'results/coverage_median.csv.gz',
        fname_upperquar = 'results/coverage_upperquartile.csv.gz'
    benchmark:
        'benchmarks/aggregate_results.benchmark.txt'
    resources:
        mem_mb = 20_000
    run:
        import pandas as pd

        df_list = []
        for fname in input.fname_list:
            tmp = pd.read_csv(fname)
            df_list.append(tmp)
        df = pd.concat(df_list, axis=1)

        # save data and basic statistics
        df.to_csv(output.fname, index=False)
        df.describe().to_csv(output.fname_stats)

        # save useful metrics
        (pd.melt(df).groupby('variable')
                    .quantile(q=.25)
                    .sort_values('value')
                    .to_csv(output.fname_lowquar))
        (pd.melt(df).groupby('variable')
                    .quantile(q=.5)
                    .sort_values('value')
                    .to_csv(output.fname_median))
        (pd.melt(df).groupby('variable')
                    .quantile(q=.75)
                    .sort_values('value')
                    .to_csv(output.fname_upperquar))


rule plot_coverage_per_locus:
    input:
        fname = 'results/coverage.csv.gz',
        fname_selection = 'results/selected_samples.csv'
    output:
        fname = report('plots/coverage_per_locus.pdf', caption='report/empty_caption.rst')
    benchmark:
        'benchmarks/plot_coverage_per_locus.benchmark.txt'
    resources:
        mem_mb = 5_000
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
        df = pd.read_csv(input.fname)
        sel = pd.read_csv(input.fname_selection, squeeze=True)
        df_sub = df[sel]

        df_genes = pd.read_csv('../../resources/references/genes.csv')
        df_primers = pd.read_csv(
            '../../resources/references/nCoV-2019.bed', sep='\t', header=None,
            names=['chrom', 'chromStart', 'chromEnd', 'name', 'foo', 'strand'])

        # compute data statistics
        df_stats = pd.DataFrame({
            'lower quartile': df_sub.quantile(q=.25, axis=1),
            'median': df_sub.quantile(q=.5, axis=1),
            'upper quartile': df_sub.quantile(q=.75, axis=1)
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

        record = GraphicRecord(sequence_length=df.shape[0], features=features)
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
        fname_lowquar = 'results/coverage_lowerquartile.csv.gz',
        fname_median = 'results/coverage_median.csv.gz',
        fname_upperquar = 'results/coverage_upperquartile.csv.gz'
    output:
        fname = report('plots/coverage_per_sample.pdf', caption='report/empty_caption.rst')
    benchmark:
        'benchmarks/plot_coverage_per_sample.benchmark.txt'
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
        df_lq = pd.read_csv(input.fname_lowquar, index_col=0)
        df_mq = pd.read_csv(input.fname_median, index_col=0)
        df_uq = pd.read_csv(input.fname_upperquar, index_col=0)

        # plot data
        accession_order = df_mq.sort_values('value').index.tolist()

        plt.figure(figsize=(10, 6))

        plt.plot(
            accession_order,
            df_mq.loc[accession_order, 'value'],
            label='median')
        plt.plot(
            accession_order,
            df_lq.loc[accession_order, 'value'],
            alpha=.5,
            label='lower quartile')
        plt.plot(
            accession_order,
            df_uq.loc[accession_order, 'value'],
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
        fname = 'results/sra_metadata.csv.gz'
    benchmark:
        'benchmarks/retrieve_sra_metadata.benchmark.txt'
    run:
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
                except KeyError:
                    print('Woopsie starting with', sub_list[0])
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
            'alignment/{accession}.cram',
            accession=accession_list),
        index_list = expand(
            'alignment/{accession}.cram.crai',
            accession=accession_list)
    output:
        fname = 'results/extra_properties.csv.gz'
    benchmark:
        'benchmarks/compute_additional_properties.benchmark.txt'
    resources:
        mem_mb = 10_000
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
                 for read in cram.fetch(until_eof=True)]
            ).astype(int)

            tmp.append({
                'accession': accession,
                'avg_read_length': read_len
            })

        pd.DataFrame(tmp).to_csv(output.fname, index=False)


rule assemble_final_dataframe:
    input:
        fname_lowquar = 'results/coverage_lowerquartile.csv.gz',
        fname_median = 'results/coverage_median.csv.gz',
        fname_upperquar = 'results/coverage_upperquartile.csv.gz',
        fname_extra = 'results/extra_properties.csv.gz',
        fname_meta = 'results/sra_metadata.csv.gz'
    output:
        fname = report('results/final_dataframe.csv.gz')
    run:
        import pandas as pd

        # read data
        df_lq = pd.read_csv(input.fname_lowquar, index_col='variable')
        df_mq = pd.read_csv(input.fname_median, index_col='variable')
        df_uq = pd.read_csv(input.fname_upperquar, index_col='variable')

        df_extra = pd.read_csv(input.fname_extra, index_col='accession')
        df_meta = pd.read_csv(input.fname_meta, index_col='run_accession')

        assert sorted(df_lq.index) == sorted(df_mq.index) == sorted(df_uq.index) == sorted(df_extra.index) == sorted(df_meta.index)

        # merge data
        df_merge = pd.concat([
            df_lq.rename(columns={'value': 'per_base_read_count_lower_quartile'}),
            df_mq.rename(columns={'value': 'per_base_read_count_median'}),
            df_uq.rename(columns={'value': 'per_base_read_count_upper_quartile'}),
            df_extra,
            df_meta
        ], axis=1)

        df_merge.to_csv(output.fname)


rule select_samples:
    input:
        fname = 'results/final_dataframe.csv.gz'
    output:
        fname = report('results/selected_samples.csv', caption='report/empty_caption.rst')
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


rule prepare_vpipe_input:
    input:
        fname_selected = 'results/selected_samples.csv',
        fname_final = 'results/final_dataframe.csv.gz'
    output:
        out_dir = directory('processed_fastq_data')
    params:
        max_batch_size = None
    run:
        import os
        import collections
        from pathlib import Path

        import pandas as pd

        from tqdm import tqdm

        # read input
        df = pd.read_csv(input.fname_selected)
        accession_list = df['accession'].tolist()

        df_final = pd.read_csv(input.fname_final, index_col=0)
        readlen_dict = df_final.to_dict()['avg_read_length']

        max_batch_size = params.max_batch_size

        name_template = 'samples_{type_}'
        if max_batch_size is not None:
            name_template += '_batch{batch:02}'
        dummy_date = '19700101'

        target_dir = Path(output.out_dir)

        # create batched output
        current_batch_id = collections.defaultdict(int)
        current_batch_size = collections.defaultdict(int)
        for accession in tqdm(accession_list):
            # gather fastq files
            available_files = list(Path('data').glob(f'{accession}*.fastq'))

            # detect type of experiment
            if len(available_files) == 1:
                # SE
                type_ = 'SE'
            elif len(available_files) == 2:
                # PE
                type_ = 'PE'
            elif len(available_files) == 3:
                # PE
                type_ = 'PE'
            else:
                raise RuntimeError(available_files)

            # assemble target names
            name = name_template.format(
                type_=type_,
                batch=current_batch_id[type_])
            target = target_dir / name / accession / dummy_date / 'raw_data'

            target.mkdir(parents=True, exist_ok=True)

            # copy/link fastq files
            for path in available_files:
                basename = path.name

                if len(available_files) == 3:
                    # if there is a varying number of reads per spot,
                    # we only consider PE reads
                    if '_1' not in basename and '_2' not in basename:
                        print('skipping', path)
                        continue

                # make V-pipe recognize PE files
                basename = basename.replace('_1', '_R1')
                basename = basename.replace('_2', '_R2')

                # create hard link
                dest = target / basename
                # print(accession, path, dest)

                os.link(path, dest)

            # add meta-info to tsv file
            with open(target_dir / f'{name}.tsv', 'a') as fd:
                fd.write(f'{accession}\t{dummy_date}\t{readlen_dict[accession]}\n')

            # handle batches
            current_batch_size[type_] += 1
            if max_batch_size is not None and current_batch_size[type_] == max_batch_size:
                current_batch_size[type_] = 0
                current_batch_id[type_] += 1

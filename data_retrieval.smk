configfile: 'data_retrieval.config.yaml'
workdir: config['workdir']

localrules: all


accession_list = config['sample_accessions']


def len_cutoff(wildcards, trim_percent_cutoff=.8):
    """Compute minimal length based on average read length."""
    import glob

    import numpy as np
    from Bio import SeqIO

    available_files = list(sorted(
        glob.glob(f'data/{wildcards.accession}*.fastq')))
    if len(available_files) not in (1, 2, 3):
        raise RuntimeError(
            'Unexpected number of FastQ files for ' +
            f'{wildcards.accession}: {len(available_files)}')

    fname = available_files[0]
    read_len = np.mean(
        [len(r.seq) for r in SeqIO.parse(fname, 'fastq')]
    ).astype(int)

    len_cutoff = int(trim_percent_cutoff * read_len)
    return len_cutoff


rule all:
    input:
        'plots/coverage_per_locus.pdf',
        'plots/coverage_per_sample.pdf',
        'results/selected_samples.csv'


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
    # group: 'data_processing'
    run:
        import os
        import glob
        import time
        import subprocess

        outdir = os.path.dirname(output.fname_marker)
        tmpdir = os.path.join(outdir, f'tmp.{wildcards.accession}')

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
            available_files = glob.glob(f'data/{wildcards.accession}*.fastq')
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
    # group: 'data_processing'
    conda:
        'envs/preprocessing.yaml'
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
        srcdir(config['reference'])
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
        '0.51.2/bio/bwa/index'


rule bwa_mem:
    input:
        fname_marker = 'trimmed/{accession}.marker',
        index = 'references/reference.amb'
    output:
        fname_bam = 'alignment/{accession}.bam'
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
    # group: 'data_processing'
    shell:
        """
        (bwa mem \
            -t {threads} \
            {params.index} \
            trimmed/{wildcards.accession}*.fastq \
            | samtools sort -o {output.fname_bam}) \
        > {log.outfile} 2> {log.errfile}
        """


rule samtools_index:
    input:
        'alignment/{accession}.bam'
    output:
        'alignment/{accession}.bam.bai'
    # group: 'data_processing'
    wrapper:
        '0.51.2/bio/samtools/index'


rule compute_coverage:
    input:
        fname = 'alignment/{accession}.bam',
        index = 'alignment/{accession}.bam.bai'
    output:
        fname = 'coverage/coverage.{accession}.csv'
    # group: 'data_processing'
    run:
        import pysam

        import numpy as np
        import pandas as pd

        bam = pysam.AlignmentFile(input.fname, 'rb')
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
            'coverage/coverage.{accession}.csv',
            accession=accession_list)
    output:
        fname = 'results/coverage.csv',
        fname_stats = report('results/statistics.csv', caption='report/empty_caption.rst'),
        fname_lowquar = 'results/coverage_lowerquartile.csv',
        fname_median = 'results/coverage_median.csv',
        fname_upperquar = 'results/coverage_upperquartile.csv'
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
        fname = 'results/coverage.csv',
        fname_selection = 'results/selected_samples.csv'
    output:
        fname = report('plots/coverage_per_locus.pdf', caption='report/empty_caption.rst')
    resources:
        mem_mb = 5_000
    run:
        import pandas as pd

        import seaborn as sns
        import matplotlib.pyplot as plt

        # fix big boom on macOS
        import matplotlib
        matplotlib.use('Agg')

        # read data
        df = pd.read_csv(input.fname)
        sel = pd.read(input.fname_selection, squeeze=True)

        df_sub = df[sel]

        # compute data statistics
        df_stats = pd.DataFrame({
            'lower quartile': df_sub.quantile(q=.25, axis=1),
            'median': df_sub.quantile(q=.5, axis=1),
            'upper quartile': df_sub.quantile(q=.75, axis=1)
        })

        # plot data
        plt.figure(figsize=(8, 6))

        plt.plot(
            df_stats.index,
            df_stats['median'],
            label='median')
        plt.plot(
            df_stats.index,
            df_stats['lower quartile'],
            alpha=.5,
            label='lower quartile')
        plt.plot(
            df_stats.index,
            df_stats['upper quartile'],
            alpha=.5,
            label='upper quartile')

        plt.xlabel('Genomic position [bp]')
        plt.ylabel('Per base read count')
        # plt.yscale('log')

        plt.legend(loc='best')

        plt.tight_layout()
        plt.savefig(output.fname)


rule plot_coverage_per_sample:
    input:
        fname_lowquar = 'results/coverage_lowerquartile.csv',
        fname_median = 'results/coverage_median.csv',
        fname_upperquar = 'results/coverage_upperquartile.csv'
    output:
        fname = report('plots/coverage_per_sample.pdf', caption='report/empty_caption.rst')
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

        plt.figure(figsize=(8, 6))

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
        fname = 'results/sra_metadata.csv'
    run:
        import io
        import subprocess

        import pandas as pd
        from tqdm import tqdm

        df_list = []
        for accession in tqdm(accession_list):
            proc = subprocess.Popen(
                ['esearch', '-db', 'sra', '-query', accession],
                stdout=subprocess.PIPE)
            proc_stdout = subprocess.check_output(
                ['efetch', '-format', 'runinfo'],
                stdin=proc.stdout)
            proc.wait()

            buf = io.StringIO(proc_stdout.decode('utf-8'))
            df_cur = pd.read_csv(buf)

            df_sel = df_cur[df_cur['Run'] == accession]
            assert df_sel.shape[0] == 1

            df_list.append(df_sel)

        df = pd.concat(df_list)
        assert df.shape[0] == len(accession_list)
        df.to_csv(output.fname, index=False)


rule compute_additional_properties:
    input:
        marker_list = expand(
            'trimmed/{accession}.marker',
            accession=accession_list)
    output:
        fname = 'results/extra_properties.csv'
    resources:
        mem_mb = 10_000
    run:
        import os
        import glob

        import numpy as np
        import pandas as pd

        from Bio import SeqIO
        from tqdm import tqdm

        tmp = []
        for marker_fname in tqdm(input.marker_list):
            accession = os.path.splitext(os.path.basename(marker_fname))[0]

            fname_glob = marker_fname.replace('.marker', '*.fastq')
            available_files = glob.glob(fname_glob)

            fname = available_files[0]
            read_len = np.mean(
                [len(r.seq) for r in SeqIO.parse(fname, 'fastq')]
            ).astype(int)

            tmp.append({
                'accession': accession,
                'avg_read_length': read_len
            })

        pd.DataFrame(tmp).to_csv(output.fname, index=False)


rule assemble_final_dataframe:
    input:
        fname_lowquar = 'results/coverage_lowerquartile.csv',
        fname_median = 'results/coverage_median.csv',
        fname_upperquar = 'results/coverage_upperquartile.csv',
        fname_extra = 'results/extra_properties.csv',
        fname_meta = 'results/sra_metadata.csv'
    output:
        fname = report('results/final_dataframe.csv')
    run:
        import pandas as pd

        # read data
        df_lq = pd.read_csv(input.fname_lowquar, index_col='variable')
        df_mq = pd.read_csv(input.fname_median, index_col='variable')
        df_uq = pd.read_csv(input.fname_upperquar, index_col='variable')

        df_extra = pd.read_csv(input.fname_extra, index_col='accession')
        df_meta = pd.read_csv(input.fname_meta, index_col='Run')

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
        fname = 'results/final_dataframe.csv'
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

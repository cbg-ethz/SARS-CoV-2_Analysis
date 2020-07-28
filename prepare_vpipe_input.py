import os
import collections
from pathlib import Path

import pandas as pd

from tqdm import tqdm


def main(root_dir, max_batch_size=None):
    df = pd.read_csv(root_dir / 'results/selected_samples.csv')
    accession_list = df['accession'].tolist()

    df_final = pd.read_csv(root_dir / 'results/final_dataframe.csv', index_col=0)
    readlen_dict = df_final.to_dict()['avg_read_length']

    name_template = 'samples_{type_}_batch{batch:02}'
    dummy_date = '19700101'

    current_batch_id = collections.defaultdict(int)
    current_batch_size = collections.defaultdict(int)
    for accession in tqdm(accession_list):
        # gather fastq files
        available_files = list(root_dir.glob(f'data/{accession}*.fastq'))

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
        name = Path(name_template.format(
            type_=type_,
            batch=current_batch_id[type_]))
        target = name / accession / dummy_date / 'raw_data'

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
        with open(f'{name}.tsv', 'a') as fd:
            fd.write(f'{accession}\t{dummy_date}\t{readlen_dict[accession]}\n')

        # handle batches
        current_batch_size[type_] += 1
        if max_batch_size is not None and current_batch_size[type_] == max_batch_size:
            current_batch_size[type_] = 0
            current_batch_id[type_] += 1


if __name__ == '__main__':
    main(Path('pipeline_run_retrieval/'))

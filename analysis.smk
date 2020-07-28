configfile: 'analysis.config.yaml'
workdir: config['workdir']

localrules: all


rule all:
    input:
        'data/vcf_data.csv'


rule gather_vcf_files:
    input:
        vcf_dir = srcdir(config['input']['vcf_directory'])
    output:
        fname = 'data/vcf_data.csv'
    run:
        import os

        import numpy as np
        import pandas as pd

        import vcf

        tmp = []
        for entry in os.scandir(input.vcf_dir):
            if not entry.name.endswith('.vcf'):
                continue

            with open(entry.path) as fd:
                vcf_reader = vcf.Reader(fd)
                for record in vcf_reader:

                    tmp.append({
                        'vcf': entry.name.split('_')[1].split('.')[0],
                        'chromosome': record.CHROM,
                        'position': record.POS,
                        'reference': record.REF,
                        'variant': [v.sequence for v in record.ALT],
                        'frequency': round(np.mean(
                            [v for k, v in record.INFO.items() if k.startswith("Freq")]
                        ), 3),
                        'posterior': round(1 - 10**(-record.QUAL / 10), 3)
                    })

        df = pd.DataFrame(tmp)
        df.to_csv(output.fname, index=False)

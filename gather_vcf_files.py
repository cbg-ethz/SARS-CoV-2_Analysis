import glob
import shutil
from pathlib import Path


def main(root_dir, target_dir):
    target_dir.mkdir(parents=True, exist_ok=True)
    vcf_iter = root_dir.glob('*/*/variants/SNVs/snvs.vcf')

    for fname in vcf_iter:
        accession = str(fname).split('/')[1]
        name_new = fname.name.replace('.vcf', f'_{accession}.vcf')

        print(fname, target_dir / name_new)
        shutil.copy(fname, target_dir / name_new)


if __name__ == '__main__':
    main(Path('samples_PE_batch00/'), Path('vcf_PE'))

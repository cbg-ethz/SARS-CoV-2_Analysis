import shutil
from pathlib import Path


# prepare V-pipe environment
vpipe_dir = Path('../rules/V-pipe/')
shutil.copy('../../config/vpipe.config', './vpipe.config')

VPIPE_BASEDIR = Path(workflow.basedir) / 'V-pipe'
include: vpipe_dir / 'vpipe.snake'


rule all_butforrealthistime:
    input:
        all_files


rule fubar:
    output:
        touch('SRR11278091/visualization/index.html')


rule gather_vcf_files:
    input:
        fname_list = [fname for fname in all_files if fname[-8:] == 'snvs.vcf']
    output:
        outdir = directory('vcf_PE')
    run:
        import glob
        import shutil
        from pathlib import Path

        # parameters
        fname_list = input.fname_list
        target_dir = Path(output.outdir)

        # copy files
        target_dir.mkdir(parents=True, exist_ok=True)

        for fname in fname_list:
            accession = fname.split('/')[4]
            name_new = Path(fname).name.replace('.vcf', f'_{accession}.vcf')

            print(fname, target_dir / name_new)
            shutil.copy(fname, target_dir / name_new)

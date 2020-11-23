import shutil
from pathlib import Path


# helper functions
def setup_config(vpipe_config_path, config_mods, target_file):
    shutil.copy(vpipe_config_path, target_file)

    # read config template
    with open(target_file) as fd:
        txt = fd.read()

    # insert template arguments
    for pattern, value in config_mods.items():
        txt = txt.replace(f'{{{pattern}}}', str(value))

    # write new config
    with open(target_file, 'w') as fd:
        fd.write(txt)


# prepare V-pipe environment
vpipe_dir = Path('../rules/V-pipe/')
setup_config(
    config['vpipe_config_path'], config['config_mods'],
    './vpipe.config'
)

VPIPE_BASEDIR = Path(workflow.basedir) / 'V-pipe'
include: vpipe_dir / 'vpipe.snake'


rule all_butforrealthistime:
    input:
        all_files


rule gather_vcf_files:
    input:
        fname_list = [fname for fname in all_files if fname[-8:] == 'snvs.vcf']
    output:
        outdir = directory('vcf_files')
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
            accession = fname.split('/')[-5]
            name_new = Path(fname).name.replace('.vcf', f'_{accession}.vcf')

            print(fname, target_dir / name_new)
            shutil.copy(fname, target_dir / name_new)

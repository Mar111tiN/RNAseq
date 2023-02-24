from yaml import CLoader as Loader, load, dump
from subprocess import run

##### setup report #####
configfile: "config/config_sarah.yml"
# set the workdir
workdir: config['work_dir']
snakedir = os.path.dirname(workflow.snakefile)
scriptdir = os.path.join(snakedir, "scripts")

# include helper functions
include: "rules/utils.smk"
include: "rules/common.smk"


samples, units = get_sample_dfs(
    sample_sheet=config['samples'],
    unit_sheet=config['units']
    )


# load the sample independent config file if things get more complex
# config = add_config(config, config_name="general")
wildcard_constraints:
    # eg sample cannot contain _ or / to prevent ambiguous wildcards
    sample = "[^_/.]+",
    unit = "[^_/.]+",
    read = "[^/.]*",
    trim = "trim|raw"

##### load rules #####
include: "rules/qc.smk"
include: "rules/align.smk"
include: "rules/trim.smk"
include: "rules/rsem.smk"

rule all:
    input:
        "/fast/groups/ag_schmueck/work/ref/hg38/star_index/hg38_ens104CX"   # ,
        # "results/counts/all.tsv",
        # # "qc/multiqc_report.html",
        # "qc/fastQC.html",
        # expand("results/RSEM/{unit.sample_name}-{unit.unit_name}.genes.results",unit=units.itertuples()) 

from yaml import CLoader as Loader, load, dump
from subprocess import run

##### setup report #####
configfile: "config/config.yaml"
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


##### load rules #####
# include: "rules/qc.smk"
include: "rules/align.smk"
include: "rules/trim.smk"

rule all:
    input:
        "results/counts/all.tsv",
        # "results/qc/multiqc_report.html",
        # "results/qc/fastQC.html"

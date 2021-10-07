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


# load the sample independent config file
# config = add_config(config, config_name="general")

# how does the report work?
# report: "report/workflow.rst"


##### setup singularity #####
# this container defines the underlying OS for each job when using the workflow
# with --use-conda --use-singularity
# container: "docker://continuumio/miniconda3"


##### load rules #####


# include: "rules/trim.smk"
include: "rules/qc.smk"
include: "rules/align.smk"
include: "rules/diffexp.smk"


rule all:
    input:
        expand(
            "results/star/{sample}-{unit}/ReadsPerGene.out.tab", 
            sample=units['sample_name'], 
            unit=units['unit_name']
        ),
        "results/counts/all.tsv",
        # get_final_output(),
        # "results/qc/multiqc_report.html",
        "results/qc/fastQC.html"

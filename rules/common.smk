import glob
import os
import re
import pandas as pd
from snakemake.utils import validate


# validate(config, schema="../schemas/config.schema.yaml")


def get_sample_dfs(sample_sheet=".", unit_sheet="."):

    # check full path or append path to snakedir
    if not sample_sheet.startswith('/'):
        sample_sheet = os.path.join(snakedir, sample_sheet)
    
    samples = (
        pd.read_csv(sample_sheet, sep="\t", dtype={"sample_name": str})
        .set_index("sample_name", drop=False)
        .sort_index()
    )
    # validate(samples, schema="../schemas/samples.schema.yaml")

    if not unit_sheet.startswith('/'):
        unit_sheet = os.path.join(snakedir, unit_sheet)

    units = (
        pd.read_csv(unit_sheet, sep="\t", dtype={"sample_name": str, "unit_name": str})
        .set_index(["sample_name", "unit_name"], drop=False)
        .sort_index()
    )

    # validate(units, schema="../schemas/units.schema.yaml")
    return (samples, units)

def static_path():
    cr = config['ref']
    return os.path.join(cr['static_path'], config['ref']['build'])


def full_path(file):

    '''
    returns the full path to a reference if file is relative to ref/hg38|19
    '''

    cr = config['ref']
    build = config['ref']['build']
    return os.path.join(static_path(), cr[build][file])


def get_gtf():
    cr = config['ref']
    build = config['ref']['build']
    release = config['ref']['release']
    gtf_file = f"{build}_ens{release}.chr.gtf"
    gtf_folder = os.path.join(static_path(), config['ref'][build]['gtf_path'])
    return os.path.join(gtf_folder, gtf_file)
    

def star_index():
    cr = config['ref']
    build = config['ref']['build']
    release = config['ref']['release'] 

    index_main_folder = os.path.join(static_path(), config['ref'][build]['index_folder'])
    star_folder = f"{build}_ens{release}"
    
    return os.path.join(index_main_folder, star_folder)


def rseqc_gene_model(_type="bed"):
    cr = config['ref']
    build = config['ref']['build']
    release = config['ref']['release'] 

    index_main_folder = os.path.join(static_path(), config['ref'][build]['rseqc_folder'])
    gene_model = f"{build}_ens{release}.{_type}"
    
    return os.path.join(index_main_folder, gene_model)


def get_cutadapt_input(wildcards):
    unit = units.loc[wildcards.sample].loc[wildcards.unit]

    if unit["fq1"].endswith("gz"):
        ending = ".gz"
    else:
        ending = ""

    if pd.isna(unit["fq2"]):
        # single end local sample
        return "pipe/cutadapt/{S}/{U}.fq1.fastq{E}".format(
            S=unit.sample_name, U=unit.unit_name, E=ending
        )
    else:
        # paired end local sample
        return expand(
            "pipe/cutadapt/{S}/{U}.{{read}}.fastq{E}".format(
                S=unit.sample_name, U=unit.unit_name, E=ending
            ),
            read=["fq1", "fq2"],
        )


def is_paired_end(sample):
    sample_units = units.loc[sample]
    fq2_null = sample_units["fastq2"].isnull()
    paired = ~fq2_null
    all_paired = paired.all()
    all_single = (~paired).all()
    assert (
        all_single or all_paired
    ), "invalid units for sample {}, must be all paired end or all single end".format(
        sample
    )
    return all_paired


fastq_path = config['input_dir']

def get_raw_fastq(w):
    # no trimming, use raw reads
    u = units.loc[(w.sample, w.unit), ["fastq1", "fastq2"]].dropna()

    if not is_paired_end(wildcards.sample):
        return {"fastq1": f"{u.fq1}"}
    else:
        return {
            "fastq1": f"{os.path.join(fastq_path, u.fastq1)}", 
            "fastq2": f"{os.path.join(fastq_path, u.fastq2)}"}


def get_star_fastq(w):
    if w.trim == "trimmed":
        # activated trimming, use trimmed data
        if is_paired_end(w.sample):
            # paired-end sample
            return dict(
                zip(
                    ["fastq1", "fastq2"],
                    expand(
                        "results/trimmed/{sample}_{unit}_{read}.fastq.gz",
                        read=["R1", "R2"],
                        **w,
                    ),
                )
            )
        # single end sample
        return {"fastq1": "trimmed/{sample}_{unit}.fastq.gz".format(**w)}
    else:
        return get_raw_fastq(w)

def get_all_fastq
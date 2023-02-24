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
        # fillna for easy checking of adapters in get_adapters
        .fillna("")
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
    gtf_file = f"{build}_{release}.chr.gtf"
    gtf_folder = os.path.join(static_path(), config['ref'][build]['gtf_path'])
    return os.path.join(gtf_folder, gtf_file)
    

def star_index():
    cr = config['ref']
    build = config['ref']['build']
    release = config['ref']['release'] 

    index_main_folder = os.path.join(static_path(), config['ref'][build]['STARindex_folder'])
    star_folder = f"{build}_{release}"
    
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


############## fastq ###############

def get_raw_fastq(w):

    fastq_path = config['input_dir']
    # no trimming, use raw reads
    u = units.loc[(w.sample, w.unit), :]

    if u['fastq2']:
    # PE case
        return {
            "fastq1": os.path.join(fastq_path, u.fastq1), 
            "fastq2": os.path.join(fastq_path, u.fastq2)
            }
    # SE case
    return {"fastq1": os.path.join(fastq_path, u.fastq1)}


def get_trim_fastq(w):
    # get the respective entry from the units_df (use tuple indexing for dual row index)
    unit = units.loc[(w.sample, w.unit), :]

    if unit['fastq2']:
    # PE case
        return {f"fastq{read}": f"results/trimmed/{w.sample}-{w.unit}_R{read}.fastq.gz" for read in [1, 2]}
    # SE case
    return {f"fastq1": f"results/trimmed/{w.sample}-{w.unit}.fastq.gz"}


############ STAR ###########################

def get_star_fastq(w):
    '''
    returns the trimmed/untrimmed mates/single fastq depending on config and existence of second fastq in units_df
    '''

    # trimming is active
    if config['trimming']['activate']:
        # activated trimming, use trimmed data

        return get_trim_fastq(w)

    # trimming is inactive --> get org fastqs from units_df  
    return get_raw_fastq(w)


############# fastQC ########################

def get_fastqc_list(_):
    '''
    returns the complete list of required fastqc files depending on trim option
    '''
    
    # create file list from the included_files tuple list
    fastqc_folder = "qc/fastqc"
    types = ["raw"]
    # add the trim folder prefix if trimming is active
    if config['trimming']['activate']:
        types.append("trim")
    reads = ["R1", "R2"]
    # create SE marker for single read units
    units['SE'] = units['fastq2'].isna() | (units['fastq2'] == "")

    PE_list = [os.path.join(fastqc_folder, f"{u[0]}-{u[1]}_{read}_{t}_fastqc.zip") for u in units[~units['SE']].index for t in types for read in reads]
    SE_list = [os.path.join(fastqc_folder, f"{u[0]}-{u[1]}_{t}_fastqc.zip") for u in units[units['SE']].index for t in types]

    return SE_list + PE_list


def get_qc_fastq(w):
    '''
    returns the trimmed/untrimmed mates/single fastq depending on trim wildcard and existence of second fastq in units_df
    '''

    # trimmed fastq as input
    if w.trim == "trim":
        # is SE or PE (read can be "" for SE
        return f"results/trimmed/{w.sample}-{w.unit}{w.read}.fastq.gz"
    # is raw fastq --> get fastq path from 
    fastq_path = config['input_dir']
    # no trimming, use raw reads
    u = units.loc[(w.sample, w.unit), :]
    # is read2
    if w.read == "_R2":
        return os.path.join(fastq_path, u.fastq2)
    # raw fastq as input
    return os.path.join(fastq_path, u.fastq1)
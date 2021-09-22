## fastqc
rule fastqc:
    input: unpack(get_fq)
    output:
        fastqc1 = "results/qc/fastqc/{sample}_{unit}_1_fastqc.zip",
        fastqc2 = "results/qc/fastqc/{sample}_{unit}_2_fastqc.zip"
    threads: 4
    log:
        "logs/fastqc/{sample}_{unit}.log"
    conda:
        f"../envs/fastQC-env.yml"
    params:
        name1 = lambda w, output: os.path.basename(output.fastqc1).replace("_fastqc.zip", ""),
        name2 = lambda w, output: os.path.basename(output.fastqc2).replace("_fastqc.zip", "")
    shell:
        "zcat {input.fq1} | fastqc stdin:{params.name1} -o results/qc/fastqc/ && "  # &>{log} "    
        "zcat {input.fq2} | fastqc stdin:{params.name2} -o results/qc/fastqc/"  # &>{log} " 


def get_fastqc_list(_):
    '''
    returns the complete list of required fastqc files depending on trim option
    '''

    # create file list from the included_files tuple list
    fastqc_list = [f"results/qc/fastqc/{s}_{u}_{r}_fastqc.zip" for s in units['sample_name'] for u in units['unit_name'] for r in [1,2]]
    return fastqc_list


rule fastq_multiQC:
    input:
        get_fastqc_list
    output:
        "results/qc/fastQC.html"
    threads: 2
    conda:
        f"../envs/fastQC-env.yml"
    shell:
        "multiqc -f -o results/qc/ -n fastQC --interactive results/qc/fastqc/; "  # interactive for big number of files
        "rm -f results/qc/fastqc/*_fastqc.html results/qc/fastq/*.sub"  # leave the zip files for accumulated multiQC of all processed samples

## RSEQC
rule rseqc_gtf2db:
    input:
        get_gtf(),
    output:
        db=rseqc_gene_model("db"),
    log:
        "logs/rseqc_gtf2db.log",
    conda:
        "../envs/gffutils.yaml"
    script:
        "../scripts/gtf2db.py"


rule rseqc_gtfdb2bed:
    input:
        db=rseqc_gene_model("db")
    output:
        bed=rseqc_gene_model()
    log:
        "logs/rseqc_gtfdb2bed.log",
    conda:
        "../envs/gffutils.yaml"
    script:
        "../scripts/gtfdb2bed.py"


rule rseqc_junction_annotation:
    input:
        bam="results/star/{sample}-{unit}/Aligned.sortedByCoord.out.bam",
        bed=rseqc_gene_model(),
    output:
        "results/qc/rseqc/{sample}-{unit}.junctionanno.junction.bed",
    priority: 1
    log:
        "logs/rseqc/rseqc_junction_annotation/{sample}-{unit}.log",
    params:
        extra=r"-q 255",  # STAR uses 255 as a score for unique mappers
        prefix=lambda w, output: output[0].replace(".junction.bed", ""),
    conda:
        "../envs/rseqc.yaml"
    shell:
        "junction_annotation.py {params.extra} -i {input.bam} -r {input.bed} -o {params.prefix} "
        # "> {log[0]} 2>&1"


rule rseqc_junction_saturation:
    input:
        bam="results/star/{sample}-{unit}/Aligned.sortedByCoord.out.bam",
        bed=rseqc_gene_model(),
    output:
        "results/qc/rseqc/{sample}-{unit}.junctionsat.junctionSaturation_plot.pdf",
    priority: 1
    log:
        "logs/rseqc/rseqc_junction_saturation/{sample}-{unit}.log",
    params:
        extra=r"-q 255",
        prefix=lambda w, output: output[0].replace(".junctionSaturation_plot.pdf", ""),
    conda:
        "../envs/rseqc.yaml"
    shell:
        "junction_saturation.py {params.extra} -i {input.bam} -r {input.bed} -o {params.prefix} "
        "> {log} 2>&1"


rule rseqc_stat:
    input:
        "results/star/{sample}-{unit}/Aligned.sortedByCoord.out.bam",
    output:
        "results/qc/rseqc/{sample}-{unit}.stats.txt",
    priority: 1
    log:
        "logs/rseqc/rseqc_stat/{sample}-{unit}.log",
    conda:
        "../envs/rseqc.yaml"
    shell:
        "bam_stat.py -i {input} > {output} 2> {log}"


rule rseqc_infer:
    input:
        bam="results/star/{sample}-{unit}/Aligned.sortedByCoord.out.bam",
        bed=rseqc_gene_model(),
    output:
        "results/qc/rseqc/{sample}-{unit}.infer_experiment.txt",
    priority: 1
    log:
        "logs/rseqc/rseqc_infer/{sample}-{unit}.log",
    conda:
        "../envs/rseqc.yaml"
    shell:
        "infer_experiment.py -r {input.bed} -i {input.bam} > {output} 2> {log}"


rule rseqc_innerdis:
    input:
        bam="results/star/{sample}-{unit}/Aligned.sortedByCoord.out.bam",
        bed=rseqc_gene_model(),
    output:
        "results/qc/rseqc/{sample}-{unit}.inner_distance_freq.inner_distance.txt",
    priority: 1
    log:
        "logs/rseqc/rseqc_innerdis/{sample}-{unit}.log",
    params:
        prefix=lambda w, output: output[0].replace(".inner_distance.txt", ""),
    conda:
        "../envs/rseqc.yaml"
    shell:
        "inner_distance.py -r {input.bed} -i {input.bam} -o {params.prefix} > {log} 2>&1"


rule rseqc_readdis:
    input:
        bam="results/star/{sample}-{unit}/Aligned.sortedByCoord.out.bam",
        bed=rseqc_gene_model(),
    output:
        "results/qc/rseqc/{sample}-{unit}.readdistribution.txt",
    priority: 1
    log:
        "logs/rseqc/rseqc_readdis/{sample}-{unit}.log",
    conda:
        "../envs/rseqc.yaml"
    shell:
        "read_distribution.py -r {input.bed} -i {input.bam} > {output} 2> {log}"


rule rseqc_readdup:
    input:
        "results/star/{sample}-{unit}/Aligned.sortedByCoord.out.bam",
    output:
        "results/qc/rseqc/{sample}-{unit}.readdup.DupRate_plot.pdf",
    priority: 1
    log:
        "logs/rseqc/rseqc_readdup/{sample}-{unit}.log",
    params:
        prefix=lambda w, output: output[0].replace(".DupRate_plot.pdf", ""),
    conda:
        "../envs/rseqc.yaml"
    shell:
        "read_duplication.py -i {input} -o {params.prefix} > {log} 2>&1"


rule rseqc_readgc:
    input:
        "results/star/{sample}-{unit}/Aligned.sortedByCoord.out.bam",
    output:
        "results/qc/rseqc/{sample}-{unit}.readgc.GC_plot.pdf",
    priority: 1
    log:
        "logs/rseqc/rseqc_readgc/{sample}-{unit}.log",
    params:
        prefix=lambda w, output: output[0].replace(".GC_plot.pdf", ""),
    conda:
        "../envs/rseqc.yaml"
    shell:
        "read_GC.py -i {input} -o {params.prefix} > {log} 2>&1"


rule multiqc:
    input:
        expand(
            "results/star/{unit.sample_name}-{unit.unit_name}/Aligned.sortedByCoord.out.bam",
            unit=units.itertuples(),
        ),
        expand(
            "results/qc/rseqc/{unit.sample_name}-{unit.unit_name}.junctionanno.junction.bed",
            unit=units.itertuples(),
        ),
        expand(
            "results/qc/rseqc/{unit.sample_name}-{unit.unit_name}.junctionsat.junctionSaturation_plot.pdf",
            unit=units.itertuples(),
        ),
        expand(
            "results/qc/rseqc/{unit.sample_name}-{unit.unit_name}.infer_experiment.txt",
            unit=units.itertuples(),
        ),
        expand(
            "results/qc/rseqc/{unit.sample_name}-{unit.unit_name}.stats.txt",
            unit=units.itertuples(),
        ),
        expand(
            "results/qc/rseqc/{unit.sample_name}-{unit.unit_name}.inner_distance_freq.inner_distance.txt",
            unit=units.itertuples(),
        ),
        expand(
            "results/qc/rseqc/{unit.sample_name}-{unit.unit_name}.readdistribution.txt",
            unit=units.itertuples(),
        ),
        expand(
            "results/qc/rseqc/{unit.sample_name}-{unit.unit_name}.readdup.DupRate_plot.pdf",
            unit=units.itertuples(),
        ),
        expand(
            "results/qc/rseqc/{unit.sample_name}-{unit.unit_name}.readgc.GC_plot.pdf",
            unit=units.itertuples(),
        ),
        expand(
            "logs/rseqc/rseqc_junction_annotation/{unit.sample_name}-{unit.unit_name}.log",
            unit=units.itertuples(),
        ),
    output:
        "results/qc/multiqc_report.html",
    log:
        "logs/multiqc.log",
    wrapper:
        "0.75.0/bio/multiqc"

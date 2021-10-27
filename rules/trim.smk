def get_adapters(w):
    '''
    create the trimming params
    '''

    is_PE = units[]
    ac = config['trimming']['adapters']
    adapter_list = [f" -g {a}" for a in ac["5prime"]] + [f"-a {a}" for a in ac["3prime"]]
    if is_paired_end(w.sample):
        adapter_list += [f" -G {a}" for a in ac["5prime"]] + [f"-A {a}" for a in ac["3prime"]]
    adapter_string = " ".join(adapter_list) + " "

    return adapter_string


rule cutadapt_pe:
    input:
        get_raw_fastq
    output:
        fastq1="results/trimmed/{sample}-{unit}_R1.fastq.gz",
        fastq2="results/trimmed/{sample}-{unit}_R2.fastq.gz",
        qc="results/trimmed/{sample}-{unit}.paired.qc.txt",
    log:
        "logs/cutadapt/{sample}-{unit}.log",
    params:
        adapters=get_adapters
    threads: 
        config['trimming']['threads']
    conda:
        "../envs/cut-env.yml"
    shell:
        "cutadapt -j {threads} {params.adapters} -o {output.fastq1} -p {output.fastq2} {input}"


rule cutadapt_se:
    input:
        get_raw_fastq
    output:
        fastq="results/trimmed/{sample}-{unit}.fastq.gz",
        qc="results/trimmed/{sample}-{unit}.qc.txt",
    log:
        "logs/cutadapt/{sample}-{unit}.log",
    params:
        adapters=get_adapters
    threads: 
        config['trimming']['threads']
    conda:
        "../envs/cut-env.yml"
    shell:
        "cutadapt -j {threads} {params.adapters} -o {output.fastq1} {input}"

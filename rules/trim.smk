def get_adapters(w):
    '''
    create the trimming params
    '''

    # select the row from units_df
    u = units.loc[(w.sample, w.unit), :]

    ac = config['trimming']['adapters']

    # get the adapters from the config and the units_df
    adapt5 = ac["5prime"]
    if (a:=u['adapters5']):
        adapt5.append(a)
    adapt3 = ac["3prime"]
    if (a:=u['adapters3']):
        adapt3.append(a)
    adapter_list = [f"-g {a}" for a in adapt5] + [f"-a {a}" for a in adapt3]
    # check for PE
    if u["fastq2"]:
        adapter_list += [f"-G {a}" for a in adapt5] + [f"-A {a}" for a in adapt3]
    adapter_string = " ".join(adapter_list) + " "

    return adapter_string


rule cutadapt_PE:
    input:
        unpack(get_raw_fastq)
    output:
        fastq1="results/trimmed/{sample}-{unit}_R1.fastq.gz",
        fastq2="results/trimmed/{sample}-{unit}_R2.fastq.gz"
    log:
        "logs/cutadapt/{sample}-{unit}.log",
    params:
        adapters=get_adapters
    threads: 
        config['trimming']['threads']
    conda:
        "../envs/cut-env.yml"
    shell:
        "cutadapt -j {threads} {params.adapters} --minimum-length 1 -q 20 "
        "-o {output.fastq1} -p {output.fastq2} {input.fastq1} {input.fastq2} > {log}"


rule cutadapt_SE:
    input:
        unpack(get_raw_fastq)
    output:
        fastq="results/trimmed/{sample}-{unit}.fastq.gz"
    log:
        "logs/cutadapt/{sample}-{unit}.log",
    params:
        adapters=get_adapters
    threads: 
        config['trimming']['threads']
    conda:
        "../envs/cut-env.yml"
    shell:
        "cutadapt -j {threads} {params.adapters} --minimum-length 1 -q 20 "
        "-o {output.fastq} {input.fastq1} > {log}"

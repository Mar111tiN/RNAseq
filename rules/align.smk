rule star_index:
    input:
        fasta=full_path("genome"),
        gtf=get_gtf(),
    output:
        directory(star_index()),
    threads: 8
    log:
        "logs/star_index_genome.log"
    params:
        overhang=config['length'],
        gtf=get_gtf()
    conda:
        "../envs/star.yml"
    shell:
        "mkdir -p {output} && "
        "STAR --runMode genomeGenerate "
        "--runThreadN {threads} "
        "--genomeDir {output} "
        "--genomeFastaFiles {input.fasta} "
        "--sjdbOverhang {params.overhang} "
        "--sjdbGTFfile {input.gtf} "
        "--outStd Log "
        "&> {log}"
        
rule star_align:
    input:
        unpack(get_fq),
        index=star_index(),
        gtf=get_gtf()
    output:
        bam = "results/star/{sample}-{unit}/Aligned.sortedByCoord.out.bam",
        tab = "results/star/{sample}-{unit}/ReadsPerGene.out.tab",
    log:
        "logs/star/{sample}-{unit}.log",
    params:
        outprefix = lambda wc, output: os.path.dirname(output.bam) + "/",
    threads: 24
    conda:
        "../envs/star.yml"
    shell:
        "STAR --runThreadN {threads} "
        "--genomeDir {input.index} "
        "--readFilesCommand zcat "
        "--readFilesIn {input.fq1} {input.fq2} "
        "--outSAMtype BAM SortedByCoordinate "
        "--quantMode GeneCounts "
        "--sjdbGTFfile {input.gtf} "
        "--outFileNamePrefix {params.outprefix} "
        "--outStd Log "
        "&> {log}; "
        "samtools index {output.bam}"

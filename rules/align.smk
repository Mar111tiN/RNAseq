rule star_index:
    input:
        fasta=full_path("genome"),
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
        "--sjdbGTFfile {params.gtf} "
        "--outStd Log "
        "&> {log}"


rule star_align:
    input:
        unpack(get_star_fastq),
        index=star_index(),
    output:
        bam = "results/star/{sample}-{unit}/Aligned.sortedByCoord.out.bam",
        tbam = "results/star/{sample}-{unit}/Aligned.toTranscriptome.out.bam",
        tab = "results/counts/{sample}-{unit}_ReadsPerGene.out.tab",
    log:
        "logs/star/{sample}-{unit}.log",
    params:
        outprefix = lambda wc, output: os.path.dirname(output.bam) + "/",
        fastqs = lambda w, input: f"{input.fastq1} {input.fastq2}" if hasattr(input, "fastq2") else input.fastq1,
        gtf = get_gtf(),
        tab = lambda w, output: output.bam.replace("Aligned.sortedByCoord.out.bam", "ReadsPerGene.out.tab")
    threads: 24
    conda:
        "../envs/star.yml"
    shell:
        "STAR --runThreadN {threads} "
        "--genomeDir {input.index} "
        "--readFilesCommand zcat "
        "--readFilesIn {params.fastqs} "
        "--outSAMtype BAM SortedByCoordinate "
        "--quantMode GeneCounts TranscriptomeSAM "
        "--sjdbGTFfile {params.gtf} "
	    "--outSAMunmapped Within "
	    "--outFilterType BySJout "
	    "--outSAMattributes NH HI AS NM MD "
	    "--alignIntronMin 20 "
	    "--alignIntronMax 1000000 "
	    "--alignMatesGapMax 1000000 "
	    "--alignSJoverhangMin 8 "
	    "--alignSJDBoverhangMin 1 " 
        "--outFileNamePrefix {params.outprefix} "
        "--outStd Log && "
        "samtools index {output.bam} && "
        # "samtools index {output.tbam} && " # samtools cannot index on transcript coords
        "mv {params.tab} {output.tab}"


rule count_matrix:
    input:
        expand(
            "results/counts/{unit.sample_name}-{unit.unit_name}_ReadsPerGene.out.tab",
            unit=units.itertuples(),
        ),
    output:
        "results/counts/all.tsv",
    log:
        "logs/count-matrix.log",
    params:
        samples=units["sample_name"].tolist()
    conda:
        "../envs/pandas.yaml"
    script:
        "../scripts/count-matrix.py"

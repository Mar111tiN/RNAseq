rule rsem:
    input:
        tbam = "results/star/{sample}-{unit}/Aligned.toTranscriptome.out.bam"
    output:
        genes = "results/RSEM/{sample}-{unit}.genes.results",
        isoforms = "results/RSEM/{sample}-{unit}.isoforms.results"
    params:
        rsem_ref = full_path('rsem_ref')
    threads: 16
    conda:
        "../envs/RSEM-env.yml"
    shell:
        "rsem-calculate-expression "
        "--bam --num-threads {threads} "
        "--paired-end --forward-prob 0.5 "
        "--no-bam-output "
        "{input.tbam} {params.rsem_ref} results/RSEM/{wildcards.sample}-{wildcards.unit}"

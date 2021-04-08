rule subset_chrM:
    input:
        "results/alignment/recal/{sample}.bqsr.bam",
    output:
        "results/sampleid/chrM/{sample}.chrM.bam",
    conda:
        "../envs/samtools.yaml",
    params:
        param='-bsh',
        chr='chrM',
    shell:
        "samtools view {params.param} {input} {chr} > {output}"

rule chrM_index:
    input:
        "results/sampleid/chrM/{sample}.chrM.bam",
    output:
        "results/sampleid/chrM/{sample}.chrM.bam.bai",
    params:
        "" # optional params string
    wrapper:
        "0.73.0/bio/samtools/index"

rule chrM_haplotype_caller:
    input:
        bam="results/sampleid/chrM/{sample}.chrM.bam",
        bai="results/sampleid/chrM/{sample}.chrM.bam.bai",
        ref=config['common']['genome']
    output:
        gvcf="results/sampleid/chrM/{sample}.chrM.vcf",
    log:
        "logs/gatk/haplotypecaller/{sample}.chrM.log"
    params:
        extra="",  # optional
        java_opts="", # optional
    resources:
        mem_mb=8192
    wrapper:
        "0.73.0/bio/gatk/haplotypecaller"

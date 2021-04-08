rule chrM_fa:
    input:
        fa=config['common']['genome'],
    output:
        "resources/chrM.fa",
    conda:
        "../envs/samtools.yaml",
    shell:
        "samtools faidx {input.fa} chrM > {output}"

rule chrM_faidx:
    input:
        "resources/chrM.fa",
    output:
        "resources/chrM.fa.fai",
    conda:
        "../envs/samtools.yaml",
    shell:
        "samtools faidx {input}"

rule chrM_dict:
    input:
        "resources/chrM.fa",
    output:
        "resources/chrM.dict",
    conda:
        "../envs/gatk.yaml",
    shell:
        "gatk CreateSequenceDictionary -R {input}"

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
        "samtools view {params.param} {input} {params.chr} > {output}"

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
        ref="resources/chrM.fa",
        refidx="resources/chrM.fa.fai",
        refdict="resources/chrM.dict",
    output:
        gvcf="results/sampleid/chrM/{sample}.chrM.gvcf",
    log:
        "logs/gatk/haplotypecaller/{sample}.gvcf.chrM.log"
    params:
        extra="",  # optional
        java_opts="", # optional
    resources:
        mem_mb=2048
    wrapper:
        "0.73.0/bio/gatk/haplotypecaller"

rule genotype_gvcfs:
    input:
        gvcf="results/sampleid/chrM/{sample}.chrM.gvcf",
        ref="resources/chrM.fa",
    output:
        vcf="results/sampleid/chrM/{sample}.chrM.vcf",
    log:
        "logs/gatk/haplotypecaller/{sample}.chrM.log",
    params:
        extra="--allow-old-rms-mapping-quality-annotation-data",  # optional
        java_opts="", # optional
    resources:
        mem_mb=2048
    wrapper:
        "0.73.0/bio/gatk/genotypegvcfs"

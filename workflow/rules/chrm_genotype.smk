rule chrM_fa:
    input:
        fa=config['common']['genome'],
    output:
        fa="resources/chrM.fa",
        idx="resources/chrM.fa.fai",
    conda:
        "../envs/samtools.yaml",
    shell:
        "samtools faidx {input.fa} chrM > {output.fa}; "
        "samtools faidx {output.fa}; "

rule chrM_dict:
    input:
        "resources/chrM.fa",
    output:
        dict="resources/chrM.dict",
        bed="resources/chrM.bed",
    conda:
        "../envs/gatk.yaml",
    shell:
        "gatk CreateSequenceDictionary -R {input} ; "
        "echo -e 'chrM\t1\t'$(grep 'chrM' resources/chrM.dict  | cut -f3 | sed 's/LN://') > {output.bed}"

rule subset_chrM:
    input:
        "results/alignment/recal/{sample}.bqsr.bam",
    output:
        bam="results/sampleid/chrM/{sample}.chrM.bam",
    conda:
        "../envs/samtools.yaml",
    params:
        param='-bsh',
        chr='chrM',
    shell:
        "samtools view {params.param} {input} {params.chr} > {output.bam} ; "

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
        ref=rules.chrM_fa.output.fa,
        refidx=rules.chrM_fa.output.fa,
        refdict=rules.chrM_dict.output.dict,
        bed=rules.chrM_dict.output.bed,
    output:
        gvcf="results/sampleid/chrM/{sample}.chrM.gvcf",
    log:
        "logs/gatk/haplotypecaller/{sample}.gvcf.chrM.log"
    params:
        extra="-L resources/chrM.bed",  # optional
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

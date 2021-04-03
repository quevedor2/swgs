rule mark_duplicates:
  input:
    "alignment/mapped_reads/{sample}.sorted.bam"
  output:
    bam="alignment/dedup/{sample}.bam",
    metrics="alignment/dedup/{sample}.metric.txt"
  log:
    "logs/picard/dedup/{sample}.log"
  threads: 8
  params:
    "REMOVE_DUPLICATES=true",
    "ASSUME_SORT_ORDER='coordinate'"
  wrapper:
    "0.73.0/bio/picard/markduplicates"

rule index_duplicates:
  input:
    "alignment/dedup/{sample}.bam"
  output:
    "alignment/dedup/{sample}.bam.bai"
  log:
    "logs/samtools/index/{sample}.log"
  wrapper:
    "0.73.0/bio/samtools/index"

rule realigner_target_creator:
    input:
        bam="alignment/dedup/{sample}.bam",
        bai=rules.index_duplicates.output,
        ref=config['common']['genome'],
        known=get_variant_paths(variant='indel'),
    output:
        intervals="alignment/realign/{sample}.intervals",
        java_temp=temp(directory("gatk3_indelrealigner/{sample}")),
    log:
        "logs/gatk/indelrealigner/{sample}.realignertargetcreator.log",
    params:
        extra="", # optional
    resources:
        mem_mb=8192,
    threads: 8
    wrapper:
        "0.73.0/bio/gatk3/realignertargetcreator"

rule indelrealigner:
    input:
        bam="alignment/dedup/{sample}.bam",
        bai="alignment/dedup/{sample}.bam.bai",
        ref=config['common']['genome'],
        known=get_variant_paths(variant='indel'),
        target_intervals="alignment/realign/{sample}.intervals"
    output:
        bam="alignment/realign/{sample}.bam",
        bai="alignment/realign/{sample}.bai",
        java_temp=temp(directory("/tmp/gatk3_indelrealigner/{sample}")),
    log:
        "logs/gatk/indelrealigner/{sample}.log"
    params:
        extra=""  # optional
    threads: 8
    resources:
        mem_mb = 8192
    wrapper:
        "0.73.0/bio/gatk3/indelrealigner"

rule baserecalibrator:
    input:
        bam="alignment/realign/{sample}.bam",
        ref=config['common']['genome'],
        known=get_variant_paths(variant='snp')
    output:
        "alignment/recal/{sample}.recal_data_table"
    log:
        "logs/gatk/bqsr/{sample}.log"
    params:
        **config["params"]["gatk"]["baserecalibrator"],
        extra=""  # optional
    resources:
        mem_mb = 8192
    threads: 8
    wrapper:
        "0.73.0/bio/gatk3/baserecalibrator"

rule printreads:
    input:
        bam="alignment/realign/{sample}.bam",
        ref=config['common']['genome'],
        recal_data="alignment/recal/{sample}.recal_data_table"
    output:
        "alignment/recal/{sample}.bqsr.bam"
    log:
        "logs/gatk/bqsr/{sample}.log"
    params:
        **config["params"]["gatk"]["printreads"],
        extra=""  # optional
    resources:
        mem_mb = 8192
    threads: 8
    wrapper:
        "0.73.0/bio/gatk3/printreads"

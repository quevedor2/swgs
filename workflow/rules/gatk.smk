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

rule realigner_target_creator:
    input:
        bam="alignment/dedup/{sample}.bam",
        ref=config['common']['genome'],
        known=get_indel_paths,
    output:
        intervals="alignment/realign/{sample}.intervals",
        java_temp=temp(directory("gatk3_indelrealigner/{sample}")),
    log:
        "logs/gatk/realigner_target_creator/{sample}.log",
    params:
        extra="", # optional
    resources:
        mem_mb=1024,
    threads: 8
    wrapper:
        "0.73.0/bio/gatk3/realignertargetcreator"

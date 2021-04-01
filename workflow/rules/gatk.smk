rule mark_duplicates:
  input:
    "alignment/mapped_reads/{sample}.bam"
  output:
    out="alignment/mapped_reads/{sample}.dup.bam",
    metrics="alignment/mapped_reads/{sample}.dup_metric.txt"
  conda:
    "../envs/gatk.yaml"
  threads: 4
  shell:
    "picard MarkDuplicates "
    "ASSUME_SORTED=TRUE "
    "INPUT={input} "
    "METRICS_FILE={output.metrics} "
    "OUTPUT={output.out} "

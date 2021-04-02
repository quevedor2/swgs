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

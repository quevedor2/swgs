rule bwa_mem:
  input:
    reads=get_fastqs
  output:
    "alignment/mapped_reads/{sample}.sorted.bam"
  log:
    "logs/bwa_mem/{sample}.log"
  threads: 8
  params:
    index=config['common']['genome'],
    extra=r"-R '@RG\tID:{sample}\tSM:{sample}'",
    sort="samtools",
    sort_order="coordinate"
  wrapper:
    "0.73.0/bio/bwa/mem"

rule samtools_index:
  input:
    "alignment/mapped_reads/{sample}.sorted.bam"
  output:
    "alignment/mapped_reads/{sample}.sorted.bam.bai"
  params:
    ""
  wrapper:
    "0.73.0/bio/samtools/index"

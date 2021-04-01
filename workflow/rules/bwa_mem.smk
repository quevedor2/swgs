rule bwa_mem:
  input:
    config['bwa']['index'],
    get_fastqs
  output:
    "alignment/mapped_reads/{sample}.bam"
  conda:
    "../envs/bwa.yaml"
  threads: 4
  shell:
    "bwa mem {input} | "
    "samtools view -bhS - | "
    "samtools sort -@4 - "

rule samtools_index:
  input:
    "alignment/mapped_reads/{sample}.bam"
  conda:
    "../envs/bwa.yaml"
  shell:
    "samtools index {input}"

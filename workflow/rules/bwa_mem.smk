rule bwa_mem:
  input:
    config['bwa']['index'],
    get_fastqs
  output:
    "alignment/mapped_reads/{sample}.bam"
  conda:
    "../envs/bwa.yaml"
  threads: 4
  params: 
    rgid=get_rgid
  shell:
    "bwa mem -R {params.rgid} {input} | "
    "samtools view -bhS - | "
    "samtools sort -@4 - > {output}"

rule samtools_index:
  input:
    "alignment/mapped_reads/{sample}.bam"
  output:
    "alignment/mapped_reads/{sample}.bam.bai"
  conda:
    "../envs/bwa.yaml"
  shell:
    "samtools index {input}"

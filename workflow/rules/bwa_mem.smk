rule bwa_mem:
  input:
    reads=get_fastqs
  output:
    "results/alignment/mapped_reads/{sample}.sorted.bam"
  log:
    "logs/bwa_mem/{sample}.log"
  threads: 8
  params:
    conda=config['env']['conda_shell'],
    env=directory(config['env']['preprocess']),
    index=config['params']['bwa']['mem']['genome'],
    extra=get_rgid,
  shell:
    """
    source {params.conda} && conda activate {params.env};
    
    bwa mem \
    -t {threads} \
    {params.extra} \
    {params.index} \
    {input.reads} | \
    samtools sort -o {output} - 2> {log}
    """

rule samtools_index:
  input:
    "results/alignment/mapped_reads/{sample}.sorted.bam"
  output:
    "results/alignment/mapped_reads/{sample}.sorted.bam.bai"
  threads: 1
  params:
    conda=config['env']['conda_shell'],
    env=directory(config['env']['preprocess']),
    extra="",
  shell:
    """
    source {params.conda} && conda activate {params.env};
    
    samtools index \
    {threads} \
    {params.extra} \
    {input} \
    {output} 2> {log}
    """

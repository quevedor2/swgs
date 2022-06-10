rule mark_duplicates:
  input:
    "results/alignment/mapped_reads/{sample}.sorted.bam"
  output:
    bam="results/alignment/dedup/{sample}.bam",
    metrics="results/alignment/dedup/{sample}.metric.txt"
  log:
    "logs/picard/dedup/{sample}.log"
  params:
    conda=config['env']['conda_shell'],
    env=directory(config['env']['preprocess']),
    mem=config['mem']['markduplicates'],
    rmduplicates="true",
    sort='coordinate',
  shell:
    """
    source {params.conda} && conda activate {params.env};
    
    java  \
    -Xmx{params.mem}G \
    -jar {params.env}/share/picard-2.27.2-0/picard.jar \
    MarkDuplicates \
    --REMOVE_DUPLICATES {params.rmduplicates} \
    --ASSUME_SORT_ORDER {params.sort} \
    -I {input} \
    -O {output.bam} \
    --METRICS_FILE {output.metrics} 2> {log}
    """
    

rule index_duplicates:
  input:
    "results/alignment/dedup/{sample}.bam"
  output:
    "results/alignment/dedup/{sample}.bam.bai"
  log:
    "logs/samtools/index/{sample}.log"
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

rule realigner_target_creator:
    input:
        bam="results/alignment/dedup/{sample}.bam",
        bai=rules.index_duplicates.output,
        ref=config['common']['genome'],
    output:
        intervals="results/alignment/realign/{sample}.intervals",
        java_temp=temp(directory("gatk3_indelrealigner/{sample}")),
    log:
        "logs/gatk/indelrealigner/{sample}.realignertargetcreator.log",
    params:
        conda=config['env']['conda_shell'],
        env=directory(config['env']['preprocess']),
        extra="", # optional (e.g. -L bedfile.bed)
        known=get_snp_paths,
    resources:
        mem_mb=8192,
    threads: 8
    shell:
        """
        source {params.conda} && conda activate {params.env};
        
        gatk3 Xmx{resources.mem_mb}M \
        -T RealignerTargetCreator \
        -nt {threads} \
        {params.extra} \
        -I {input.bam} \
        -R {input.ref} \
        -known {params.known} \
        -o {output.intervals} 2> {log}
        """

rule indelrealigner:
    input:
        bam="results/alignment/dedup/{sample}.bam",
        bai="results/alignment/dedup/{sample}.bam.bai",
        ref=config['common']['genome'],
        target_intervals="results/alignment/realign/{sample}.intervals"
    output:
        bam="results/alignment/realign/{sample}.bam",
        bai="results/alignment/realign/{sample}.bai",
#        java_temp=temp(directory("/tmp/gatk3_indelrealigner/{sample}")),
    log:
        "logs/gatk/indelrealigner/{sample}.log"
    params:
        conda=config['env']['conda_shell'],
        env=directory(config['env']['preprocess']),
        extra="",  # optional
        known=get_indel_paths,
    threads: 8
    resources:
        mem_mb = 8192
    shell:
        """
        source {params.conda} && conda activate {params.env};
        
        gatk3 -Xmx{resources.mem_mb}M \
        -T IndelRealigner \
        {params.extra} \
        -I {input.bam} \
        -R {input.ref} \
        -known {params.known} \
        --targetIntervals {input.target_intervals} \
        -o {output.bam} 2> {log}
        """

rule baserecalibrator:
    input:
        #bam="results/alignment/realign/{sample}.bam",
        bam="results/alignment/dedup/{sample}.bam",
        ref=config['common']['genome'],
    output:
        "results/alignment/recal/{sample}.recal_data_table"
    log:
        "logs/gatk/bqsr/{sample}.recal.log",
    params:
        conda=config['env']['conda_shell'],
        env=directory(config['env']['preprocess']),
        extra=combine_args(config["params"]["gatk"]["baserecalibrator"]),
        known=get_indel_paths,
    resources:
        mem_mb = 8192
    threads: 8
    shell:
        """
        source {params.conda} && conda activate {params.env};
        
        gatk \
        BaseRecalibrator \
        {params.extra} \
        -I {input.bam} \
        -R {input.ref} \
        --knownSites {params.known} \
        -O {output} 2> {log}
        """
        

rule applybqsr:
    input:
        bam="results/alignment/realign/{sample}.bam",
        ref=config['common']['genome'],
        recal_data="results/alignment/recal/{sample}.recal_data_table"
    output:
        "results/alignment/recal/{sample}.bqsr.bam"
    log:
        "logs/gatk/bqsr/{sample}.print.log"
    params:
        conda=config['env']['conda_shell'],
        env=directory(config['env']['preprocess']),
        extra=combine_args(config["params"]["gatk"]["printreads"]),
    resources:
        mem_mb = 8192
    threads: 8
    shell:
        """
        source {params.conda} && conda activate {params.env};
        
        gatk \
        ApplyBQSR \
        {params.extra} \
        -I {input.bam} \
        -R {input.ref} \
        --bqsr-recal-file {input.recal_data} \
        -o {output} 2> {log}
        """

rule symlink_bai:
    input:
        "results/alignment/recal/{sample}.bqsr.bam",
    params:
        "results/alignment/recal/{sample}.bqsr.bai",
    output:
        "results/alignment/recal/{sample}.bqsr.bam.bai",
    shell:
        "ln -s $(readlink -f {params}) {output}"

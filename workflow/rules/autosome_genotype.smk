rule collectAllelicCounts:
  input:
    "results/alignment/recal/{sample}.bqsr.bam",
  output:
    temp("results/zygosity/counts/{sample}.allelicCounts.tsv"),
  conda:
    "../envs/gatk.yaml",
  log:
    "logs/gatk/collectalleliccounts/{sample}.log"
  params:
    ref=config['common']['genome'],
    target=config['params']['gatk']['collectalleliccounts']['target'],
  shell:
    "gatk --java-options '-Xmx8G -XX:ParallelGCThreads=10' CollectAllelicCounts "
    "-I {input} "
    "-R {params.ref} "
    "-L {params.target} "
    "-O {output} > {log} 2>&1"

rule categorizeAD_GATK:
  input:
    "results/zygosity/counts/{sample}.allelicCounts.tsv",
  output:
    intermediate=temp("results/zygosity/counts/{sample}_out.tmp"),
    simple=temp("results/zygosity/counts/{sample}_out.tsv"),
  params:
    ref=2,
    alt=3,
  conda:
    "../envs/perl.yaml",
  shell:
    "grep -v '^@' {input} | tail -n +2 > {output.intermediate}; "
    "perl workflow/scripts/allelic_count_helper.pl categorize "
    "{output.intermediate} {params.ref} {params.alt} > {output.simple}"

rule getHetSNPs:
  input:
    "results/zygosity/counts/{sample}_out.tsv"
  output:
    ""
  conda:
    ""
  params:
    ""
  shell:
    

rule aggregateAD:
  input:
    expand("results/zygosity/counts/{sample}_out.tsv", sample=samples.index),
  output:
    "results/zygosity/AD/aggregate.csv",
  params:
    ""
  conda:
    "../envs/perl.yaml",
  shell:
    "paste -d',' {input} > {output}; "

rule getADlines:
  input:
    "results/zygosity/AD/aggregate.csv",
  output:
    "results/zygosity/AD/aggregate_lines.txt",
  params:
    min_n=2,
  conda:
    "../envs/perl.yaml",
  shell:
    "perl workflow/scripts/allelic_count_helper.pl setlines "
    "{input} {params.min_n} > {output}; "

rule filtADraw:
  input:
    filt_lines="results/zygosity/AD/aggregate_lines.txt",
    aggregate_raw="results/zygosity/AD/aggregate.csv",
  output:
    "results/zygosity/AD/aggregate_filt.csv",
  params:
    "",
  conda:
    "../envs/perl.yaml",
  shell:
    "perl workflow/scripts/allelic_count_helper.pl getlines "
    "{input.filt_lines} {input.aggregate_raw} > {output}; "

rule filtADdbsnp:
  input:
    filt_lines="results/zygosity/AD/aggregate_lines.txt",
  output:
    "results/zygosity/AD/aggregate_pos.txt",
  params:
    target=config['params']['gatk']['collectalleliccounts']['target'],
  conda:
    "../envs/perl.yaml",
  shell:
    "perl workflow/scripts/allelic_count_helper.pl getlines "
    "{input.filt_lines} {params.target} > {output}; "

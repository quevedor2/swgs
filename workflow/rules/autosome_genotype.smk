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
    "gatk CollectAllelicCounts "
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
    "perl workflow/scripts/categorizeAD_GATK.pl "
    "{output.intermediate} {params.ref} {params.alt} > {output.simple}"

rule aggregateAD:
  input:
    expand("results/zygosity/counts/{sample}_out.tsv", sample=samples.index)
  output:
    aggregate_raw="results/zygosity/AD/aggregate.tsv",
    aggregate_filt="results/zygosity/AD/aggregate_filt.tsv",
    filt_lines="results/zygosity/AD/aggregate_lines.txt",
    filt_pos="results/zygosity/AD/aggregate_pos.txt",
  params:
    target=config['params']['gatk']['collectalleliccounts']['target'],
  conda:
    "../envs/perl.yaml",
  shell:
    "paste -d',' {input} > {output.aggregate_raw}; "
    "awk -F\',\' \'{ for(i=1; i<=NF;i++) if ($i!=0) j++; if (j >= {params.n}) print NR; j=0 }\' {output.aggregate_raw} > {output.filt_lines}; "
    "perl workflow/scripts/getLinesFromFile.py "
    "{output.filt_lines} {output.aggregate_raw} > {output.aggregate_filt}; "
    "perl workflow/scripts/getLinesFromFile.py "
    "{output.filt_lines} {params.target} > {output.filt_pos}; "

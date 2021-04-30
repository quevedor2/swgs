rule collectAllelicCounts:
  input:
    "results/alignment/recal/{sample}.bqsr.bam",
  output:
    temp("results/sampleid/autosome/{sample}.allelicCounts.tsv"),
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
    "results/sampleid/autosome/{sample}.allelicCounts.tsv",
  output:
    intermediate=temp("results/sampleid/autosome/simple/{sample}_out.tmp"),
    simple="results/sampleid/autosome/simple/{sample}_out.tsv",
  params:
    ref=2,
    alt=3,
  conda:
    "../envs/perl.yaml",
  shell:
    "grep -v '^@' {input} | tail -n +2 > {output.intermediate}; "
    "perl workflow/scripts/categorizeAD_GATK.pl "
    "{output.intermediate} {params.ref} {params.alt} > {output.simple}"

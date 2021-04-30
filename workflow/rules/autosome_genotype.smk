rule collectAllelicCounts:
  input:
    "results/alignment/recal/{sample}.bqsr.bam",
  output:
    temp("results/sampleid/autosome/{sample}.allelicCounts.tsv"),
  conda:
    "../envs/gatk.yaml",
  params:
    ref=params['common']['genome'],
    target=params['params']['gatk']['collectalleliccounts']['target'],
  shell:
    "gatk CollectAllelicCounts "
    "-I {input} "
    "-R {params.ref} "
    "-L {params.target} "
    "-O {output}"

rule categorizeAD_GATK:
  input:
    "results/sampleid/autosome/{sample}.allelicCounts.tsv",
  output:
    intermediate=temp("results/sampleid/autosome/simple/{sample}_out.tmp"),
    simple=temp("results/sampleid/autosome/simple/{sample}_out.tsv"),
  params:
    ref=2,
    alt=3,
  conda:
    "../envs/perl.yaml",
  shell:
    "grep -v '^@' {input} | tail -n +2 > {output.intermediate}"
    "perl workflow/scripts/categorizeAD_GATK.pl "
    "{output.intermediate} {params.ref} {params.alt} > {output.simple}"

rule isolate_AD:
    input:
        fa=config['common']['genome'],
    output:
        fa="resources/chrM.fa",
        idx="resources/chrM.fa.fai",
    conda:
        "../envs/samtools.yaml",
    shell:
        "samtools faidx {input.fa} chrM > {output.fa}; "
        "samtools faidx {output.fa}; "

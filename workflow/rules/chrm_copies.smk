rule chrM_make_interval:
    input:
        bed=rules.chrM_dict.output.bed,
        dict=rules.chrM_dict.output.dict,
    output:
        bed="resources/chrM.headered.bed",
    conda:
        "../envs/picard.yaml"
    shell:
        "picard BedToIntervalList "
        "I={input.bed} "
        "O={output.bed} "
        "SD={input.dict} "

rule collect_chr_metrics:
    input:
        bam=rules.subset_chrM.output.bam,
        ref=rules.chrM_fa.output.fa,
        bed=rules.chrM_make_interval.output.bed,
    output:
        "results/chrM/copies/{sample}.wgs_metrics"
    log:
        "logs/picard/wgs_metrics/{sample}_chrM.log"
    conda:
        "../envs/picard.yaml"
    shell:
        "picard CollectWgsMetrics "
        "I={input.bam} "
        "O={output} "
        "R={input.ref} "
        "INTERVALS={input.bed} 2> {log}"

'''
rule estimate_mt_copies:
    input:
        chrm=rules.estimate_mt_copies.output,
        genome=rules.estimate_mt_copies.output,
    output:
        "results/chrM/copies/{sample}.wgs_metrics"
    log:
        "logs/picard/wgs_metrics/{sample}_chrM.log"
    conda:
        "../envs/picard.yaml"
    shell:
        "picard CollectWgsMetrics "
        "I={input.bam} "
        "O={output} "
        "R={input.ref} "
        "INTERVALS={input.bed} 2> {log}"
'''

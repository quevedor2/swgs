rule get_chromosomes:
    input:
        "results/alignment/recal/{sample}.bqsr.bam",
    output:
        "results/cnv/ichorcna/{sample}.chrs"
    params:
        read=1000,
        chrpattern="\"^chr[0-9XY]*$\"",
    conda:
        "../envs/samtools.yaml"
    shell:
        "samtools idxstats {input} | "  # returns list of chr and reads
        "awk '$3 > {params.read}' - | " # remove chromosomes with <1000 reads
        "cut -f1 | "                    # selects chr column
        "grep {params.chrpattern} | "   # select chrs 1-22,X,Y
        "paste -s -d, - > {output}"     # comma-separated concatenation of chrs

rule readcounter:
    input:
        bam="results/alignment/recal/{sample}.bqsr.bam",
        bai="results/alignment/recal/{sample}.bqsr.bam.bai",
        chrs="results/cnv/ichorcna/{sample}.chrs",
    output:
        "results/cnv/ichorcna/{sample}.wig",
    params:
        window=config['params']['readcounter']['window'],
        quality=config['params']['readcounter']['quality'],
    log:
        "logs/cnv/ichorcna/{sample}_readcounter.log"
    conda:
        "../envs/r.yaml"
    shell:
        "readCounter "
        "--window {params.window} "
        "--quality {params.quality} "
        "--chromosome $(cat {input.chrs}) "
        "{input.bam} | "
        "sed \"s/chrom=chr/chrom=/\" "
        "> {output} 2> {log}"

'''
rule ichorcna:
    input:
        bam="results/alignment/recal/{sample}.bqsr.bam",
        ref=config['common']['genome'],
    output:
        "results/cnv/{sample}.ichor"
    log:
        "logs/picard/wgs_metrics/{sample}.log"
    conda:
        "../envs/r.yaml"
    shell:
        "picard CollectWgsMetrics "
        "I={input.bam} "
        "O={output} "
        "R={input.ref} 2> {log}"
'''

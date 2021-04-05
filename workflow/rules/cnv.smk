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
        "runIchorCNA.R "
        "--WIG {input.wig} "
        "--NORMWIG {params.normalWig} "
        "--gcWig {params.gc_wig} "
        "--mapWig {params.mapWig} "
        "--normalPanel {params.normal_panel} "
        "--exons.bed {params.exons_bed} "
        "--id {params.outputFileNamePrefix} "
        "--centromere {params.centromere} "
        "--minMapScore {params.min_map_score} "
        "--rmCentromereFlankLength {params.rmCentromereFlankLength} "
        "--normal {params.normal} "
        "--scStates {params.scStates} "
        "--coverage {params.coverage} "
        "--lambda {params.lambda} "
        "--lambdaScaleHyperParam {params.lambdaScaleHyperParam} "
        "--ploidy {params.ploidy} "
        "--maxCN {params.maxCN} "
        "true="--estimateNormal True" false="--estimateNormal False" estimateNormal} \
        "true="--estimateScPrevalence True" false="--estimateScPrevalence  False" estimateScPrevalence} \
        "true="--estimatePloidy True" false="--estimatePloidy False" estimatePloidy} \
        "--maxFracCNASubclone {params.maxFracCNASubclone} "
        "--maxFracGenomeSubclone {params.maxFracGenomeSubclone} "
        "--minSegmentBins {params.minSegmentBins} "
        "--altFracThreshold {params.altFracThreshold} "
        "--chrNormalize {params.chrNormalize} "
        "--chrTrain {params.chrTrain} "
        "--chrs "c(~{sep="," chrs})" "
        "--genomeBuild {params.genomeBuild} "
        "--genomeStyle {params.genomeStyle} "
        "true="--normalizeMaleX True" false="--normalizeMaleX False" normalizeMaleX} \
        "--fracReadsInChrYForMale " + fracReadsInChrYForMale} \
        "true="--includeHOMD True" false="--includeHOMD False" includeHOMD} \
        "--txnE {params.txnE} "
        "--txnStrength {params.txnStrength} "
        "--plotFileType {params.plotFileType} "
        "--plotYLim {params.plotYLim} "
        "--libdir {params.libdir} "
        "--outDir {output} "

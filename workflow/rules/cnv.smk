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

rule get_RlibPath:
    output:
        "results/cnv/ichorcna/libpath"
    conda:
        "../envs/r.yaml"
    shell:
        "Rscript -e \"cat(.libPaths(), '\n')\" > {output}"

rule ichorcna:
    input:
        bam="results/alignment/recal/{sample}.bqsr.bam",
        rlib="results/cnv/ichorcna/libpath",
        chrs="results/cnv/ichorcna/{sample}.chrs",
        ref=config['common']['genome'],
    output:
        dir="results/cnv/ichorcna/{sample}",
        file="results/cnv/ichorcna/{sample}/{sample}_tpdf.pdf",
    log:
        "logs/cnv/ichorcna/{sample}.log"
    conda:
        "../envs/r.yaml"
    params:
        wig: "None",
        ploidy: config['params']['ichorcna']['ploidy'],
        normal: config['params']['ichorcna']['normal'],
        maxCN: config['params']['ichorcna']['maxCN'],
        gc_wig: get_ichorPath({input.rlib})['gc'],
        map_wig: get_ichorPath({input.rlib})['map'],
        centromere: get_ichorPath({input.rlib})['cen'],
        normal_panel: get_ichorPath({input.rlib})['norm'],
        include_HOMD: config['params']['ichorcna']['include_HOMD'],
        chrs: get_ichorChrs({input.chrs})['all'],
        chr_train: get_ichorChrs({input.chrs})['train'],
        estimateNormal: config['params']['ichorcna']['estimateNormal'],
        estimatePloidy: config['params']['ichorcna']['estimatePloidy'],
        estimateScPrevalence: config['params']['ichorcna']['estimateScPrevalence'],
        sc_states: config['params']['ichorcna']['sc_states'],
        txnE: config['params']['ichorcna']['txnE'],
        txn_strength: config['params']['ichorcna']['txn_strength'],
        sc_states: config['params']['ichorcna']['sc_states'],
    shell:
        "runIchorCNA.R "
        "--id {wildcards.sample} "
        "--WIG {input.wig} "
        "--ploidy '{params.ploidy}' "
        "--normal '{params.normal}' "
        "--maxCN {params.maxCN} "
        "--gcWig '{params.gc_wig}' "
        "--mapWig '{params.map_wig}' "
        "--centromere '{params.centromere}' "
        "--normalPanel '{params.normal_panel}' "
        "--includeHOMD {params.HOMD} "
        "--chrs '{params.chrs}' "
        "--chrTrain '{params.chrTrain}' "
        "--estimateNormal {params.estimateNormal} "
        "--estimatePloidy {params.estimatePloidy} "
        "--estimateScPrevalence {params.estimateScPrevalence} "
        "--scStates '{params.sc_states}' "
        "--txnE {params.txnE} "
        "--txnStrength {params.txn_strength} "
        "--outDir '{output}'"

'''
        "--NORMWIG {params.normalWig} "
        "--exons.bed {params.exons_bed} "
        "--minMapScore {params.min_map_score} "
        "--rmCentromereFlankLength {params.rmCentromereFlankLength} "
        "--coverage {params.coverage} "
        "--lambda {params.lambda} "
        "--lambdaScaleHyperParam {params.lambdaScaleHyperParam} "
        "--maxFracCNASubclone {params.maxFracCNASubclone} "
        "--maxFracGenomeSubclone {params.maxFracGenomeSubclone} "
        "--minSegmentBins {params.minSegmentBins} "
        "--altFracThreshold {params.altFracThreshold} "
        "--chrNormalize {params.chrNormalize} "
        "--genomeBuild {params.genomeBuild} "
        "--genomeStyle {params.genomeStyle} "
        "--normalizeMaleX {params.normalizeMalX} "
        "--fracReadsInChrYForMale {params.fracReadsInChrYForMale} "
        "--plotFileType {params.plotFileType} "
        "--plotYLim {params.plotYLim} "
        "--libdir {params.libdir} "
'''

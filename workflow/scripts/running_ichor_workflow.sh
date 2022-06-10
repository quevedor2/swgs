cd /cluster/projects/pughlab/projects/KOMBAT/test/data
outdir='/cluster/projects/pughlab/projects/KOMBAT/test/swgs'
sample='1428_RES_S11'

#bwa_mem
index='/cluster/projects/mcgahalab/ref/genomes/human/GRCh38/BWAIndex/genome.fasta'
extra="-R '@RG\\tID:${sample}\\tSM:${sample}\\tPL:ILLUMINA\\tPU:L001\\tLB:ATCACG' "
output="${outdir}/${sample}.sorted.bam"
threads=8

bwa mem \
-t ${threads} \
-R '@RG\tID:1428_RES_S11\tSM:1428_RES_S11\tPL:illumina\tPU:L001\tLB:lib' \
${index} \
${sample}_R1_001.fastq.gz ${sample}_R2_001.fastq.gz | \
samtools sort -o ${output} -

cd ${outdir}
output="${sample}.sorted.bam"
samtools index \
${output} \
${output}.bai


# gatk
outdir='/cluster/projects/pughlab/projects/KOMBAT/test/swgs'
cd ${outdir}

sample='1428_RES_S11'
mem=10
rmduplicates='true'
sort='coordinate'
input="${sample}.sorted.bam"
output="${sample}.markdup.bam"
metrics='${sample}.metric.txt'

java  \
-Xmx8G  \
-jar /cluster/home/quever/miniconda3/envs/preprocess/share/picard-2.27.2-0/picard.jar \
MarkDuplicates \
--REMOVE_DUPLICATES true \
--ASSUME_SORT_ORDER coordinate \
-I 1428_RES_S11.sorted.bam \
-O 1428_RES_S11.markdup.bam \
--METRICS_FILE ${sample}.metric.txt

samtools index \
${output} \
${output}.bai


mem=10
threads=8
input="${sample}.markdup.bam"
ref="/cluster/projects/mcgahalab/ref/genomes/human/GRCh38/genome.fasta"
recal="${sample}.recal_data_table"
output="${sample}.bqsr.bam"
knownsite='/cluster/projects/mcgahalab/ref/gatk/GRCh38/resources_broad_hg38_v0_Mills_and_1000G_gold_standard.indels.GRCh38.vcf.gz'

gatk \
BaseRecalibrator \
-I ${input} \
-R ${ref} \
--known-sites ${knownsite} \
-O ${recal}

gatk \
ApplyBQSR \
-I ${input} \
-R ${ref} \
--bqsr-recal-file ${recal} \
-O ${output}


# cnv
cd $outdir
input="${sample}.bqsr.bam"
output="${sample}.chrs"
read=1000
chrpattern="^[0-9XY]*$"
    
samtools idxstats ${input} | \
awk '$3>1000' - | \
cut -f1 | \
grep ${chrpattern} | \
paste -s -d, - > ${output}



cd /cluster/projects/pughlab/projects/KOMBAT/test/data
outdir='/cluster/projects/pughlab/projects/KOMBAT/test/swgs'
cd $outdir
conda activate r-ichor
sample='1428_RES_S11'

bam="${sample}.bqsr.bam"
chrs="${sample}.chrs"
output="${sample}.wig"
window=1000000
quality=20

readCounter \
--window ${window} \
--quality ${quality} \
--chromosome $(cat ${chrs}) \
${bam} | \
sed "s/chrom=chr/chrom=/" \
> ${output}


sample='1428_RES_S11'
wig="${sample}.wig"
chrs="${sample}.chrs"
echo "/cluster/home/quever/miniconda3/envs/r-ichor/lib/R/library/" > rlib_path
rlib='rlib_path'
ploidy="c(2,3)"
normal="c(0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9)"
genome_style='NCBI'
maxCN=5
include_HOMD='TRUE'
estimateNormal='TRUE'
estimatePloidy='TRUE'
estimateScPrevalence='TRUE'
sc_states="c(1,3)"
txnE=0.9999
txn_strength=10000

Rscript runIchorCNA.R \
--id ${sample} \
--WIG ${wig} \
--ploidy ${ploidy} \
--normal "${normal}" \
--maxCN ${maxCN} \
--gcWig "$(python parse_paths.py -i ${rlib} -r 'gc')"  \
--mapWig "$(python parse_paths.py -i ${rlib} -r 'map')" \
--centromere "$(python parse_paths.py -i ${rlib} -r 'cen')" \
--normalPanel "$(python parse_paths.py -i ${rlib} -r 'norm')" \
--includeHOMD "${HOMD}" \
--chrs "$(python parse_paths.py -i ${chrs} -r 'all')" \
--chrTrain "$(python parse_paths.py -i ${chrs} -r 'train')" \
--estimateNormal "${estimateNormal}" \
--estimatePloidy "${estimatePloidy}" \
--estimateScPrevalence "${estimateScPrevalence}" \
--scStates "${sc_states}" \
--txnE ${txnE} \
--txnStrength ${txn_strength} \
--genomeStyle ${genome_style} \
--outDir ${sample}

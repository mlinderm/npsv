#!/usr/bin/env bash

# Run NPSV "end-to-end" on GIAB data (NGS data and SVs). Note that aspects of this script
# are specific to the local computing infrastructure.

#SBATCH --job-name=example
#SBATCH --output=example-%j.out
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=36
#SBATCH --mem=88G                              
#SBATCH --partition=long
#SBATCH --time=7-00:00:00

set -euo pipefail

REFERENCE=/storage/mlinderman/ngs/resources/gatk/b37/human_g1k_v37.fasta
NPSV_ROOT=/home/mlinderman/research/npsv

## Setup environment
export TMPDIR=$(mktemp -d --tmpdir="$SCRATCH") || exit 1
trap "rm -rf ${TMPDIR};" 0

## Download NGS data for HG002 (corresponding to approximately 25x coverage)
# Adapted from: https://github.com/brentp/duphold/blob/master/paper/evaluation.sh
ALIGN_FILE=hg002.bam
wget -N -r ftp://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/data/AshkenazimTrio/HG002_NA24385_son/NIST_HiSeq_HG002_Homogeneity-10953946/HG002_HiSeq300x_fastq/140528_D00360_0018_AH8VC6ADXX/

for f in $(find . -type f -name "*_R1_*.fastq.gz" | sort); do zcat $f; done | bgzip -@ 8 -c > hg002_R1.fastq.gz &
for f in $(find . -type f -name "*_R2_*.fastq.gz" | sort); do zcat $f; done | bgzip -@ 8 -c > hg002_R2.fastq.gz &
wait

## Align NGS data using simplified pipeline
bwa mem -c 250 -M -v 1 -t 31 \
    -R '@RG\tID:HG002\tSM:HG002\tPL:illumina\tPU:HG002\tLB:HG002' \
    $REFERENCE \
    hg002_R1.fastq.gz hg002_R2.fastq.gz | \
    samblaster -q -M --addMateTags | \
    samtools sort -T "${TMPDIR}/samtools" -@ 4 -m 16G --output-fmt BAM --reference $REFERENCE -o $ALIGN_FILE
samtools index $ALIGN_FILE

## Download GIAB variants, BED files and filter callset
wget -N ftp://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/release//AshkenazimTrio/HG002_NA24385_son/NIST_SV_v0.6/HG002_SVs_Tier1_v0.6.vcf.gz
wget -N ftp://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/release//AshkenazimTrio/HG002_NA24385_son/NIST_SV_v0.6/HG002_SVs_Tier1_v0.6.vcf.gz.tbi
wget -N ftp://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/release//AshkenazimTrio/HG002_NA24385_son/NIST_SV_v0.6/HG002_SVs_Tier1_v0.6.bed
wget -N ftp://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/release//AshkenazimTrio/HG002_NA24385_son/NIST_SV_v0.6/HG002_SVs_Tier1plusTier2_v0.6.1.bed

GIAB_VCF=HG002_SVs_Tier1_v0.6.vcf.gz
GIAB_BED=HG002_SVs_Tier1_v0.6.bed
GIAB_ALL_TIERS_BED=HG002_SVs_Tier1and2_v0.6.bed
GIAB_FILTERED_VCF="${GIAB_VCF%.vcf.gz}.genotyped.passing.tier1and2.vcf.gz"

cat HG002_SVs_Tier1_v0.6.bed HG002_SVs_Tier1plusTier2_v0.6.1.bed | sort -k1,1 -k2,2n > "${TMPDIR}/tmp.bed"
bedtools merge -i "${TMPDIR}/tmp.bed" > $GIAB_ALL_TIERS_BED

# Alternately filter with $GIAB_ALL_TIERS_BED
bcftools view -i '(INFO/sizecat != "20to49")' -g ^miss -f 'PASS,LongReadHomRef' -R $GIAB_BED $GIAB_VCF | \
    bcftools sort -Oz -o $GIAB_FILTERED_VCF
bcftools index -t $GIAB_FILTERED_VCF

## Genotype SVs

# Preprocess alignment file prior to genotyping to obtain coverage, insert size distribution, etc.
# As a "shortcut" this step can use pre-computed Picard metrics, or be bypassed directly and the relevant
# statistics provided as arguments to the NPSV genotyper.
npsvg preprocess \
    -r $REFERENCE \
    -b $ALIGN_FILE \
    --genome $NPSV_ROOT/etc/human_g1k_v37.genome \
    -o stats.json

# Load the BWA index into shared memory so it can be reused during simulation and across threads
bwa shm /storage/mlinderman/ngs/resources/gatk/b37/human_g1k_v37.fasta
trap "rm -rf ${TMPDIR}; bwa shm -d;" 0

NPSV_PREFIX="$(basename -s ".vcf.gz" $GIAB_FILTERED_VCF)"
npsv \
    --tempdir $TMPDIR \
    --threads 36 \
    -r $REFERENCE \
    --genome $NPSV_ROOT/etc/human_g1k_v37.genome \
    --gaps $NPSV_ROOT/etc/human_g1k_v37.gaps.bed.gz \
    --stats-path stats.json \
    --profile HS25 \
    --filter-bed $GIAB_BED \
    -i $GIAB_FILTERED_VCF \
    -b $ALIGN_FILE \
    -o $PWD \
    --prefix $NPSV_PREFIX


## Concordance analysis with Truvari (incorporating LongReadHomRef calls)

GIAB_FILTERED_DEL_LONGREADHOMREF_VCF="${GIAB_FILTERED_VCF%.vcf.gz}.DEL.longreadhomref.vcf.gz"
GIAB_FILTERED_INS_LONGREADHOMREF_VCF="${GIAB_FILTERED_VCF%.vcf.gz}.INS.longreadhomref.vcf.gz"

GENOTYPED_DEL_LONGREADHOMREF_VCF="${NPSV_PREFIX}.npsv.DEL.longreadhomref.vcf.gz"
GENOTYPED_INS_LONGREADHOMREF_VCF="${NPSV_PREFIX}.npsv.INS.longreadhomref.vcf.gz"

bcftools annotate -i '(SVTYPE ~ "^DEL")' -x "FILTER/LongReadHomRef" -Oz -o $GIAB_FILTERED_DEL_LONGREADHOMREF_VCF $GIAB_FILTERED_VCF
bcftools index -t $GIAB_FILTERED_DEL_LONGREADHOMREF_VCF

bcftools annotate -i '(SVTYPE ~ "^INS")' -x "FILTER/LongReadHomRef" -Oz -o $GIAB_FILTERED_INS_LONGREADHOMREF_VCF $GIAB_FILTERED_VCF
bcftools index -t $GIAB_FILTERED_INS_LONGREADHOMREF_VCF

bcftools annotate -i '(SVTYPE ~ "^DEL")' -x "FILTER/LongReadHomRef" -Oz -o $GENOTYPED_DEL_LONGREADHOMREF_VCF "${NPSV_PREFIX}.npsv.vcf" 
bcftools index -t $GENOTYPED_DEL_LONGREADHOMREF_VCF

bcftools annotate -i '(SVTYPE ~ "^INS")' -x "FILTER/LongReadHomRef" -Oz -o $GENOTYPED_INS_LONGREADHOMREF_VCF "${NPSV_PREFIX}.npsv.vcf" 
bcftools index -t $GENOTYPED_INS_LONGREADHOMREF_VCF

module load truvari/genotype

run_truvari() {
    rm -rf "truvari_${1}"
    truvari bench \
        -f $REFERENCE \
        -b $3 \
        -c $2 \
        -o "truvari_${1}" \
        --includebed $GIAB_BED \
        --sizemax 15000000 -s 50 -S 30 --pctsim=0 -r 20 -O 0.6 \
        --passonly --bSample HG002 --cSample HG002
}

run_truvari DEL $GENOTYPED_DEL_LONGREADHOMREF_VCF $GIAB_FILTERED_DEL_LONGREADHOMREF_VCF
run_truvari INS $GENOTYPED_INS_LONGREADHOMREF_VCF $GIAB_FILTERED_INS_LONGREADHOMREF_VCF

#!/usr/bin/env bash

# Benchmark SV genotypers. Note that aspects of this script are specific to the local computing infrastructure.

#SBATCH --job-name=benchmarking
#SBATCH --output=benchmarking-%j.out
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=36
#SBATCH --mem=88G                              
#SBATCH --partition=long
#SBATCH --time=7-00:00:00

set -euo pipefail

# Metadata
THREADS=$SLURM_CPUS_PER_TASK
TIME_FORMAT="%E,%M,%U,%S,${THREADS},${SLURM_JOB_ID}"
OUTPUT_DIR=

# Data paths
REFERENCE=/storage/mlinderman/ngs/resources/gatk/b37/human_g1k_v37.fasta
SNV_FILE=/storage/mlinderman/ngs/resources/ashkenazi-trio/final/2019-06-13_project/ashkenazi-gatk-haplotype-annotated.vcf.gz

# NPSV paths
NPSV_ROOT=/home/mlinderman/research/npsv
NPSV_STATS_FILE=HG002-ready.stats.json

# GIAB inputs
GIAB_BED=HG002_SVs_Tier1_v0.6.bed
GIAB_VCF=HG002_SVs_Tier1_v0.6.vcf.gz

usage()
{
    cat << EOF
usage: $(basename "$0") [options] GENOTYPER BAM 

Benchmark SV gneotypers

Options:
  -h            Print this message
EOF
}

while getopts "hR:S:c:m:s:t:l:p:z:f:i:p:n:g:j:" Option
do
    case $Option in
        h)
            usage
            exit 0
            ;;
        ?)
            usage
            exit 85
            ;;
    esac
done

shift $((OPTIND-1))
if [[ $# -ne 2 ]]; then
    >&2 echo "Error: Missing positional arguments"
    >&2 usage
    exit 1
fi

GENOTYPER=$1
BAM_FILE=$2

# Setup environment
export TMPDIR=$(mktemp -d --tmpdir="$SCRATCH") || exit 1
trap "rm -rf ${TMPDIR};" 0

if [[ -z $OUTPUT_DIR ]]; then
    OUTPUT_DIR="${TMPDIR}/output"
fi
mkdir -p $OUTPUT_DIR

module load bcftools
module load htslib

# Prepare inputs from GIAB-provided files
GIAB_ALL_TIERS_BED="${TMPDIR}/HG002_SVs_Tier1and2_v0.6.bed"
GIAB_FILTERED_VCF="${TMPDIR}/${GIAB_VCF%.vcf.gz}.genotyped.passing.tier1and2.vcf.gz"
GIAB_FILTERED_DEL_VCF="${GIAB_FILTERED_VCF%.vcf.gz}.DEL.vcf.gz"

# Combine GIAB BED files
cat $GIAB_BED HG002_SVs_Tier1plusTier2_v0.6.1.bed | sort -k1,1 -k2,2n > "${TMPDIR}/tmp.bed"
bedtools merge -i "${TMPDIR}/tmp.bed" > $GIAB_ALL_TIERS_BED

# Filter GIAB VCF
bcftools view -i '(INFO/sizecat != "20to49")' -g ^miss -f 'PASS,LongReadHomRef' -R $GIAB_ALL_TIERS_BED $GIAB_VCF | \
    bcftools sort -Oz -o $GIAB_FILTERED_VCF
bcftools index -t $GIAB_FILTERED_VCF

# Create deletion-only input for GenomeSTRiP
bcftools view -i '(SVTYPE ~ "^DEL")' -Oz -o $GIAB_FILTERED_DEL_VCF $GIAB_FILTERED_VCF
bcftools index -t $GIAB_FILTERED_DEL_VCF

# Which filters to reset with BCFTools prior to Truvari analysis
BCFTOOLS_FILTER="FILTER/LongReadHomRef"
MODE="NA"

# Call set sizes
echo "Variants,ALL,$(zgrep -c -v "^#" $GIAB_FILTERED_VCF)"
echo "Variants,DEL,$(zgrep -c -v "^#" $GIAB_FILTERED_DEL_VCF)"

if [[ $GENOTYPER =~ ^npsv ]]; then
    # Load reference genome into shared memory
    bwa shm $REFERENCE
    trap "rm -rf ${TMPDIR}; bwa shm -d;" 0
    
    # Preprocessing step
    /bin/time -f "Timing,npsv,${MODE},preprocessing,${TIME_FORMAT}" \
    npsvg preprocess \
        -r $REFERENCE \
        -b $BAM_FILE \
        --genome $NPSV_ROOT/etc/human_g1k_v37.genome \
        -o "${OUTPUT_DIR}/stats.json"

fi

if [[ $GENOTYPER == "npsv" ]]; then
    MODE="default"
    NPSV_PREFIX="$(basename -s ".vcf.gz" $GIAB_FILTERED_VCF)"

    /bin/time -f "Timing,npsv,${MODE},genotyping,${TIME_FORMAT}" \
    npsv \
        --tempdir $TMPDIR \
        --threads $SLURM_CPUS_PER_TASK \
        -r $REFERENCE \
        --genome $NPSV_ROOT/etc/human_g1k_v37.genome \
        --gaps $NPSV_ROOT/etc/human_g1k_v37.gaps.bed.gz \
        --stats-path $NPSV_STATS_FILE \
        --filter-bed $GIAB_BED \
        -i $GIAB_FILTERED_VCF \
        -b $BAM_FILE \
        -o $OUTPUT_DIR \
        --prefix $NPSV_PREFIX
    
    GENOTYPED_VCF="${OUTPUT_DIR}/${NPSV_PREFIX}.npsv.vcf"
    bgzip $GENOTYPED_VCF && tabix "${GENOTYPED_VCF}.gz"
    GENOTYPED_VCF="${GENOTYPED_VCF}.gz"    

elif [[ $GENOTYPER =~ npsv_(single|variant|hybrid) ]]; then
    if [[ $GENOTYPER == "npsv_single" ]]; then
        N=1
        KLASS="svm"
        MODE="single"
    elif [[ $GENOTYPER == "npsv_variant" ]]; then
        N=100
        KLASS="rf"
        MODE="variant"
    elif [[ $GENOTYPER == "npsv_hybrid" ]]; then
        N=100
        KLASS="rf"
        MODE="hybrid"
    fi
    
    NPSV_PREFIX="$(basename -s ".vcf.gz" $GIAB_FILTERED_VCF)"

    # Simulation and genotyping
    /bin/time -f "Timing,npsv,${MODE},genotyping,${TIME_FORMAT}" \
    npsv \
        --tempdir $TMPDIR \
        --threads $THREADS \
        -r $REFERENCE \
        --genome $NPSV_ROOT/etc/human_g1k_v37.genome \
        --gaps $NPSV_ROOT/etc/human_g1k_v37.gaps.bed.gz \
        --stats-path $NPSV_STATS_FILE \
        --filter-bed $GIAB_BED \
        -i $GIAB_FILTERED_VCF \
        -b $BAM_FILE \
        -o $OUTPUT_DIR \
        --prefix $NPSV_PREFIX \
        --DEL-gt-mode $MODE --DEL-n $N --DEL-classifier $KLASS --DEL-hybrid-threshold 1000 \
        --INS-gt-mode $MODE --INS-n $N --DEL-classifier $KLASS --INS-hybrid-threshold 1000

    GENOTYPED_VCF="${OUTPUT_DIR}/${NPSV_PREFIX}.npsv.vcf"
    bgzip $GENOTYPED_VCF && tabix "${GENOTYPED_VCF}.gz"
    GENOTYPED_VCF="${GENOTYPED_VCF}.gz" 

elif [[ $GENOTYPER == "svviz2_mapq" ]]; then
    module load svviz2
    
    MODE="mapq"
    REPORT_DIR="${OUTPUT_DIR}/svviz2_reports"
    mkdir -p $REPORT_DIR

    /bin/time -f "Timing,svviz2,${MODE},genotyping,${TIME_FORMAT}" \
    svviz2 \
        --report-only \
        --ref $REFERENCE \
        -o $REPORT_DIR \
		--variants <(bcftools annotate --set-id +'%CHROM\_%POS\_%INFO/END' $GIAB_FILTERED_VCF) \
        $BAM_FILE

    GENOTYPED_VCF="${OUTPUT_DIR}/HG002_SVs_Tier1_v0.6.genotyped.passing.tier1and2.svviz2_mapq.vcf.gz"
    
    /bin/time -f "Timing,svviz2,${MODE},postprocessing,${TIME_FORMAT}" \
    svviz22vcf \
        --model mapq \
        -i $GIAB_FILTERED_VCF \
        -r $REPORT_DIR \
		-s $(basename -s .bam $BAM_FILE) \
		-o /dev/stdout | \
		bcftools reheader -s <(echo -e "HG002") - | \
		sed 's/\bnan\b/./g' | \
		bgzip > $GENOTYPED_VCF
	bcftools index -t $GENOTYPED_VCF

elif [[ $GENOTYPER == "svtyper" ]]; then
    module load svtyper
    
    GENOTYPED_VCF="${OUTPUT_DIR}/HG002_SVs_Tier1_v0.6.genotyped.passing.tier1and2.svtyper.vcf.gz"
    /bin/time -f "Timing,svtyper,${MODE},genotyping,${TIME_FORMAT}" \
    svtyper \
        -i <(bcftools view -G $GIAB_FILTERED_VCF) \
        -B $BAM_FILE | \
		sed 's/=""/="/' | \
		bgzip -c > $GENOTYPED_VCF
	bcftools index -t $GENOTYPED_VCF

elif [[ $GENOTYPER == "paragraph" ]]; then   
    PARAGRAPH_INPUT_VCF="${TMPDIR}/paragraph.vcf.gz"
    padAlleles -r $REFERENCE -i $GIAB_FILTERED_VCF | bcftools sort -Oz -o $PARAGRAPH_INPUT_VCF
    bcftools index -t $PARAGRAPH_INPUT_VCF

    module load paragraph/2.4a
    
    PARAGRAPH_SAMPLE_FILE="${TMPDIR}/paragraph.sample.txt"
    cat <<EOF > $PARAGRAPH_SAMPLE_FILE
id	path	depth	read length	sex
HG002	$BAM_FILE	25	148	male
EOF

	/bin/time -f "Timing,paragraph,${MODE},genotyping,${TIME_FORMAT}" \
    multigrmpy.py \
		--scratch-dir $TMPDIR -t $SLURM_CPUS_PER_TASK \
		-i $PARAGRAPH_INPUT_VCF \
		-m $PARAGRAPH_SAMPLE_FILE \
		-r $REFERENCE \
		-o $OUTPUT_DIR
	bcftools index -t "${OUTPUT_DIR}/genotypes.vcf.gz"   
    
    GENOTYPED_VCF="${OUTPUT_DIR}/HG002_SVs_Tier1_v0.6.genotyped.passing.tier1and2.paragraph.vcf.gz"
	
    # Need to set ploidy to 2 for truvari comparison with GIAB and we remove the filters
    bcftools +fixploidy -Oz -o $GENOTYPED_VCF "${OUTPUT_DIR}/genotypes.vcf.gz" -- -f 2    
	bcftools index -t $GENOTYPED_VCF

    BCFTOOLS_FILTER="FILTER"

elif [[ $GENOTYPER == "sv2" ]]; then
    module load sv2/1.5

    SV2_PED_FILE="${TMPDIR}/sv2.ped"
    cat <<EOF > $SV2_PED_FILE
8392	HG002	HG003	HG004	1	0       
EOF

    GENOTYPED_VCF="${OUTPUT_DIR}/sv2_genotypes/HG002_SVs_Tier1_v0.6.genotyped.passing.tier1and2.sv2.vcf"
	/bin/time -f "Timing,sv2,${MODE},genotyping,${TIME_FORMAT}" sv2 \
        -M \
        -i $BAM_FILE \
        -v $GIAB_FILTERED_VCF \
        -snv $SNV_FILE \
        -p $SV2_PED_FILE \
        -O $OUTPUT_DIR \
        -o HG002_SVs_Tier1_v0.6.genotyped.passing.tier1and2.sv2    

    sed '/^#/b; s/NA/./g' $GENOTYPED_VCF | \
        sed 's/\(^##FILTER=<ID=GAP\)>/\1,Description="Missing description">/' | \
        bcftools +fixploidy -- -f 2 - | \
        bcftools norm -f $REFERENCE -c s - | \
        bcftools annotate -x INFO/GENES -Oz -o "${GENOTYPED_VCF}.gz" -
	
    GENOTYPED_VCF="${GENOTYPED_VCF}.gz"
    bcftools index -t $GENOTYPED_VCF

    BCFTOOLS_FILTER="FILTER"

elif [[ $GENOTYPER == "delly" ]]; then
    module load delly/0.8.3

    DELLY_INPUT_BCF="${TMPDIR}/delly.vcf"
    bcftools view -Ob -o $DELLY_INPUT_BCF $GIAB_FILTERED_VCF
    bcftools index -c $DELLY_INPUT_BCF

    GENOTYPED_BCF="${OUTPUT_DIR}/HG002_SVs_Tier1_v0.6.genotyped.passing.tier1and2.delly.bcf"
    /bin/time -f "Timing,delly,${MODE},genotyping,${TIME_FORMAT}" \
    delly \
        call \
        -g $REFERENCE -x "${DELLY_HOME}/excludeTemplates/human.hg19.excl.tsv" \
        -v $DELLY_INPUT_BCF -o $GENOTYPED_BCF $BAM_FILE

    GENOTYPED_VCF="${OUTPUT_DIR}/HG002_SVs_Tier1_v0.6.genotyped.passing.tier1and2.delly.vcf.gz"
    bcftools view -Oz -o $GENOTYPED_VCF $GENOTYPED_BCF
    bcftools index -t $GENOTYPED_VCF

    BCFTOOLS_FILTER="FILTER"

elif [[ $GENOTYPER == "genomestrip" ]]; then
    module load GenomeSTRiP/2.00.1958 

    SV_METADATA_DIR=/storage/mlinderman/ngs/resources/genomestrip/1000G_phase1
    SAMPLE_METADATA_DIR="${OUTPUT_DIR}/genomestrip_preprocessing"
    mkdir -p $SAMPLE_METADATA_DIR

    /bin/time -f "Timing,genomestrip,${MODE},preprocessing,${TIME_FORMAT}" \
    java -Xmx4g -cp $CLASSPATH \
		-Djava.io.tmpdir=$TMPDIR \
		org.broadinstitute.gatk.queue.QCommandLine \
		-S $SV_DIR/qscript/SVPreprocess.q \
		-S $SV_DIR/qscript/SVQScript.q \
		-cp $CLASSPATH \
		-gatk $SV_DIR/lib/gatk/GenomeAnalysisTK.jar \
		-configFile $SV_DIR/conf/genstrip_parameters.txt \
		-ploidyMapFile $SV_METADATA_DIR/human_g1k_v37.ploidymap.txt \
		-R $REFERENCE \
		-I $BAM_FILE \
		-md $SAMPLE_METADATA_DIR \
		-bamFilesAreDisjoint true \
		-jobLogDir "${OUTPUT_DIR}/genomestrip_log" \
		-jobRunner ParallelShell -maxConcurrentRun $SLURM_CPUS_PER_TASK \
		-run

    GENOMESTRIP_GENDER_MAP="${TMPDIR}/genomestrip.gender.txt"
    echo -e "HG002\tMale" > $GENOMESTRIP_GENDER_MAP

    GENOTYPED_VCF="${OUTPUT_DIR}/$(basename -s ".vcf.gz" $GIAB_FILTERED_DEL_VCF).genomestrip.vcf.gz"

    /bin/time -f "Timing,genomestrip,${MODE},genotyping,${TIME_FORMAT}" \
    java -Xmx4g -cp $CLASSPATH \
		-Djava.io.tmpdir=$TMPDIR \
    	org.broadinstitute.gatk.queue.QCommandLine \
    	-S $SV_DIR/qscript/SVGenotyper.q \
		-S $SV_DIR/qscript/SVQScript.q \
		-cp $CLASSPATH \
		-gatk $SV_DIR/lib/gatk/GenomeAnalysisTK.jar \
		-configFile $SV_DIR/conf/genstrip_parameters.txt \
		-rmd $SV_METADATA_DIR \
		-R $REFERENCE \
		-genomeMaskFile $SV_METADATA_DIR/human_g1k_v37.svmask.fasta \
		-genderMapFile $GENOMESTRIP_GENDER_MAP \
		-ploidyMapFile $SV_METADATA_DIR/human_g1k_v37.ploidymap.txt \
		-md $SAMPLE_METADATA_DIR \
		-runDirectory $OUTPUT_DIR \
		-vcf $GIAB_FILTERED_DEL_VCF \
		-I $BAM_FILE \
		-O $GENOTYPED_VCF \
		-bamFilesAreDisjoint true \
		-jobLogDir "${OUTPUT_DIR}/genomestrip_log" \
		-jobRunner ParallelShell -maxConcurrentRun $SLURM_CPUS_PER_TASK \
		-run

    # Remove FORMAT field that seems to be causing bcftools parsing errors
    gzip -dc $GENOTYPED_VCF | \
        sed 's/:NA:/:.:/g' | \
        bcftools annotate -Oz -o "${GENOTYPED_VCF%.vcf.gz}.clean.vcf.gz" -x FORMAT/CNF -
    GENOTYPED_VCF="${GENOTYPED_VCF%.vcf.gz}.clean.vcf.gz"
    bcftools index -t $GENOTYPED_VCF

elif [[ $GENOTYPER == "graphtyper" ]]; then
    module load graphtyper/2.5.1

    GRAPHTYPER_BAMS="${OUTPUT_DIR}/graphtyper.bams.list"
    echo "$BAM_FILE" > $GRAPHTYPER_BAMS

    GRAPHTYPER_OUTPUT_DIR="${OUTPUT_DIR}/graphtyper"
    /bin/time -f "Timing,graphtyper,${MODE},genotyping,${TIME_FORMAT}" \
    graphtyper \
        genotype_sv $REFERENCE $GIAB_FILTERED_VCF \
		--output=$GRAPHTYPER_OUTPUT_DIR \
		--sams=$GRAPHTYPER_BAMS \
		--region_file=<(awk '{ print $1 ":" $2+1 "-" $3; }' $GIAB_ALL_TIERS_BED) \
		--force_no_copy_reference --threads=1
	
    GENOTYPED_VCF="${OUTPUT_DIR}/HG002_SVs_Tier1_v0.6.genotyped.passing.tier1and2.graphtyper.vcf.gz"
    bcftools concat --naive $GRAPHTYPER_OUTPUT_DIR/**/*.vcf.gz | \
        bcftools view -i '((SVTYPE ~ "^DEL" || SVTYPE ~ "^INS" || SVTYPE ~ "^DUP") && INFO/SVMODEL=="AGGREGATED")' | \
        sed '/^#/!s/DUP/INS/g' | \
        bcftools sort -Oz -o $GENOTYPED_VCF -
	bcftools index -t $GENOTYPED_VCF

    BCFTOOLS_FILTER="FILTER"
fi

# Deactivate NPSV virtual environment to switch to Truvari environemnt
module load truvari/genotype

GIAB_FILTERED_LONGREADHOMREF_VCF="${OUTPUT_DIR}/HG002_SVs_Tier1_v0.6.genotyped.passing.tier1and2.longreadhomref.vcf.gz"
GIAB_FILTERED_DEL_LONGREADHOMREF_VCF="${OUTPUT_DIR}/HG002_SVs_Tier1_v0.6.genotyped.passing.tier1and2.DEL.longreadhomref.vcf.gz"
GIAB_FILTERED_INS_LONGREADHOMREF_VCF="${OUTPUT_DIR}/HG002_SVs_Tier1_v0.6.genotyped.passing.tier1and2.INS.longreadhomref.vcf.gz"

GENOTYPED_LONGREADHOMREF_VCF="${GENOTYPED_VCF%.vcf.gz}.longreadhomref.vcf.gz"
GENOTYPED_DEL_LONGREADHOMREF_VCF="${GENOTYPED_VCF%.vcf.gz}.DEL.longreadhomref.vcf.gz"
GENOTYPED_INS_LONGREADHOMREF_VCF="${GENOTYPED_VCF%.vcf.gz}.INS.longreadhomref.vcf.gz"

# Prepare GIAB files
bcftools annotate -x "FILTER/LongReadHomRef" -Oz -o $GIAB_FILTERED_LONGREADHOMREF_VCF $GIAB_FILTERED_VCF
bcftools index -t $GIAB_FILTERED_LONGREADHOMREF_VCF

bcftools annotate -i '(SVTYPE ~ "^DEL")' -x "FILTER/LongReadHomRef" -Oz -o $GIAB_FILTERED_DEL_LONGREADHOMREF_VCF $GIAB_FILTERED_VCF
bcftools index -t $GIAB_FILTERED_DEL_LONGREADHOMREF_VCF

bcftools annotate -i '(SVTYPE ~ "^INS")' -x "FILTER/LongReadHomRef" -Oz -o $GIAB_FILTERED_INS_LONGREADHOMREF_VCF $GIAB_FILTERED_VCF
bcftools index -t $GIAB_FILTERED_INS_LONGREADHOMREF_VCF

# Prepare output from the tools
bcftools annotate -x $BCFTOOLS_FILTER -Oz -o $GENOTYPED_LONGREADHOMREF_VCF $GENOTYPED_VCF 
bcftools index -t $GENOTYPED_LONGREADHOMREF_VCF

bcftools annotate -i '(SVTYPE ~ "^DEL")' -x $BCFTOOLS_FILTER -Oz -o $GENOTYPED_DEL_LONGREADHOMREF_VCF $GENOTYPED_VCF 
bcftools index -t $GENOTYPED_DEL_LONGREADHOMREF_VCF

bcftools annotate -i '(SVTYPE ~ "^INS")' -x $BCFTOOLS_FILTER -Oz -o $GENOTYPED_INS_LONGREADHOMREF_VCF $GENOTYPED_VCF 
bcftools index -t $GENOTYPED_INS_LONGREADHOMREF_VCF

run_truvari() {  
    # sv2 converts complex alleles to symbolic so pctsize can prevent variant matching
    rm -rf "$TMPDIR/truvari"
    truvari bench \
        -f $REFERENCE \
        -b $3 \
        -c $2 \
        -o $TMPDIR/truvari \
        --includebed $GIAB_BED \
        --sizemax 15000000 -s 50 -S 30 --pctsim=0 -r 20 -O 0.6 \
        $([ "$GENOTYPER" == "sv2" ] && echo -n "--pctsize=0") \
        --passonly --bSample HG002 --cSample HG002

    echo -n "Truvari,${1},"
    python3 -c "import json; conc = json.load(open(f'${TMPDIR}/truvari/summary.txt')); print(conc['gt_concordance'], conc['gt_nonref_concordance'], sep=',');"
}

run_truvari "${GENOTYPER},${MODE},ALL" $GENOTYPED_LONGREADHOMREF_VCF $GIAB_FILTERED_LONGREADHOMREF_VCF
run_truvari "${GENOTYPER},${MODE},DEL" $GENOTYPED_DEL_LONGREADHOMREF_VCF $GIAB_FILTERED_DEL_LONGREADHOMREF_VCF
run_truvari "${GENOTYPER},${MODE},INS" $GENOTYPED_INS_LONGREADHOMREF_VCF $GIAB_FILTERED_INS_LONGREADHOMREF_VCF


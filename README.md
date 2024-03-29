# NPSV: Non-parametric Structural Variant Genotyper

NPSV is a Python-based tool for stand-alone genotyping of previously detected/reported deletion and insertion structural variants (SVs) in short-read whole genome sequencing (WGS) data. NPSV implements a machine learning-based approach for SV genotyping that employs NGS simulation to model the combined effects of the genomic region, sequencer and alignment pipeline. 

NPSV is described in the following publication:

Linderman MD, Paudyal C, Shakeel M, Kelley W, Bashir A, Gelb BD. [NPSV: A simulation-driven approach to genotyping structural variants in whole-genome sequencing data](https://doi.org/10.1093/gigascience/giab046). Gigascience. 2021;10(7).

## Installation

When cloning NPSV, make sure to recursively clone all of the submodules, i.e. `git clone --recursive git@github.com:mlinderm/npsv.git`.

NPSV requires Python 3.6+ and a suite of command-line genomics tools. For convenience, a Docker file is provided that installs all of the dependencies. To build that image:
```
docker build -t npsv .
```

### Manual installation

To manually install and run NPSV from the source, you will need the following dependencies:

* ART (NGS simulator)
* bwa
* bedtools
* bcftools
* goleft
* htslib (i.e., tabix and bgzip)
* samblaster
* sambamba
* samtools

along with standard command-line utilities, such as bc, gawk, etc., CMake and a C++14 compiler. After installing the dependencies listed above, install the Python dependencies, and then NPSV itself via:
```
python3 -m pip install -r requirements.txt
python3 setup.py install
```

## Running NPSV

NPSV requires basic information about the aligned reads (i.e. sequencer model, coverage, insert size distribution). These data can be provided as a command line parameters enabling you to immediately start genotyping. An optional preprocessing step (the typical workflow) will collect that data and more from the BAM file (to inform both simulation and feature extraction) into a stats file that can be used with the genotyper.

### Running the NPSV tools with Docker

Given the multi-step workflow, the typical approach when using the Docker image is to run NPSV from a shell. The following command will start a Bash session in the Docker container (replace `/path/to/reference/directory` with the path to directory containing the reference genome and associated BWA indices). NPSV is most efficient when the BWA indices are loaded into shared memory. To load BWA indices into shared memory you will need to configure the Docker container with at least 10G of memory and set the shared memory size to 6G or more.

```
docker run --entrypoint /bin/bash \
    --shm-size=6g \
    -v `pwd`:/opt/npsv \
    -v /path/to/reference/directory:/data \
    -w /opt/npsv \
    -it \
    npsv
```

## NPSV Genotyping

The NPSV package installs two key executables, `npsv`, the main entry point for the genotyper, and `npsvg`, which implements multiple sub-commands for preprocessing and individual steps in the different NPSV workflows.

### Prerequisites

NPSV requires the reference genome and these examples, in particular, require the "b37" reference. To obtain and index those files from within the Docker container:
```
cd /data
curl ftp://ftp.ncbi.nlm.nih.gov/1000genomes/ftp/technical/reference/human_g1k_v37.fasta.gz -o human_g1k_v37.fasta.gz
gunzip human_g1k_v37.fasta.gz
bwa index human_g1k_v37.fasta
samtools faidx human_g1k_v37.fasta
```

### Basic Workflow

The minimal NPSV workflow requires the aligned reads, the putative SV(s) as a VCF file and basic sequencing statistics (the sequencer model, read length, the mean and SD of the insert size, and depth). A minimal example follows.

*Creating the simulated replicates is more efficient when the BWA indices are loaded into shared memory prior to running NPSV (and thus doesn't need to re-loaded for each replicate). To load the BWA indices into shared memory:*
```
bwa shm /data/human_g1k_v37.fasta
```

Run NPSV genotyping:
```
mkdir -p tests/results
npsv \
    -r /data/human_g1k_v37.fasta \
    --genome etc/human_g1k_v37.genome \
    --gaps etc/human_g1k_v37.gaps.bed.gz \
    -i tests/data/12_22129565_22130387_DEL.vcf.gz \
    -b tests/data/12_22125565_22134387.bam \
    -o tests/results \
    --prefix 12_22129565_22130387_DEL.result \
    --read-length 148 --fragment-mean 573 --fragment-sd 164 --depth 25 --profile HS25 \
    --sim-ref \
    --DEL-n 50
```

This will produce a VCF file `tests/results/12_22129565_22130387_DEL.result.npsv.vcf` (determined by the output directory and prefix) with the genotypes, along with TSV files with the real and simulated features. The input variant is derived from the Genome-in-a-Bottle SV dataset; NPSV successfully genotypes this variant as homozygous alternate.

By default, NPSV uses "hybrid" mode for deletions (i.e., build per-variant classifiers trained on multiple simulated replicates of each zygosity for smaller variants and a single classifier for larger variants) and "single" mode for insertions (i.e., build just a single classifier using 1 replicate per variant per zygosity as the training data). Since this variant is a deletion < 1 kbp in length, NPSV will create a variant-specific classifier. The genotyping mode (single, variant, hybrid), classifier type, replicates, and threshold for choosing between per-variant and single classifiers in hybrid mode, are configurable for each variant type. To speed up this example we reduced the number of replicates to 50 (`--DEL-n 50`) from the default of 100.

The `--profile` argument specifies the sequencer model and thus the profile to use with the [ART NGS simulator](https://www.niehs.nih.gov/research/resources/software/biostatistics/art/index.cfm). Currently available profiles in ART are:
```
GA1 - GenomeAnalyzer I (36bp,44bp), GA2 - GenomeAnalyzer II (50bp, 75bp)
HS10 - HiSeq 1000 (100bp),          HS20 - HiSeq 2000 (100bp),      HS25 - HiSeq 2500 (125bp, 150bp)
HSXn - HiSeqX PCR free (150bp),     HSXt - HiSeqX TruSeq (150bp),   MinS - MiniSeq TruSeq (50bp)
MSv1 - MiSeq v1 (250bp),            MSv3 - MiSeq v3 (250bp),        NS50 - NextSeq500 v2 (75bp)
```

The `--sim-ref` argument is used here because the test data (`-b`) only includes a small set of the data. By default `npsv` samples random size-matched SVs from the genome to serve as the "null model" with homozygous reference genotypes, but that requires sequencing data from the whole genome. `--sim-ref` will use simulation to generate homozygous reference data.

The `--genome` file is used to determine chromosome sizes for various operations and the `--gaps` file contains regions that should not be sampled when generating random variants to model homozygous reference genotypes. The relevant files are distributed with this package for the `human_g1k_v37.fasta` and `Homo_sapiens_assembly38.genome` reference genomes.

### Preprocessing to create a "stats" file

NPSV can utilize more information about the aligned reads to improve simulation and feature extraction. The preprocessing step, run with the `preprocess` sub-command for `npsv`, will create a JSON file with the relevant stats. It can compute those stats directly, or extract many of them from the Picard metrics that may already have been generated as part of many pipelines. For example, the following command would construct the stats file from a combination of BAM file analysis with goleft and previously computed Picard metrics. Note that since this example BAM file only includes reads in a small region on chromosome 12, the results for this example command will not be meaningful.

```
npsvg preprocess \
    -r /data/human_g1k_v37.fasta \
    -b tests/data/12_22125565_22134387.bam \
    --picard-gc tests/data/gc_bias.detail_metrics \
    --picard-insert tests/data/insert_size_metrics \
    --picard-wgs tests/data/wgs_metrics \
    -o tests/results/stats.json
```

The stats file can then be used with `npsv` via the `--stats-path` option in lieu of directly specifying the sequencing statistics as shown below (here with a stats file previously generated from the entire HG002 genome).

```
npsv \
    -r /data/human_g1k_v37.fasta \
    --genome etc/human_g1k_v37.genome \
    --gaps etc/human_g1k_v37.gaps.bed.gz \
    -i tests/data/12_22129565_22130387_DEL.vcf.gz \
    -b tests/data/12_22125565_22134387.bam \
    -o tests/results \
    --prefix 12_22129565_22130387_DEL.result \
    --stats-path tests/data/stats.json --profile HS25 \
    --sim-ref \
    --DEL-n 50
```

If the Picard metrics are not available, the `preprocess` sub-command can compute the necessary metrics directly, as shown below. Note that since this example BAM file only includes reads in a small region on chromosome 12, the results for this example command will not be meaningful.

```
npsvg preprocess \
    -r /data/human_g1k_v37.fasta \
    --genome etc/human_g1k_v37.genome \
    -b tests/data/12_22125565_22134387.bam \
    -o tests/results/stats.json
```

### "End-to-end" example

The `paper` directory includes an `example.sh` script that downloads the HG002 short-read sequencing data and the GIAB SV calls, aligns the reads with BWA and then genotypes those SVs with NPSV using a representative workflow. The `benchmark.sh` script genotypes the GIAB SVs in previously aligned HG002 NGS data with NPSV along with a set of other SV genotypers. 

Aspects of both scripts are specific to the local computing infrastructure (e.g., directory paths, number of cores, executable paths) and so will need to be modified prior to use. Both scripts assume you have a customized version of [Truvari](https://github.com/mlinderm/truvari/tree/genotype_stats) installed.

## Proposing alternate SV representations

NPSV includes experimental support for automatically identifying "better" SV representations during genotyping using the simulated data. This workflow is implemented with a pre-processing step that generates possible alternate SV representations and a post-processing step that updates the genotype for the original SV from the alternate SV whose simulated data is most similar to the real data.

### Prerequisites

Variant proposal requires a BED file (`--simple-repeats-bed`) derived from the UCSC Genome Browser [simpleRepeats.txt.gz](http://hgdownload.soe.ucsc.edu/goldenPath/hg19/database/simpleRepeat.txt.gz) table dump that contains the standard BED columns plus the repeat period, number of copies and consensus repeat sequence. Alternative representations will only be generated for variants that overlap regions in this file. For convenience `simple_repeats.b37.bed.gz` and the index file (along with the hg38 version `simple_repeats.hg38.bed.gz`) are available at <http://skylight.middlebury.edu/~mlinderman/data/simple_repeats.b37.bed.gz>. To download these files in the Docker container:
```
curl -k https://www.cs.middlebury.edu/~mlinderman/data/simple_repeats.b37.bed.gz -o /data/simple_repeats.b37.bed.gz
curl -k https://www.cs.middlebury.edu/~mlinderman/data/simple_repeats.b37.bed.gz.tbi -o /data/simple_repeats.b37.bed.gz.tbi 
```

### Workflow

To generate possible alternate representations, use the `propose` sub-command for `npsvg`, e.g.

```
npsvg \
    propose \
    -r /data/human_g1k_v37.fasta \
    --simple-repeats-bed /data/simple_repeats.b37.bed.gz \
    -i tests/data/1_1865644_1866241_DEL.vcf.gz \
    -o tests/results/1_1865644_1866241_DEL.propose.vcf
```

The `tests/results/1_1865644_1866241_DEL.propose.vcf` file contains the original variant along with the proposed alternative descriptions (linked by the "INFO/ORIGINAL" field).

Then genotype the expanded set of putative variant. Note that the refinement workflow requires "variant" mode and the `--dm2` option to compute the Mahalanobis distance between the real and simulated data. Since this commands will genotype tens of putative variants, using multiple cores (if available) is recommended (see FAQ below).

```
npsv \
    -r /data/human_g1k_v37.fasta \
    --genome etc/human_g1k_v37.genome \
    --gaps etc/human_g1k_v37.gaps.bed.gz \
    -i tests/results/1_1865644_1866241_DEL.propose.vcf \
    -b tests/data/1_1861644_1871561.bam \
    --stats-path tests/data/stats.json --profile HS25 \
    -o tests/results \
    --prefix 1_1865644_1866241_DEL.propose \
    --sim-ref \
    --DEL-gt-mode variant \
    --dm2 \
    --DEL-n 50
```

Then select the best of the proposed representations with the `refine` sub-command. Refinement will update the original VCF with genotypes for the best representation.

```
npsvg \
    refine \
    --select dm2 \
    -i tests/results/1_1865644_1866241_DEL.propose.npsv.vcf \
    -o tests/results/1_1865644_1866241_DEL.propose.refined.vcf
```

When reviewing the pileup, the GIAB SV description appears to be "left shifted" from the true location as estimated from long-read sequencing data (approximately 1:1866429-1867023). NPSV (and other tools) incorrectly genotype the original SV description as homozygous reference. The NPSV proposal algorithm selects the alternative description where the actual data is most similar to simulated data for non-reference genotypes. The VCF calls produced by `refine` (shown below for this example) contain the alternate and original genotypes and PLs, the alternate and original Mahalanobis distance (smaller is better) and the alternate SV description. For this variant, `refine` selects 1:1866388-1867000 as the best SV description. The minimum non-reference distance for that SV description is 8.0, compared to 570.1 for the original description. The alternate SV description is correctly genotyped as heterozygous. 

```
GT:PL:DM:AD:CL:OGT:OPL:ODM:OAD	0/1:99,0,99:7348.1,8.0,21366.2:31,21:1_1866388_1867000_DEL:0/0:0,99,99:17.4,570.1,85419.0:67,0
```

Note that due to the random simulations the distances will differ between runs.

## FAQ

### Parallelization

`npsv` can simulate and extract variant evidence in parallel (controlled via the `--threads` option), before performing the genotyping in serial.  In "variant" mode, each variant can be genotyped independently. When employing that genotyping mode, a typical approach is to partition the input VCF file into chunks that are analyzed concurrently.

### Data Availability

The `example.sh` script in `paper` directory includes an example of downloading and preparing both the HG002 short-read sequencing data and the GIAB SV calls for use with the NPSV genotyper. Similar NGS data is available for the parental [HG003](ftp://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/data/AshkenazimTrio/HG003_NA24149_father/NIST_HiSeq_HG003_Homogeneity-12389378/HG003_HiSeq300x_fastq/140721_D00360_0044_AHA66RADXX) and [HG004](ftp://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/data//AshkenazimTrio/HG004_NA24143_mother/NIST_HiSeq_HG004_Homogeneity-14572558/HG004_HiSeq300x_fastq/140818_D00360_0046_AHA5R5ADXX) samples. The NA12878 "Platinum Genomes" NGS data is available in the European Nucleotide Archive under project [PRJEB3381](https://www.ebi.ac.uk/ena/browser/view/PRJEB3381). The Polaris SV call set is available at <https://github.com/Illumina/Polaris> and the SV-plaudit call set is available via the [supplemental materials](http://gigadb.org/dataset/100450) for the describing [publication](https://doi.org/10.1093/gigascience/giy064).

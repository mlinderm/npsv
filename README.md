# NPSV: Non-parametric Structural Variant Genotyper

## Installation

When cloning NPSV, make sure to recursively clone all of the submodules, i.e. `git clone --recursive git@github.com:mlinderm/npsv.git`.

NPSV requires Python 3.6+ and a suite of command-line genomics tools. For convenience, a Docker file is provided that installs all of the dependencies. To build that image:
```
docker build -t npsv .
```

### Manual installation

To manually install and run NPSV from the source, you will need the following dependencies:
* bwa
* bedtools
* bcftools
* samtools
* htslib (i.e. tabix and bgzip)
* samblaster
* sambamba
* goleft
along with standard command-line utilities, such as bc, gawk, etc. and a C++14 compiler. After installing the dependencies listed above, install the Python dependencies, and then NPSV itself via:
```
python3 -m pip install -r requirements.txt
python3 setup.py install
```

## Running NPSV

The full NPSV workflow has at least two steps, preprocessing and genoptyping, and is most efficient when the BWA indices are loaded into shared memory. Thus the typical approach, even when using the Docker image, is to run NPSV from a shell.

### Running the NPSV tools with Docker

The following command will start a Bash session in the Docker container (replace `/path/to/reference/directory` with the path to directory containing the reference genome and associated BWA indices). To load BWA indices into shared memory you will need to configure the Docker container with at least 10G of memory and set the shared memory size to 6G or more.

```
docker run --entrypoint /bin/bash \
    --shm-size=6g \
    -v `pwd`:/opt/npsv \
    -v /path/to/reference/directory:/data \
    -w /opt/npsv \
    -it \
    npsv
```

## NPSV Workflow

The NPSV package installs two key executables, `npsv`, the main entry point for the genotyper, and `npsvg`, which implements multiple sub-commands for preprocessing and individual steps in the different NPSV workflows.

### Basic Workflow

The minimal NPSV workflow requires the aligned reads, the putative SV(s) as a VCF file and basic sequencing statistics (the read length, the mean and SD of the insert size, and depth). A minimal example:

```
npsv \
    -r /data/human_g1k_v37.fasta \
    --genome etc/human_g1k_v37.genome \
    --gaps etc/human_g1k_v37.gaps.bed.gz \
    -i tests/data/1_1598414_1598580_DEL.vcf.gz \
    -b tests/data/1_1598414_1598580_DEL.bam \
    -o tests/results \
    --prefix 1_1598414_1598580_DEL.result \
    --read-length 148 --fragment-mean 573 --fragment-sd 164 --depth 25 \
    --n 5 \
    --sim-ref \
    --gt-mode variant \
    --classifier rf
```

This will produce a VCF file `test/results/1_1598414_1598580_DEL.result.npsv.vcf` (determined by the output directory and prefix) with the genotypes, along with TSV files with the real and simulated features.

The above will genotype the single deletion by training a unique classifier for that specific variant (`--gt-mode variant`) using 5 simulated replicates of each zygosity (`--n 5`) as the training data. Five replicates is chosen for speeding up the examples, 50-100 replicates is more typical for effective genotyping.

Creating the simulated replicates is more efficient when the BWA indices are loaded into shared memory prior to running NPSV (and thus doesn't need to re-loaded for each replicate). To do so run the following before launching NPSV:
```
bwa shm /data/human_g1k_v37.fasta
```

The `--sim-ref` argument is used here because the test data (`-b`) only includes a small set of the data. By default `npsv` samples random size-matched SVs from the genome to serve as the "null model" with homozygous reference genotypes, but that requires a whole genome sequencing data. `sim-ref` will use simulation to generate homozygous reference data.

### Preprocessing to create a "stats" file

NPSV can utilize more information about the aligned reads to improve simulation and feature extraction. The preprocessing step, run with the `preprocess` sub-command for `npsv`, will create a JSON file with the relevant stats. It can compute those stats directly, or extract them from the Picard metrics that may already have been generated as part of many pipelines. For example, the following command would 

```
npsvg preprocess \
    -r /data/human_g1k_v37.fasta \
    -b HG002-ready.bam \
    --picard-gc HG002-ready.gc_bias.detail_metrics \
    --picard-insert HG002-ready.insert_size_metrics \
    --picard-wgs HG002-ready.wgs_metrics
    -o HG002-ready.stats.json
```

construct the stats file form previously computed Picard metrics.

The stats file can then be used with `npsv` via the `--stats-path` option in lieu of directly specifying the sequencing statistics, e.g.

```
npsv \
    -r /data/human_g1k_v37.fasta \
    --genome etc/human_g1k_v37.genome \
    --gaps etc/human_g1k_v37.gaps.bed.gz \
    -i tests/data/1_1598414_1598580_DEL.vcf.gz \
    -b tests/data/1_1598414_1598580_DEL.bam \
    -o tests/results \
    --prefix 1_1598414_1598580_DEL.result \
    --stats-path tests/data/stats.json \
    --n 5 \
    --sim-ref \
    --gt-mode variant \
    --classifier rf
```

### Proposing alternate SV representations

To generate possible alternate representations, use the `propose` sub-command for `npsvg`, e.g.

```
npsvg \
    propose \
    -r /data/human_g1k_v37.fasta \
    --simple-repeats-bed /data/simple_repeats.b37.bed.gz \
    -i tests/data/1_1598414_1598580_DEL.vcf.gz \
    -o tests/results/1_1598414_1598580_DEL.propose.vcf
```

then genotype the expanded set of putative variant. Note that the refinement workflow requires the "per-variant" classification and workflow and the `--dm2` option to compute the Mahalanobis distance between the real and simulated data.

```
npsv \
    -r /data/human_g1k_v37.fasta \
    --genome etc/human_g1k_v37.genome \
    --gaps etc/human_g1k_v37.gaps.bed.gz \
    -i tests/results/1_1598414_1598580_DEL.propose.vcf \
    -b tests/data/1_1598414_1598580_DEL.bam \
    --stats-path tests/data/stats.json \
    -o tests/results \
    --prefix 1_1598414_1598580_DEL.propose \
    --n 5 \
    --sim-ref \
    --gt-mode variant \
    --classifier rf \
    --dm2
```

Then select the best of the proposed representations with the `refine` sub-command. Refinement will update the original VCF with genotypes for the best representation.

```
npsvg \
    refine \
    --select dm2 \
    -i tests/results/1_1598414_1598580_DEL.propose.npsv.vcf \
    -o tests/results/1_1598414_1598580_DEL.propose.refined.vcf
```

## FAQ

### Parallelization

`npsv` can simulate and extract variant evidence in parallel (controlled via the `--threads` option), before performing the genotyping in serial.  In "variant" mode, each variant can be genotyped independently. When employing that genotyping mode, a typical approach is to partition the input VCF file into chunks that are analyzed concurrently.
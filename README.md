
# NPSV: Non-parametric Structural Variant Genotyper

## Installation

NPSV requires Python 3.6+ and a suite of command-line genomics tools. To facilitate use, a Docker file is provided that installs all dependencies.

### Build the Docker Image

```
docker build -t npsv .
```

### Running the NPSV tools with Docker

The following will start a Bash session in the Docker container (replace `/path/to/reference/directory` with the path to directory containing the reference genome and associated BWA indices). NPSV loads the BWA index into shared memory and so you need to configure the Docker container with at least 10G of memory and set the shared memory size to 6G or more.

```
docker run --entrypoint /bin/bash \
    --shm-size=6g \
    -v `pwd`:/opt/npsv \
    -v /path/to/reference/directory:/data \
    -w /opt/npsv \
    -it \
    npsv
```

Then build the package for development:
```
python3 setup.py develop
```

and run the genotyper on the available test data:
```
mkdir -p tests/results
npsv \
    -r /data/human_g1k_v37.fasta \
    --genome etc/human_g1k_v37.genome \
    --gaps etc/human_g1k_v37.gaps.bed.gz \
    -i tests/data/1_1598414_1598580_DEL.vcf.gz \
    -b tests/data/1_1598414_1598580_DEL.bam \
    --stats-path tests/data/stats.json \
    -o tests/results \
    --prefix 1_1598414_1598580_DEL.result \
    --n 5 --reuse --sim-ref --local -c rf
```

The above command runs the genotyper on the original variants, to generate possible alternate representations, use the `propose` sub-command for `npsvg`, e.g.

```
mkdir -p tests/results
npsvg \
    propose \
    -r /data/human_g1k_v37.fasta \
    --simple-repeats-bed /data/simple_repeats.b37.bed.gz \
    -i tests/data/1_1598414_1598580_DEL.vcf.gz \
    -o tests/results/1_1598414_1598580_DEL.propose.vcf
```

then genotype the expanded set of putative variants (note the addition the `--dm2` option to compute the Mahalanobis distance between the real and simulated data):

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
    --n 5 --reuse --sim-ref --local -c rf --dm2
```

and refine the expanded set by updated the genotypes of the original variant representations:

```
npsvg \
    refine \
    -i tests/results/1_1598414_1598580_DEL.propose.npsv.vcf \
    -o tests/results/1_1598414_1598580_DEL.propose.refined.vcf
```


```
npsv \
    -r /data/human_g1k_v37.fasta \
    --genome etc/human_g1k_v37.genome \
    --gaps etc/human_g1k_v37.gaps.bed.gz \
    -i tests/data/1_931634_931634_INS.vcf.gz \
    -b tests/data/1_1598414_1598580_DEL.bam \
    --stats-path tests/data/stats.json \
    -o tests/results \
    --prefix 1_931634_931634_INS.result \
    --n 5 --reuse --sim-ref --local -c rf
```

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
python3 scripts/npsv \
    --vcf2paragraph /opt/paragraph/bin/vcf2paragraph.py \
    --grmpy /opt/paragraph/bin/grmpy \
    -r /data/human_g1k_v37.fasta \
    --genome etc/human_g1k_v37.genome \
    --gaps etc/human_g1k_v37.gaps.bed.gz \
    -i tests/data/1_1598414_1598580_DEL.vcf.gz \
    -b tests/data/1_1598414_1598580_DEL.bam \
    --stats-path tests/data/stats.json \
    -o tests/results \
    --n 5 \
    --reuse \
    --prefix 1_1598414_1598580_DEL.result \
    --sim-ref
```


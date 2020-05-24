FROM ubuntu:18.04

RUN apt-get -qq update && apt-get install --no-install-recommends -yq \
  art-nextgen-simulation-tools \
  bcftools \
  bedtools \
  build-essential \
  bwa \
  cmake \
  curl \
  gawk \
  git \
  libbz2-dev \
  liblzma-dev \
  python3-dev \
  python3-pip \
  python3-pkgconfig \
  python3-setuptools \
  samtools \
  unzip \
  zlib1g-dev \
  && \
  apt-get clean -y && \
  rm -rf /var/lib/apt/lists/*

RUN pip3 install --upgrade pip wheel

RUN curl -SL https://github.com/Illumina/paragraph/releases/download/v2.4a/paragraph-v2.4a-binary.zip \
    -o /opt/paragraph-v2.4a-binary.zip \
    && unzip -d /opt/paragraph /opt/paragraph-v2.4a-binary.zip \
    && rm -rf /opt/paragraph-v2.4a-binary.zip

RUN mkdir -p /opt/samblaster \
    && curl -SL https://github.com/GregoryFaust/samblaster/releases/download/v.0.1.25/samblaster-v.0.1.25.tar.gz \
    | tar -xzC /opt/samblaster --strip-components=1 \
    && make -C /opt/samblaster \
    && cp /opt/samblaster/samblaster /usr/local/bin/.

RUN curl -SL https://github.com/biod/sambamba/releases/download/v0.7.1/sambamba-0.7.1-linux-static.gz \
    | gzip -dc > /usr/local/bin/sambamba \
    && chmod +x /usr/local/bin/sambamba

RUN curl -SL https://github.com/brentp/goleft/releases/download/v0.2.3/goleft_linux64 \
    -o /usr/local/bin/goleft \
    && chmod +x /usr/local/bin/goleft

ADD . /opt/npsv

# Install npsv and paragraph dependencies
RUN pip3 install -r /opt/npsv/requirements.txt pysam intervaltree jsonschema

ENV PATH="/opt/npsv/scripts:/opt/paragraph/bin:${PATH}"

ENTRYPOINT ["/usr/bin/python3"]
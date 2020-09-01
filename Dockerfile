FROM ubuntu:18.04

RUN apt-get -qq update && apt-get install --no-install-recommends -yq \
  art-nextgen-simulation-tools \
  bc \
  bcftools \
  bedtools \
  build-essential \
  bwa \
  cmake \
  curl \
  gawk \
  libbz2-dev \
  liblzma-dev \
  python3-dev \
  python3-pip \
  python3-pkgconfig \
  python3-setuptools \
  samtools \
  tabix \
  unzip \
  zlib1g-dev \
  && \
  apt-get clean -y && \
  rm -rf /var/lib/apt/lists/*

RUN pip3 install --upgrade pip wheel

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

# Install npsv dependencies
RUN pip3 install -r /opt/npsv/requirements.txt

# Install npsv
WORKDIR /opt/npsv
RUN python3 setup.py install
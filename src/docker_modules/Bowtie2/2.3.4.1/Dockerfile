FROM sambamba:0.6.7
MAINTAINER Laurent Modolo

ENV BOWTIE2_VERSION=2.3.4.1
ENV SAMTOOLS_VERSION=1.7
ENV PACKAGES bowtie2=${BOWTIE2_VERSION}* \
    samtools=${SAMTOOLS_VERSION}* \
    perl=5.26.1*

RUN apt-get update && \
    apt-get install -y --no-install-recommends ${PACKAGES} && \
    apt-get clean

ENV PACKAGES curl=7.58.0* \
    ca-certificates=20180409 \
    build-essential=12.4* \
    zlib1g-dev=1:1.2.11*

RUN apt-get update && \
    apt-get install -y --no-install-recommends ${PACKAGES} && \
    apt-get clean

ENV SAMBLASTER_VERSION=0.1.24
RUN curl -k -L https://github.com/GregoryFaust/samblaster/releases/download/v.${SAMBLASTER_VERSION}/samblaster-v.${SAMBLASTER_VERSION}.tar.gz -o samblaster-v.${SAMBLASTER_VERSION}.tar.gz && \
tar xvf samblaster-v.${SAMBLASTER_VERSION}.tar.gz && \
cd samblaster-v.${SAMBLASTER_VERSION}/ && \
make && \
cp samblaster /usr/bin && \
cd .. && \
rm -R samblaster-v.${SAMBLASTER_VERSION}/

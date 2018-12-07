FROM ubuntu:18.04
MAINTAINER Laurent Modolo

ENV CUTADAPT_VERSION=1.14
ENV PACKAGES build-essential=12.4* \
             python3-pip=9.0.1* \
             python3-setuptools=39.0.1* \
             python3-dev=3.6.5* \
             python3-wheel=0.30.0*

RUN apt-get update && \
    apt-get install -y --no-install-recommends ${PACKAGES} && \
    apt-get clean

RUN pip3 install cutadapt==${CUTADAPT_VERSION}

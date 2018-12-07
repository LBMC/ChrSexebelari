FROM ubuntu:18.04
MAINTAINER Laurent Modolo

ENV DEEPTOOLS_VERSION=3.0.2
ENV PACKAGES build-essential=12.4* \
             python3-pip=9.0.1* \
             python3-setuptools=39.0.1* \
             python3-dev=3.6.5* \
             python3-wheel=0.30.0* \
             zlib1g-dev=1:1.2.11* \
             libcurl4-gnutls-dev=7.58.0*

RUN apt-get update && \
    apt-get install -y --no-install-recommends ${PACKAGES} && \
    apt-get clean

RUN pip3 install deeptools==${DEEPTOOLS_VERSION}

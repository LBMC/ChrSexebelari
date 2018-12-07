FROM ubuntu:18.04
MAINTAINER Laurent Modolo

ENV MACS_VERSION=2.1.1.20160309
ENV PACKAGES git=1:2.17* \
   build-essential=12.4* \
   python-pip=9.0.1* \
   ca-certificates=20180409 \
   python-setuptools=39.0.1* \
   python-dev=2.7* \
   python-numpy=1:1.13* \
   python-wheel=0.30.0* \
   zlib1g-dev=1:1.2.11*

RUN apt-get update && \
    apt-get install -y --no-install-recommends ${PACKAGES} && \
    apt-get clean

RUN pip install MACS2==${MACS_VERSION}


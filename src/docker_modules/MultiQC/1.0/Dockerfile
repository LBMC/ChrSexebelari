FROM ubuntu:18.04
MAINTAINER Laurent Modolo

ENV MULTIQC_VERSION=1.0
ENV PACKAGES build-essential=12.4* \
             python3-pip=9.0.1* \
             python3-setuptools=39.0.1* \
             python3-dev=3.6.5* \
             python3-wheel=0.30.0* \
             locales

RUN apt-get update && \
    apt-get install -y --no-install-recommends ${PACKAGES} && \
    apt-get clean

RUN locale-gen en_US.UTF-8
ENV LC_ALL=en_US.utf-8
ENV LANG=en_US.utf-8

RUN pip3 install multiqc==${MULTIQC_VERSION}


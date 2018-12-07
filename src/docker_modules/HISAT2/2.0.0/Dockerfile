FROM ubuntu:18.04
MAINTAINER Nicolas Fontrodona

ENV HISAT2_VERSION=2.0.0
ENV PACKAGES unzip=6.0* \
   gcc=4:7.3.0* \
   g++=4:7.3.0* \
   make=4.1* \
   curl=7.58.0* \
   ca-certificates=20180409 

RUN apt-get update && \
    apt-get install -y --no-install-recommends ${PACKAGES} && \
    apt-get clean

RUN curl -k -L http://ccb.jhu.edu/software/hisat2/downloads/hisat2-${HISAT2_VERSION}-beta-source.zip -o hisat2_linux-v${HISAT2_VERSION}.zip && \
unzip hisat2_linux-v${HISAT2_VERSION}.zip && \
cd hisat2-${HISAT2_VERSION}-beta && \
make && \
cp hisat2 /usr/bin && \
cp hisat2-* /usr/bin && \
rm -Rf hisat2-${HISAT2_VERSION}-beta 


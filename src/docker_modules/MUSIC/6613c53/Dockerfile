FROM samtools:1.7
MAINTAINER Laurent Modolo

ENV PACKAGES git=1:2.17* \
   build-essential=12.4* \
   ca-certificates=20180409 \
   zlib1g-dev=1:1.2.11*

RUN apt-get update && \
    apt-get install -y --no-install-recommends ${PACKAGES} && \
    apt-get clean

RUN git clone https://github.com/gersteinlab/MUSIC.git && \
  cd MUSIC && \
  git checkout ${MUSIC_VERSION} && \
  make clean && \
  make  && \
  cd .. && \
  mv MUSIC/bin/MUSIC /usr/bin/ && \
  mv MUSIC/bin/generate_multimappability_signal.csh /usr/bin/ && \
  mv MUSIC/bin/run_MUSIC.csh /usr/bin/ && \
  rm -Rf MUSIC

RUN chmod +x /usr/bin/*


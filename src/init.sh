#!/bin/sh

################################## Docker ######################################
# provide Docker with file_handle
cp src/file_handle/src/file_handle.py src/pipe/quality_control/
# build docker image for quality_control.nf
docker build src/pipe/quality_control -t 'quality_control:0.0.1'

# provide Docker with file_handle
cp src/file_handle/src/file_handle.py src/pipe/mapping/
# build docker image for mapping.nf
docker build src/pipe/mapping -t 'mapping:0.0.1'

# install nextflow
cd bin
wget -qO- get.nextflow.io | bash
pip install argparse
pip install docker
pip install numpy
pip2 install multiqc
pip2 install click
pip2 install jinja2
sudo pip install cutadapt
sudo pip2 install htseq
apt-get install pigz

cd /tmp
git clone https://github.com/COMBINE-lab/salmon.git
cd /tmp/salmon
mkdir build
cmake -DFETCH_BOOST=TRUE -DCMAKE_INSTALL_PREFIX=/usr/local/
make
sudo make install
rm -Rf /tmp/salmon

cd /tmp
git clone https://github.com/pachterlab/kallisto.git
cd /tmp/kallisto
mkdir build
cd build
cmake ..
make
sudo make install
rm -Rf /tmp/kallisto

cd /tmp
wget https://sourceforge.net/projects/bowtie-bio/files/bowtie2/2.3.2/bowtie2-2.3.2-linux-x86_64.zip/download
unzip download
cd bowtie2-2.3.2
sudo cp bowtie2* /usr/local/bin/
rm -Rf /tmp/bowtie2* /tmp/download

mkdir /tmp/rsem
cd /tmp/rsem
git clone https://github.com/deweylab/RSEM.git
cd RSEM
make -j 8
sudo make install
rm -Rf /tmp/rsem
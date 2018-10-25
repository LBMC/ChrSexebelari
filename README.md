# nextflow pipeline

This repository is a template and a library repository to help you build nextflow pipeline.
You can fork this repository to build your own pipeline.
To get the last commits from this repository into your fork use the following commands:

```sh
git remote add upstream gitlab_lbmc:pipelines/nextflow.git
git pull upstream master
```

## Getting Started

These instructions will get you a copy of the project up and running on your local machine for development and testing purposes. See deployment for notes on how to deploy the project on a live system.

### Prerequisites

To run nextflow on you computer you need to have java (>= 1.8) installed.

```sh
java --version
```

To be able to easily test tools already implemented for nextflow on your computer (`src/nf_modules/` to see their list). You need to have docker installed.

```sh
docker run hello-world
```

### Installing

To install nextflow on you computer simply run the following command:

```sh
src/install_nextflow.sh
```

Then to initialize a given tools run the following command:

```sh
src/docker_modules/<tool_name>/<tool_version>/docker_init.sh
```

For example to initialize `file_handle` version `0.1.1`, run:

```sh
src/docker_modules/file_handle/0.1.1/docker_init.sh
```

To initialize all the tools:
```sh
find src/docker_modules/ -name "docker_init.sh" | awk '{system($0)}'
```

## Running the tests

To run tests we first need to get a training set
```sh
cd data
git clone -c http.sslVerify=false https://gitlab.biologie.ens-lyon.fr/LBMC/tiny_dataset.git
cp tiny_dataset/fastq/tiny_R1.fastq tiny_dataset/fastq/tiny2_R1.fastq
cp tiny_dataset/fastq/tiny_R2.fastq tiny_dataset/fastq/tiny2_R2.fastq
cp tiny_dataset/fastq/tiny_S.fastq tiny_dataset/fastq/tiny2_S.fastq
cd ..
```

Then to run the tests for a given tools run the following command:

```sh
src/nf_modules/<tool_name>/<tool_version>/tests.sh
```

For example to run the tests on `Bowtie2` run:

```sh
src/nf_modules/Bowtie2/tests.sh
```

## Available tools

| tool | nf module | docker module | psmn module |
|------|:---------:|:-------------:|:----------:|
BEDtools | ok | ok | ok
BFCtools |**no**  | ok | ok
bioawk |**no**  | ok | ok
Bowtie | ok | ok | **no**
Bowtie2 | ok | ok | ok
BWA | ok | ok | ok
canu | ok | ok | ok
cutadapt | ok | ok | ok
deepTools | **no** | ok | ok
FastQC | ok | ok | ok
file_handle | **no** | ok | ok
GATK | **no** | ok | ok
HISAT2 | **no** | ok | **no**
HTSeq | ok | ok | ok
Kallisto | ok | ok | ok
MACS2 | **no** | ok | ok
MultiQC | ok | ok | ok
MUSIC | ok | ok | ok
picard | **no** | ok | ok
pigz | **no** | ok | ok
RSEM | ok | ok | ok
sambamba | ok | ok | ok
samblaster | ok | ok | ok
SAMtools | ok | ok | ok
SRAtoolkit | ok | ok | ok
Salmon | **no** | ok | ok
TopHat | **no** | ok | ok
Trimmomatic | **no** | ok | ok
UrQt | ok | ok | ok


## Contributing

Please read [CONTRIBUTING.md](CONTRIBUTING.md) for details on our code of conduct, and the process for submitting pull requests to us.

## Versioning

We use [SemVer](http://semver.org/) for versioning. For the versions available, see the [tags on this repository](https://gitlab.biologie.ens-lyon.fr/pipelines/nextflow/tags). 

## Authors

* **Laurent Modolo** - *Initial work*

See also the list of [contributors](https://gitlab.biologie.ens-lyon.fr/pipelines/nextflow/graphs/master) who participated in this project.

## License

This project is licensed under the CeCiLL License- see the [LICENSE](LICENSE) file for details


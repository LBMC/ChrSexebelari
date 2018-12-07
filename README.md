# SNP calling pipeline

SNP calling pipeline to find homozygote SNPs present in JU2817 strain but not in JU2859

## Getting Started

These instructions will get you a working version of the pipeline the SNP calling pipeline.

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

Then to initialize the necessary Docker tools with the following command:

```sh
src/docker_modules/cutadapt/1.14/docker_init.sh
src/docker_modules/UrQt/d62c1f8/docker_init.sh
src/docker_modules/bioawk/1.0/docker_init.sh
src/docker_modules/Bowtie2/2.3.4.1/docker_init.sh
src/docker_modules/sambamba/0.6.7/docker_init.sh
src/docker_modules/sambamba/0.6.7/docker_init.sh
src/docker_modules/GATK/4.0.8.1/docker_init.sh
src/docker_modules/SAMtools/1.7/docker_init.sh
src/docker_modules/bcftools/1.7/docker_init.sh
```

Necessary R packages

```sh
R -e 'install.packages(c("tidyverse", "seqinr"), repos = "https://cloud.r-project.org")'
```


### Running

To launch the analysis, you can execute the content of the script `src/1_JU28_59vs17_SNP_calling.sh`.
There, is a first section to run the pipeline locally with Docker on a training set (after generating the training set), a second to run it with Docker on the full data set, and a last seciton to run it on the PSMN.

After running the `src/SNP_calling.nf` pipeline, the `src/intersect_SNP.R` R scripts will format the `.vcf` files into `.csv` table.
The final output is filtered to keep only SNP matching a list of enzymes and SNP that are homozygote in one strain and not present in the other.

## Contributing

Please read [CONTRIBUTING.md](CONTRIBUTING.md) for details on our code of conduct, and the process for submitting pull requests to us.

## Versioning

We use [SemVer](http://semver.org/) for versioning. For the versions available, see the [tags on this repository](https://gitlab.biologie.ens-lyon.fr/pipelines/nextflow/tags). 

## Authors

* **Laurent Modolo** - *Initial work*

See also the list of [contributors](https://gitlab.biologie.ens-lyon.fr/pipelines/nextflow/graphs/master) who participated in this project.

## License

This project is licensed under the CeCiLL License- see the [LICENSE](LICENSE) file for details


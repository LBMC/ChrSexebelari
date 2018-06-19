---
title: "TP for computational biologists"
author: Laurent Modolo [laurent.modolo@ens-lyon.fr](mailto:laurent.modolo@ens-lyon.fr)
date: 20 Jun 2018
output:
pdf_document:
toc: true
toc_depth: 3
    number_sections: true
highlight: tango
    latex_engine: xelatex
---

The goal of this practical is to learn how to *wrap* tools in [Docker](https://www.docker.com/what-docker) or [Environment Module](http://www.ens-lyon.fr/PSMN/doku.php?id=documentation:tools:modules) to make them available to nextflow on a personal computer or at the [PSMN](http://www.ens-lyon.fr/PSMN/doku.php).

Here we assume that you followed the [TP for experimental biologists](./TP_experimental_biologists.md), and that you know the basics on [Docker containers](https://www.docker.com/what-container) and [Environment Module](http://www.ens-lyon.fr/PSMN/doku.php?id=documentation:tools:modules) usage. We are also going to assume that you know how to build and use a nextflow pipeline from the template [pipelines/nextflow](https://gitlab.biologie.ens-lyon.fr/pipelines/nextflow).

For the practical you can either work with the WebIDE of Gitlab, or locally as described in [git : the basis formation](https://gitlab.biologie.ens-lyon.fr/formations/git_basis).

# Docker

To run a tool within a [Docker container](https://www.docker.com/what-container) you need to write a `Dockerfile`.

[`Dockerfile`](./src/docker_modules/Kallisto/0.44.0/Dockerfile) are found in the [pipelines/nextflow](https://gitlab.biologie.ens-lyon.fr/pipelines/nextflow) project under `src/docker_modules/`. Each [`Dockerfile`](./src/docker_modules/Kallisto/0.44.0/Dockerfile) are paired with a [`docker_init.sh`](./src/docker_modules/Kallisto/0.44.0/docker_init.sh) file like following the example for `Kallisto` version `0.43.1`:

```sh
$ ls -l src/docker_modules/Kallisto/0.43.1/
total 16K                                                                        
drwxr-xr-x 2 laurent users 4.0K Jun 5 19:06 ./                                  
drwxr-xr-x 3 laurent users 4.0K Jun 6 09:49 ../                                 
-rw-r--r-- 1 laurent users  587 Jun  5 19:06 Dockerfile                          
-rwxr-xr-x 1 laurent users 79 Jun 5 19:06 docker_init.sh*                     
```

## [`docker_init.sh`](./src/docker_modules/Kallisto/0.44.0/docker_init.sh)
The [`docker_init.sh`](./src/docker_modules/Kallisto/0.44.0/docker_init.sh) is a simple sh script with the executable right (`chmod +x`).By executing this script, the user creates the [Docker container](https://www.docker.com/what-container) for the tools in a specific version. You can check the [`docker_init.sh`](./src/docker_modules/Kallisto/0.44.0/docker_init.sh) file of any implemented tools as a template. Remember that the name of the [container](https://www.docker.com/what-container) must be in lower case.

## [`Dockerfile`](./src/docker_modules/Kallisto/0.44.0/Dockerfile)

The recipe to wrap your tool in a [Docker container](https://www.docker.com/what-container) is written in a [`Dockerfile`](./src/docker_modules/Kallisto/0.44.0/Dockerfile) file.

For `Kallisto` version `0.44.0` the header of the `Dockerfile` is :

```Docker
FROM ubuntu:18.04
MAINTAINER Laurent Modolo

ENV KALLISTO_VERSION=0.44.0
```

This means that we initialize the [container](https://www.docker.com/what-container) from a bare installation of Ubuntu 18.04. You can check the ubuntu available versions [here](https://hub.docker.com/_/ubuntu/) or others operating systems like [debian](https://hub.docker.com/_/debian/) or [worst](https://hub.docker.com/r/microsoft/windowsservercore/).

Then we declare the *maintainer* of the container. Before declaring a environment variable for the container named `KALLISTO_VERSION` which contains the version of the tools wrapped. This means that this bash variable will be declared within the [container](https://www.docker.com/what-container).

You should always declare a variable `TOOLSNAME_VERSION` that contains the version number of commit number of the tools you wrap. Therefore in simple case you just have to modify this line to create a new `Dockerfile` for another version of the tool.

The following of the [`Dockerfile`](./src/docker_modules/Kallisto/0.44.0/Dockerfile) is a succession of `bash` commands executed as the **root** user within the container.
When you build your [`Dockerfile`](./src/docker_modules/Kallisto/0.44.0/Dockerfile), instead of launching many times the [`docker_init.sh`](./src/docker_modules/Kallisto/0.44.0/docker_init.sh) script you can connect to a base container in interactive mode to launch tests your commands.

```sh
docker run -it ubuntu:18.04 bash
KALLISTO_VERSION=0.44.0
```

Each `RUN` block is run sequentially by `Docker`. If there is an error or modifications in a `RUN` block, only this block and the following `RUN` will be executed.

You can learn more about the building of Docker containers [here](https://docs.docker.com/engine/reference/builder/#usage).

# SGE

To run easily tools on the PSMN, you need to build your own [Environment Module](http://www.ens-lyon.fr/PSMN/doku.php?id=documentation:tools:modules).

You can read the Contributing guide of the [PMSN/modules](https://gitlab.biologie.ens-lyon.fr/PSMN/modules) [here](https://gitlab.biologie.ens-lyon.fr/PSMN/modules/blob/master/CONTRIBUTING.md)

# Nextflow

The last step to wrap your tool, is to make it available in nextflow. For this you need to create at least 4 files, like the following for Kallisto version `0.44.0`:

```sh
ls -lR src/nf_modules/Kallisto
src/nf_modules/Kallisto/:
total 12
-rw-r--r-- 1 laurent users 866 Jun 18 17:13 kallisto.config
-rw-r--r-- 1 laurent users 2711 Jun 18 17:13 kallisto.nf
drwxr-xr-x 2 laurent users 4096 Jun 18 17:14 tests/

src/nf_modules/Kallisto/tests:
total 16
-rw-r--r-- 1 laurent users 551 Jun 18 17:14 index.nf
-rw-r--r-- 1 laurent users 901 Jun 18 17:14 mapping_paired.nf
-rw-r--r-- 1 laurent users 1037 Jun 18 17:14 mapping_single.nf
-rwxr-xr-x 1 laurent users 627 Jun 18 17:14 tests.sh*
```

The [`kallisto.config`](./src/nf_modules/Kallisto/kallisto.config) file contains intruction for two profiles : `sge` and `docker`.
The [`kallisto.nf`](./src/nf_modules/Kallisto/kallisto.nf) file contains nextflow processes to use `Kallisto`.

The [`tests/tests.sh`](./src/nf_modules/Kallisto/tests/tests.sh) script, contains a serie of nextflow calls on the other `.nf` files of the [`tests/`](./src/nf_modules/kallisto/tests/) folder. Those tests correspond to execution of the processes present in the [`kallisto.nf`](./src/nf_modules/Kallisto/kallisto.nf) file on the [LBMC/tiny_dataset](https://gitlab.biologie.ens-lyon.fr/LBMC/tiny_dataset) dataset. You can read the *Running the tests* section of the [README.md](https://gitlab.biologie.ens-lyon.fr/pipelines/nextflow/blob/master/README.md).


# 2023-amblyomma-americanum-txome-assembly

## Getting started with this repository

### Running the pipeline on an AWS EC2 instance

We used the Canonical, Ubuntu, 22.04 LTS, amd64 jammy image build on 2023-05-16 with 64 bit architecture and AMI ID ami-0f8e81a3da6e2510a.
We initially launched an m5a.large instance, and after configuration ran the pipeline on m5a.2xlarge instance type.
To set up the instance to run the pipeline, we installed and configured miniconda with the commands below.

```
curl -JLO https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh # download the miniconda installation script
bash Miniconda3-latest-Linux-x86_64.sh # run the miniconda installation script. Accept the license and follow the defaults.
source ~/.bashrc # source the .bashrc for miniconda to be available in the environment

# configure miniconda channel order
conda config --add channels defaults 
conda config --add channels bioconda
conda config --add channels conda-forge
conda config --set channel_priority strict
conda install mamba # install mamba for faster software installation.
```

### Running the snakemake pipeline

This repository uses snakemake to run the pipeline and conda to manage software environments and installations.
You can find operating system-specific instructions for installing miniconda [here](https://docs.conda.io/en/latest/miniconda.html) (see above for linux).
After installing conda and [mamba](https://mamba.readthedocs.io/en/latest/), run the following command to create the pipeline run environment.
We installaed Miniconda3 version `py311_23.5.2-0` and mamba version `1.4.9`.

```
mamba env create -n amam --file environment.yml
conda activate amam
```

To start the pipeline, run:

```
snakemake --use-conda -j 2
```

# 2023-amblyomma-americanum-txome-assembly

## Getting started with this repository

### Running the pipeline on an AWS EC2 instance

We used the Canonical, Ubuntu, 22.04 LTS, amd64 jammy image build on 2023-05-16 with 64 bit architecture and AMI ID ami-0f8e81a3da6e2510a.
We initially launched an m5a.large instance, and after configuration ran the pipeline on XXX instance type.
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

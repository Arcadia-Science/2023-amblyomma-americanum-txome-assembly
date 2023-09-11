# Starting with an older version of Ubuntu 
FROM ubuntu:18.04

# Set environment variable to allow non-interactive installations (avoiding user prompts)
ENV DEBIAN_FRONTEND=noninteractive

RUN apt-get update \
&& apt-get -y install wget curl build-essential coreutils

# If mkdir is not found in /usr/bin/, create a symbolic link to it
RUN if [ ! -f /usr/bin/mkdir ]; then ln -s $(which mkdir) /usr/bin/mkdir; fi

# Install Miniconda
RUN wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O /tmp/miniconda.sh \
&& bash /tmp/miniconda.sh -b -p /opt/conda \
&& rm /tmp/miniconda.sh

# Create the Conda environment
RUN /opt/conda/bin/conda create -y -n transfuse ruby=2.7.2 transrate-tools=1.0.0 make=4.3 gxx_linux-64=13.1.0 -c conda-forge -c bioconda -c defaults

# Explicitly set PATH for subsequent RUN commands
ENV PATH=/opt/conda/envs/transfuse/bin:$PATH

# Install the gems
RUN gem install trollop -v 2.1.3 \
&& gem install transfuse

# Install deps and catch exit code 1 which isn't real 
RUN /opt/conda/bin/conda run -n transfuse transfuse --install || true
# mv a executable that causes salmon to error out
RUN mv /opt/conda/envs/transfuse/share/rubygems/bin/librt.so.1 /opt/conda/envs/transfuse/share/rubygems/bin/so.1

# Activate the conda environment
RUN echo "source /opt/conda/etc/profile.d/conda.sh" >> ~/.bashrc && \
    echo "conda activate transfuse" >> ~/.bashrc

LABEL tool=transfuse

WORKDIR /data

# Set the default command to launch an interactive Ubuntu shell
CMD ["/bin/bash", "-l"]

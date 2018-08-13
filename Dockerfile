FROM ubuntu:18.04
MAINTAINER Rob<rob@linkage.io>
LABEL Description "https://github.com/UMN-EGGL/BuildNiceEquCab3Fasta"

ENV DEBIAN_FRONTEND=noninteractive


# Install the necessary packages ontop of base ubuntu installation 
RUN apt-get -y update && apt-get install -y \
    curl \
    lsb-release \ 
    wget \
    git \
    gcc \
    vim \
    build-essential \
    apt-transport-https \
    python3 \
    python3-dev \
    python3-pip \
    s3cmd \
    zlib1g-dev

RUN cd /root

# Install STAR
RUN mkdir -p /src/ && \
    cd /src && \
    git clone https://github.com/alexdobin/STAR.git && \ 
    cd STAR/source && \
    make STAR && \
    cp STAR /usr/local/bin

# Install miniconda
RUN wget https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh -O miniconda.sh
RUN bash miniconda.sh -b -f
ENV PATH=/root/miniconda3/bin:${PATH}
RUN conda update -n base conda

# Create a python environment
RUN conda create -y -n default python=3 
RUN /bin/bash -c "source activate default"

#RUN wget https://raw.github.com/pypa/pip/master/contrib/get-pip.py
#RUN python3 get-pip.py
RUN pip install ipython
RUN pip install locuspocus
RUN pip install snakemake
RUN pip install boto3
RUN pip install ftputil

# Install AdapterREmoval
RUN conda install -c maxibor adapterremoval2

COPY . /root/EquCab3Nice 

RUN cd /root/EquCab3Nice

ENTRYPOINT ["snakemake", "-s", "/root/EquCab3Nice/Snakefile"]


# -----------------------------
#  Build Instructions
# -----------------------------


# Build the Container with:
# $ docker build -t hga:latest .

# Run the Container passing through a port to the host
# $ docker run -p 4000:4000 -it hga


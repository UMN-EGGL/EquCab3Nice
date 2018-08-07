FROM ubuntu:18.04
MAINTAINER Rob<rob@linkage.io>
LABEL Description "https://github.com/UMN-EGGL/BuildNiceEquCab3Fasta"

# Install the necessary packages ontop of base ubuntu installation 
RUN apt-get -y update && apt-get install -y \
    curl \
	lsb-release \ 
    wget \
    git \
    gcc \
    build-essential \
	apt-transport-https \
    python3 \
    python3-dev \
    python3-pip \
    zlib1g-dev

#RUN wget https://raw.github.com/pypa/pip/master/contrib/get-pip.py
#RUN python3 get-pip.py
RUN pip3 install ipython
RUN pip3 install locuspocus

# Build the Container with:
# $ docker build -t hga:latest .

# Run the Container passing through a port to the host
# $ docker run -p 4000:4000 -it hga

# Inside the container
# $ cd BuildNiceEquCab3Fasta

# Use an official Ubuntu base image
FROM ubuntu:20.04

# Set environment variables
ENV DEBIAN_FRONTEND=noninteractive

# Update packages and install necessary dependencies
RUN apt-get update && apt-get install -y \
    wget \
    build-essential \
    automake \
    autoconf \ 
    libtool \
    pkg-config \
    git \
    ffmpeg \
    libsm6 \
    libxext6 \
    ghostscript \
    libx11-dev \
    libqt5gui5 \
    curl \
    latexmk \
    cmake \
    openmpi-bin

# Install MultiMarkdown from the official GitHub repository
RUN git clone https://github.com/fletcher/MultiMarkdown-6.git /opt/MultiMarkdown-6 \
    && cd /opt/MultiMarkdown-6 \
    && make release \
    && cd build \
    && make \
    && cp multimarkdown /usr/local/bin/

# Install Anaconda
RUN wget https://repo.continuum.io/archive/Anaconda3-2024.02-1-Linux-x86_64.sh \
    && bash Anaconda3-2024.02-1-Linux-x86_64.sh -b \
    && rm Anaconda3-2024.02-1-Linux-x86_64.sh

# Add Anaconda binaries to the PATH environment variable
ENV PATH /root/anaconda3/bin:$PATH

# Update Anaconda packages from conda-forge and bioconda channels
RUN conda config --set channel_priority false \
    && conda config --add channels conda-forge \
    && conda config --add channels bioconda \
    && conda config --remove channels defaults \
    && conda config --set channel_priority true

# Install additional packages with conda from specific channels
RUN conda install -c conda-forge cartopy -y \
    && conda install -c conda-forge jupyterlab -y

# Install greasy
WORKDIR /opt/greasy

# Clone the greasy repository
RUN git clone https://gitlab.bsc.es/support/greasy.git /opt/greasy

# Ensure m4 directory exists
RUN mkdir -p /opt/greasy/m4

# Install greasy
WORKDIR /opt/greasy

# Build and install greasy
RUN aclocal \
    && autoconf \
    && automake --add-missing \
    && ./configure --disable-doc --enable-mpi-engine

RUN make clean \
    && VERBOSE=1 make \
    && make install

# Clean up unnecessary packages
RUN apt-get autoremove -y \
    && apt-get clean -y \
    && rm -rf /var/lib/apt/lists/*

WORKDIR /tmp

# Set entrypoint or CMD to run greasy or provide any other instructions
CMD ["/bin/bash"]


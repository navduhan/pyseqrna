# Use the official Python image as base
FROM --platform=linux/amd64 ubuntu:20.04
FROM python:3.8

# Author: Naveen Duhan
# Title: Dockerfile for pySeqRNA

# Install Miniconda
RUN apt-get update && apt-get install -y wget && \
    wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O miniconda.sh && \
    bash miniconda.sh -b -p /usr/local/miniconda && \
    rm miniconda.sh

# Add Miniconda to PATH
ENV PATH="/usr/local/miniconda/bin:${PATH}"

# Initialize conda
RUN conda init bash

# Install Git
RUN apt-get install -y git

# Clone pySeqRNA repository
RUN git clone https://github.com/navduhan/pyseqrna.git

# Navigate to the pySeqRNA directory
WORKDIR /pyseqrna

# Create and activate a Conda environment
RUN conda env create -f pyseqrna_environment.yml && \
    echo "conda activate pyseqrna-0.2" >> ~/.bashrc

# Install pySeqRNA
RUN bash -c "source /usr/local/miniconda/etc/profile.d/conda.sh && \
             conda activate pyseqrna-0.2 && \
             pip install ."

# Expose port if needed
# EXPOSE 80

RUN cd ..

# Create data and output directories
RUN mkdir /data /output

# Command to run when the container starts
CMD ["bash"]

# docker build -t pyseqrna .
#docker run -it --rm pyseqrna

##################################################################
##                                                              ##
##  Written by Naveen Duhan (naveen.duhan@outlook.com)          ##
##  Kaundal Bioinformatics Lab, Utah State University           ##
##  Released under the terms of GNU General Public Licence v3   ## 
##                                                              ##
##################################################################

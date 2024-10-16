# Use the official Miniconda base image
FROM ubuntu:24.04
FROM python:3.10
FROM continuumio/miniconda3:latest

# Author: Naveen Duhan
# Title: Dockerfile for pySeqRNA

# Install system dependencies
RUN apt-get update && \
    apt-get install -y git bzip2 wget curl nano htop

# Clone the pySeqRNA repository
RUN git clone https://github.com/navduhan/pyseqrna.git /pyseqrna

# Set the working directory
WORKDIR /pyseqrna

# Install dependencies from the YAML file
RUN conda env create -f pyseqrna_environment.yml

# Activate the conda environment
# Note: The activation command is generally not needed in RUN instructions, so we directly use conda.
RUN echo "source activate pyseqrna-0.2" > ~/.bashrc

# Install pySeqRNA package
RUN /opt/conda/bin/conda run -n pyseqrna-0.2 pip install .

WORKDIR /home

# Create data and output directories
RUN mkdir -p /data /output

# Create a startup script to check architecture and append the STAR script at runtime

# Append the necessary logic to the STAR script
RUN STAR_PATH="/opt/conda/envs/pyseqrna-0.2/bin/STAR" && \
    if [ -f "$STAR_PATH" ]; then \
        echo "Appending commands to the existing STAR script..." && \
        echo 'if [ -x "${BASE}-plain" ]; then' >> "$STAR_PATH" && \
        echo '    cmd="${BASE}-plain"' >> "$STAR_PATH" && \
        echo '    "${cmd}" "$@"' >> "$STAR_PATH" && \
        echo 'else' >> "$STAR_PATH" && \
        echo '    echo "No suitable command found for this CPU. Exiting."' >> "$STAR_PATH" && \
        echo '    exit 1' >> "$STAR_PATH" && \
        echo 'fi' >> "$STAR_PATH" && \
        chmod +x "$STAR_PATH"; \
    else \
        echo "STAR script not found at $STAR_PATH"; \
        exit 1; \
    fi
# Set the default command to run the startup script, remove it, and then start bash
CMD ["bash"]

##################################################################
##                                                              ##
##  Written by Naveen Duhan (naveen.duhan@outlook.com)          ##
##  Kaundal Bioinformatics Lab, Utah State University           ##
##  Released under the terms of GNU General Public Licence v3   ## 
##                                                              ##
##################################################################
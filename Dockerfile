# Use Python base image
FROM python:3.10-slim

# Install system dependencies
RUN apt-get update && \
    apt-get install -y git build-essential python3-dev bedtools && \
    rm -rf /var/lib/apt/lists/*

# Install uv
COPY --from=ghcr.io/astral-sh/uv:latest /uv /usr/local/bin/uv

# Clone the Minda repository
RUN git clone https://github.com/shahcompbio/minda.git && \
    cd minda && git checkout 80cb5cc && \
    git rev-parse --short HEAD > /opt/minda_version.txt && \
    git rev-parse HEAD > /opt/minda_version_full.txt

# Set working directory
WORKDIR /minda

# Install Python dependencies using uv
RUN uv pip install --system pandas>=2.1.1 && \
    uv pip install --system numpy>=1.26.0 && \
    uv pip install --system pybedtools>=0.9.1 && \
    uv pip install --system intervaltree
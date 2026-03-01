# Use Python base image
FROM community.wave.seqera.io/library/uv:0.9.22--2326eebc25ee1806

# Install system dependencies
RUN apt-get update && \
    apt-get install -y git build-essential python3-dev bedtools && \
    rm -rf /var/lib/apt/lists/*

# Clone the Minda repository
RUN git clone https://github.com/shahcompbio/minda.git /opt/minda && \
    cd /opt/minda && git checkout 8b6d81c && \
    git rev-parse --short HEAD > /opt/minda_version.txt && \
    git rev-parse HEAD > /opt/minda_version_full.txt && \
    chmod +x minda.py && \
    ln -s /opt/minda/minda.py /usr/local/bin/minda

# Set working directory
WORKDIR /opt/minda

# Install Python dependencies using uv
RUN uv pip install --system --break-system-packages pandas>=2.1.1 && \
    uv pip install --system --break-system-packages numpy>=1.26.0 && \
    uv pip install --system --break-system-packages pybedtools>=0.9.1 && \
    uv pip install --system --break-system-packages intervaltree
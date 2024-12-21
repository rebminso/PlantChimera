# Base image with Python 3.12.4
FROM python:3.12.4-slim

# Set environment variables
ENV LANG=C.UTF-8 LC_ALL=C.UTF-8

# Install required dependencies and tools
RUN apt-get update && apt-get install -y \
    build-essential \
    zlib1g-dev \
    libbz2-dev \
    liblzma-dev \
    libcurl4-openssl-dev \
    libncurses5-dev \
    libncursesw5-dev \
    libssl-dev \
    wget \
    git \
    r-base \
    tzdata \
    bedtools \
    ncbi-blast+ \
    parallel \
    && apt-get clean

# Create a virtual environment and set environment path
RUN python -m venv /opt/venv
ENV PATH="/opt/venv/bin:$PATH"

# Install Python libraries
RUN pip install --upgrade pip \
    pandas==2.2.2 \
    numpy==2.1.1 \
    argparse==1.4.0 \
    pysam==0.22.1 \
    tqdm==4.66.5 \
    biopython==1.84 \
    memory-profiler==0.61.0 \
    joblib==1.3.2

# Install bwa
RUN git clone https://github.com/lh3/bwa.git \
    && cd bwa \
    && make CFLAGS="-w" \
    && cp bwa /usr/local/bin/

# Install HTSlib
RUN git clone https://github.com/samtools/htslib.git \
    && cd htslib \
    && git checkout 1.13 \
    && git submodule update --init --recursive \
    && make \
    && make install

# Install Samtools 1.13 from source
RUN git clone https://github.com/samtools/samtools.git \
    && cd samtools \
    && git checkout 1.13 \
    && make \
    && make install

# Install yq from source
RUN wget https://github.com/mikefarah/yq/releases/download/v4.18.1/yq_linux_amd64 \
    && chmod +x yq_linux_amd64 \
    && mv yq_linux_amd64 /usr/local/bin/yq

# Install R packages (GenomicFeatures, biomaRt)
RUN R -e "if (!requireNamespace('BiocManager', quietly = TRUE)) install.packages('BiocManager')" \
    && R -e "BiocManager::install(c('GenomicFeatures', 'biomaRt'))"

# Create necessary directories
RUN mkdir -p /app/Plantchimera/data /app/PlantChimera

# Set the working directory
WORKDIR /app/PlantChimera/data

# Copy the entire PlantChimera folder into the container
COPY PlantChimera /app/PlantChimera

# Make all scripts in the PlantChimera folder executable
RUN chmod -R +x /app/PlantChimera
# Copy the script into a directory in the PATH
RUN ln -s /app/PlantChimera/PlantChimera.sh /usr/local/bin/PlantChimera.sh

# Make the script executable
RUN chmod +x /usr/local/bin/PlantChimera.sh

# Set default command
CMD ["bash"]

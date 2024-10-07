# Base image with Python 3.12.4
FROM python:3.12.4-slim

# Set environment variables
ENV LANG=C.UTF-8 LC_ALL=C.UTF-8

# Install required dependencies and tools
# Install required dependencies for Samtools, HTSlib, and related tools
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
RUN pip install --upgrade pip\
    pandas==2.2.2 \
    numpy==2.1.1 \
    argparse==1.4.0 \
    pysam==0.22.1 \
    tqdm==4.66.5 \
    biopython==1.84

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

# Set working directory
WORKDIR /pipeline
RUN chmod -R 755 /pipeline
# Copy your scripts into the container
COPY ./scripts /pipeline/scripts
COPY ./PlantChimera.sh /pipeline/
COPY ./config.yaml /pipeline/

# Make the PlantChimera.sh script executable
RUN chmod +x /pipeline/PlantChimera.sh
 
# Copy the script to /usr/local/bin and rename it
RUN cp /pipeline/PlantChimera.sh /usr/local/bin/PlantChimera

# Make it executable
RUN chmod +x /usr/local/bin/plantchimera

# (Optional) Create a data directory if needed
RUN mkdir /data

# Set default command (this can be changed depending on how you want to run your pipeline)
CMD ["bash"]

FROM continuumio/miniconda3:latest

LABEL maintainer="Brown Beckley <brownbeckley94@gmail.com>"
LABEL description="EcoliTyper - Complete E. coli typing pipeline"

# Install system dependencies (only essential)
RUN apt-get update && apt-get install -y \
    procps \
    && rm -rf /var/lib/apt/lists/*

# Copy project files
WORKDIR /opt/ecolityper
COPY . /opt/ecolityper/

# Create the Conda environment from environment.yml
RUN conda env create -f environment.yml && \
    conda clean -afy

# Make the environment the default for RUN commands
SHELL ["conda", "run", "-n", "ecolityper", "/bin/bash", "-c"]

# Run abricate database setup (one-time)
RUN abricate --setupdb

# Set entrypoint
ENTRYPOINT ["conda", "run", "--no-capture-output", "-n", "ecolityper", "ecolityper"]
CMD ["-h"]

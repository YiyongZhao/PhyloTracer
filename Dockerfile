# =============================================================================
# Dockerfile for PhyloTracer
# A toolkit for gene tree rooting, gene duplication identification, ortholog
# retrieval, phylogenetic noise elimination, species hybridization detection,
# and visualization.
# =============================================================================

FROM continuumio/miniconda3

LABEL maintainer="Yiyong Zhao <yiyong.zhao@yale.edu>"
LABEL description="PhyloTracer: a toolkit for gene tree rooting, gene duplication identification, ortholog retrieval, phylogenetic noise elimination, species hybridization detection, and visualization."
LABEL version="1.0.0"
LABEL license="MIT"
LABEL url="https://github.com/YiyongZhao/PhyloTracer"

# Avoid interactive prompts during package installation
ENV DEBIAN_FRONTEND=noninteractive

# Install minimal system dependencies required by Qt/PyQt5 for rendering
RUN apt-get update && \
    apt-get install -y --no-install-recommends \
        libgl1-mesa-glx \
        libglib2.0-0 \
        libxrender1 \
        libxext6 \
        libsm6 \
        gcc \
        g++ && \
    apt-get clean && \
    rm -rf /var/lib/apt/lists/*

# Copy environment definition and install conda dependencies
COPY environment.yml /tmp/environment.yml
RUN conda env update -n base --file /tmp/environment.yml && \
    conda clean -afy && \
    rm /tmp/environment.yml

# Copy the project source and install via pip
COPY . /opt/PhyloTracer
RUN pip install --no-cache-dir /opt/PhyloTracer && \
    rm -rf /opt/PhyloTracer

# Headless rendering — required by ete3/PyQt5 when no display is available
ENV QT_QPA_PLATFORM=offscreen

# Working directory for user data (mount volumes here)
WORKDIR /data

ENTRYPOINT ["PhyloTracer"]

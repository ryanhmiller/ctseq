FROM continuumio/miniconda3

RUN conda create -n ctseqEnv python=3.7 && \
    echo "source activate ctseqEnv" > ~/.bashrc && \
    conda install -c bioconda cutadapt=1.18 && \
    conda install -c bioconda bismark=0.22.3 && \
    conda install -c bioconda samtools=1.9 && \
    conda install -c bioconda umi_tools=1.0.1 && \
    conda install -c bioconda simplesam=0.1.3.1 && \
    conda install -c conda-forge ncurses=6.1 && \
    conda install -c conda-forge openssl=1.0.2p && \
    conda install -c conda-forge r-base=3.5.1 && \
    conda install -c r r-ggplot2=3.0.0 && \
    conda install -c r r-reshape=0.8.7 && \
    conda install -c conda-forge r-pheatmap=1.0.12 && \
    apt-get install git && \
    git clone https://github.com/ryanhmiller/ctseq.git


WORKDIR /ctseq
RUN python setup.py install

WORKDIR /

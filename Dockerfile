FROM continuumio/miniconda3

RUN conda create -n ctseqEnv python=3.7 && \
    echo "source activate ctseqEnv" > ~/.bashrc && \
    conda install -c bioconda cutadapt && \
    conda install -c bioconda bismark && \
    conda install -c bioconda samtools=1.9 && \
    conda install -c bioconda umi_tools && \
    conda install -c bioconda simplesam && \
    conda install -c conda-forge ncurses=6.1=he6710b0_1 && \
    conda install -c conda-forge openssl=1.0.2p=h14c3975_1002 && \
    conda install -c conda-forge r-base=3.5.1 && \
    conda install -c r r-ggplot2 && \
    conda install -c r r-reshape && \
    conda install -c conda-forge r-pheatmap && \
    apt-get install git && \
    git clone https://github.com/ryanhmiller/ctseq.git && \
    python ./ctseq/setup.py install && \
    rm -r ctseq

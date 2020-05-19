FROM continuumio/miniconda3
RUN conda create -n ctseqEnv python=3.7
RUN echo "source activate ctseqEnv" > ~/.bashrc

RUN conda install -c bioconda cutadapt
RUN conda install -c bioconda bismark
RUN conda install -c bioconda samtools
RUN conda install -c bioconda umi_tools
RUN conda install -c bioconda simplesam
RUN conda install -c conda-forge ncurses=6.1=he6710b0_1
RUN conda install -c conda-forge openssl=1.0.2p=h14c3975_1002

RUN conda install -c conda-forge r-base=3.5.1
RUN conda install -c r r-ggplot2
RUN conda install -c r r-reshape

RUN apt-get install git

RUN git clone https://github.com/ryanhmiller/ctseq.git

WORKDIR /ctseq
RUN python setup.py install

WORKDIR /

RUN rm -r ctseq

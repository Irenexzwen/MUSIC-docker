FROM continuumio/miniconda3
#FROM mambaorg/micromamba

ENV LANG en_US.UTF-8

RUN apt-get update && \
    apt-get install -y --no-install-recommends \
        coreutils libfontconfig1-dev libfreetype-dev libfreetype6 libfreetype6-dev libarchive13\
        vim \
        less 

RUN conda config --add channels conda-forge \
    && conda config --add channels bioconda \
	&& conda install -y mamba>=0.22.1

RUN mamba install -y python=3.8 pandas openpyxl numpy gzip \
    pybktree h5py time psutil biopython Levenshtein pathlib pysam \
    samtools snakemake bowtie2 bwa cutadapt fastp pandoc

RUN mamba install -c r r-base r-stringi r-stringr r-systemfonts r-tidyverse r-data.table r-magrittr r-seqinr r-gdtools r-flextable r-knitr r-kableextra r-janitor r-markdown bioconductor-genomicranges

RUN R -e "install.packages(c('stringi'), repos = 'http://cran.us.r-project.org')" 
						   
ARG TAG=unknown
RUN mkdir -p /snakemake /res /input
COPY src/ /snakemake

WORKDIR /snakemake

ENTRYPOINT ["python", "/snakemake/run_pipeline.py"]





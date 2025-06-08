# Use the base image specified in Apptainer
FROM ubuntu:20.04

# Set non-interactive mode for apt-get
ENV DEBIAN_FRONTEND=noninteractive

# Install dependencies and clean up
RUN apt-get update && \
    apt-get install -y --no-install-recommends \
        build-essential \
        gdb \
        lcov \
        pkg-config \
        libbz2-dev \
        libffi-dev \
        libgdbm-dev \
        libgdbm-compat-dev \
        liblzma-dev \
        libncurses5-dev \
        libreadline6-dev \
        libsqlite3-dev \
        libssl-dev \
        lzma \
        lzma-dev \
        tk-dev \
        uuid-dev \
        zlib1g-dev \
        libnss3-dev \
        wget \
        dbus \
        software-properties-common \
        autoconf \
        automake \
        flex \
        bison \
        g++-11 && \
    apt-get clean && \
    rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*

# Install Python 3.11.9 from source
RUN cd /tmp && \
    wget https://www.python.org/ftp/python/3.11.9/Python-3.11.9.tgz && \
    tar -xf Python-3.11.9.tgz && \
    cd Python-3.11.9 && \
    ./configure --enable-optimizations && \
    make install && \
    cd && \
    rm -rf /tmp/*Python* /var/tmp/*

# Installing python packages
RUN pip3 install --upgrade pip --no-cache-dir && \
    pip3 install genewalk --no-cache-dir 

ENV R_VERSION=4.3.3

# Update and install dependencies for R
RUN sed -i.bak "/^#.*deb-src.*universe$/s/^# //g" /etc/apt/sources.list && \
    apt-get update && \
    apt-get install -y \
        curl \
	cmake \
        locales \
        build-essential \
        imagemagick \
        libpoppler-cpp-dev \
        libcurl4-gnutls-dev \
	gfortran \
	libpcre2-dev \
        libxml2-dev \
	libharfbuzz-dev \ 
	libfribidi-dev \
	libfreetype6-dev \ 
	libpng-dev \
	libtiff5-dev \ 
	libjpeg-dev \
	libmagick++-dev \
	libatlas3-base \ 
	liblapack3 && \
    locale-gen en_US.UTF-8 && \
    dpkg-reconfigure locales

# Install R
RUN curl -O https://cran.rstudio.com/src/base/R-4/R-${R_VERSION}.tar.gz && \
    tar -xzvf R-${R_VERSION}.tar.gz && \
    cd R-${R_VERSION} && \
    ./configure \
      --prefix=/opt/R/${R_VERSION} \
      --enable-R-shlib \
      --enable-memory-profiling \
      --with-blas \
      --with-lapack && \
    make && \
    make install && \
    ln -s /opt/R/${R_VERSION}/bin/R /usr/local/bin/R && \
    ln -s /opt/R/${R_VERSION}/bin/Rscript /usr/local/bin/Rscript

# Install R packages
RUN R -e 'install.packages("rlang", repos="https://cloud.r-project.org/")' && \
    R -e 'install.packages("ragg", repos="https://cloud.r-project.org/")' && \
    R -e 'install.packages("remotes", repos="https://cloud.r-project.org/")' && \
    R -e 'install.packages("R.utils", repos="https://cloud.r-project.org/")' && \
    R -e 'install.packages("dplyr", repos="https://cloud.r-project.org/")' && \
    R -e 'install.packages("ggplot2", repos="https://cloud.r-project.org/")' && \
    R -e 'install.packages("ggpubr", repos="https://cloud.r-project.org/")' && \
    R -e 'install.packages("ggrepel", repos="https://cloud.r-project.org/")' && \
    R -e 'install.packages("tidyverse", repos="https://cloud.r-project.org/")' && \
    R -e 'install.packages("rmarkdown", repos="https://cloud.r-project.org/")' && \
    R -e 'install.packages("pdftools", repos="https://cloud.r-project.org/")' && \
    R -e 'install.packages("magick", repos="https://cloud.r-project.org/")' && \
    R -e 'install.packages("pheatmap", repos="https://cloud.r-project.org/")' && \
    R -e 'install.packages("argparse", repos="https://cloud.r-project.org/")' && \
    R -e 'install.packages("patchwork", repos="https://cloud.r-project.org/")' && \
    R -e 'install.packages("devtools", repos="https://cloud.r-project.org/")' && \
    R -e 'install.packages("enrichR", repos="https://cloud.r-project.org/")' && \
    R -e 'devtools::install_version("BiocManager", version = "1.30.22", repos="https://cloud.r-project.org/")' && \
    R -e 'BiocManager::install(version = "3.18")' && \
    R -e 'BiocManager::install("DESeq2", version = "3.18")' && \
    R -e 'BiocManager::install("biomaRt", version = "3.18")'
    
    
# Debugging sleuth installation
RUN apt-get install -y libhdf5-dev && \
    R -e 'Sys.setenv(GITHUB_PAT = "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")' && \
    R -e 'install.packages("ellipsis", repos="https://cloud.r-project.org/")' && \
    R -e 'install.packages("ps", repos="https://cloud.r-project.org/")' && \
    R -e 'BiocManager::install(version = "3.18")' && \
    R -e 'BiocManager::install("rhdf5", version = "3.18")' && \
    R -e 'devtools::install_github("pachterlab/sleuth", auth_token="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")'

# Other packages
RUN R -e 'BiocManager::install("tximport", version = "3.18")' && \
    R -e 'BiocManager::install("tximportData", version = "3.18")' && \
    R -e 'BiocManager::install("fgsea", version = "3.18")' && \
    R -e 'BiocManager::install("EnsDb.Mmusculus.v79", force = TRUE, version = "3.18")' && \
    R -e 'BiocManager::install("EnsDb.Hsapiens.v86", version = "3.18")' && \
    R -e 'devtools::install_version("msigdbr", version = "7.5.1", repos="https://cloud.r-project.org/")' && \
    R -e 'install.packages("RColorBrewer", repos="https://cloud.r-project.org/")' && \
    R -e 'BiocManager::install("ComplexHeatmap", version = "3.18")' && \
    R -e 'install.packages("UpSetR", repos="https://cloud.r-project.org/")'
        
# Installing MAJIQ
RUN apt-get update && \
    	apt-get install -y \
	libhts-dev \
	git-all && \
	cd opt && \
	git clone https://bitbucket.org/biociphers/majiq_academic.git && \
	cd majiq_academic && \
	pip install .

# Adding MAJIQ license
COPY majiq_license_academic_official.lic majiq_license_academic_official.lic
ENV MAJIQ_LICENSE_FILE=/majiq_license_academic_official.lic

# More misc packages
RUN R -e 'BiocManager::install("EnhancedVolcano", version = "3.18")' && \
	R -e 'BiocManager::install("sva", version = "3.18")'
	
# troubleshooting fgsea
RUN R -e 'remotes::install_version("BH", version = "1.84.0-0", repos = "https://cloud.r-project.org")' && \
	R -e 'BiocManager::install("fgsea", version = "3.18")' 

# Clean up
RUN apt-get clean && \
    rm -rf /var/lib/apt/lists/* /root/.cache

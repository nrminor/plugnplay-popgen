# This Dockerfile was written with ChatGPT Feb 13 Version. Free Research Preview 
# https://help.openai.com/en/articles/6825453-chatgpt-release-notes
# Start with the latest miniconda3 image
FROM continuumio/miniconda3:latest

# Update and install necessary packages
RUN apt-get update && apt-get install -y \
	build-essential \
	curl \
	default-jre \
	wget \
	r-base \
	bcftools \
	vcftools \
	procps \
	&& apt-get clean \
	&& rm -rf /var/lib/apt/lists/*

# Install R-tidyverse
RUN Rscript -e 'install.packages("tidyverse")'

# Install necessary packages with conda
RUN conda install -y -c conda-forge -c anaconda -c bioconda \
	numpy \
	pandas \
	scipy \
	plink \
	angsd \
	&& conda clean -ya

# Install Nextflow
RUN curl -s https://get.nextflow.io | bash

# Add Nextflow to the PATH
ENV PATH="/nextflow:${PATH}"

# Install Dsuite
RUN git clone https://github.com/millanek/Dsuite.git && \
	cd Dsuite && \
	make

# Install ADMIXTURE
RUN wget https://dalexander.github.io/admixture/binaries/admixture_linux-1.3.0.tar.gz -O /tmp/admixture.tar.gz && \
	tar xvf /tmp/admixture.tar.gz -C /usr/local/bin && \
	rm /tmp/admixture.tar.gz

# Install STRUCTURE
RUN wget https://web.stanford.edu/group/pritchardlab/structure_software/release_versions/v2.3.4/release/structure_linux_console.tar.gz -O /tmp/structure.tar.gz && \
	tar xvf /tmp/structure.tar.gz -C /usr/local/bin && \
	rm /tmp/structure.tar.gz

# Install EEMS
RUN wget https://github.com/dipetkov/eems/archive/master.tar.gz -O /tmp/eems.tar.gz && \
	tar xvf /tmp/eems.tar.gz -C /usr/local/bin && \
	mv /usr/local/bin/eems-master /usr/local/bin/eems && \
	rm /tmp/eems.tar.gz

# Install stairway-plot-v2
RUN wget https://github.com/xiaoming-liu/stairway-plot-v2/blob/master/stairway_plot_v2.1.1.zip?raw=true /tmp && \
	mv /tmp/stairway_plot_v2.1.1.zip?raw=true /tmp/stairway.zip && \ 
	unzip /tmp/stairway.zip -d /usr/local/bin/stairway_plot_v2.1.1 && \
	rm /tmp/stairway.zip

# Install easySFS
RUN git clone https://github.com/isaacovercast/easySFS.git -O /usr/local/bin/easysfs && \
	chmod 777 /usr/local/bin/easysfs/easySFS.py

# Set up the working directory
WORKDIR /app

# Copy the current directory into the container
COPY . .

# Set up the entry point
ENTRYPOINT ["bash"]


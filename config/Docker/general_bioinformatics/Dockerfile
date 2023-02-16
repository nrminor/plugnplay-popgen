FROM continuumio/miniconda3

LABEL AUTHOR="Nicholas R. Minor" \
	  MAINTAINER="Nicholas R. Minor <nrminor@wisc.edu>" \
	  VERSION="1.0.0" \
	  DESCRIPTION="Docker image for bioinformatics and genomic analysis tools"

# Install basic utilities and system dependencies
RUN apt-get update && apt-get install -y build-essential wget curl unzip vim git make bzip2 zlib1g-dev liblzma-dev libncurses5-dev libbz2-dev libcurl4-openssl-dev libssl-dev libgsl-dev libboost-all-dev libhts-dev

# Install some bioinformatics and genomic analysis tools
RUN conda update -n base -c defaults conda && \
	conda install -c conda-forge mamba && \
	conda install -y python=3.9 -c conda-forge -c bioconda -c terhorst -c anaconda \
	pysam biopython samtools bcftools vcftools plink admixture structure pandas scipy numpy r-base r-tidyverse r-biocmanager bwa bowtie2 minimap2 muscle bbmap hisat2 star cufflinks stringtie kallisto salmon fastqc multiqc spades quast trimmomatic picard gatk4 freebayes ivar varscan snpeff snpsift bedtools blast sra-tools ncbi-datasets pangolin angsd dsuite smcpp nf-core nextflow && \
	git clone https://github.com/ANGSD/angsd.git && \
	cd angsd && \
	make HTSLIB_PATH=/usr/local && \
	git clone https://github.com/dipetkov/eems.git && \
	cd eems/src && \
	make && \
	wget https://sites.google.com/site/jpopgen/stairway-plot-v2/stairway_plot_v2.1.1.zip -O /tmp/stairway_plot_v2.zip && \
	unzip /tmp/stairway_plot_v2.zip -d /usr/local/bin && \
	chmod +x /usr/local/bin/stairway_plot_v2 && \
	git clone https://github.com/isaacovercast/easySFS.git && \
	git clone --recursive https://github.com/petrelharp/local_pca.git && \
	git clone --recursive https://github.com/cscherrer/souporcell.git && \
	git clone --recursive https://github.com/stephenturner/qqman.git && \
	git clone --recursive https://github.com/DReichLab/AdmixTools.git && \
	wget https://github.com/millanek/sgAdm/archive/v0.80.tar.gz && tar -xvf v0.80.tar.gz && \
	git clone https://github.com/millanek/IBD.git && \
	wget https://github.com/stschiff/msmc/archive/v1.1.0.tar.gz && tar -xvf v1.1.0.tar.gz && \
	git clone --recursive https://github.com/dportik/locuszoom.git && \
	git clone https://github.com/odelaneau/shapeit4.git && \
	wget https://github.com/igordot/genotype_imputation/raw/master/shapeit4.2.r407.linux.x64.dynamic.zip && unzip shapeit4.2.r407.linux.x64.dynamic.zip && \
	wget https://github.com/stschiff/easyDadi/archive/refs/tags/v1.3.tar.gz && tar -xvf v1.3.tar.gz && \
	https://github.com/dportik/moments/archive/v0.11.4.tar.gz && tar -xvf v0.11.4.tar.gz

RUN apt-get clean && \
		rm -rf /var/lib/apt/lists/* && \
		conda clean -ya
	
ENV PATH=/root/miniconda/bin:$PATH
	
CMD ["/bin/bash"]

FROM r-base:3.6.3

RUN apt-get update
RUN apt-get -y install cmake  python3-pip python3-dev
RUN pip install cget

RUN R -e "install.packages('data.table', version = '1.11.8', dependencies=TRUE, repos='http://cran.rstudio.com/')"

RUN mkdir /app-tools
WORKDIR /app-tools
RUN cget install --prefix /app-tools statgen/Minimac4

RUN mkdir /app
WORKDIR /app

RUN apt-get -y install bzip2
WORKDIR /app-tools
COPY required_tools/samtools-1.11.tar.bz2 ./samtools-1.11.tar.bz2
RUN bzip2 -d samtools-1.11.tar.bz2
RUN tar -xvf samtools-1.11.tar
RUN ls -l
WORKDIR /app-tools/samtools-1.11
RUN ./configure --prefix=/app-tools/
RUN make
RUN make install
ENV PATH /app-tools/bin:$PATH

WORKDIR /app-tools
COPY required_tools/bcftools-1.11.tar.bz2 ./bcftools-1.11.tar.bz2
RUN bzip2 -d bcftools-1.11.tar.bz2
RUN tar -xvf bcftools-1.11.tar
RUN ls -l
WORKDIR /app-tools/bcftools-1.11
RUN ./configure --prefix=/app-tools/
RUN make
RUN make install

WORKDIR /app-tools
COPY required_tools/plink2 bin/plink2
COPY required_tools/plink bin/plink

RUN apt-get -y update

RUN apt-get -y update
RUN apt-get -y install tabix parallel

WORKDIR /app

COPY 1_lift_vcfs_to_GRCh37.slurm.sh 1_lift_vcfs_to_GRCh37.slurm.sh
COPY 2_Genotype_Harmonizer_QC1.slurm.sh 2_Genotype_Harmonizer_QC1.slurm.sh
COPY 3_ancestry_analysis.slurm.sh 3_ancestry_analysis.slurm.sh
COPY 4_split_QC2.slurm.sh 4_split_QC2.slurm.sh
COPY 5_phase.slurm.sh 5_phase.slurm.sh
COPY 6_impute.slurm.sh 6_impute.slurm.sh

COPY required_tools/chainfiles required_tools/chainfiles
COPY required_tools/lift required_tools/lift

RUN apt-get -y install python2
RUN cp /usr/bin/python2 /usr/bin/python

ENV PATH /app/required_tools/lift:$PATH

RUN apt-get -y install default-jre

COPY required_tools/GenotypeHarmonizer required_tools/GenotypeHarmonizer

WORKDIR /app-tools
# install admixture
COPY required_tools/admixture_linux-1.3.0.tar.gz .
RUN tar -zxvf admixture_linux-1.3.0.tar.gz
ENV PATH /app-tools/admixture_linux-1.3.0:$PATH

WORKDIR /app

COPY required_tools required_tools
COPY 0_check_vcf_build.slurm.sh 0_check_vcf_build.slurm.sh

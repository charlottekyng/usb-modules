ifndef CONFIG_INC
REF ?= b37
SEQ_PLATFORM ?= illumina
SAMPLE_SET_FILE ?= sample_sets.txt
SAMPLE_FILE ?= samples.txt
SAMPLE_SPLIT_FILE ?= samples.split.txt # generate using scripts/prepareFastq.sh [remove underscores from sample names]

##Interpreters
#ANACONDA_PYTHON ?= $(HOME)/share/usr/anaconda/bin/python
#PYTHON ?= $(ANACONDA_PYTHON)
PERL_MODULE ?= Perl/5.22.2-goolf-1.7.20
PERL ?= "module load $(PERL_MODULE); perl"

R_MODULE ?= R/3.2.4-goolf-1.7.20
R ?= "module load $(R_MODULE); R"
RSCRIPT ?= "module load $(R_MODULE); Rscript"


JAVA7_MODULE ?= Java/1.7.0_75
JAVA8_MODULE ?= Java/1.8.0_74
JAVA7 ?= "module load JAVA7_MODULE; java"
JAVA8 ?= "module load JAVA8_MODULE; java"


#JAVA_BIN ?= /usr/bin/java
#JAVA_TMP := $(HOME)/tmp
#JAVA_BIN ?=  /opt/common/CentOS_6/java/jdk1.6.0_45/bin/java
#JAVA6_BIN ?= /opt/common/CentOS_6/java/jdk1.6.0_45/bin/java
#JAVA7_BIN ?= /opt/common/CentOS_6/java/jdk1.7.0_75/bin/java
#JAVA8_BIN ?= $(HOME)/usr/jdk1.8.0_74/bin/java
#JAVA ?= /usr/bin/java $(JAVA_ARGS)


# General python 2.7 environment
#ANACONDA_27_ENV ?= /home/debruiji/share/usr/anaconda-envs/anaconda-2.7
# SUFAM python environment
#SUFAM_ENV ?= /home/debruiji/share/usr/anaconda-envs/sufam-0.4.2

JARDIR ?= /scicore/home/terracci/ngkiuy/usr_nobackup/local/lib/java
# jrflab modules environment
#JRFLAB_MODULES_ENV ?= $(HOME)/share/usr/anaconda-envs/jrflab-modules-0.0.2

### Applications
UNZIP ?= /usr/bin/unzip

FASTQC_MODULE ?= FastQC/0.11.4-Java-1.7.0_21
FASTQC ?= $(PERL) fastqc ####### BROKEN

MUTECT_JAR ?= $(JARDIR)/muTect-1.1.4.jar ##### NOT UPDATED

BOWTIE_MODULE ?= Bowtie2/2.2.6-goolf-1.4.10
BOWTIE ?= "module load $(BOWTIE_MODULE); bowtie2" 

BWA_MODULE ?= BWA/0.7.12-goolf-1.4.10
BWA ?= "module load $(BWA_MODULE); bwa"

SAMTOOLS_MODULE ?= SAMtools/1.2-goolf-1.4.10
SAMTOOLS ?= module load $(SAMTOOLS_MODULE); samtools

VCFUTILS ?= $(HOME)/share/usr/bin/vcfutils.pl ##### NOT UPDATED

BCFTOOLS_MODULE ?= BCFtools/1.2-goolf-1.4.10
BCFTOOLS ?= "module load $(BCFTOOLS_MODULE); bcftools"

BEDTOOLS_MODULE ?= BEDTools/2.25.0-goolf-1.4.10
BEDTOOLS ?= "module load $(BEDTOOLS_MODULE); bedtools"

BGZIP ?= $(HOME)/share/usr/bin/bgzip #### NOT UPDATED

IGVTOOLS_MODULE ?= IGVTools/2.3.60-Java-1.7.0_21
IGVTOOLS ?= igvtools

VCFTOOLS_MODULE ?= VCFtools/0.1.14-goolf-1.4.10-Perl-5.16.3
VCFTOOLS ?= "module load $(VCFTOOLS_MODULE); vcftools"

VCF_SORT ?= $(PERL) $(HOME)/share/usr/bin/vcfsorter.pl ##### NOT UPDATED

SNP_EFF_JAR ?= $(JARDIR)/snpEff-4.1.jar
SNP_SIFT_JAR ?= $(JARDIR)/snpSift-4.1.jar
SNP_EFF_CONFIG ?= $(HOME)/share/usr/snpEff-4.1/snpEff.config
#DB_NSFP ?= $(HOME)/share/reference/snpEff/dbNSFP2.4.txt.gz
DB_NSFP ?= $(HOME)/share/reference/snpEff-4.1/dbNSFP3.0b1a.hg19.txt.gz
NSFP_FIELDS ?= Uniprot_acc_Polyphen2 Uniprot_id_Polyphen2 Polyphen2_HVAR_score Polyphen2_HVAR_pred 1000Gp3_AF ESP6500_AA_AF ESP6500_EA_AF MutationAssessor_pred MutationAssessor_score MutationTaster_pred MutationTaster_score PROVEAN_pred ExAC_Adj_AF clinvar_rs clinvar_clnsig Interpro_domain

CUFFLINKS ?= $(HOME)/share/usr/bin/cufflinks-2.2.1
CUFFCMP ?= $(HOME)/share/usr/bin/cuffcompare-2.2.1
#SNVMIX ?= $(HOME)/share/usr/bin/SNVMix2
#BAM2FASTQ ?= $(HOME)/share/usr/bin/bam2fastq
#QUALIMAP ?= $(HOME)/share/usr/qualimap/qualimap
TOPHAT ?= $(HOME)/share/usr/bin/tophat2-2.0.13
DEFUSE ?= $(PERL) $(HOME)/share/usr/defuse-0.6.1/scripts/defuse.pl

#VARSCAN
ONCOFUSE_JAR ?= $(HOME)/share/usr/oncofuse-1.0.9b2/Oncofuse.jar
VARSCAN_JAR ?= $(JARDIR)/VarScan.v2.4.1.jar

# PICARD tools
PICARD_DIR ?= $(JARDIR)/picard-1.92
PICARD_JAR ?= $(JARDIR)/picard-tools-1.141/picard.jar

# scripts
SCRIPTS_DIR ?= usb-modules/scripts
MERGE ?= $(SCRIPTS_DIR)/merge.R
VCF_TO_TABLE ?= $(SCRIPTS_DIR)/vcfToTable.R
INTRON_POSN_LOOKUP ?= $(SCRIPTS_DIR)/posnGeneLookup.pl
RBIND ?= $(SCRIPTS_DIR)/rbind.R
NORMAL_FILTER ?= $(PERL) $(SCRIPTS_DIR)/normalFilterVCF.pl
SOMATIC_FILTER_VCF ?= $(PERL) $(SCRIPTS_DIR)/somaticFilterVCF.pl
FACETS ?= $(RSCRIPT) modules/copy_number/runFacets.R

#JAVA_TMP := $(HOME)/tmp
JAVA_BIN ?=  /opt/common/CentOS_6/java/jdk1.6.0_45/bin/java
JAVA6_BIN ?= /opt/common/CentOS_6/java/jdk1.6.0_45/bin/java
JAVA7_BIN ?= /opt/common/CentOS_6/java/jdk1.7.0_75/bin/java
JAVA8_BIN ?= $(HOME)/usr/jdk1.8.0_74/bin/java
#JAVA ?= /usr/bin/java $(JAVA_ARGS)

GET_INSERT_SIZE ?= $(HOME)/share/usr/bin/getInsertSize.py

#GATK
GATK_JAR ?= $(JARDIR)/GenomeAnalysisTK.jar
GATK_JAR2 ?= $(JARDIR)/GenomeAnalysisTK-3.3-0.jar

VT ?= $(HOME)/share/usr/bin/vt
TABIX ?= $(HOME)/share/usr/bin/tabix

EMBL_MYSQLDB_HOST ?= 10.0.200.71
EMBL_MYSQLDB_PORT ?= 38493
EMBL_MYSQLDB_USER ?= embl
EMBL_MYSQLDB_PW ?= embl
EMBL_MYSQLDB_DB ?= homo_sapiens_core_75_37

INTEGRATE ?= $(HOME)/share/usr/Integrate_0.2.0/bin/Integrate

VCF2MAF ?= $(PERL) $(HOME)/share/usr/vcf2maf.pl
MAF2MAF ?= $(PERL) $(HOME)/share/usr/maf2maf.pl
MAF2VCF ?= $(PERL) $(HOME)/share/usr/maf2vcf.pl
VEP_PATH ?= /opt/common/CentOS_6/vep/v82
VEP_DATA ?= /opt/common/CentOS_6/vep/v82

HISAT ?= $(HOME)/share/usr/hisat2-2.0.2-beta/hisat2



endif
CONFIG_INC = true

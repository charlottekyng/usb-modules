ifndef CONFIG_INC

##Interpreters
SHELL = /bin/bash

#ANACONDA_PYTHON ?= $(HOME)/share/usr/anaconda/bin/python
#PYTHON ?= $(ANACONDA_PYTHON)
PERL_MODULE ?= Perl/5.22.2-goolf-1.7.20
PERL ?= module load $(PERL_MODULE); perl

R_MODULE ?= R/3.3.0-goolf-1.7.20
LOAD_R_MODULE ?= module load $(R_MODULE)
R ?= R
RSCRIPT ?= Rscript

JAVA_ARGS = -Djava.io.tmpdir=$(TMPDIR)
JAVA7_MODULE ?= Java/1.7.0_75
JAVA8_MODULE ?= Java/1.8.0_74
LOAD_JAVA7_MODULE ?= module load $(JAVA7_MODULE)
LOAD_JAVA8_MODULE ?= module load $(JAVA8_MODULE) 
JAVA7 ?= java $(JAVA_ARGS)
JAVA8 ?= java $(JAVA_ARGS)
JARDIR ?= /scicore/home/terracci/ngkiuy/usr_nobackup/local/lib/java

# General python 2.7 environment
#ANACONDA_27_ENV ?= /home/debruiji/share/usr/anaconda-envs/anaconda-2.7
# SUFAM python environment
#SUFAM_ENV ?= /home/debruiji/share/usr/anaconda-envs/sufam-0.4.2

# jrflab modules environment
#JRFLAB_MODULES_ENV ?= $(HOME)/share/usr/anaconda-envs/jrflab-modules-0.0.2

### Applications
UNZIP ?= /usr/bin/unzip

FASTQC_MODULE ?= FastQC/0.11.4-Java-1.7.0_21
LOAD_FASTQC_MODULE ?= module load $FASTQC_MODULE)
FASTQC ?= $(PERL) fastqc ####### BROKEN

MUTECT_JAR ?= $(JARDIR)/muTect-1.1.4.jar ##### NOT UPDATED

BOWTIE2_MODULE ?= Bowtie2/2.2.6-goolf-1.4.10
LOAD_BOWTIE2_MODULE ?= module load $(BOWTIE_MODULE)
BOWTIE2 ?= bowtie2 

BWA_MODULE ?= BWA/0.7.12-goolf-1.4.10
LOAD_BWA_MODULE ?= module load $(BWA_MODULE)
BWA ?= bwa

SAMTOOLS_MODULE ?= SAMtools/1.2-goolf-1.4.10
LOAD_SAMTOOLS_MODULE ?= module load $(SAMTOOLS_MODULE)
SAMTOOLS ?= samtools

VCFUTILS ?= $(HOME)/share/usr/bin/vcfutils.pl ##### NOT UPDATED

BCFTOOLS_MODULE ?= BCFtools/1.2-goolf-1.4.10
LOAD_BCFTOOLS_MODULE ?= module load $(BCFTOOLS_MODULE)
BCFTOOLS ?= bcftools

BEDTOOLS_MODULE ?= BEDTools/2.25.0-goolf-1.4.10
LOAD_BEDTOOLS_MODULE ?= module load $(BEDTOOLS_MODULE)
BEDTOOLS ?= bedtools

BGZIP ?= $(HOME)/share/usr/bin/bgzip #### NOT UPDATED

IGVTOOLS_MODULE ?= IGVTools/2.3.60-Java-1.7.0_21
LOAD_IGVTOOLS_MODULE ?= module load $(IGVTOOLS_MODULE)
IGVTOOLS ?= igvtools

VCFTOOLS_MODULE ?= VCFtools/0.1.14-goolf-1.4.10-Perl-5.16.3
LOAD_VCFTOOLS_MODULE ?= module load $(VCFTOOLS_MODULE)
VCFTOOLS ?= vcftools

VCF_SORT ?= $(PERL) $(HOME)/share/usr/bin/vcfsorter.pl ##### NOT UPDATED

SNP_EFF_JAR ?= $(JARDIR)/snpEff-4.1.jar
SNP_EFF = $(JAVA7) -Xmx$(DEFAULT_JAVA_MEM) -jar $(SNP_EFF_JAR)
SNP_SIFT_JAR ?= $(JARDIR)/snpSift-4.1.jar
SNP_SIFT = $(JAVA) -Xmx$(DEFAULT_JAVA_MEM) -jar $(SNP_SIFT_JAR)
SNP_EFF_CONFIG ?= $(HOME)/share/usr/snpEff-4.1/snpEff.config
#DB_NSFP ?= $(HOME)/share/reference/snpEff/dbNSFP2.4.txt.gz
DB_NSFP ?= $(HOME)/share/reference/snpEff-4.1/dbNSFP3.0b1a.hg19.txt.gz
NSFP_FIELDS ?= Uniprot_acc_Polyphen2 Uniprot_id_Polyphen2 Polyphen2_HVAR_score Polyphen2_HVAR_pred 1000Gp3_AF ESP6500_AA_AF ESP6500_EA_AF MutationAssessor_pred MutationAssessor_score MutationTaster_pred MutationTaster_score PROVEAN_pred ExAC_Adj_AF clinvar_rs clinvar_clnsig Interpro_domain

VCF_EFF_ONE_PER_LINE = $(HOME)/share/usr/snpEff-4.1/scripts/vcfEffOnePerLine.pl
VCF_JOIN_EFF = usb-modules/scripts/joinEff.pl

CUFFLINKS ?= $(HOME)/share/usr/bin/cufflinks-2.2.1
CUFFCMP ?= $(HOME)/share/usr/bin/cuffcompare-2.2.1
TOPHAT ?= $(HOME)/share/usr/bin/tophat2-2.0.13
DEFUSE ?= $(PERL) $(HOME)/share/usr/defuse-0.6.1/scripts/defuse.pl

#VARSCAN
ONCOFUSE_JAR ?= $(JARDIR)/Oncofuse.jar
VARSCAN_JAR ?= $(JARDIR)/varscan/VarScan.v2.4.2.jar
VARSCAN = $(JAVA8) -Xmx6G -jar $(VARSCAN_JAR)

# PICARD tools
PICARD_MEM = 10G
PICARD_OPTS = VALIDATION_STRINGENCY=LENIENT MAX_RECORDS_IN_RAM=4000000 TMP_DIR=$(TMPDIR) REFERENCE_SEQUENCE=$(REF_FASTA) CREATE_INDEX=true
PICARD_JAR ?= $(JARDIR)/picard-tools-2.4.1/picard.jar
PICARD = $(JAVA8) -Xmx$(PICARD_MEM) -jar $(PICARD_JAR) $(1) $(PICARD_OPTS)

SORT_SAM = $(call PICARD, SortSam)
REORDER_SAM = $(call PICARD, ReorderSam)
MARK_DUP = $(call PICARD, MarkDuplicates)
MARK_DUP_CIGAR = $(call PICARD, MarkDuplicatesWithMateCigar)
MERGE_SAMS = $(call PICARD, MergeSamFiles)
FIX_MATE = $(call PICARD, FixMateInformation)
SAM_TO_FASTQ = $(call PICARD, SamToFastq)
CLEANBAM = $(call PICARD, CleanSam)
ADD_RG = $(call PICARD, AddOrReplaceReadGroups)
VALIDATE_SAM = $(call PICARD, ValidateSamFile)

INTERVAL_LIST_TOOL = $(call PICARD, IntervalListTools)
CREATE_SEQ_DICT = $(call PICARD, CreateSequenceDictionary)

COLLECT_HS_METRICS = $(call PICARD, CollectHsMetrics)
COLLECT_MULT_METRICS = $(call PICARD, CollectMultipleMetrics)
COLLECT_TARGETED_METRICS = $(call PICARD, CollectTargetedPcrMetrics)
COLLECT_ALIGNMENT_METRICS = $(call PICARD, CollectAlignmentSummaryMetrics)
COLLECT_ILMN_LANE_METRICS = $(call PICARD, CollectIlluminaLaneMetrics)
COLLECT_OXOG_METRICS = $(call PICARD, CollectOxoGMetrics)
COLLECT_RNA_METRICS = $(call PICARD, CollectRnaSeqMetrics)
COLLECT_SEQ_ARTIFACT_METRICS = $(call PICARD, CollectSequencingArtifactMetrics)
COLLECT_WGS_METRICS = $(call PICARD, CollectWgsMetrics)

# scripts
SCRIPTS_DIR ?= usb-modules/scripts
MERGE ?= $(SCRIPTS_DIR)/merge.R
VCF_TO_TABLE ?= $(SCRIPTS_DIR)/vcfToTable.R
INTRON_POSN_LOOKUP ?= $(SCRIPTS_DIR)/posnGeneLookup.pl
RBIND ?= $(SCRIPTS_DIR)/rbind.R
NORMAL_FILTER ?= $(PERL) $(SCRIPTS_DIR)/normalFilterVCF.pl
SOMATIC_FILTER_VCF ?= $(PERL) $(SCRIPTS_DIR)/somaticFilterVCF.pl

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

ifndef CONFIG_INC

USR_DIR ?= /scicore/home/terracci/GROUP/usr_nobackup
JAR_DIR ?= $(USR_DIR)/local/lib/java
SCRIPTS_DIR ?= $(USR_DIR)/local/scripts
BIN_DIR ?= $(USR_DIR)/local/bin
MODULE_SCRIPTS_DIR ?= usb-modules/scripts

############################################################
########## OPTIONS MOST LIKELY NEED TO CUSTOMIZED ##########
############################################################

# default reference:
# possible values: mm10, hg19, b37, hg19_ionref, b37_hbv_hcv
REF ?= b37
# possible values: ILLUMINA, IONTORRENT
SEQ_PLATFORM ?= ILLUMINA
# possible values: NONE, BAITS, PCR, RNA
CAPTURE_METHOD ?= NONE
# possible values: NONE, AGILENT_CLINICAL_EXOME, IMPACTv5, CCP, CHPv2, IMPACT410, IMPACT341, AGILENT_ALLEXON_MOUSE, HCC, NOVARTIS_PC2
PANEL ?= NONE
# SE/PE
PAIRED_END ?= true
# for RNA-seq, possible values NONE, FIRST_READ_TRANSCRIPTION_STRAND, SECOND_READ_TRANSCRIPTION_STRAND
STRAND_SPECIFICITY ?= NONE

# possible values: CANCER, PSYCH
ANALYSIS_TYPE ?= CANCER 

## process FASTQs before alignment
TRIM_LENGTH ?= 150
# trim is fixed length trim, cutadapt
ifeq ($(strip $(CAPTURE_METHOD)),RNA)
FASTQ_FILTER ?= cutadapt
else
FASTQ_FILTER ?=
endif

CLIP_FASTQ_R1 ?=
CLIP_FASTQ_R2 ?=


MERGE_SPLIT_BAMS ?= false  # merge processed split bams
BAM_REPROCESS ?= false

# not primary alignment 256
# read fails platform/vendor quality checks 512
BAM_FILTER_FLAGS ?= 768
SAMTOOLS_DOWNSAMPLE_FACTOR ?= 0.2

ifeq ($(CAPTURE_METHOD),RNA)
BAM_DUP_TYPE ?= none
BAM_NO_FILTER ?= true
BAM_SPLIT_N_TRIM ?= true
BAM_NO_RECAL ?= false
BAM_NO_REALN ?= true
else
ifeq ($(CAPTURE_METHOD),PCR)
BAM_DUP_TYPE ?= none
BAM_NO_FILTER ?= true
BAM_SPLIT_N_TRIM ?= false
BAM_NO_RECAL ?= true
BAM_NO_REALN ?= true
else
BAM_DUP_TYPE ?= markdup
BAM_NO_FILTER ?= false
BAM_SPLIT_N_TRIM ?= false
BAM_NO_RECAL ?= false
BAM_NO_REALN ?= false
endif
BAM_DUP_TYPE ?= markdup
BAM_NO_FILTER ?= false
BAM_SPLIT_N_TRIM ?= false
BAM_NO_RECAL ?= false
BAM_NO_REALN ?= false
endif

INCLUDE_CHR_Y ?= true

### VARIANT FILTERS
ifeq ($(findstring ILLUMINA,$(SEQ_PLATFORM)),ILLUMINA)
MIN_NORMAL_DEPTH ?= 5
MIN_TUMOR_DEPTH ?= 10
MIN_TUMOR_AD ?= 3
MIN_TN_AD_RATIO ?= 5.0
endif

ifeq ($(findstring IONTORRENT,$(SEQ_PLATFORM)),IONTORRENT)
MIN_NORMAL_DEPTH ?= 10
MIN_TUMOR_DEPTH ?= 10
MIN_TUMOR_AD ?= 8
MIN_TN_AD_RATIO ?= 10.0
endif

MIN_MQ ?= 10
MIN_AF_SNP ?= 0.01
MIN_AF_INDEL ?= 0.02

HRUN ?= false
HRUN_FILTER ?= 4

USE_SUFAM ?= true

# default to medium for exomes, other values are "low" for targeted panels and "high" for WGS
include usb-modules/genome_inc/$(REF).inc
RESOURCE_REQ ?= medium
ifeq ($(findstring low,$(RESOURCE_REQ)),low)
RESOURCE_REQ_VSHORT ?= 00:29:59
RESOURCE_REQ_SHORT ?= 00:29:59
RESOURCE_REQ_MEDIUM ?= 05:59:59
RESOURCE_REQ_LONG ?= 05:59:59
RESOURCE_REQ_VLONG ?= 05:59:59
RESOURCE_REQ_LOWMEM ?= 4G
RESOURCE_REQ_MEDIUM_MEM ?= 8G
RESOURCE_REQ_HIGHMEM ?= 16G
RESOURCE_REQ_VHIGHMEM ?= 16G
RESOURCE_REQ_VVHIGHMEM ?= 32G
endif
ifeq ($(findstring medium,$(RESOURCE_REQ)),medium)
RESOURCE_REQ_VSHORT ?= 00:29:59
RESOURCE_REQ_SHORT ?= 05:29:59
RESOURCE_REQ_MEDIUM ?= 05:59:59
RESOURCE_REQ_LONG ?= 23:59:59
RESOURCE_REQ_VLONG ?= 23:59:59
RESOURCE_REQ_LOWMEM ?= 4G
RESOURCE_REQ_MEDIUM_MEM ?= 8G
RESOURCE_REQ_HIGHMEM ?= 16G
RESOURCE_REQ_VHIGHMEM ?= 32G
RESOURCE_REQ_VVHIGHMEM ?= 64G
endif
ifeq ($(findstring high,$(RESOURCE_REQ)),high)
RESOURCE_REQ_VSHORT ?= 00:29:59
RESOURCE_REQ_SHORT ?= 05:59:59
RESOURCE_REQ_MEDIUM ?= 23:59:59
RESOURCE_REQ_LONG ?= 23:59:59
RESOURCE_REQ_VLONG ?= 167:59:59
RESOURCE_REQ_LOWMEM ?= 8G
RESOURCE_REQ_MEDIUM_MEM ?= 16G
RESOURCE_REQ_HIGHMEM ?= 32G
RESOURCE_REQ_VHIGHMEM ?= 64G
RESOURCE_REQ_VVHIGHMEM ?= 128G
$(info HERE)
endif

##################################
####### Environmental var ########
##################################

export TMPDIR := /scicore/scratch
PATH "= $(USR_DIR):$(PATH)
export PATH
R_LIBS = $(HOME)/R/x86_64-pc-linux-gnu-library/3.3

#LD_LIBRARY_PATH := $(HOME)/share/usr/mysql/lib:$(LD_LIBRARY_PATH)
#export LD_LIBRARY_PATH

#PERL5LIB := $(HOME)/share/usr/perl5/lib/perl5:$(HOME)/share/usr/ActivePerl-5.18/lib
#PERL5LIB := $(HOME)/share/usr/src/bioperl-live:$(PERL5LIB)
#PERL5LIB := $(PERL5LIB):$(HOME)/share/usr/src/ensembl/modules
#PERL5LIB := $(PERL5LIB):$(HOME)/share/usr/src/ensembl-compara/modules
#PERL5LIB := $(PERL5LIB):$(HOME)/share/usr/src/ensembl-variation/modules
#PERL5LIB := $(PERL5LIB):$(HOME)/share/usr/src/ensembl-funcgen/modules
PERL5LIB := $(HOME)/perl5/lib
export PERL5LIB

ifneq ("$(wildcard config.inc)", "")
include config.inc
endif
ifneq ("$(wildcard project_config.inc)", "")
include project_config.inc
endif

##################################
########## Interpreters ##########
##################################

SHELL = /bin/bash

PERL_MODULE ?= Perl/5.22.2-goolf-1.7.20
LOAD_PERL_MODULE ?= module load $(PERL_MODULE)
PERL ?= perl

R_MODULE ?= R/3.3.0-goolf-1.7.20
LOAD_R_MODULE ?= module load $(R_MODULE)
R ?= R
RSCRIPT ?= Rscript

SWEAVE = $(RSCRIPT) $(MODULE_SCRIPTS_DIR)/Sweave.R

JAVA_ARGS = -Djava.io.tmpdir=$(TMPDIR)
JAVA6_MODULE ?= Java/1.6.0_24
JAVA7_MODULE ?= Java/1.7.0_80
JAVA8_MODULE ?= Java/1.8.0_92
LOAD_JAVA6_MODULE ?= module load $(JAVA6_MODULE)
LOAD_JAVA7_MODULE ?= module load $(JAVA7_MODULE)
LOAD_JAVA8_MODULE ?= module load $(JAVA8_MODULE) 
JAVA ?= java $(JAVA_ARGS)

PYTHON_MODULE ?= Python/3.5.1-goolf-1.4.10
LOAD_PYTHON_MODULE ?= module load $(PYTHON_MODULE)

#ANACONDA_PYTHON ?= $(HOME)/share/usr/anaconda/bin/python
#ANACONDA_27_ENV ?= /home/debruiji/share/usr/anaconda-envs/anaconda-2.7
#PYTHON ?= $(ANACONDA_PYTHON)

# SUFAM python environment
#SUFAM_ENV ?= /home/debruiji/share/usr/anaconda-envs/sufam-0.4.2

# jrflab modules environment
#JRFLAB_MODULES_ENV ?= $(HOME)/share/usr/anaconda-envs/jrflab-modules-0.0.2

##################################
#### General system tools ########
##################################

UNZIP ?= /usr/bin/unzip

##################################
#### BAM/FASTQ OPTIONS ###########
##################################

BAM_PHRED64 ?= false
SPLIT_CHR ?= true
SPLIT_FASTQ ?= false
BAM_NO_SORT ?= false
BAM_FIX_RG ?= false # fix the readgroups so that it corresponds with the bam filename

MERGE_SPLIT_FASTQ ?= true


##################################
#### General MPS tools ###########
##################################

SAMTOOLS_MODULE ?= SAMtools/1.2-goolf-1.4.10
LOAD_SAMTOOLS_MODULE ?= module load $(SAMTOOLS_MODULE)
SAMTOOLS ?= samtools

SAMTOOLS_SORT_MEM = 2000000000

#VCFUTILS ?= $(HOME)/share/usr/bin/vcfutils.pl ##### NOT UPDATED

BCFTOOLS_MODULE ?= BCFtools/1.2-goolf-1.4.10
LOAD_BCFTOOLS_MODULE ?= module load $(BCFTOOLS_MODULE)
BCFTOOLS ?= bcftools

BEDTOOLS_MODULE ?= BEDTools/2.25.0-goolf-1.4.10
LOAD_BEDTOOLS_MODULE ?= module load $(BEDTOOLS_MODULE)
BEDTOOLS ?= bedtools

IGVTOOLS_MODULE ?= IGVTools/2.3.60-Java-1.7.0_21
LOAD_IGVTOOLS_MODULE ?= module load $(IGVTOOLS_MODULE)
IGVTOOLS ?= igvtools

VCFTOOLS_MODULE ?= VCFtools/0.1.14-goolf-1.4.10-Perl-5.16.3
LOAD_VCFTOOLS_MODULE ?= module load $(VCFTOOLS_MODULE)
VCFTOOLS ?= vcftools
VCFTOOLS_MERGE ?= vcf-merge

VCF_SORT ?= $(PERL) $(SCRIPTS_DIR)/vcfsorter.pl 

##VCF2MAF ?= $(PERL) $(SCRIPTS_DIR)/vcf2maf.pl
#MAF2MAF ?= $(PERL) $(SCRIPTS_DIR)/maf2maf.pl
#MAF2VCF ?= $(PERL) $(SCRIPTS_DIR)/maf2vcf.pl
#VEP_PATH ?= /opt/common/CentOS_6/vep/v82
#VEP_DATA ?= /opt/common/CentOS_6/vep/v82

BCL2FASTQ_MODULE ?= bcl2fastq/2.17.1.14-goolf-1.7.20
LOAD_BCL2FASTQ_MODULE ?= module load $(BCL2FASTQ_MODULE)
BCL2FASTQ ?= bcl2fastq

#BAM2FASTX = $(HOME)/usr/bin/bam2fastx
EXTRACT_TOOL ?= PICARD
#BAM_TO_FASTQ ?= xxx

FASTQ_TRIMMER = $(PERL) usb-modules/fastq_tools/trimFastq.pl

TRIM_GALORE_MODULE = Trim_Galore/0.4.1-goolf-1.4.10-Perl-5.16.3
LOAD_TRIM_GALORE_MODULE ?= module load $(TRIM_GALORE_MODULE)
TRIM_GALORE ?= trim_galore

MERGE ?= $(MODULE_SCRIPTS_DIR)/merge.R
VCF_TO_TABLE ?= $(MODULE_SCRIPTS_DIR)/vcfToTable.R
INTRON_POSN_LOOKUP ?= $(MODULE_SCRIPTS_DIR)/posnGeneLookup.pl
RBIND ?= $(MODULE_SCRIPTS_DIR)/rbind.R

#GET_INSERT_SIZE ?= $(HOME)/share/usr/bin/getInsertSize.py

TABIX_MODULE ?= tabix/0.2.6-goolf-1.4.10
LOAD_TABIX_MODULE ?= module load $(TABIX_MODULE)
BGZIP ?= bgzip
TABIX ?= tabix -p vcf

BAM_READCOUNT ?= $(BIN_DIR)/bam-readcount

##################################
########## QC/metrics ############
##################################

FASTQC_MODULE ?= FastQC/0.11.4-Java-1.7.0_21
LOAD_FASTQC_MODULE ?= module load $(FASTQC_MODULE)
FASTQC ?= fastqc
FASTQC_SUMMARY_PLOT = $(RSCRIPT) usb-modules/qc/fastqcSummaryPlot.R

PLOT_RNASEQ_METRICS = $(RSCRIPT) usb-modules/qc/plotRnaseqMetrics.R
#NON_REF_FREQ = usb-modules/qc/nonRefFreqFromPileup.pl
#NON_REF_FREQ_BIN_SIZE = 0.01


##################################
########## Aligners ##############
##################################

BOWTIE2_MODULE ?= Bowtie2/2.2.6-goolf-1.4.10
LOAD_BOWTIE2_MODULE ?= module load $(BOWTIE_MODULE)
BOWTIE2 ?= bowtie2

BWA_MODULE ?= BWA/0.7.12-goolf-1.4.10
LOAD_BWA_MODULE ?= module load $(BWA_MODULE)
BWA_ALN ?= bwa aln
BWA_SAMPE ?= bwa sampe
BWA_MEM ?= bwa mem

BWA_NUM_CORES ?= 8
BWA_ALN_OPTS ?= -M

FASTQ_CHUNKS := 10
FASTQ_CHUNK_SEQ := $(shell seq 1 $(FASTQ_CHUNKS))
FASTQUTILS = $(HOME)/share/usr/ngsutils/bin/fastqutils

#TMAP = $(HOME)/share/usr/bin/tmap
#TMAP_MODE ?= map3
#TMAP_OPTS =

TOPHAT_MODULE ?= TopHat/2.1.1-goolf-1.7.20
LOAD_TOPHAT_MODULE ?= module load $(TOPHAT_MODULE)
TOPHAT ?= tophat
TOPHAT2 ?= tophat2
TOPHAT_NUM_CORES ?= 4

NO_NOVEL_SPLICING ?= false

TOPHAT_OPTS = --keep-fasta-order --no-sort-bam -G $(GENES_GTF) -p ${TOPHAT_NUM_CORES} --tmp-dir $(TMPDIR)/$*
ifeq ($(BAM_PHRED64),true)
        TOPHAT_OPTS += --solexa1.3-quals
endif
ifeq ($(NO_NOVEL_SPLICING),true)
        TOPHAT_OPTS += --no-novel-juncs
endif

STAR_MODULE ?= STAR/2.5.2a-goolf-1.7.20
LOAD_STAR_MODULE ?= module load $(STAR_MODULE)
STAR ?= STAR

STAR_FUSION_MODULE ?= STAR-Fusion/0.6.0-goolf-1.4.10-Perl-5.16.3
LOAD_STAR_FUSION_MODULE ?= module load $(STAR_FUSION_MODULE)
STAR_FUSION ?= STAR-Fusion

STAR_FUSION_MAX_PROMISCUITY ?= 3
STAR_FUSION_MIN_JUNC_SUPP ?= 2
STAR_FUSION_MIN_ALT_PCT_JUNC ?= 5

HISAT2_MODULE ?= HISAT2/2.0.4-goolf-1.7.20
LOAD_HISAT2_MODULE ?= module load $(HISAT2_MODULE)
HISAT2 ?= hisat2
HISAT_NUM_CORES ?= 8
HISAT_OPTS = -q -p $(HISAT_NUM_CORES) -x $(HISAT_REF)
ifeq ($(SEQ_PLATFORM),true)
        HISAT_OPTS += --fr
endif
ifeq ($(PHRED64),true)
        HISAT_OPTS += --phred64
endif



##################################
########## Picard/GATK ###########
##################################

PICARD_OPTS = VALIDATION_STRINGENCY=LENIENT MAX_RECORDS_IN_RAM=4000000 TMP_DIR=$(TMPDIR) REFERENCE_SEQUENCE=$(REF_FASTA) CREATE_INDEX=true
PICARD_JAR ?= $(JAR_DIR)/picard-tools-2.4.1/picard.jar
PICARD = $(JAVA) -Xmx$$(2) -jar $(PICARD_JAR) $$(1) $(PICARD_OPTS)

GATK_JAR ?= $(JAR_DIR)/GenomeAnalysisTK-3.6/GenomeAnalysisTK-3.6.jar
GATK ?= $(JAVA) -Xmx$$(2) -jar $(GATK_JAR) -S LENIENT -T $$(1) -filterNoBases

GENOTYPE_WITH_CHR2 = ?= false
BAM_CHR2_BASE_RECAL ?= false
BAM_BASE_RECAL_OPTS = -knownSites $(DBSNP) \
	$(if $(findstring true,$(BAM_CHR2_BASE_RECAL)),-L $(word 2,$(CHROMOSOMES))) \
	$(if $(TARGETS_FILE_INTERVALS),-L $(TARGETS_FILE_INTERVALS))

ifneq ($(KNOWN_INDELS),)
BAM_REALN_OPTS = --knownAlleles $(KNOWN_INDELS)
BAM_REALN_TARGET_OPTS = --known $(KNOWN_INDELS)
endif

CONTEST_OPTS ?= -T Contamination -BTI genotypes
CONTEST += -R $(REF_FASTA) -B:pop$(,)vcf $(CONTEST_REF_VCF) $(CONTEST_OPTS)

COVARIATES = -cov ReadGroupCovariate -cov QualityScoreCovariate -cov DinucCovariate -cov CycleCovariate

REPORT_STRATIFICATION := Filter

##################################
########## Variant callers #######
##################################

VALIDATION ?= false
ANNOTATE_VARIANTS ?= true

FIX_GATK_VCF = $(PERL) usb-modules/variant_callers/fix_vcf_scripts/fix_gatk_vcf.pl

MUTECT2 ?= $(call GATK,MuTect2,$$(1))

MUTECT_JAR ?= $(JAR_DIR)/muTect-1.1.4.jar
MUTECT ?= $(JAVA) -Xmx$$(1) -jar $(MUTECT_JAR) --analysis_type MuTect

MUTECT_MAX_ALT_IN_NORMAL ?= 500
MUTECT_MAX_ALT_IN_NORMAL_FRACTION ?= 0.05
MUTECT_GT_MAX_ALT_IN_NORMAL_FRACTION ?= 0.5

#MUT_FREQ_REPORT = modules/variant_callers/somatic/mutectReport.Rmd

VARSCAN_JAR ?= $(JAR_DIR)/varscan/VarScan.v2.4.3.jar
VARSCAN = $(JAVA) -Xmx6G -jar $(VARSCAN_JAR)
VARSCAN_OPTS=--output-vcf 1 --mpileup 1 --min-var-freq $(MIN_AF_SNP)
FIX_VARSCAN_VCF = $(PERL) usb-modules/variant_callers/fix_vcf_scripts/fix_varscan_vcf.pl

CONFIGURE_STRELKA = $(PERL) $(USR_DIR)/local/strelka_workflow-1.0.15/bin/configureStrelkaWorkflow.pl
STRELKA_CONFIG = usb-modules/variant_callers/somatic/strelka_config_bwa_default.ini
FIX_STRELKA_VCF = $(PERL) usb-modules/variant_callers/fix_vcf_scripts/fix_strelka_vcf.pl

TVC_ROOT_DIR ?= $(USR_DIR)/local/tvc-5.0.3-CentOS_6.5_x86_64-binary
TVC = $(TVC_ROOT_DIR)/bin/variant_caller_pipeline.py
TVC_SOMATIC_JSON = usb-modules/variant_callers/somatic/TVC_somatic_parameters.json
TVC_SENSITIVE_JSON = usb-modules/variant_callers/somatic/TVC_somatic_parameters_sensitive.json
TVC_MOTIF ?= $(TVC_ROOT_DIR)/share/TVC/sse/ampliseqexome_germline_p1_hiq_motifset.txt
USE_FPFILTER_FOR_TVC ?= true
FIX_TVC_VCF = $(PERL) usb-modules/variant_callers/fix_vcf_scripts/fix_tvc_vcf.pl

#SOMATIC_FILTER_VCF ?= $(PERL) $(SCRIPTS_DIR)/somaticFilterVCF.pl
#NORMAL_FILTER ?= $(PERL) $(SCRIPTS_DIR)/normalFilterVCF.pl

GATK_HARD_FILTER_SNPS ?= true
GATK_POOL_SNP_RECAL ?= false

VARIANT_CALL_THRESHOLD = 50
VARIANT_EMIT_THRESHOLD = 10
HAPLOTYPE_CALLER_OPTS = --dbsnp $(DBSNP) -stand_call_conf $(VARIANT_CALL_THRESHOLD) -stand_emit_conf $(VARIANT_EMIT_THRESHOLD) -R $(REF_FASTA)

VARIANT_RECAL_TRUTH_SENSITIVITY_LEVEL = 99.0
VARIANT_RECAL_ANNOTATIONS = QD MQRankSum FS MQ ReadPosRankSum

HAPCALL_SNP_MQ_THRESHOLD ?= 40.0
HAPCALL_SNP_QD_THRESHOLD ?= 2.0
HAPCALL_SNP_FS_THRESHOLD ?= 60.0
HAPCALL_SNP_HAP_SCORE_THRESHOLD ?= 13.0
HAPCALL_SNP_MQ_RANKSUM_THRESHOLD ?= -12.5
HAPCALL_SNP_READPOS_RANKSUM_THRESHOLD ?= -8.0
HAPCALL_SNP_FILTERS := --filterName 'QD' --filterExpression 'QD < $(HAPCALL_SNP_QD_THRESHOLD)' \
--filterName 'MQ' --filterExpression 'MQ < $(HAPCALL_SNP_MQ_THRESHOLD)' \
--filterName 'FS' --filterExpression 'FS > $(HAPCALL_SNP_FS_THRESHOLD)' \
--filterName 'HapScore' --filterExpression 'HaplotypeScore > $(HAPCALL_SNP_HAP_SCORE_THRESHOLD)' \
--filterName 'MQRankSum' --filterExpression 'MQRankSum < $(HAPCALL_SNP_MQ_RANKSUM_THRESHOLD)' \
--filterName 'ReadPosRankSum' --filterExpression 'ReadPosRankSum < $(HAPCALL_SNP_READPOS_RANKSUM_THRESHOLD)' \
--filterName 'Depth' --filterExpression 'DP < 5'

HAPCALL_INDEL_INBREED_COEFF_THRESHOLD ?= -0.8
HAPCALL_INDEL_QD_THRESHOLD ?= 2.0
HAPCALL_INDEL_FS_THRESHOLD ?= 200.0
HAPCALL_INDEL_HAP_SCORE_THRESHOLD ?= 13.0
HAPCALL_INDEL_MQ_RANKSUM_THRESHOLD ?= -20.0
HAPCALL_INDEL_READPOS_RANKSUM_THRESHOLD ?= -8.0

INDEL_FILTERS = --filterName 'QD' --filterExpression 'QD < $(HAPCALL_INDEL_QD_THRESHOLD)' \
--filterName 'ReadPosRankSum' --filterExpression 'ReadPosRankSum < ($HAPCALL_INDEL_MQ_RANKSUM_THRESHOLD)' \
--filterName 'InbreedingCoeff' --filterExpression 'InbreedingCoeff < $(HAPCALL_INDEL_INBREED_COEFF_THRESHOLD)'  \
--filterName 'FS' --filterExpression 'FS > $(HAPCALL_INDEL_FS_THRESHOLD)' \
--filterName 'DP' --filterExpression 'DP < 5'

VARIANT_EVAL_GATK_REPORT = $(RSCRIPT) usb-modules/variant_callers/variantEvalGatkReport.R

##################################
########## Annotation ############
##################################

SNP_EFF_MEM ?= 8g
SNP_EFF_MODULE ?= snpEff/4.1l-Java-1.7.0_75
LOAD_SNP_EFF_MODULE ?= module load $(SNP_EFF_MODULE)
SNP_EFF_DIR ?= /scicore/soft/apps/$(SNP_EFF_MODULE)/

SNP_EFF_JAR ?= $(SNP_EFF_DIR)/snpEff.jar
SNP_EFF = $(JAVA) -Xmx$(SNP_EFF_MEM) -jar $(SNP_EFF_JAR)
SNP_SIFT_JAR ?= $(SNP_EFF_DIR)/SnpSift.jar
SNP_SIFT = $(JAVA) -Xmx$(SNP_EFF_MEM) -jar $(SNP_SIFT_JAR)

SNP_EFF_CONFIG ?= usb-modules/snpEff.config
SNP_EFF_FLAGS ?= -canon -no NEXT_PROT -no-intergenic
SNP_EFF_OPTS = -c $(SNP_EFF_CONFIG) -i vcf -o vcf $(SNP_EFF_FLAGS)
SNP_SIFT_OPTS = -c $(SNP_EFF_CONFIG)

EXAC_INFO_FIELDS ?= AC_AFR AC_AMR AC_Adj AC_EAS AC_FIN AC_Hemi AC_Het AC_Hom AC_NFE AC_OTH AC_SAS AN_AFR AN_AMR AN_Adj AN_EAS AN_FIN AN_NFE AN_OTH AN_SAS Hemi_AFR Hemi_AMR Hemi_EAS Hemi_FIN Hemi_NFE Hemi_OTH Hemi_SAS Het_AFR Het_AMR Het_EAS Het_FIN Het_NFE Het_OTH Het_SAS Hom_AFR Hom_AMR Hom_EAS Hom_FIN Hom_NFE Hom_OTH Hom_SAS AC_MALE AC_FEMALE AN_MALE AN_FEMALE
VCF_FIELDS = CHROM POS ID REF ALT FILTER
#NSFP_FIELDS ?= ExAC_Adj_AF 1000Gp3_AF ESP6500_AA_AF ESP6500_EA_AF \
#	clinvar_rs clinvar_clnsig \
#	Uniprot_acc_Polyphen2 Uniprot_id_Polyphen2 Polyphen2_HVAR_pred \
#	MutationAssessor_pred MutationTaster_pred \
#	PROVEAN_pred Interpro_domain
ANN_FIELDS = $(addprefix ANN[*].,ALLELE EFFECT IMPACT GENE GENEID FEATURE FEATUREID BIOTYPE RANK \
	HGVS_C HGVS_P CDNA_POS CDNA_LEN CDS_POS CDS_LEN AA_POS AA_LEN DISTANCE ERRORS)

VCF_EFF_ONE_PER_LINE = $(PERL) $(SNP_EFF_DIR)/scripts/vcfEffOnePerLine.pl
VCF_JOIN_EFF = $(PERL) $(MODULE_SCRIPTS_DIR)/joinEff.pl

#ANNOVAR = $(PERL) $(HOME)/share/usr/annovar/table_annovar.pl
#ANNOVAR_PROTOCOL ?= refGene$(,)cytoBand$(,)genomicSuperDups$(,)esp6500siv2_all$(,)1000g2014oct_all$(,)1000g2014oct_afr$(,)1000g2014oct_eas$(,)1000g2014oct_eur$(,)snp138$(,)ljb26_all
#ANNOVAR_OPERATION ?= g$(,)r$(,)r$(,)f$(,)f$(,)f$(,)f$(,)f$(,)f$(,)f
#ANNOVAR_OPTS = --dot2underline -remove -protocol $(ANNOVAR_PROTOCOL) -operation $(ANNOVAR_OPERATION) -nastring . -vcfinput -buildver $(ANNOVAR_REF)


#CHASM = $(RSCRIPT) modules/vcf_tools/chasmVcf.R
#CHASM_DIR = /ifs/opt/common/CHASM/CHASMDL.1.0.7
#CHASM_DIR = $(HOME)/share/usr/CHASM
#CHASM_PYTHON_ENV = $(HOME)/share/usr/anaconda-envs/pyenv27-chasm
#CHASM_CLASSIFIER ?= Breast

#FATHMM = $(MY_RSCRIPT) modules/vcf_tools/fathmmVcf.R
#FATHMM_DIR = $(HOME)/share/usr/fathmm
#FATHMM_PYTHON = $(HOME)/share/usr/bin/python
#FATHMM_PYTHONPATH = $(HOME)/share/usr/lib/python:$(HOME)/share/usr/lib/python2.7
#FATHMM_OPTS = --genome $(REF) --ensemblTxdb $(ENSEMBL_TXDB) --ref $(REF_FASTA) --fathmmDir $(FATHMM_DIR) --python $(FATHMM_PYTHON) \
#	--mysqlHost $(EMBL_MYSQLDB_HOST) --mysqlPort $(EMBL_MYSQLDB_PORT) --mysqlUser $(EMBL_MYSQLDB_USER) \
#	--mysqlPassword $(EMBL_MYSQLDB_PW) --mysqlDb (EMBL_MYSQLDB_DB)
#TRANSFIC = $(RSCRIPT) modules/vcf_tools/transficVcf.R
#TRANSFIC_PERL_SCRIPT = $(HOME)/share/usr/transfic/bin/transf_scores.pl

#PROVEAN = $(RSCRIPT) modules/vcf_tools/proveanVcf.R
#PROVEAN_OPTS = --genome $(REF) --aaTable $(AA_TABLE) --ensemblTxdb $(ENSEMBL_TXDB) --mysqlHost $(EMBL_MYSQLDB_HOST) \
#			   --mysqlPort $(EMBL_MYSQLDB_PORT) --mysqlUser $(EMBL_MYSQLDB_USER) --mysqlPassword $(EMBL_MYSQLDB_PW) \
#			   --mysqlDb $(EMBL_MYSQLDB_DB) --numThreads 8 --memPerThread 1G --queue $(QUEUE) --qsubPriority $(QSUB_PRIORITY)

#AA_TABLE = $(HOME)/share/reference/aa_table.tsv

#MUT_ASS = $(RSCRIPT) usb-modules/vcf_tools/mutAssVcf.R
#MUTATION_TASTER = $(PYTHON) usb-modules/vcf_tools/mutation_taster_vcf.py

ADD_GENE_LIST_ANNOTATION = $(RSCRIPT) usb-modules/vcf_tools/addGeneListAnnotationToVcf.R

#NON_SILENT_EFF = START_GAINED SPLICE_SITE_ACCEPTOR SPLICE_SITE_DONOR START_LOST NON_SYNONYMOUS_CODING FRAME_SHIFT CODON_CHANGE CODON_INSERTION CODON_CHANGE_PLUS_CODON_INSERTION CODON_DELETION CODON_CHANGE_PLUS_CODON_DELETION STOP_GAINED STOP_LOST NON_SYNONYMOUS_START
#NON_SILENT_CODING_EFF = START_GAINED START_LOST NON_SYNONYMOUS_CODING FRAME_SHIFT CODON_CHANGE CODON_INSERTION CODON_CHANGE_PLUS_CODON_INSERTION CODON_DELETION CODON_CHANGE_PLUS_CODON_DELETION STOP_GAINED STOP_LOST NON_SYNONYMOUS_START

#EMBL_MYSQLDB_HOST ?= 10.0.200.71
#EMBL_MYSQLDB_PORT ?= 38493
#EMBL_MYSQLDB_USER ?= embl
#EMBL_MYSQLDB_PW ?= embl
#EMBL_MYSQLDB_DB ?= homo_sapiens_core_75_37

ANNOTATE_FACETS_VCF = $(RSCRIPT) usb-modules/copy_number/annotateFacets2Vcf.R

MUTATION_SUMMARY_FORMAT ?= EXCEL
INCLUDE_LINCRNA_IN_SUMMARY ?= false

#NORMAL_VCF ?= $(HOME)/share/reference/spowellnormal.gatk_variants.vcf
#CANCER_HOTSPOT_ANNOTATION_SCRIPT = python modules/vcf_tools/hotspot_vcf.py
#CLASSIFY_PATHOGENICITY = $(PYTHON) modules/vcf_tools/classify_pathogenicity_vcf.py
#COMMON_FILTER_VCF = $(PYTHON) modules/vcf_tools/common_filter_vcf.py
##################################
########## Copy number ###########
##################################

FACETS_GATK_VARIANTS ?= false
FACETS_SNP_PILEUP = $(BIN_DIR)/snp-pileup
FACETS_MINGC ?= 0
FACETS_MAXGC ?= 1
ifeq ($(findstring ILLUMINA,$(SEQ_PLATFORM)),ILLUMINA)
FACETS_SNP_PILEUP_MIN_DEPTH ?= 25
FACETS_SNP_PILEUP_MAX_DEPTH ?= 1000
FACETS_WINDOW_SIZE ?= 200
endif

ifeq ($(findstring IONTORRENT,$(SEQ_PLATFORM)),IONTORRENT)
FACETS_SNP_PILEUP_MIN_DEPTH ?= 50
FACETS_SNP_PILEUP_MAX_DEPTH ?= 3000
FACETS_WINDOW_SIZE ?= 10
endif

FACETS_SNP_PILEUP_MINMAPQ ?= 1
FACETS_SNP_PILEUP_MINBASEQ ?= 13

FACETS_PRE_CVAL ?=
FACETS_CVAL1 ?= 75 100 200 300 400
FACETS_CVAL2 ?=
FACETS_MIN_NHET ?= 15
FACETS = $(RSCRIPT) usb-modules/copy_number/runFacets.R

FACETS_GENE_CN = $(RSCRIPT) usb-modules/copy_number/facetsGeneCN.R
FACETS_FILL_GENE_CN = $(RSCRIPT) usb-modules/copy_number/facetsFillGeneCN.R
FACETS_GENE_CN_OPTS = $(if $(TARGETS_FILE_GENES),--genesFile $(TARGETS_FILE_GENES))
                                     #     --mysqlHost $(EMBL_MYSQLDB_HOST) --mysqlPort $(EMBL_MYSQLDB_PORT) \
                                     #     --mysqlUser $(EMBL_MYSQLDB_USER) --mysqlPassword $(EMBL_MYSQLDB_PW) \
                                     #     --mysqlDb $(EMBL_MYSQLDB_DB)

FACETS_PLOT_GENE_CN = $(RSCRIPT) usb-modules/copy_number/facetsGeneCNPlot.R
FACETS_PLOT_GENE_CN_OPTS = --sampleColumnPostFix '_EM'

ANN_FACETS ?= false

CBS_SEGMENTCNV = $(RSCRIPT) usb-modules/copy_number/segmentVarscanCNV.R

CBS_SEG_SD ?= 2
CBS_SEG_SMOOTH ?= 10
CBS_SEG_ALPHA ?= 0.000001

CBS_MULTIPARAM_SEGMENT ?= false
CBS_SEG_SDS ?= 1.5 2 2.5 3
CBS_SEG_SMOOTHS ?= 5 10
CBS_SEG_ALPHAS ?= 0.01 0.0001 0.000001 0.0000000001

VARSCAN_GENE_CN = $(RSCRIPT) usb-modules/copy_number/varscanCNVGeneCN.R
VARSCAN_GENE_CN_OPTS = $(if $(GENES_FILE),--genesFile $(GENES_FILE))
CBS_MPILEUP_OPTS = -q 1 $(if $(TARGETS_FILE_INTERVALS),-l $(TARGETS_FILE_INTERVALS))

#################################
########## MutSigCV #############


MUTSIGCV_DIR ?= $(USR_DIR)/local/MutSigCV_1.4
MUTSIGCV ?= $(MUTSIGCV_DIR)/run_MutSigCV.sh
MCR ?= $(USR_DIR)/local/MCR_R2013a_glnxa64/v81/

##################################
############# GISTIC ############
GISTIC ?= $(USR_DIR)/local/gistic_2.0.23/gp_gistic2_from_seg
#GISTIC = $(HOME)/share/usr/gistic_2_0_21/gp_gistic2_from_seg
GISTIC_THRESHOLD ?= 0.3
GISTIC_JS ?= 15
GISTIC_OPTS = -genegistic 0 -smallmem 1 -maxseg 5000 -savegene 1 -saveseg 1 -savedata 0 -v 30 -ta $(GISTIC_THRESHOLD) -td $(GISTIC_THRESHOLD) -js $(GISTIC_JS) -qvt 0.25 -conf 0.99 -broad 1 -brlen 0.5 -rx 0


##################################
########## Gene expression #######
##################################

CUFFLINKS_MODULE ?= Cufflinks/2.2.1-Linux_x86_64
LOAD_CUFFLINKS_MODULE ?= module load $(CUFFLINKS_MODULE)
CUFFLINKS ?= cufflinks
CUFFCMP ?= currcompare

STRINGTIE_MODULE ?= StringTie/1.3.3-goolf-1.7.20
LOAD_STRINGTIE_MODULE ?= module load $(STRINGTIE_MODULE)

DESEQ_CONDITION ?= condition
DESEQ_REF_CONDITION ?= ref
DESEQ_ALT_CONDITION ?= alt
DESEQ_PHENO_FILE ?= pheno.txt

#CUFFLINKS ?= $(HOME)/share/usr/bin/cufflinks-2.2.1
#CUFFCMP ?= $(HOME)/share/usr/bin/cuffcompare-2.2.1

##################################
########## Fusions ###############
##################################

#DEFUSE ?= $(PERL) $(HOME)/share/usr/defuse-0.6.1/scripts/defuse.pl
#INTEGRATE ?= $(HOME)/share/usr/Integrate_0.2.0/bin/Integrate
#ONCOFUSE_JAR ?= $(JARDIR)/Oncofuse.jar
endif
CONFIG_INC = true

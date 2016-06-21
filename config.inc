ifndef CONFIG_INC


JARDIR ?= $(HOME)/usr_nobackup/local/lib/java
SCRIPTSDIR ?= $(HOME)/usr_nobackup/local/scripts
BINDIR ?= $(HOME)/usr_nobackup/local/bin

############################################################
########## OPTIONS MOST LIKELY NEED TO CUSTOMIZED ##########
############################################################

# default reference:
# possible values: mm10, hg19, b37
REF ?= b37
# possible values: ILLUMINA, IONTORRENT
SEQ_PLATFORM ?= ILLUMINA
# possible values: NONE, BAITS, PCR
CAPTURE_METHOD ?= NONE
# possible values: NONE, EXOME, IMPACTv5, CCP, CHPv2, IMPACT410, IMPACT341
PANEL ?= NONE


MERGE_SPLIT_BAMS ?= false  # merge processed split bams
BAM_REPROCESS ?= false

TRIM_READS ?= false
FASTQ_FILTER := trim
TRIM_LENGTH ?= 150
TRIM_OPTS ?= -l $(TRIM_LENGTH)

BAM_DUP_TYPE ?= markdup
BAM_NO_FILTER ?= false
BAM_NO_RECAL ?= false
BAM_NO_REALN ?= false

INCLUDE_CHR_Y ?= false

### VARIANT FILTERS
MIN_NORMAL_DEPTH ?= 10
MIN_TUMOR_DEPTH ?= 20
MIN_TUMOR_AD ?= 5
MIN_TN_AD_RATIO ?= 5.0

HRUN ?= false
HRUN_FILTER ?= 4

ANN_FACETS ?= true ## useful mostly for exomes

include usb-modules/genome_inc/$(REF).inc

##################################
####### Environmental var ########
##################################

export TMPDIR := /scicore/scratch
PATH "= $(HOME)/usr_nobackup:$(PATH)
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

JAVA_ARGS = -Djava.io.tmpdir=$(TMPDIR)
JAVA6_MODULE ?= Java/1.6.0_24
JAVA7_MODULE ?= Java/1.7.0_80
JAVA8_MODULE ?= Java/1.8.0_92
LOAD_JAVA6_MODULE ?= module load $(JAVA6_MODULE)
LOAD_JAVA7_MODULE ?= module load $(JAVA7_MODULE)
LOAD_JAVA8_MODULE ?= module load $(JAVA8_MODULE) 
JAVA ?= java $(JAVA_ARGS)
JARDIR ?= /scicore/home/terracci/ngkiuy/usr_nobackup/local/lib/java

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

BAM_SUFFIX := $(subst $( ),.,$(strip \
        $(if $(findstring false,$(BAM_NO_SORT)),sorted)\
        $(if $(findstring false,$(BAM_NO_FILTER)),filtered)\
        $(if $(findstring true,$(BAM_FIX_RG)),rg)\
        $(if $(findstring false,$(BAM_NO_REALN)),realn)\
        $(if $(findstring rmdup,$(BAM_DUP_TYPE)),rmdup)\
        $(if $(findstring markdup,$(BAM_DUP_TYPE)),markdup)\
        $(if $(findstring false,$(BAM_NO_RECAL)),recal)\
        bam))

MERGE_SPLIT_FASTQ ?= true

##################################
#### General MPS tools ###########
##################################

SAMTOOLS_MODULE ?= SAMtools/1.2-goolf-1.4.10
LOAD_SAMTOOLS_MODULE ?= module load $(SAMTOOLS_MODULE)
SAMTOOLS ?= samtools

SAMTOOLS_SORT_MEM = 2000000000

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

VCF_SORT ?= $(PERL) $(SCRIPTSDIR)/vcfsorter.pl 

VCF2MAF ?= $(PERL) $(SCRIPTSDIR)/vcf2maf.pl
MAF2MAF ?= $(PERL) $(SCRIPTSDIR)/maf2maf.pl
MAF2VCF ?= $(PERL) $(SCRIPTSDIR)/maf2vcf.pl
VEP_PATH ?= /opt/common/CentOS_6/vep/v82
VEP_DATA ?= /opt/common/CentOS_6/vep/v82

BAM2FASTX = $(HOME)/usr/bin/bam2fastx
EXTRACT_TOOL ?= picard

FASTQ_TRIMMER = $(PERL) usb-modules/fastq_tools/trimFastq.pl

MODULE_SCRIPTS_DIR ?= usb-modules/scripts
MERGE ?= $(MODULE_SCRIPTS_DIR)/merge.R
VCF_TO_TABLE ?= $(MODULE_SCRIPTS_DIR)/vcfToTable.R
INTRON_POSN_LOOKUP ?= $(MODULE_SCRIPTS_DIR)/posnGeneLookup.pl
RBIND ?= $(MODULE_SCRIPTS_DIR)/rbind.R

GET_INSERT_SIZE ?= $(HOME)/share/usr/bin/getInsertSize.py

VT ?= $(HOME)/share/usr/bin/vt
TABIX ?= $(HOME)/share/usr/bin/tabix


##################################
########## QC/metrics ############
##################################

FASTQC_MODULE ?= FastQC/0.11.4-Java-1.7.0_21
LOAD_FASTQC_MODULE ?= module load $(FASTQC_MODULE)
FASTQC ?= $(PERL) fastqc ####### BROKEN

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

TMAP = $(HOME)/share/usr/bin/tmap
TMAP_MODE ?= map3
TMAP_OPTS =

TOPHAT ?= $(HOME)/share/usr/bin/tophat2-2.0.13
TOPHAT_NUM_CORES ?= 4

NO_NOVEL_SPLICING ?= false

TOPHAT_OPTS = --keep-fasta-order --no-sort-bam -G $(GENES_GTF) -p ${TOPHAT_NUM_CORES} --tmp-dir $(TMPDIR)/$*
ifeq ($(BAM_PHRED64),true)
        TOPHAT_OPTS += --solexa1.3-quals
endif
ifeq ($(NO_NOVEL_SPLICING),true)
        TOPHAT_OPTS += --no-novel-juncs
endif

HISAT ?= $(HOME)/share/usr/hisat2-2.0.2-beta/hisat2
HISAT_NUM_CORES ?= 8
HISAT_OPTS = -q -p $(HISAT_NUM_CORES) -x $(HISAT_REF)
ifeq ($(SEQ_PLATFORM),true)
        HISAT_OPTS += --fr
endif
ifeq ($(PHRED64),true)
        HISAT_OPTS += --phred64
endif



##################################
########## Picard ################
##################################

PICARD_OPTS = VALIDATION_STRINGENCY=LENIENT MAX_RECORDS_IN_RAM=4000000 TMP_DIR=$(TMPDIR) REFERENCE_SEQUENCE=$(REF_FASTA) CREATE_INDEX=true
PICARD_JAR ?= $(JARDIR)/picard-tools-2.4.1/picard.jar
PICARD = $(JAVA) -Xmx$($2) -jar $(PICARD_JAR) $(1) $(PICARD_OPTS)

SORT_SAM = $(call PICARD,SortSam,19G)
REORDER_SAM = $(call PICARD,ReorderSam,20G)
MARK_DUP = $(call PICARD,MarkDuplicates,12G)
MARK_DUP_CIGAR = $(call PICARD,MarkDuplicatesWithMateCigar,20G)
MERGE_SAMS = $(call PICARD,MergeSamFiles,10G)
FIX_MATE = $(call PICARD,FixMateInformation,8G)
SAM_TO_FASTQ = $(call PICARD,SamToFastq,9G)
CLEANBAM = $(call PICARD,CleanSam,5G)
ADD_RG = $(call PICARD,AddOrReplaceReadGroups,12G)
VALIDATE_SAM = $(call PICARD,ValidateSamFile,10G)

INTERVAL_LIST_TOOL = $(call PICARD,IntervalListTools,6G)
CREATE_SEQ_DICT = $(call PICARD,CreateSequenceDictionary,4G)

COLLECT_HS_METRICS = $(call PICARD,CollectHsMetrics,8G)
COLLECT_MULT_METRICS = $(call PICARD,CollectMultipleMetrics,8G)
COLLECT_TARGETED_METRICS = $(call PICARD,CollectTargetedPcrMetrics,8G)
COLLECT_ALIGNMENT_METRICS = $(call PICARD,CollectAlignmentSummaryMetrics,8G)
COLLECT_ILMN_LANE_METRICS = $(call PICARD,CollectIlluminaLaneMetrics,8G)
COLLECT_OXOG_METRICS = $(call PICARD,CollectOxoGMetrics,8G)
COLLECT_RNA_METRICS = $(call PICARD,CollectRnaSeqMetrics,8G)
COLLECT_SEQ_ARTIFACT_METRICS = $(call PICARD,CollectSequencingArtifactMetrics,8G)
COLLECT_WGS_METRICS = $(call PICARD,CollectWgsMetrics,10G)


##################################
########## GATK ### ##############
##################################

GATK_JAR ?= $(JARDIR)/GenomeAnalysisTK-3.6/GenomeAnalysisTK-3.6.jar
GATK ?= $(JAVA) -Xmx$(2) -jar $(GATK_JAR) -S LENIENT -T $(1)

REALIGN_TARGET_CREATOR ?= $(call GATK,RealignerTargetCreator,11G)
INDEL_REALIGN ?= $(call GATK,IndelRealigner,8G)

ANALYZE_COVARIATES ?= $(call GATK,AnalyzeCovariates,18G)
BASE_RECALIBRATOR ?= $(call GATK,BaseRecalibrator,10G)

LEFT_ALIGN_VCF ?= $(call GATK,LeftAlignAndTrimVariants,8G)
VARIANT_FILTRATION ?= $(call GATK,VariantFiltration,6G)
VARIANT_ANNOTATOR ?= $(call GATK,VariantAnnotator,8G)
COMBINE_VARIANTS ?= $(call GATK,CombineVariants,8G)

CONTEST ?= $(call GATK,ContEst,6G)
ASSIGN_SOMATIC_STATUS ?= $(call GATK,AssignSomaticStatus,8G)
HAPLOTYPE_CALLER ?= $(call GATK,HaplotypeCaller,18G)
MUTECT2 ?= $(call GATK,MuTect2,18G)
PRINT_READS ?= $(call GATK,PrintReads,10G)
APPLY_RECALIBRATION ?= $(call GATK,ApplyRecalibration,8G)
VARIANT_RECALIBRATOR ?= $(call GATK,VariantRecalibrator,18G)

BAM_CHR1_BASE_RECAL ?= false
BAM_BASE_RECAL_OPTS = -knownSites $(DBSNP) $(if $(findstring true,$(BAM_CHR1_BASE_RECAL)),-L $(word 1,$(CHROMOSOMES)))

ifneq ($(KNOWN_INDELS),)
BAM_REALN_OPTS = --knownAlleles $(KNOWN_INDELS)
BAM_REALN_TARGET_OPTS = --known $(KNOWN_INDELS)
endif

CONTEST_OPTS ?= -T Contamination -BTI genotypes
CONTEST += -R $(REF_FASTA) -B:pop$(,)vcf $(CONTEST_REF_VCF) $(CONTEST_OPTS)

COVARIATES = -cov ReadGroupCovariate -cov QualityScoreCovariate -cov DinucCovariate -cov CycleCovariate

VARIANT_CALL_THRESHOLD = 50
VARIANT_EMIT_THRESHOLD = 10
VARIANT_RECAL_TRUTH_SENSITIVITY_LEVEL = 99.0
VARIANT_RECAL_ANNOTATIONS = QD MQRankSum FS MQ ReadPosRankSum
HAPLOTYPE_CALLER_OPTS = --dbsnp $(DBSNP) -rf BadCigar -stand_call_conf $(VARIANT_CALL_THRESHOLD) -stand_emit_conf $(VARIANT_EMIT_THRESHOLD) $(REF_FASTA)

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

REPORT_STRATIFICATION := Filter

##################################
########## Variant callers #######
##################################

VALIDATION ?= false

MUTECT_JAR ?= $(JARDIR)/muTect-1.1.4.jar
MUTECT_MEM ?= 11g
MUTECT_MAX_ALT_IN_NORMAL ?= 500
MUTECT_MAX_ALT_IN_NORMAL_FRACTION ?= 0.05
MUTECT_OPTS = --max_alt_alleles_in_normal_count $(MUTECT_MAX_ALT_IN_NORMAL) --max_alt_allele_in_normal_fraction $(MUTECT_MAX_ALT_IN_NORMAL_FRACTION)
#MUT_FREQ_REPORT = modules/variant_callers/somatic/mutectReport.Rmd
MUTECT = $(JAVA) -Xmx$(MUTECT_MEM) -jar $(MUTECT_JAR) --analysis_type MuTect

VARSCAN_JAR ?= $(JARDIR)/varscan/VarScan.v2.4.2.jar
VARSCAN = $(JAVA8) -Xmx6G -jar $(VARSCAN_JAR)

CONFIGURE_STRELKA = $(PERL) /scicore/home/terracci/ngkiuy/usr_nobackup/local/strelka_workflow-1.0.15/bin/configureStrelkaWorkflow.pl
STRELKA_CONFIG = /scicore/home/terracci/ngkiuy/usr_nobackup/local/strelka_workflow-1.0.15/etc/strelka_config_bwa_default.ini

#SOMATIC_FILTER_VCF ?= $(PERL) $(SCRIPTS_DIR)/somaticFilterVCF.pl
#NORMAL_FILTER ?= $(PERL) $(SCRIPTS_DIR)/normalFilterVCF.pl


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
SNP_EFF_FLAGS ?= -canon # -ud 0  -no-intron -no-intergenic -no-utr
SNP_EFF_OPTS = -c $(SNP_EFF_CONFIG) -i vcf -o vcf $(SNP_EFF_FLAGS)
SNP_SIFT_OPTS = -c $(SNP_EFF_CONFIG)

VCF_FIELDS = CHROM POS ID REF ALT FILTER
#NSFP_FIELDS ?= ExAC_Adj_AF 1000Gp3_AF ESP6500_AA_AF ESP6500_EA_AF \
#	clinvar_rs clinvar_clnsig \
#	Uniprot_acc_Polyphen2 Uniprot_id_Polyphen2 Polyphen2_HVAR_pred \
#	MutationAssessor_pred MutationTaster_pred \
#	PROVEAN_pred Interpro_domain
ANN_FIELDS = $(addprefix ANN[*].,ALLELE EFFECT IMPACT GENE GENEID FEATURE FEATUREID BIOTYPE RANK \
	HGVS_C HGVS_P CDNA_POS CDNA_LEN CDS_POS CDS_LEN AA_POS AA_LEN DISTANCE ERRORS)

VCF_EFF_ONE_PER_LINE = $(PERL) $(SNP_EFF_DIR)/scripts/vcfEffOnePerLine.pl
VCF_JOIN_EFF = $(PERL) usb-modules/scripts/joinEff.pl

ANNOVAR = $(PERL) $(HOME)/share/usr/annovar/table_annovar.pl
ANNOVAR_PROTOCOL ?= refGene$(,)cytoBand$(,)genomicSuperDups$(,)esp6500siv2_all$(,)1000g2014oct_all$(,)1000g2014oct_afr$(,)1000g2014oct_eas$(,)1000g2014oct_eur$(,)snp138$(,)ljb26_all
ANNOVAR_OPERATION ?= g$(,)r$(,)r$(,)f$(,)f$(,)f$(,)f$(,)f$(,)f$(,)f
ANNOVAR_OPTS = --dot2underline -remove -protocol $(ANNOVAR_PROTOCOL) -operation $(ANNOVAR_OPERATION) -nastring . -vcfinput -buildver $(ANNOVAR_REF)


CHASM = $(RSCRIPT) modules/vcf_tools/chasmVcf.R
#CHASM_DIR = /ifs/opt/common/CHASM/CHASMDL.1.0.7
CHASM_DIR = $(HOME)/share/usr/CHASM
CHASM_PYTHON_ENV = $(HOME)/share/usr/anaconda-envs/pyenv27-chasm
CHASM_CLASSIFIER ?= Breast

FATHMM = $(MY_RSCRIPT) modules/vcf_tools/fathmmVcf.R
FATHMM_DIR = $(HOME)/share/usr/fathmm
FATHMM_PYTHON = $(HOME)/share/usr/bin/python
FATHMM_PYTHONPATH = $(HOME)/share/usr/lib/python:$(HOME)/share/usr/lib/python2.7
FATHMM_OPTS = --genome $(REF) --ensemblTxdb $(ENSEMBL_TXDB) --ref $(REF_FASTA) --fathmmDir $(FATHMM_DIR) --python $(FATHMM_PYTHON) \
	--mysqlHost $(EMBL_MYSQLDB_HOST) --mysqlPort $(EMBL_MYSQLDB_PORT) --mysqlUser $(EMBL_MYSQLDB_USER) \
	--mysqlPassword $(EMBL_MYSQLDB_PW) --mysqlDb (EMBL_MYSQLDB_DB)
TRANSFIC = $(RSCRIPT) modules/vcf_tools/transficVcf.R
TRANSFIC_PERL_SCRIPT = $(HOME)/share/usr/transfic/bin/transf_scores.pl

PROVEAN = $(RSCRIPT) modules/vcf_tools/proveanVcf.R
PROVEAN_OPTS = --genome $(REF) --aaTable $(AA_TABLE) --ensemblTxdb $(ENSEMBL_TXDB) --mysqlHost $(EMBL_MYSQLDB_HOST) \
			   --mysqlPort $(EMBL_MYSQLDB_PORT) --mysqlUser $(EMBL_MYSQLDB_USER) --mysqlPassword $(EMBL_MYSQLDB_PW) \
			   --mysqlDb $(EMBL_MYSQLDB_DB) --numThreads 8 --memPerThread 1G --queue $(QUEUE) --qsubPriority $(QSUB_PRIORITY)

AA_TABLE = $(HOME)/share/reference/aa_table.tsv

MUT_ASS = $(RSCRIPT) usb-modules/vcf_tools/mutAssVcf.R
MUTATION_TASTER = $(PYTHON) usb-modules/vcf_tools/mutation_taster_vcf.py

ADD_GENE_LIST_ANNOTATION = $(RSCRIPT) usb-modules/vcf_tools/addGeneListAnnotationToVcf.R

NON_SILENT_EFF = START_GAINED SPLICE_SITE_ACCEPTOR SPLICE_SITE_DONOR START_LOST NON_SYNONYMOUS_CODING FRAME_SHIFT CODON_CHANGE CODON_INSERTION CODON_CHANGE_PLUS_CODON_INSERTION CODON_DELETION CODON_CHANGE_PLUS_CODON_DELETION STOP_GAINED STOP_LOST NON_SYNONYMOUS_START
NON_SILENT_CODING_EFF = START_GAINED START_LOST NON_SYNONYMOUS_CODING FRAME_SHIFT CODON_CHANGE CODON_INSERTION CODON_CHANGE_PLUS_CODON_INSERTION CODON_DELETION CODON_CHANGE_PLUS_CODON_DELETION STOP_GAINED STOP_LOST NON_SYNONYMOUS_START

EMBL_MYSQLDB_HOST ?= 10.0.200.71
EMBL_MYSQLDB_PORT ?= 38493
EMBL_MYSQLDB_USER ?= embl
EMBL_MYSQLDB_PW ?= embl
EMBL_MYSQLDB_DB ?= homo_sapiens_core_75_37

ANNOTATE_FACETS_VCF = $(RSCRIPT) usb-modules/copy_number/annotateFacets2Vcf.R

#NORMAL_VCF ?= $(HOME)/share/reference/spowellnormal.gatk_variants.vcf

CANCER_HOTSPOT_ANNOTATION_SCRIPT = python modules/vcf_tools/hotspot_vcf.py

CLASSIFY_PATHOGENICITY = $(PYTHON) modules/vcf_tools/classify_pathogenicity_vcf.py
#COMMON_FILTER_VCF = $(PYTHON) modules/vcf_tools/common_filter_vcf.py
##################################
########## Copy number ###########
##################################

GET_BASE_COUNTS = $(PERL) usb-modules/copy_number/facets-snp-pileup.pl
GET_BASE_COUNTS_PARAMS = $(REF_FASTA) $(TARGETS_FILE_INTERVALS) 25 1000

FACETS_PRE_CVAL ?= 50
FACETS_CVAL1 ?= 150
FACETS_CVAL2 ?= 50
FACETS_MIN_NHET ?= 25
FACETS_GATK_VARIANTS ?= false
FACETS_OPTS = --cval2 $(FACETS_CVAL2) --cval1 $(FACETS_CVAL1) --genome $(REF) --min_nhet $(FACETS_MIN_NHET) --pre_cval $(FACETS_PRE_CVAL)
FACETS = $(RSCRIPT) usb-modules/copy_number/runFacets.R

FACETS_UNION_GATK_DBSNP ?= false
ifeq ($(FACETS_UNION_GATK_DBSNP),true)
FACETS_SNP_VCF = facets/vcf/dbsnp_het_gatk.snps.vcf.gz
else
FACETS_SNP_VCF = $(FACETS_DBSNP)
endif

FACETS_GENE_CN = $(RSCRIPT) usb-modules/copy_number/facetsGeneCN.R
FACETS_FILL_GENE_CN = $(RSCRIPT) usb-modules/copy_number/facetsFillGeneCN.R
FACETS_GENE_CN_OPTS = $(if $(GENES_FILE),--genesFile $(GENES_FILE)) \
                                          --mysqlHost $(EMBL_MYSQLDB_HOST) --mysqlPort $(EMBL_MYSQLDB_PORT) \
                                          --mysqlUser $(EMBL_MYSQLDB_USER) --mysqlPassword $(EMBL_MYSQLDB_PW) \
                                          --mysqlDb $(EMBL_MYSQLDB_DB)

FACETS_PLOT_GENE_CN = $(RSCRIPT) usb-modules/copy_number/facetsGeneCNPlot.R
FACETS_PLOT_GENE_CN_OPTS = --sampleColumnPostFix '_EM'

ANN_FACETS ?= false

##################################
########## Gene expression #######
##################################

CUFFLINKS ?= $(HOME)/share/usr/bin/cufflinks-2.2.1
CUFFCMP ?= $(HOME)/share/usr/bin/cuffcompare-2.2.1

##################################
########## Fusions ###############
##################################

DEFUSE ?= $(PERL) $(HOME)/share/usr/defuse-0.6.1/scripts/defuse.pl
INTEGRATE ?= $(HOME)/share/usr/Integrate_0.2.0/bin/Integrate
ONCOFUSE_JAR ?= $(JARDIR)/Oncofuse.jar


endif
CONFIG_INC = true
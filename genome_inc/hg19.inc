ifndef HG19_INC
REF = hg19
REF_DIR = $(HOME)/ref_nobackup

ifeq ($(SEQ_PLATFORM),IONTORRENT)
REF_FASTA ?= $(REF_DIR)/ion_reference/hg19.fasta
ifeq ($(PANEL),CHPv2)
TARGETS_FILE_INTERVALS = $(REF_DIR)/target_beds/AmpliSeq_CancerHotspotPanelv2.106003.designed.bed
TARGETS_FILE_MUTATIONS = $(REF_DIR)/target_beds/AmpliSeq_CancerHotspotPanelv2.106003.mutations.bed
CAPTURE_METHOD = PCR
endif
ifeq ($(PANEL),CCP)
TARGETS_FILE_INTERVALS = $(REF_DIR)/target_beds/AmpliSeq.CCP.20131001.designed.bed
TARGETS_FILE_MUTATIONS = $(REF_DIR)/target_beds/AmpliSeq.CCP.20131001.hotspots.bed
CAPTURE_METHOD = PCR
endif
ifeq ($(PANEL),EXOME)
TARGETS_FILE_INTERVALS = $(REF_DIR)/target_beds/AmpliSeqExome.20131001.designed.bed
TARGETS_FILE_MUTATIONS = $(REF_DIR)/target_beds/AmpliSeqExome.20131001.designed.bed
CAPTURE_METHOD = PCR
endif
else
GSA_REF_DIR = $(REF_DIR)/gsa_resource_bundle/2.8/hg19
REF_FASTA ?= $(GSA_REF_DIR)/human_g1k_v37.fasta.gz
endif

#TOPHAT_REF = GRCh37
#GENES_GTF = $(REF_DIR)/Homo_sapiens/Ensembl/$(TOPHAT_REF)/Annotation/Genes/genes.gtf

#ANNOVAR_REF = hg19

CLINVAR = $(REF_DIR)/clinvar_20151102.vcf.gz
EXAC_NONTCGA = $(REF_DIR)/annotation_data/ExAC_nonTCGA.r0.3.1.sites.vep.vcf.gz

#UCSC_REFSEQ_GTF = $(REF_DIR)s/genomes/ucsc.refseq.gtf
#REFSEQ_GTF := $(UCSC_REFSEQ_GTF)

#CHR_LEN = $(REF_DIR)/GATK_bundle/2.3/b37.len

SNP_TXT = $(REF_DIR)/b37_snp137.SingleDiNucl.1based.nochr.txt

CHIMERASCAN_REF = $(REF_DIR)/chimerascan_index

RIBOSOMAL_INTERVALS = $(REF_DIR)/rrna_intervals.txt
GENE_REF_FLAT = $(REF_DIR)/b37_genes.refFlat.txt

CENTROMERE_TABLE = $(REF_DIR)/annotation_data/centromeres_b37_hg19.txt

ENSEMBL_TXDB = $(REF_DIR)/hsapiens_ensembl_biomart.sqlite

MAPSPLICE_REF_DIR = $(REF_DIR)/GATK_bundle/2.3/mapsplice_ref/chr
MAPSPLICE_REF_BASENAME = $(REF_DIR)/GATK_bundle/2.3/mapsplice_ref/hg19

GISTIC_REF ?= $(HOME)/share/usr/gistic_2_0_21/refgenefiles/hg19.mat

#REF_2BIT ?= $(REF_DIR)/GATK_bundle/2.3/human_g1k_v37.2bit
DBSNP ?= $(GSA_REF_DIR)/dbsnp_138.b37.excluding_sites_after_129.vcf.gz
KNOWN_INDELS ?= $(GSA_REF_DIR)/Mills_and_1000G_gold_standard.indels.b37.vcf.gz
OMNI ?= $(GSA_REF_DIR)/1000G_omni2.5.b37.vcf.gz
#HAPMAP ?= $(REF_DIR)/hapmap_3.3.b37.vcf.gz

COSMIC ?= $(REF_DIR)/CosmicCodingMuts_v69.vcf
BOWTIE_REF = $(REF_DIR)/GATK_bundle/2.3/human_g1k_v37
ILMN_EXOME_BED ?= $(REF_DIR)/SureSelect_50MB_S02972011_Regions_nochr_noAnn.bed
ILMN_EXOME_BED_NOHEADER ?= $(REF_DIR)/SureSelect_50MB_S02972011_Regions_nochr_noAnn_noheader.bed
MAP_BIGWIG ?= $(REF_DIR)/wgEncodeCrgMapabilityAlign100mer.nochr.bigWig

EXOME ?= false
ifeq ($(EXOME),true)
TARGETS_FILE = $(EXOME_BED_NOHEADER)
QUALIMAP_TARGETS_FILE = $(REF_DIR)/SureSelect_50MB_S02972011_Regions_nochr.bed
EMU_TARGETS_FILE = $(REF_DIR)/SureSelect_50MB_S02972011_Regions_nochr_noAnn_noheader_numchr_spacedelim.bed
endif

INCLUDE_CHR_Y ?= false
ifneq ($(and $(TARGETS_FILE),$(findstring false,$(EXOME))),)
CHROMOSOMES := $(shell grep -v '@' $(TARGETS_FILE) | cut -f1 | sort | uniq)
else
CHROMOSOMES ?= $(shell seq 1 22) X $(if $(findstring true,$(INCLUDE_CHR_Y)),Y) MT
endif

KNOWN_CNVS = $(REF_DIR)/GRCh37_hg19_variants_2013-07-23.txt


REF_NAME := 1000GenomesPilot-NCBI37
REF_DATE := 20100517
REF_DICT = $(REF_FASTA:.fasta=.dict)
SNP_EFF_GENOME = hg19

EMU_REF_DIR ?= $(REF_DIR)/GATK_bundle/2.3/emu_ref
EMU_REF_MUTATIONS ?= $(REF_DIR)/emu_ref/exome_ref/Breast_clean_somatic_mutations_for_signature_analysis.emu.mutations.txt
EMU_REF_CNV ?= $(REF_DIR)/emu_ref/exome_ref/Breast_clean_somatic_mutations_for_signature_analysis.emu.copynum.txt

HMMCOPY_GC_WIG = $(REF_DIR)/GATK_bundle/2.3/gc_hg19.wig
HMMCOPY_MAP_WIG = $(REF_DIR)/map_hg19.wig

REF_HOUSEKEEPING_GENE_BED = $(HOME)/share/reference/hg19.HouseKeepingGenes.bed

INTEGRATE_ANN = $(REF_DIR)/annot.ucsc.txt
INTEGRATE_BWTS = $(REF_DIR)/GATK_bundle/2.3/bwts
ONCOTATOR_DB ?= $(REF_DIR)/oncotator_v1_ds_Jan262014

ANNOVAR_DB = $(HOME)/share/usr/annovar/humandb

HISAT_REF = $(HOME)/share/reference/hisat/grch37_snp_tran/genome_snp_tran

endif
HG19_INC = true


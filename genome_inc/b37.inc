ifndef B37_INC
REF = b37
REF_DIR = $(HOME)/ref_nobackup
GSA_REF_DIR = $(REF_DIR)/gsa_resource_bundle/2.8/b37

ifeq ($(SEQ_PLATFORM),IONTORRENT)
ifeq ($(PANEL),CHPv2)
TARGETS_FILE_INTERVALS = $(REF_DIR)/target_beds/AmpliSeq_CancerHotspotPanelv2.106003.designed.bed
TARGETS_FILE_MUTATIONS = $(REF_DIR)/target_beds/AmpliSeq_CancerHotspotPanelv2.106003.mutations.bed
endif
ifeq ($(PANEL),CCP)
TARGETS_FILE_INTERVALS = $(REF_DIR)/target_beds/AmpliSeq.CCP.20131001.designed.nochr.bed
TARGETS_FILE_MUTATIONS = $(REF_DIR)/target_beds/AmpliSeq.CCP.20131001.hotspots.bed
endif
ifeq ($(PANEL),EXOME)
TARGETS_FILE_INTERVALS = $(REF_DIR)/target_beds/AmpliSeqExome.20131001.designed.bed
TARGETS_FILE_MUTATIONS = $(REF_DIR)/target_beds/AmpliSeqExome.20131001.designed.bed
endif
endif

ifeq ($(SEQ_PLATFORM),ILLUMINA)
ifeq ($(PANEL),EXOME)
TARGETS_FILE_INTERVALS = $(REF_DIR)/target_beds/
TARGETS_FILE_MUTATIONS = $(REF_DIR)/target_beds/
endif
ifeq ($(PANEL),IMPACT410)
TARGETS_FILE_INTERVALS = $(REF_DIR)/target_beds/impact410_b37.bed
TARGETS_FILE_MUTATIONS = $(REF_DIR)/target_beds/impact410_b37.bed
endif
ifeq ($(PANEL),IMPACT341)
TARGETS_FILE_INTERVALS = $(REF_DIR)/target_beds/impact341_b37.bed
TARGETS_FILE_MUTATIONS = $(REF_DIR)/target_beds/impact341_b37.bed
endif
endif

REF_FASTA ?= $(GSA_REF_DIR)/human_g1k_v37.fasta

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
DBSNP ?= $(GSA_REF_DIR)/dbsnp_138.b37.excluding_sites_after_129.vcf
KNOWN_INDELS ?= $(GSA_REF_DIR)/Mills_and_1000G_gold_standard.indels.b37.vcf
OMNI ?= $(GSA_REF_DIR)/1000G_omni2.5.b37.vcf
#HAPMAP ?= $(REF_DIR)/hapmap_3.3.b37.vcf.gz

COSMIC ?= $(REF_DIR)/CosmicCodingMuts_v69.vcf
BOWTIE_REF = $(REF_DIR)/GATK_bundle/2.3/human_g1k_v37
ILMN_EXOME_BED ?= $(REF_DIR)/SureSelect_50MB_S02972011_Regions_nochr_noAnn.bed
ILMN_EXOME_BED_NOHEADER ?= $(REF_DIR)/SureSelect_50MB_S02972011_Regions_nochr_noAnn_noheader.bed
MAP_BIGWIG ?= $(REF_DIR)/wgEncodeCrgMapabilityAlign100mer.nochr.bigWig

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
SNP_EFF_GENOME = GRCh37.75

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
B37_INC = true


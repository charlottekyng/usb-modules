ifndef HG19_INC
REF = hg19
REF_DIR = $(HOME)/ref_nobackup

ifeq ($(findstring IONTORRENT,$(SEQ_PLATFORM)),IONTORRENT)
ifeq ($(PANEL),CHPv2)
TARGETS_FILE_INTERVALS = $(REF_DIR)/target_beds/hg19/AmpliSeq_CancerHotspotPanelv2.106003.designed.bed
TARGETS_FILE_MUTATIONS = $(REF_DIR)/target_beds/hg19/AmpliSeq_CancerHotspotPanelv2.106003.mutations.bed
endif
ifeq ($(PANEL),CCP)
TARGETS_FILE_INTERVALS = $(REF_DIR)/target_beds/hg19/AmpliSeq.CCP.20131001.designed.bed
TARGETS_FILE_MUTATIONS = $(REF_DIR)/target_beds/hg19/AmpliSeq.CCP.20131001.hotspots.bed
endif
ifeq ($(PANEL),EXOME)
TARGETS_FILE_INTERVALS = $(REF_DIR)/target_beds/hg19/AmpliSeqExome.20131001.designed.bed
TARGETS_FILE_MUTATIONS = $(REF_DIR)/target_beds/hg19/AmpliSeqExome.20131001.designed.bed
endif
endif

ifeq ($(findstring ILLUMINA,$(SEQ_PLATFORM)),ILLUMINA)
ifeq ($(PANEL),AGILENT_CLINICAL_EXOME)
TARGETS_FILE_INTERVALS = $(REF_DIR)/target_beds/hg19/S06588914_Agilent_Clinical_Exome_Padded_noAnn.bed
TARGETS_FILE_MUTATIONS = $(REF_DIR)/target_beds/hg19/S06588914_Agilent_Clinical_Exome_Padded_noAnn.bed
endif
#ifeq ($(findstring IMPACT410,$(PANEL)),IMPACT410)
#TARGETS_FILE_INTERVALS = $(REF_DIR)/target_beds/hg19/impact410_b37.bed
#TARGETS_FILE_MUTATIONS = $(REF_DIR)/target_beds/hg19/impact410_b37.bed
#endif
#ifeq ($(PANEL),IMPACT341)
#TARGETS_FILE_INTERVALS = $(REF_DIR)/target_beds/hg19/impact341_b37.bed
#TARGETS_FILE_MUTATIONS = $(REF_DIR)/target_beds/hg19/impact341_b37.bed
#endif
endif

ifndef $(TARGETS_FILE_INTERVALS)
CHROMOSOMES ?= chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22 chrX $(if $(findstring true,$(INCLUDE_CHR_Y)),chrY) chrM
else
CHROMOSOMES := $(shell grep -v '@' $(TARGETS_FILE_INTERVALS) | cut -f1 | sort | uniq)
endif

### GENERAL REFERENCES ###

GSA_REF_DIR = $(REF_DIR)/gsa_resource_bundle/2.8/hg19

ifeq ($(findstring ILLUMINA,$(SEQ_PLATFORM)),ILLUMINA)
	REF_FASTA ?= $(GSA_REF_DIR)/ucsc.hg19.fasta
endif
ifeq ($(findstring IONTORRENT,$(SEQ_PLATFORM)),IONTORRENT)
	REF_FASTA ?= $(REF_DIR)/ion_reference/hg19.fasta
endif

### ANNOTATION SETS ####

ANN_DIR ?= $(REF_DIR)/annotation_data/hg19

CONTEST_REF_VCF ?= $(ANN_DIR)/hg19_population_stratified_af_hapmap_3.3.vcf

DBSNP ?= $(GSA_REF_DIR)/dbsnp_138.hg19.vcf
KNOWN_INDELS ?= $(GSA_REF_DIR)/Mills_and_1000G_gold_standard.indels.hg19.sites.vcf
OMNI ?= $(GSA_REF_DIR)/1000G_omni2.5.hg19.sites.vcf
HAPMAP ?= $(GSA_REF_DIR)/hapmap_3.3.hg19.sites.vcf

CLINVAR = $(ANN_DIR)/clinvar.vcf
EXAC_NONTCGA = $(ANN_DIR)/ExAC_nonTCGA.r0.3.1.sites.vep.vcf
COSMIC ?= $(ANN_DIR)/CosmicCodingMuts.sorted.vcf

HAPLOTYPE_INSUF_BED = $(ANN_DIR)/haploinsufficient.genes.bed
CANCER_GENE_CENSUS_BED = $(ANN_DIR)/cancer_gene_census.bed
KANDOTH_BED = $(ANN_DIR)/kandoth127smg.bed
LAWRENCE_BED = $(ANN_DIR)/lawrence_cancer5000s.bed
DGD_BED = $(ANN_DIR)/dgd_Hsa_all_v71.bed
CANCER_HOTSPOT_VCF = $(ANN_DIR)/hotspot.vcf

#FALSE_POSITIVE_BED = $(HOME)/share/reference/fuentes_blacklist.include_cosmic.hg19.bed
#ENCODE_BED = $(HOME)/share/reference/wgEncodeDacMapabilityConsensusExcludable.include_cosmic.bed

CENTROMERE_TABLE = $(ANN_DIR)/centromeres_hg19.txt

KNOWN_CNVS = $(ANN_DIR)/GRCh37_hg19_variants_2013-07-23.txt

REF_HOUSEKEEPING_GENE_BED = $(ANN_DIR)/share/reference/hg19.HouseKeepingGenes.bed

ENSEMBL_TXDB = $(ANN_DIR)/hsapiens_ensembl_biomart.sqlite

RIBOSOMAL_INTERVALS = $(ANN_DIR)/rrna_intervals.txt
GENE_REF_FLAT = $(ANN_DIR)/b37_genes.refFlat.txt

SNP_TXT = $(ANN_DIR)/b37_snp137.SingleDiNucl.1based.nochr.txt

### TOOL-SPECIFIC REF ###

GISTIC_REF ?= $(HOME)/share/usr/gistic_2_0_21/refgenefiles/hg19.mat
CHIMERASCAN_REF ?= $(REF_DIR)/chimerascan_index
BOWTIE_REF = $(REF_DIR)/GATK_bundle/2.3/human_g1k_v37

MAPSPLICE_REF_DIR = $(REF_DIR)/GATK_bundle/2.3/mapsplice_ref/chr
MAPSPLICE_REF_BASENAME = $(REF_DIR)/GATK_bundle/2.3/mapsplice_ref/hg19

DB_NSFP ?= $(REF_DIR)/snpeff/GRCh37/dbNSFP-v2.9.txt.gz

#TOPHAT_REF = GRCh37
#GENES_GTF = $(REF_DIR)/Homo_sapiens/Ensembl/$(TOPHAT_REF)/Annotation/Genes/genes.gtf

#ANNOVAR_REF = hg19

#UCSC_REFSEQ_GTF = $(REF_DIR)s/genomes/ucsc.refseq.gtf
#REFSEQ_GTF := $(UCSC_REFSEQ_GTF)

#CHR_LEN = $(REF_DIR)/GATK_bundle/2.3/b37.len

#REF_2BIT ?= $(REF_DIR)/GATK_bundle/2.3/human_g1k_v37.2bit

MAP_BIGWIG ?= $(REF_DIR)/wgEncodeCrgMapabilityAlign100mer.nochr.bigWig

REF_NAME := 1000GenomesPilot-NCBI37
REF_DATE := 20100517
REF_DICT = $(REF_FASTA:.fasta=.dict)
SNP_EFF_GENOME = GRCh37.75

#EMU_REF_DIR ?= $(REF_DIR)/GATK_bundle/2.3/emu_ref
#EMU_REF_MUTATIONS ?= $(REF_DIR)/emu_ref/exome_ref/Breast_clean_somatic_mutations_for_signature_analysis.emu.mutations.txt
#EMU_REF_CNV ?= $(REF_DIR)/emu_ref/exome_ref/Breast_clean_somatic_mutations_for_signature_analysis.emu.copynum.txt

#HMMCOPY_GC_WIG = $(REF_DIR)/GATK_bundle/2.3/gc_hg19.wig
#HMMCOPY_MAP_WIG = $(REF_DIR)/map_hg19.wig


INTEGRATE_ANN = $(REF_DIR)/annot.ucsc.txt
INTEGRATE_BWTS = $(REF_DIR)/GATK_bundle/2.3/bwts
ONCOTATOR_DB ?= $(REF_DIR)/oncotator_v1_ds_Jan262014

ANNOVAR_DB = $(HOME)/share/usr/annovar/humandb

HISAT_REF = $(HOME)/share/reference/hisat/grch37_snp_tran/genome_snp_tran

endif
HG19_INC = true


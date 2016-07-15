ifndef B37_INC
REF = b37
REF_DIR = $(HOME)/ref_nobackup

ifeq ($(SEQ_PLATFORM),IONTORRENT)
ifeq ($(PANEL),CHPv2)
TARGETS_FILE_INTERVALS = $(REF_DIR)/target_beds/b37/AmpliSeq_CancerHotspotPanelv2.106003.designed.bed
TARGETS_FILE_MUTATIONS = $(REF_DIR)/target_beds/b37/AmpliSeq_CancerHotspotPanelv2.106003.mutations.bed
endif
ifeq ($(PANEL),CCP)
TARGETS_FILE_INTERVALS = $(REF_DIR)/target_beds/b37/AmpliSeq.CCP.20131001.designed.nochr.bed
TARGETS_FILE_MUTATIONS = $(REF_DIR)/target_beds/b37/AmpliSeq.CCP.20131001.hotspots.nochr.bed
endif
ifeq ($(PANEL),EXOME)
TARGETS_FILE_INTERVALS = $(REF_DIR)/target_beds/b37/AmpliSeqExome.20131001.designed.nochr.bed
TARGETS_FILE_MUTATIONS = $(REF_DIR)/target_beds/b37/AmpliSeqExome.20131001.designed.nochr.bed
endif
endif

ifeq ($(findstring ILLUMINA,$(SEQ_PLATFORM)),ILLUMINA)
ifeq ($(PANEL),AGILENT_CLINICAL_EXOME)
TARGETS_FILE_COVERED_INTERVALS = $(REF_DIR)/target_beds/b37/S06588914_Agilent_Clinical_Exome_Covered_noHeader_noAnn.bed
TARGETS_FILE_INTERVALS = $(REF_DIR)/target_beds/b37/S06588914_Agilent_Clinical_Exome_Padded_noHeader_noAnn.bed
TARGETS_FILE_MUTATIONS = $(REF_DIR)/target_beds/b37/S06588914_Agilent_Clinical_Exome_Padded_noHeader_noAnn.bed
endif
ifeq ($(findstring IMPACT410,$(PANEL)),IMPACT410)
TARGETS_FILE_INTERVALS = $(REF_DIR)/target_beds/b37/impact410_b37.bed
TARGETS_FILE_MUTATIONS = $(REF_DIR)/target_beds/b37/impact410_b37.bed
endif
ifeq ($(PANEL),IMPACT341)
TARGETS_FILE_INTERVALS = $(REF_DIR)/target_beds/b37/impact341_b37.bed
TARGETS_FILE_MUTATIONS = $(REF_DIR)/target_beds/b37/impact341_b37.bed
endif
endif

ifeq ($(findstring bed,$(TARGETS_FILE_INTERVALS)),bed)
CHROMOSOMES := $(shell grep -v '@' $(TARGETS_FILE_INTERVALS) | cut -f1 | sort | uniq)
else
CHROMOSOMES ?= $(shell seq 1 22) X $(if $(findstring true,$(INCLUDE_CHR_Y)),Y) MT
endif
### GENERAL REFERENCES ###

GSA_REF_DIR = $(REF_DIR)/gsa_resource_bundle/2.8/b37
REF_FASTA ?= $(GSA_REF_DIR)/human_g1k_v37_decoy.fasta


### ANNOTATION SETS ####

ANN_DIR ?= $(REF_DIR)/annotation_data/b37

CONTEST_REF_VCF ?= $(ANN_DIR)/hg19_population_stratified_af_hapmap_3.3.vcf

DBSNP ?= $(GSA_REF_DIR)/dbsnp_138.b37.vcf
KNOWN_INDELS ?= $(GSA_REF_DIR)/Mills_and_1000G_gold_standard.indels.b37.vcf
OMNI ?= $(GSA_REF_DIR)/1000G_omni2.5.b37.vcf
HAPMAP ?= $(GSA_REF_DIR)/hapmap_3.3.b37.vcf

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

CENTROMERE_TABLE = $(ANN_DIR)/centromeres_b37.txt

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
B37_INC = true


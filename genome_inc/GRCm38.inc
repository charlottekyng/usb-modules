# vim: set ft=make:
ifndef GRCM38_INC
REF = GRCm38

REF_DIR = /scicore/home/terracci/GROUP/ref_nobackup
ANN_DIR ?= $(REF_DIR)/annotation_data/$(REF)
BED_DIR ?= $(REF_DIR)/target_beds/$(REF)

ifeq ($(findstring ILLUMINA,$(SEQ_PLATFORM)),ILLUMINA)
ifeq ($(PANEL),AGILENT_ALLEXON_MOUSE)
TARGETS_FILE_COVERED_INTERVALS ?= $(BED_DIR)/S0276129_Covered_noAnn_noHeader_GRCm37liftover.bed
TARGETS_FILE_INTERVALS ?= $(BED_DIR)/S0276129_Padded_noAnn_noHeader_GRCm37liftover.bed
TARGETS_FILE_MUTATIONS ?= $(BED_DIR)/S0276129_Padded_noAnn_noHeader_GRCm37liftover.bed
TARGETS_FILE_GENES ?= $(BED_DIR)/S0276129_Targets.txt
RESOURCE_REQ ?= medium
endif
endif

ifeq ($(findstring bed,$(TARGETS_FILE_INTERVALS)),bed)
CHROMOSOMES := $(shell grep -v '@' $(TARGETS_FILE_INTERVALS) | cut -f1 | sort | uniq)
else
CHROMOSOMES ?= $(shell seq 1 19) X $(if $(findstring true,$(INCLUDE_CHR_Y)),Y) MT
endif

DBSNP := $(ANN_DIR)/mgp.v5.merged.snps_all.indels.dbSNP142_noGT.vcf
KNOWN_INDELS := $(ANN_DIR)/mgp.v5.merged.indels.dbSNP142.normed.noGT.vcf

FACETS_TARGETS_INTERVALS ?= $(DBSNP)

### GENERAL REFERENCES ###

REF_FASTA := $(REF_DIR)/genomes/GRCm38_68.resorted.fa
REF_DICT := $(REF_DIR)/genomes/GRCm38_68.resorted.dict

SNP_EFF_GENOME = GRCm38.79

GENCODE_GTF = $(ANN_DIR)/gencode.vM11.annotation_nochr.gtf

GENE_REF_FLAT = $(ANN_DIR)/refFlat.txt
RIBOSOMAL_INTERVALS = $(ANN_DIR)/ribosomalRNA.interval_list


#### TOOL SPECIFIC OPTIONS ####

STAR_GENOME_DIR = $(REF_DIR)/genome_index/star_GRCm38
RSEM_INDEX = $(REF_DIR)/genome_index/rsem_GRCm38/rsem_index


#INTEGRATE_ANN = $(REF_DIR)/mm10_annot.ucsc.txt
#INTEGRATE_BWTS =$(REF_DIR)/Mus_musculus_GRCm38/bwts
#GENES_GTF = $(REF_DIR)/Mus_musculus/NCBI/$(TOPHAT_REF)/Annotation/Genes/genes.gtf
#CHR_LEN = $(REF_DIR)/Mus_musculus_GRCm38/mm10.len
#BOWTIE_REF = $(REF_DIR)/Mus_musculus_GRCm38/Mus_musculus.GRCm38.71.dna.chromosome.genome
#ENSEMBL_TXDB = $(REF_DIR)/mus_musculus_ensembl_biomart.2014-04-28.sqlite
#TXDB = $(HOME)/ensmusg70.08032013.sqlite

endif
GRCM38_INC = true

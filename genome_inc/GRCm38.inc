# vim: set ft=make:
ifndef GRCM38_INC
REF = GRCm38

REF_DIR = $(HOME)/ref_nobackup

ifeq ($(findstring ILLUMINA,$(SEQ_PLATFORM)),ILLUMINA)
ifeq ($(PANEL),AGILENT_ALLEXON)
TARGETS_FILE_COVERED_INTERVALS = $(REF_DIR)/target_beds/GRCm38/S0276129_Covered_noAnn_noHeader_GRCm37liftover.bed
TARGETS_FILE_INTERVALS = $(REF_DIR)/target_beds/GRCm38/S0276129_Padded_noAnn_noHeader_GRCm37liftover.bed
TARGETS_FILE_MUTATIONS = $(REF_DIR)/target_beds/GRCm38/S0276129_Padded_noAnn_noHeader_GRCm37liftover.bed
TARGETS_FILE_GENES = $(REF_DIR)/target_beds/GRCm38/S0276129_Targets.txt
endif
endif

ifeq ($(findstring bed,$(TARGETS_FILE_INTERVALS)),bed)
CHROMOSOMES := $(shell grep -v '@' $(TARGETS_FILE_INTERVALS) | cut -f1 | sort | uniq)
else
CHROMOSOMES ?= $(shell seq 1 19) X $(if $(findstring true,$(INCLUDE_CHR_Y)),Y) MT
endif

### GENERAL REFERENCES ###

REF_FASTA := $(REF_DIR)/genomes/GRCm38_68.resorted.fa
REF_DICT := $(REF_DIR)/genomes/GRCm38_68.fa.resorted.dict

SNP_EFF_GENOME = GRCm38.79
DBSNP := $(REF_DIR)/annotation_data/GRCm38/mgp.v5.merged.snps_all.indels.dbSNP142_noGT.vcf
KNOWN_INDELS := $(REF_DIR)/annotation_data/GRCm38/mgp.v5.merged.indels.dbSNP142.normed.noGT.vcf

GET_BASE_COUNTS_POS ?= $(REF_DIR)/facets/pileup_GRCm38/snplocs_chr

#INTEGRATE_ANN = $(REF_DIR)/mm10_annot.ucsc.txt
#INTEGRATE_BWTS =$(REF_DIR)/Mus_musculus_GRCm38/bwts
#GENES_GTF = $(REF_DIR)/Mus_musculus/NCBI/$(TOPHAT_REF)/Annotation/Genes/genes.gtf
#CHR_LEN = $(REF_DIR)/Mus_musculus_GRCm38/mm10.len
#BOWTIE_REF = $(REF_DIR)/Mus_musculus_GRCm38/Mus_musculus.GRCm38.71.dna.chromosome.genome
#ENSEMBL_TXDB = $(REF_DIR)/mus_musculus_ensembl_biomart.2014-04-28.sqlite
#TXDB = $(HOME)/ensmusg70.08032013.sqlite

#RIBOSOMAL_INTERVALS = $(REF_DIR)/mm10_rrna_intervals.txt
#GENE_REF_FLAT = $(REF_DIR)/mm10_genes.refFlat.txt

endif
GRCM38_INC = true

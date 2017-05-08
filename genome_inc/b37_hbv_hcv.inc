ifndef B37_HBV_HCV_INC
include usb-modules/genomes/b37.inc
REF = b37_hbv_hcv
REF_DIR = /scicore/home/terracci/GROUP/ref_nobackup


ANN_DIR ?= $(REF_DIR)/annotation_data/b37_hbv_hcv

REF_FASTA ?= $(REF_DIR)/genomes/human_g1k_v37_decoy_HBV_NC_003977.2_HCV.fasta

STAR_GENOME_DIR = $(REF_DIR)/genome_index/star_b37_hbv_hcv
STAR_CTAT_DIR = $(ANN_DIR)/FusionFilter/

### TOOL-SPECIFIC REF ###

#CHIMERASCAN_REF ?= $(REF_DIR)/chimerascan_index
#BOWTIE_REF = $(REF_DIR)/GATK_bundle/2.3/human_g1k_v37

#MAPSPLICE_REF_DIR = $(REF_DIR)/GATK_bundle/2.3/mapsplice_ref/chr
#MAPSPLICE_REF_BASENAME = $(REF_DIR)/GATK_bundle/2.3/mapsplice_ref/hg19

endif
B37_HBV_HCV_INC = true


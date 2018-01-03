ifndef B37_HBV_HCV_GRCM38_INC
include usb-modules/genome_inc/b37.inc
REF = b37_hbv_hcv_GRCm38
REF_DIR = /scicore/home/terracci/GROUP/ref_nobackup

#ANN_DIR := $(REF_DIR)/annotation_data/b37_hbv_hcv
REF_FASTA := $(REF_DIR)/genomes/human_g1k_v37_decoy_HBV_NC_003977.2_HCV_GRCm38.fasta
REF_DICT = $(REF_FASTA:.fasta=.dict)

STAR_GENOME_DIR := $(REF_DIR)/genome_index/star_b37_hbv_hcv_GRCm38
endif
B37_HBV_HCV_GRCM38_INC = true


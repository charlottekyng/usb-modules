ifndef B37_INC
REF = b37
REF_DIR = /scicore/home/terracci/GROUP/ref_nobackup

GSA_REF_DIR ?= $(REF_DIR)/gsa_resource_bundle/2.8/b37
ANN_DIR ?= $(REF_DIR)/annotation_data/b37
BED_DIR ?= $(REF_DIR)/target_beds/b37

REF_FASTA ?= $(GSA_REF_DIR)/human_g1k_v37_decoy.fasta

DBSNP ?= $(GSA_REF_DIR)/dbsnp_138.b37.vcf
DBSNP_TARGETS_INTERVALS ?= $(DBSNP)
DBSNP_TARGETS_INTERVALS_COMMON ?= $(DBSNP_TARGETS_INTERVALS)


ifeq ($(SEQ_PLATFORM),IONTORRENT)
ifeq ($(PANEL),CHPv2)
TARGETS_FILE_INTERVALS = $(BED_DIR)/AmpliSeq_CancerHotspotPanelv2.106003.designed.bed
TARGETS_FILE_MUTATIONS = $(BED_DIR)/AmpliSeq_CancerHotspotPanelv2.106003.mutations.bed
endif
ifeq ($(PANEL),CCP)
TARGETS_FILE_INTERVALS = $(BED_DIR)/AmpliSeq.CCP.20131001.designed.nochr.bed
TARGETS_FILE_MUTATIONS = $(BED_DIR)/AmpliSeq.CCP.20131001.hotspots.nochr.bed
endif
ifeq ($(PANEL),EXOME)
TARGETS_FILE_INTERVALS = $(BED_DIR)/AmpliSeqExome.20131001.designed.nochr.bed
TARGETS_FILE_MUTATIONS = $(BED_DIR)/AmpliSeqExome.20131001.designed.nochr.bed
endif
endif

ifeq ($(findstring ILLUMINA,$(SEQ_PLATFORM)),ILLUMINA)
ifeq ($(PANEL),AGILENT_CLINICAL_EXOME)
TARGETS_FILE_COVERED_INTERVALS = $(BED_DIR)/S06588914_Agilent_Clinical_Exome_Covered_noHeader_noAnn.bed
TARGETS_FILE_INTERVALS = $(BED_DIR)/S06588914_Agilent_Clinical_Exome_Padded_noHeader_noAnn.bed
TARGETS_FILE_MUTATIONS = $(BED_DIR)/S06588914_Agilent_Clinical_Exome_Padded_noHeader_noAnn.bed
TARGETS_FILE_GENES = $(BED_DIR)/S06588914_Agilent_Clinical_Exome_TargetGeneSymbols.txt
DBSNP_TARGETS_INTERVALS = $(BED_DIR)/S06588914_Agilent_Clinical_Exome.dbsnp_138.vcf
DBSNP_TARGETS_INTERVALS_COMMON = $(BED_DIR)/S06588914_Agilent_Clinical_Exome.dbsnp_138.common.vcf
RESOURCE_REQ ?= medium
endif
ifeq ($(PANEL),AGILENT_ALLEXON_V4)
TARGETS_FILE_COVERED_INTERVALS = $(BED_DIR)/S03723314_Agilent_AllExon_V4_Covered_noHeader_noAnn.bed
TARGETS_FILE_INTERVALS = $(BED_DIR)/S03723314_Agilent_AllExon_V4_Padded_noHeader_noAnn.bed
TARGETS_FILE_MUTATIONS = $(BED_DIR)/S03723314_Agilent_AllExon_V4_Padded_noHeader_noAnn.bed
TARGETS_FILE_GENES = $(BED_DIR)/S03723314_Agilent_AllExon_V4_TargetGeneSymbols.txt
DBSNP_TARGETS_INTERVALS = $(BED_DIR)/S03723314_Agilent_AllExon_V4.dbsnp_138.vcf
DBSNP_TARGETS_INTERVALS_COMMON = $(BED_DIR)/S03723314_Agilent_AllExon_V4.dbsnp_138.common.vcf
RESOURCE_REQ ?= medium
endif
ifeq ($(PANEL),AGILENT_ALLEXON_V6)
TARGETS_FILE_COVERED_INTERVALS = $(BED_DIR)/S07604514_Agilent_AllExon_V6_Covered_noHeader_noAnn.bed
TARGETS_FILE_INTERVALS = $(BED_DIR)/S07604514_Agilent_AllExon_V6_Padded_noHeader_noAnn.bed
TARGETS_FILE_MUTATIONS = $(BED_DIR)/S07604514_Agilent_AllExon_V6_Padded_noHeader_noAnn.bed
TARGETS_FILE_GENES = $(BED_DIR)/S07604514_Agilent_AllExon_V6_TargetGeneSymbols.txt
DBSNP_TARGETS_INTERVALS = $(BED_DIR)/S07604514_Agilent_AllExon_V6.dbsnp_138.vcf
DBSNP_TARGETS_INTERVALS_COMMON = $(BED_DIR)/S07604514_Agilent_AllExon_V6.dbsnp_138.common.vcf
RESOURCE_REQ ?= medium
endif
ifeq ($(PANEL),NOVARTIS_PC2)
TARGETS_FILE_COVERED_INTERVALS = $(BED_DIR)/pc2_novartis_Covered_noHeader_noAnn.bed
TARGETS_FILE_INTERVALS = $(BED_DIR)/pc2_novartis_Padded_noHeader_noAnn.bed
TARGETS_FILE_MUTATIONS = $(BED_DIR)/pc2_novartis_Padded_noHeader_noAnn.bed
DBSNP_TARGETS_INTERVALS = $(BED_DIR)/pc2_novartis.dbsnp_138.vcf
RESOURCE_REQ ?= low
endif
ifeq ($(PANEL),IMPACT410)
TARGETS_FILE_INTERVALS = $(BED_DIR)/impact410_b37.bed
TARGETS_FILE_MUTATIONS = $(BED_DIR)/impact410_b37.bed
RESOURCE_REQ ?= low
endif
ifeq ($(PANEL),IMPACT341)
TARGETS_FILE_INTERVALS = $(BED_DIR)/impact341_b37.bed
TARGETS_FILE_MUTATIONS = $(BED_DIR)/impact341_b37.bed
RESOURCE_REQ ?= low
endif
ifeq ($(PANEL),NONE)
RESOURCE_REQ ?= high
BAM_CHR2_BASE_RECAL = true
endif

endif

$(info RESOURCE_REQ is $(RESOURCE_REQ))

FACETS_TARGETS_INTERVALS ?= $(DBSNP_TARGETS_INTERVALS)


ifeq ($(findstring bed,$(TARGETS_FILE_INTERVALS)),bed)
CHROMOSOMES := $(shell grep -v '@' $(TARGETS_FILE_INTERVALS) | cut -f1 | sort | uniq)
else
CHROMOSOMES ?= $(shell seq 1 22) X $(if $(findstring true,$(INCLUDE_CHR_Y)),Y) MT
endif


### ANNOTATION SETS ####

CONTEST_REF_VCF ?= $(ANN_DIR)/hg19_population_stratified_af_hapmap_3.3.vcf

KNOWN_INDELS ?= $(GSA_REF_DIR)/Mills_and_1000G_gold_standard.indels.b37.vcf
OMNI ?= $(GSA_REF_DIR)/1000G_omni2.5.b37.vcf
HAPMAP ?= $(GSA_REF_DIR)/hapmap_3.3.b37.vcf

CLINVAR = $(ANN_DIR)/clinvar.vcf
EXAC_NONTCGA = $(ANN_DIR)/ExAC_nonTCGA.r0.3.1.sites.vep.vcf
COSMIC ?= $(ANN_DIR)/CosmicCodingMuts.sorted.vcf
EXAC_NONPSYCH = $(ANN_DIR)/ExAC.r0.3.nonpsych.sites.edited.vcf

HAPLOTYPE_INSUF_BED = $(ANN_DIR)/haploinsufficient.genes.bed
CANCER_GENE_CENSUS_BED = $(ANN_DIR)/cancer_gene_census.bed
KANDOTH_BED = $(ANN_DIR)/kandoth127smg.bed
LAWRENCE_BED = $(ANN_DIR)/lawrence_cancer5000s.bed
DGD_BED = $(ANN_DIR)/dgd_Hsa_all_v71.bed
CANCER_HOTSPOT_VCF = $(ANN_DIR)/publication_hotspots_inclnoncoding_incl3d.vcf

CENTROMERE_TABLE = $(ANN_DIR)/centromeres_b37.txt

DGV_FILE = $(ANN_DIR)/DGV_GRCh37_hg19_variants_2016-05-15.txt
#KNOWN_CNVS = $(ANN_DIR)/GRCh37_hg19_variants_2013-07-23.txt
#REF_HOUSEKEEPING_GENE_BED = $(ANN_DIR)/share/reference/hg19.HouseKeepingGenes.bed
#ENSEMBL_TXDB = $(ANN_DIR)/hsapiens_ensembl_biomart.sqlite

STAR_GENOME_DIR = $(REF_DIR)/genome_index/star_b37
STAR_CTAT_DIR = $(ANN_DIR)/GRCh37_gencode_v19_CTAT_lib

RIBOSOMAL_INTERVALS ?= $(ANN_DIR)/ribosomalRNA.interval_list
GENE_REF_FLAT = $(ANN_DIR)/refFlat.txt

#SNP_TXT = $(ANN_DIR)/b37_snp137.SingleDiNucl.1based.nochr.txt

### TOOL-SPECIFIC REF ###

GISTIC_REF ?= /scicore/home/terracci/GROUP/usr_nobackup/local/gistic_2.0.23/refgenefiles/hg19.UCSC.add_miR.140312.refgene.mat

#CHIMERASCAN_REF ?= $(REF_DIR)/chimerascan_index
#BOWTIE_REF = $(REF_DIR)/GATK_bundle/2.3/human_g1k_v37

#MAPSPLICE_REF_DIR = $(REF_DIR)/GATK_bundle/2.3/mapsplice_ref/chr
#MAPSPLICE_REF_BASENAME = $(REF_DIR)/GATK_bundle/2.3/mapsplice_ref/hg19

DB_NSFP ?= $(REF_DIR)/snpeff/GRCh37/dbNSFP-v2.9.txt.gz

#TOPHAT_REF = GRCh37
#GENES_GTF = $(REF_DIR)/Homo_sapiens/Ensembl/$(TOPHAT_REF)/Annotation/Genes/genes.gtf

#ANNOVAR_REF = hg19

#UCSC_REFSEQ_GTF = $(REF_DIR)s/genomes/ucsc.refseq.gtf
#REFSEQ_GTF := $(UCSC_REFSEQ_GTF)

#CHR_LEN = $(REF_DIR)/GATK_bundle/2.3/b37.len

#REF_2BIT ?= $(REF_DIR)/GATK_bundle/2.3/human_g1k_v37.2bit

#MAP_BIGWIG ?= $(REF_DIR)/wgEncodeCrgMapabilityAlign100mer.nochr.bigWig

REF_NAME := 1000GenomesPilot-NCBI37
REF_DATE := 20100517
REF_DICT = $(REF_FASTA:.fasta=.dict)
SNP_EFF_GENOME = GRCh37.75

MUTSIGCV_REF_DIR ?= $(REF_DIR)/MutSigCV
MUTSIGCV_COVERAGE_REF ?= $(MUTSIGCV_REF_DIR)/coverage.exome_full192.txt
MUTSIGCV_COV_REF ?= $(MUTSIGCV_REF_DIR)/gene.covariates.txt
MUTSIGCV_DICT_REF ?= $(MUTSIGCV_REF_DIR)/mutation_type_dictionary_file.txt
MUTSIGCV_SEQ_REF_DIR ?= $(MUTSIGCV_REF_DIR)/chr_files_hg19/

#EMU_REF_DIR ?= $(REF_DIR)/GATK_bundle/2.3/emu_ref
#EMU_REF_MUTATIONS ?= $(REF_DIR)/emu_ref/exome_ref/Breast_clean_somatic_mutations_for_signature_analysis.emu.mutations.txt
#EMU_REF_CNV ?= $(REF_DIR)/emu_ref/exome_ref/Breast_clean_somatic_mutations_for_signature_analysis.emu.copynum.txt

#HMMCOPY_GC_WIG = $(REF_DIR)/GATK_bundle/2.3/gc_hg19.wig
#HMMCOPY_MAP_WIG = $(REF_DIR)/map_hg19.wig


#INTEGRATE_ANN = $(REF_DIR)/annot.ucsc.txt
#INTEGRATE_BWTS = $(REF_DIR)/GATK_bundle/2.3/bwts
#ONCOTATOR_DB ?= $(REF_DIR)/oncotator_v1_ds_Jan262014

#ANNOVAR_DB = $(HOME)/share/usr/annovar/humandb

#HISAT_REF = $(HOME)/share/reference/hisat/grch37_snp_tran/genome_snp_tran
endif
B37_INC = true


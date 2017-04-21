include usb-modules/Makefile.inc
include usb-modules/config.inc
include usb-modules/variant_callers/somatic/somaticVariantCaller.inc

LOGDIR = log/summary.$(NOW)

#.DELETE_ON_ERROR:
.SECONDARY: 
PHONY : mutation_summary

ifeq ($(findstring ILLUMINA,$(SEQ_PLATFORM)),ILLUMINA)
CALLER_PREFIX ?= mutect strelka_indel
endif
ifeq ($(findstring IONTORRENT,$(SEQ_PLATFORM)),IONTORRENT)
CALLER_PREFIX ?= tvc_snps tvc_indels varscan_snps varscan_indels
endif

# Add optional absolute results to excel
# the $(wildcard x) syntax is used to check for existence of file
ABSOLUTE_SOMATIC_TXTS ?= $(wildcard $(foreach set,$(SAMPLE_SETS),absolute/tables/$(set).somatic.txt))
ABSOLUTE_SEGMENTS ?= $(wildcard $(foreach set,$(SAMPLE_SETS),absolute/reviewed/SEG_MAF/$(set)_ABS_MAF.txt))
ifneq ($(ABSOLUTE_SOMATIC_TXTS)$(ABSOLUTE_SEGMENTS),)
EXCEL_ABSOLUTE_PARAMS = --absolute_segments $(subst $(space),$(,),$(strip $(ABSOLUTE_SEGMENTS))) --absolute_somatic_txts $(subst $(space),$(,),$(strip $(ABSOLUTE_SOMATIC_TXTS)))
endif

# Add optional annotations file to excel, used to add any other desired columns
# to the excel file inner merge on TUMOR_SAMPLE CHROM POS REF ALT
EXCEL_ANNOTATION ?= $(wildcard summary/annotation.tsv)
ifneq ($(EXCEL_ANNOTATION),)
EXCEL_ANNOTATION_PARAMS = --annotation_tsv $(EXCEL_ANNOTATION)
endif

# Set max for ExAC_AF column in MUTATION_SUMMARY sheet and tsv
EXCEL_MAX_EXAC_AF ?= 1


ifeq ($(findstring EXCEL,$(MUTATION_SUMMARY_FORMAT)),EXCEL)
mutation_summary : $(shell rm -f summary/mutation_summary.xlsx) summary/mutation_summary.xlsx

ALLTABLES_NS_SYNON_HS_LNC = $(foreach prefix,$(CALLER_PREFIX),alltables/allTN.$(call SOMATIC_VCF_SUFFIXES,$(prefix)).tab.nonsynonymous_synonymous_hotspot_lincRNA.txt)
ALLTABLES_NS_SYNON_HS = $(foreach prefix,$(CALLER_PREFIX),alltables/allTN.$(call SOMATIC_VCF_SUFFIXES,$(prefix)).tab.nonsynonymous_synonymous_hotspot.txt)

summary/mutation_summary.xlsx : $(if $(findstring false,$(INCLUDE_LINCRNA_IN_SUMMARY)),$(ALLTABLES_NS_SYNON_HS),$(ALLTABLES_NS_SYNON_HS_LNC))
	$(call LSCRIPT_CHECK_MEM,6G,00:29:59,"$(LOAD_R_MODULE); $(RSCRIPT) usb-modules/summary/mutation_summary_excel.R \
	--outFile $@ $^")
endif

ifeq ($(findstring TXT,$(MUTATION_SUMMARY_FORMAT)),TXT)
mutation_summary : $(foreach prefix,$(CALLER_PREFIX),$(shell rm -f summary/mutation_summary_$(prefix).txt) summary/mutation_summary_$(prefix).txt)

define mut_sum_txt
summary/mutation_summary_$1.txt : $$(if $$(findstring false,$$(INCLUDE_LINCRNA_IN_SUMMARY)),alltables/allTN.$$(call SOMATIC_VCF_SUFFIXES,$1).tab.nonsynonymous_synonymous_hotspot.txt,alltables/allTN.$$(call SOMATIC_VCF_SUFFIXES,$1).tab.nonsynonymous_synonymous_hotspot_lincRNA.txt)
	$$(call LSCRIPT_CHECK_MEM,4G,00:29:29,"$$(LOAD_R_MODULE); $$(RSCRIPT) usb-modules/summary/mutation_summary_excel.R \
	--outFile $$@ --outputFormat TXT $$<")
endef
$(foreach prefix,$(CALLER_PREFIX),$(eval $(call mut_sum_txt,$(prefix))))
endif



#
#summary/mutation_summary_snps.txt : $(if $(findstring false,$(INCLUDE_LINCRNA_IN_SUMMARY)),$(ALLTABLES_NONSYNONYMOUS_SYNONYMOUS_HOTSPOT_SNPS),$(ALLTABLES_COMPLETE_SNPS))
#	$(call LSCRIPT_CHECK_MEM,4G,00:29:29,"$(LOAD_R_MODULE); $(RSCRIPT) usb-modules/summary/mutation_summary_excel.R \
#	--outFile $@ --outputFormat TXT $<")
#summary/mutation_summary_indels.txt : $(if $(findstring false,$(INCLUDE_LINCRNA_IN_SUMMARY)),$(ALLTABLES_NONSYNONYMOUS_SYNONYMOUS_HOTSPOT_INDELS),$(ALLTABLES_COMPLETE_INDELS))
#	$(call LSCRIPT_CHECK_MEM,4G,00:29:29,"$(LOAD_R_MODULE); $(RSCRIPT) usb-modules/summary/mutation_summary_excel.R \
#	--outFile $@ --outputFormat TXT $<")

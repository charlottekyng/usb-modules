include usb-modules/Makefile.inc
include usb-modules/config.inc
include usb-modules/variant_callers/somatic/somaticVariantCaller.inc

LOGDIR = log/summary.$(NOW)

#.DELETE_ON_ERROR:
.SECONDARY: 
PHONY : mutation_summary

ifeq ($(findstring ILLUMINA,$(SEQ_PLATFORM)),ILLUMINA)
ALLTABLES_COMPLETE_SNPS = alltables/allTN.$(call SOMATIC_VCF_SUFFIXES,mutect).tab.nonsynonymous_synonymous_hotspot_lincRNA.txt
ALLTABLES_COMPLETE_INDELS = alltables/allTN.$(call SOMATIC_VCF_SUFFIXES,strelka_indels).tab.nonsynonymous_synonymous_hotspot_lincRNA.txt
#ALLTABLES_LOW_MODIFIER_SNPS = alltables/allTN.$(call SOMATIC_VCF_SUFFIXES,mutect).tab.low_modifier.txt
#ALLTABLES_SYNONYMOUS_SNPS = alltables/allTN.$(call SOMATIC_VCF_SUFFIXES,mutect).tab.synonymous.txt
#ALLTABLES_NONSYNONYMOUS_SNPS = alltables/allTN.$(call SOMATIC_VCF_SUFFIXES,mutect).tab.nonsynonymous.txt
#ALLTABLES_LOW_MODIFIER_INDELS = alltables/allTN.$(call SOMATIC_VCF_SUFFIXES,strelka_indels).tab.low_modifier.txt
#ALLTABLES_NONSYNONYMOUS_INDELS = alltables/allTN.$(call SOMATIC_VCF_SUFFIXES,strelka_indels).tab.nonsynonymous.txt
endif
ifeq ($(findstring IONTORRENT,$(SEQ_PLATFORM)),IONTORRENT)
ALLTABLES_COMPLETE_SNPS = alltables/allTN.$(call SOMATIC_VCF_SUFFIXES,tvc_snps).tab.nonsynonymous_synonymous_hotspot_lincRNA.txt
ALLTABLES_COMPLETE_INDELS = alltables/allTN.$(call SOMATIC_VCF_SUFFIXES,tvc_indels).tab.nonsynonymous_synonymous_hotspot_lincRNA.txt
#ALLTABLES_LOW_MODIFIER_SNPS = alltables/allTN.$(call SOMATIC_VCF_SUFFIXES,tvc_snps).tab.low_modifier.txt
#ALLTABLES_SYNONYMOUS_SNPS = alltables/allTN.$(call SOMATIC_VCF_SUFFIXES,tvc_snps).tab.synonymous.txt
#ALLTABLES_NONSYNONYMOUS_SNPS = alltables/allTN.$(call SOMATIC_VCF_SUFFIXES,tvc_snps).tab.nonsynonymous.txt
#ALLTABLES_LOW_MODIFIER_INDELS = alltables/allTN.$(call SOMATIC_VCF_SUFFIXES,tvc_indels).tab.low_modifier.txt
#ALLTABLES_NONSYNONYMOUS_INDELS = alltables/allTN.$(call SOMATIC_VCF_SUFFIXES,tvc_indels).tab.nonsynonymous.txt
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

summary/mutation_summary.xlsx : $(ALLTABLES_COMPLETE_SNPS) $(ALLTABLES_COMPLETE_INDELS)
	$(call LSCRIPT_CHECK_MEM,6G,00:29:59,"$(LOAD_R_MODULE); $(RSCRIPT) usb-modules/summary/mutation_summary_excel.R \
	--outFile $@ $(wordlist 1,2,$^)")
endif

ifeq ($(findstring TXT,$(MUTATION_SUMMARY_FORMAT)),TXT)
mutation_summary : $(shell rm -f summary/mutation_summary_snps.txt) summary/mutation_summary_snps.txt $(shell rm -f summary/mutation_summary_indels.txt) summary/mutation_summary_indels.txt

summary/mutation_summary_snps.txt : $(ALLTABLES_COMPLETE_SNPS)
	$(call LSCRIPT_CHECK_MEM,4G,00:29:29,"$(LOAD_R_MODULE); $(RSCRIPT) usb-modules/summary/mutation_summary_excel.R \
	--outFile $@ --outputFormat TXT $<")
summary/mutation_summary_indels.txt : $(ALLTABLES_COMPLETE_INDELS)
	$(call LSCRIPT_CHECK_MEM,4G,00:29:29,"$(LOAD_R_MODULE); $(RSCRIPT) usb-modules/summary/mutation_summary_excel.R \
	--outFile $@ --outputFormat TXT $<")
endif

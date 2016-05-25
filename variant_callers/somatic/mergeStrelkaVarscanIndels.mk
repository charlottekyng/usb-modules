# Merge strelka and varscan indel results
LOGDIR ?= log/merge_strelka_varscan_indels.$(NOW)

include modules/Makefile.inc
include modules/config.inc
include modules/variant_callers/gatk.inc
include modules/variant_callers/somatic/somaticVariantCaller.inc

.DELETE_ON_ERROR:
.SECONDARY: 
.PHONY : strelka_varscan_merge strelka_varscan_merge_vcfs strelka_varscan_merge_tables strelka_varscan_all_tables

strelka_varscan_merge : strelka_varscan_merge_vcfs strelka_varscan_merge_tables strelka_varscan_all_tables
strelka_varscan_merge_vcfs : $(foreach pair,$(SAMPLE_PAIRS),vcf/$(pair).strelka_varscan_indels.vcf)
strelka_varscan_merge_tables : $(foreach pair,$(SAMPLE_PAIRS),\
	$(foreach ext,$(SOMATIC_TABLE_EXTENSIONS),tables/$(pair).strelka_varscan_indels.$(ext).txt))
strelka_varscan_all_tables : $(foreach ext,$(SOMATIC_TABLE_EXTENSIONS),alltables/allTN.strelka_varscan_indels.$(ext).txt)

vcf/%.strelka_varscan_indels.vcf : vcf/%.$(call SOMATIC_VCF_SUFFIXES,varscan_indels).vcf vcf/%.$(call SOMATIC_VCF_SUFFIXES,strelka_indels).vcf
	$(call LSCRIPT_MEM,9G,12G,"grep -P '^#' $< > $@ && $(BEDTOOLS) intersect -a $< -b $(<<) >> $@")

include modules/vcf_tools/vcftools.mk

include usb-modules/Makefile.inc
include usb-modules/config.inc
include usb-modules/variant_callers/somatic/somaticVariantCaller.inc

LOGDIR ?= log/tvc_somtic.$(NOW)

VPATH ?= bam
VARIANT_TYPES ?= tvc_snps tvc_indels
PHONY += tvc_somatic

tvc_somatic : tvc_vcfs tvc_tables
tvc_vcfs : $(foreach type,$(VARIANT_TYPES),$(call SOMATIC_VCFS,$(type)) $(addsuffix .idx,$(call SOMATIC_VCFS,$(type))))
tvc_tables : $(foreach type,$(VARIANT_TYPES),$(call SOMATIC_TABLES,$(type)))

define tvc-somatic-vcf
tvc/vcf/$1_$2/TSVC_variants.vcf.gz : bam/$1.bam bam/$1.bam.bai bam/$2.bam bam/$2.bam.bai
	$$(call LSCRIPT_PARALLEL_MEN,8,3G,00:59:59,"$$(TVC) -i $$< -n $$(word 3,$$^) -r $$(REF_FASTA) -o $$(@D) -N 8 \
	$$(if $$(TARGETS_FILE_INTERVALS),-b $$(TARGETS_FILE_INTERVALS)) -p $$(TVC_SOMATIC_JSON) -m $$(TVC_MOTIF) \
	-t $$(TVC_ROOT_DIR) --primer-trim-bed $$(PRIMER_TRIM_BED)")
endef
$(foreach pair,$(SAMPLE_PAIRS), \
	$(eval $(call tvc-somatic-vcf,$(tumor.$(pair)),$(normal.$(pair)))))

tvc/vcf/%/TSVC_variants.snps.vcf : tvc/vcf/%/TSVC_variants.vcf.gz
	$(call LSCRIPT_CHECK_MEM,5G,00:29:59,"$(LOAD_JAVA8_MODULE); $(call SELECT_VARIANTS,7G) \
	-R $(REF_FASTA) --variant $< -o $@ -selectType SNP")

tvc/vcf/%/TSVC_variants.indels.vcf : tvc/vcf/%/TSVC_variants.vcf.gz
	$(call LSCRIPT_CHECK_MEM,5G,00:29:59,"$(LOAD_JAVA8_MODULE); $(call SELECT_VARIANTS,7G) \
	-R $(REF_FASTA) --variant $< -o $@ -selectType INDEL")

vcf/%.tvc_snps.vcf : tvc/vcf/%/TSVC_variants.snps.vcf
	$(INIT) ln -f $< $@

vcf/%.tvc_indels.vcf : tvc/vcf/%/TSVC_variants.indels.vcf
	$(INIT) ln -f $< $@

.DELETE_ON_ERROR:
.SECONDARY:
.PHONY : $(PHONY)

include usb-modules/vcf_tools/vcftools.mk

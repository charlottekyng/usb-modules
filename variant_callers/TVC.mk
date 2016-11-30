include usb-modules/Makefile.inc
include usb-modules/config.inc
include usb-modules/variant_callers/variantCaller.inc

LOGDIR ?= log/tvc.$(NOW)

VPATH ?= bam
VARIANT_TYPES ?= tvc_snps tvc_indels

tvc : tvc_vcfs tvc_tables
tvc_vcfs : $(foreach type,$(VARIANT_TYPES),$(call VCFS,$(type)) $(addsuffix .idx,$(call VCFS,$(type))))
tvc_tables : $(foreach type,$(VARIANT_TYPES),$(call TABLES,$(type)))

tvc/dbsnp/%/TSVC_variants.vcf : bam/%.bam
	$(call LSCRIPT_PARALLEL_MEM,4,10G,05:59:59,"$(TVC) -s $(DBSNP_TARGETS_INTERVALS) -i $< -r $(REF_FASTA) -o $(@D) -N 4 \
	$(if $(TARGETS_FILE_INTERVALS),-b $(TARGETS_FILE_INTERVALS)) -m $(TVC_MOTIF) \
	-t $(TVC_ROOT_DIR) --primer-trim-bed $(PRIMER_TRIM_BED)")


define tvc-vcf
tvc/vcf/$1/TSVC_variants.vcf.gz : bam/$1.bam bam/$1.bai
	$$(call LSCRIPT_PARALLEL_MEN,4,8G,00:29:59,"$$(TVC) -i $$< -r $$(REF_FASTA) -o $$(@D) -N 4 \
	$$(if $(TARGETS_FILE_INTERVALS),-b $$(TARGETS_FILE_INTERVALS)) -m $$(TVC_MOTIF) \
	-t $$(TVC_ROOT_DIR) --primer-trim-bed $$(PRIMER_TRIM_BED)")
endef
$(foreach sample,$(SAMPLES), \
	$(eval $(call tvc-vcf,$(sample))))

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

vcf/%.tvc_snps_sufam.vcf : tvc/sufam/%.tvc_snps_sufam.vcf
	$(INIT) ln -f $< $@

vcf/%.tvc_indels_sufam.vcf : tvc/sufam/%.tvc_indels_sufam.vcf
	$(INIT) ln -f $< $@

.DELETE_ON_ERROR:
.SECONDARY:
.PHONY : $(PHONY)

include usb-modules/vcf_tools/vcftools.mk

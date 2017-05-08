include usb-modules/Makefile.inc
include usb-modules/config.inc
include usb-modules/variant_callers/variantCaller.inc

LOGDIR ?= log/tvc.$(NOW)

VPATH ?= bam
VARIANT_TYPES ?= tvc_snps tvc_indels

tvc : tvc_vcfs tvc_tables
tvc_vcfs : $(foreach type,$(VARIANT_TYPES),$(call VCFS,$(type)) $(addsuffix .idx,$(call VCFS,$(type))))
tvc_tables : $(foreach type,$(VARIANT_TYPES),$(call TABLES,$(type)))

MUT_CALLER = tvc

tvc/dbsnp/%/TSVC_variants.vcf : bam/%.bam
	$(call LSCRIPT_PARALLEL_MEM,4,10G,05:59:59,"$(TVC) -s $(DBSNP_TARGETS_INTERVALS) -i $< -r $(REF_FASTA) -o $(@D) -N 4 \
	$(if $(TARGETS_FILE_INTERVALS),-b $(TARGETS_FILE_INTERVALS)) -m $(TVC_MOTIF) \
	-t $(TVC_ROOT_DIR) --primer-trim-bed $(PRIMER_TRIM_BED)")


define tvc-vcf
tvc/vcf/$1/TSVC_variants_final.vcf : bam/$1.bam bam/$1.bai
	$$(call LSCRIPT_PARALLEL_MEN,4,8G,00:29:59,"$$(LOAD_JAVA8_MODULE);  $$(TVC) -i $$< -r $$(REF_FASTA) -o $$(@D) -N 4 \
	$$(if $(TARGETS_FILE_INTERVALS),-b $$(TARGETS_FILE_INTERVALS)) -m $$(TVC_MOTIF) \
	-t $$(TVC_ROOT_DIR) --primer-trim-bed $$(PRIMER_TRIM_BED) && \
	$$(call LEFT_ALIGN_VCF,6G) -R $$(REF_FASTA) --variant $$(@D)/TSVC_variants.vcf -o $$(@D)/TSVC_variants.vcf.tmp --splitMultiallelics && \
	$$(FIX_TVC_VCF) $$(@D)/TSVC_variants.vcf.tmp > $$@ && $$(RM) $$@.gz")
endef
$(foreach sample,$(SAMPLES), \
	$(eval $(call tvc-vcf,$(sample))))

tvc/vcf/%/TSVC_variants.snps.vcf : tvc/vcf/%/TSVC_variants_final.vcf
	$(call LSCRIPT_CHECK_MEM,5G,00:29:59,"$(LOAD_VCFTOOLS_MODULE); $(VCFTOOLS) --vcf $< --remove-indels --recode --recode-INFO-all --out $@ && mv $@.recode.vcf $@")

tvc/vcf/%/TSVC_variants.indels.vcf : tvc/vcf/%/TSVC_variants_final.vcf
	$(call LSCRIPT_CHECK_MEM,5G,00:29:59,"$(LOAD_VCFTOOLS_MODULE); $(VCFTOOLS) --vcf $< --keep-only-indels --recode --recode-INFO-all --out $@ && mv $@.recode.vcf $@")

vcf/%.tvc_snps.vcf : tvc/vcf/%/TSVC_variants.snps.vcf
	$(INIT) $(VCF_SORT) $(REF_DICT) $< | $(FIX_TVC_VCF) > $@

vcf/%.tvc_indels.vcf : tvc/vcf/%/TSVC_variants.indels.vcf
	$(INIT) $(VCF_SORT) $(REF_DICT) $< | $(FIX_TVC_VCF) > $@

.DELETE_ON_ERROR:
.SECONDARY:
.PHONY : $(PHONY)

include usb-modules/vcf_tools/vcftools.mk

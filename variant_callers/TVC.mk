include usb-modules/Makefile.inc
include usb-modules/config.inc
include usb-modules/variant_callers/variantCaller.inc

LOGDIR ?= log/tvc.$(NOW)

VPATH ?= bam
VARIANT_TYPES ?= tvc_snps tvc_indels

PHONY += tvc tvc_vcfs tvc_tables
tvc : tvc_vcfs tvc_tables
tvc_vcfs : $(foreach type,$(VARIANT_TYPES),$(call VCFS,$(type)) $(addsuffix .idx,$(call VCFS,$(type))))
tvc_tables : $(foreach type,$(VARIANT_TYPES),$(call TABLES,$(type)))

MUT_CALLER = tvc

tvc/dbsnp/%/TSVC_variants.vcf : bam/%.bam
	$(call LSCRIPT_PARALLEL_MEM,4,$(RESOURCE_REQ_MEDIUM_MEM),$(RESOURCE_REQ_LONG),"$(TVC) \
	-s $(DBSNP_TARGETS_INTERVALS) -i $< -r $(REF_FASTA) -o $(@D) -N 4 -p $(TVC_SOMATIC_JSON) \
	$(if $(TARGETS_FILE_INTERVALS),-b $(TARGETS_FILE_INTERVALS)) -m $(TVC_MOTIF) \
	-t $(TVC_ROOT_DIR) --primer-trim-bed $(PRIMER_TRIM_BED)")


define tvc-vcf
tvc/vcf/$1/TSVC_variants.vcf.gz : bam/$1.bam bam/$1.bai
	$$(call LSCRIPT_PARALLEL_MEN,4,$$(RESOURCE_REQ_MEDIUM_MEM),$$(RESOURCE_REQ_VSHORT),"\
	$$(TVC) -i $$< -r $$(REF_FASTA) -o $$(@D) -N 4 \
	$$(if $(TARGETS_FILE_INTERVALS),-b $$(TARGETS_FILE_INTERVALS)) -m $$(TVC_MOTIF) \
	-t $$(TVC_ROOT_DIR) --primer-trim-bed $$(PRIMER_TRIM_BED)")

tvc/vcf/$1/isec/0000.vcf : tvc/vcf/$1/TSVC_variants.multiallelic_ft.norm.left_align.vcf tvc/vcf/$1/TSVC_variants.biallelic_ft.vcf tvc/vcf/$1/TSVC_variants.multiallelic_ft.norm.left_align.vcf.gz tvc/vcf/$1/TSVC_variants.biallelic_ft.vcf.gz tvc/vcf/$1/TSVC_variants.multiallelic_ft.norm.left_align.vcf.gz.tbi tvc/vcf/$1/TSVC_variants.biallelic_ft.vcf.gz.tbi
	$$(call LSCRIPT_MEM,$$(RESOURCE_REQ_LOWMEM),$$(RESOURCE_REQ_VSHORT),"$$(LOAD_BCFTOOLS_MODULE); \
	$$(BCFTOOLS) isec -O v -p $$(dir $$@) $$(word 3,$$^) $$(word 4,$$^)")

tvc/vcf/$1/isec/0001.vcf : tvc/vcf/$1/isec/0000.vcf
	
tvc/vcf/$1/isec/0003.vcf : tvc/vcf/$1/isec/0000.vcf
	

tvc/vcf/$1/TSVC_variants_final.vcf : tvc/vcf/$1/isec/0000.post_bcftools.vcf tvc/vcf/$1/isec/0001.post_bcftools.vcf tvc/vcf/$1/isec/0003.post_bcftools.vcf
	$$(call LSCRIPT_MEM,$$(RESOURCE_REQ_MEDIUM_MEM),$$(RESOURCE_REQ_VSHORT),"$$(LOAD_JAVA8_MODULE); \
	$$(call GATK,CombineVariants,$$(RESOURCE_REQ_MEDIUM_MEM)) -R $$(REF_FASTA) \
	$$(foreach vcf,$$^,-V $$(vcf)) -o $$@ --assumeIdenticalSamples")

endef
$(foreach sample,$(SAMPLES), \
	$(eval $(call tvc-vcf,$(sample))))

tvc/vcf/%/TSVC_variants.snps.vcf : tvc/vcf/%/TSVC_variants_final.vcf
	$(call LSCRIPT_CHECK_MEM,$(RESOURCE_REQ_LOWMEM),$(RESOURCE_REQ_VSHORT),"$(LOAD_VCFTOOLS_MODULE); \
	$(VCFTOOLS) --vcf $< --remove-indels --recode --recode-INFO-all --out $@ && mv $@.recode.vcf $@")

tvc/vcf/%/TSVC_variants.indels.vcf : tvc/vcf/%/TSVC_variants_final.vcf
	$(call LSCRIPT_CHECK_MEM,$(RESOURCE_REQ_LOWMEM),$(RESOURCE_REQ_VSHORT),"$(LOAD_VCFTOOLS_MODULE); \
	$(VCFTOOLS) --vcf $< --keep-only-indels --recode --recode-INFO-all --out $@ && mv $@.recode.vcf $@")

vcf/%.tvc_snps.vcf : tvc/vcf/%/TSVC_variants.snps.vcf
	$(INIT) $(VCF_SORT) $(REF_DICT) $< > $@
#	$(INIT) $(VCF_SORT) $(REF_DICT) $< | $(FIX_TVC_VCF) > $@

vcf/%.tvc_indels.vcf : tvc/vcf/%/TSVC_variants.indels.vcf
	$(INIT) $(VCF_SORT) $(REF_DICT) $< > $@
#	$(INIT) $(VCF_SORT) $(REF_DICT) $< | $(FIX_TVC_VCF) > $@

.DELETE_ON_ERROR:
.SECONDARY:
.PHONY : $(PHONY)

include usb-modules/vcf_tools/vcftools.mk

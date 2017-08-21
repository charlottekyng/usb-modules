include usb-modules/Makefile.inc
include usb-modules/config.inc
include usb-modules/variant_callers/somatic/somaticVariantCaller.inc

LOGDIR ?= log/tvc_somatic.$(NOW)

PHONY += tvc_somatic tvc_somatic_vcfs tvc_somatic_tables 
tvc_somatic : tvc_somatic_vcfs tvc_somatic_tables 

VARIANT_TYPES ?= tvc_snps tvc_indels
tvc_somatic_vcfs : $(foreach type,$(VARIANT_TYPES),$(call SOMATIC_VCFS,$(type)) $(addsuffix .idx,$(call SOMATIC_VCFS,$(type))))
tvc_somatic_tables : $(foreach type,$(VARIANT_TYPES),$(call SOMATIC_TABLES,$(type)))

.DELETE_ON_ERROR:
.SECONDARY:
.PHONY : $(PHONY)

MUT_CALLER = tvc
define tvc-somatic-vcf
tvc/vcf/$1_$2/TSVC_variants.vcf.gz : bam/$1.bam bam/$1.bam.bai bam/$2.bam bam/$2.bam.bai
	$$(call LSCRIPT_PARALLEL_MEM,8,$$(RESOURCE_REQ_HIGHMEM),$$(RESOURCE_REQ_LONG),"\
	$$(TVC) -i $$< -n $$(word 3,$$^) -r $$(REF_FASTA) -o $$(@D) -N 8 \
	$$(if $$(TARGETS_FILE_INTERVALS),-b $$(TARGETS_FILE_INTERVALS)) -p $$(TVC_SOMATIC_JSON) -m $$(TVC_MOTIF) \
	-t $$(TVC_ROOT_DIR) --primer-trim-bed $$(PRIMER_TRIM_BED) -g $$(basename $$(notdir $$<))")


tvc/vcf/$1_$2/isec/0000.vcf : tvc/vcf/$1_$2/TSVC_variants.multiallelic_ft.norm.left_align.vcf tvc/vcf/$1_$2/TSVC_variants.biallelic_ft.vcf tvc/vcf/$1_$2/TSVC_variants.multiallelic_ft.norm.left_align.vcf.gz tvc/vcf/$1_$2/TSVC_variants.biallelic_ft.vcf.gz tvc/vcf/$1_$2/TSVC_variants.multiallelic_ft.norm.left_align.vcf.gz.tbi tvc/vcf/$1_$2/TSVC_variants.biallelic_ft.vcf.gz.tbi
	$$(call LSCRIPT_MEM,$$(RESOURCE_REQ_LOWMEM),$$(RESOURCE_REQ_VSHORT),"$$(LOAD_BCFTOOLS_MODULE); \
	$$(BCFTOOLS) isec -O v -p $$(dir $$@) $$(word 3,$$^) $$(word 4,$$^)")

tvc/vcf/$1_$2/isec/0001.vcf : tvc/vcf/$1_$2/isec/0000.vcf
tvc/vcf/$1_$2/isec/0003.vcf : tvc/vcf/$1_$2/isec/0000.vcf

tvc/vcf/$1_$2/TSVC_variants_final.vcf : tvc/vcf/$1_$2/isec/0000.post_bcftools.vcf tvc/vcf/$1_$2/isec/0001.post_bcftools.vcf tvc/vcf/$1_$2/isec/0003.post_bcftools.vcf
	$$(call LSCRIPT_MEM,$$(RESOURCE_REQ_MEDIUM_MEM),$$(RESOURCE_REQ_VSHORT),"$$(LOAD_JAVA8_MODULE); \
	$$(call GATK,CombineVariants,$$(RESOURCE_REQ_MEDIUM_MEM)) -R $$(REF_FASTA) \
	$$(foreach vcf,$$^,-V $$(vcf)) -o $$@ --assumeIdenticalSamples")

endef
$(foreach pair,$(SAMPLE_PAIRS), \
	$(eval $(call tvc-somatic-vcf,$(tumor.$(pair)),$(normal.$(pair)))))

include usb-modules/vcf_tools/vcftools.mk
include usb-modules/variant_callers/TVC.mk
include usb-modules/variant_callers/somatic/pon.mk

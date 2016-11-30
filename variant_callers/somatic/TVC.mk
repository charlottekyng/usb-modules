include usb-modules/Makefile.inc
include usb-modules/config.inc
include usb-modules/variant_callers/somatic/somaticVariantCaller.inc

LOGDIR ?= log/tvc_somtic.$(NOW)

PHONY += tvc_somatic tvc_somatic_vcfs tvc_somatic_tables
tvc_somatic : tvc_somatic_vcfs tvc_somatic_tables

VARIANT_TYPES ?= tvc_snps tvc_indels
tvc_somatic_vcfs : $(foreach type,$(VARIANT_TYPES),$(call SOMATIC_VCFS,$(type)) $(addsuffix .idx,$(call SOMATIC_VCFS,$(type))))
tvc_somatic_tables : $(foreach type,$(VARIANT_TYPES),$(call SOMATIC_TABLES,$(type)))

.DELETE_ON_ERROR:
.SECONDARY:
.PHONY : $(PHONY)

define tvc-somatic-vcf
tvc/vcf/$1_$2/TSVC_variants.vcf.gz : bam/$1.bam bam/$1.bam.bai bam/$2.bam bam/$2.bam.bai
	$$(call LSCRIPT_PARALLEL_MEM,4,10G,05:59:59,"$$(TVC) -i $$< -n $$(word 3,$$^) -r $$(REF_FASTA) -o $$(@D) -N 4 \
	$$(if $$(TARGETS_FILE_INTERVALS),-b $$(TARGETS_FILE_INTERVALS)) -p $$(TVC_SOMATIC_JSON) -m $$(TVC_MOTIF) \
	-t $$(TVC_ROOT_DIR) --primer-trim-bed $$(PRIMER_TRIM_BED) -g $$(basename $$(notdir $$<))")
endef
$(foreach pair,$(SAMPLE_PAIRS), \
	$(eval $(call tvc-somatic-vcf,$(tumor.$(pair)),$(normal.$(pair)))))



#define tvc-somatic-sufam-get-positions
#tvc/sufam/$1/sufampos.vcf : $2
#	$$(call LSCRIPT_MEM,22G,03:59:59,"$$(LOAD_JAVA8_MODULE); $$(call COMBINE_VARIANTS,21G) \
#		$$(foreach vcf,$$^,--variant $$(vcf) ) -o $$@ --genotypemergeoption UNSORTED -R $$(REF_FASTA)")
#endif
#$(foreach set,$(SAMPLE_SETS),\
#	$(eval $(call tvc-somatic-sufam-get-positions,$(set),\
#		$(foreach tumor,$(tumor.$(set)),vcf/$(tumor)_$(normal.$(set)).$(SOMATIC_SUFFIX).vcf))))

#define tvc-somatic-sufam
#tvc/sufam/$2/$1.vcf : $1.bam tvc/sufam/$2/sufampos.vcf
#	$$(call LSCRIPT_PARALLEL_MEM,4,10G,05:59:59,"$$(TVC) -s $$(<<) -i $$(<<) -r $$(REF_FASTA) -o $$(@D) -N 4 \
#		-m $$(TVC_MOTIF) -t $$(TVC_ROOT_DIR) --primer-trim-bed $$(PRIMER_TRIM_BED)")
#endef
#$(foreach pair,$(SAMPLE_PAIRS),\
#	$(eval $(call tvc-somatic-sufam,$(tumor.$(pair)),$(normal.$(pair)))))
#$(foreach pair,$(SAMPLE_PAIRS),\
#	$(eval $(call tvc-somatic-sufam,$(normal.$(pair)),$(normal.$(pair)))))

#define tvc-somatic-sufam-merge
#tvc/sufam/$1.sufam.vcf : $(foreach tumor,$2,tvc/sufam/$1/$2.vcf) tvc/sufam/$1.vcf
#	$$(call LSCRIPT_MEM,22G,03:59:59,"$$(LOAD_JAVA8_MODULE); $$(call COMBINE_VARIANTS,21G) \
#		$$(foreach vcf,$$^,--variant $$(vcf) ) -o $$@ --genotypemergeoption UNSORTED -R $$(REF_FASTA)")
#endef
#$(foreach set,$(SAMPLE_SETS),\
#	$(eval $(call tvc-somatic-sufam-merge,$(normal.$(set)),$(foreach tumor,$(tumor.$(set)),$(tumor))))

include usb-modules/vcf_tools/vcftools.mk
include usb-modules/variant_callers/TVC.mk

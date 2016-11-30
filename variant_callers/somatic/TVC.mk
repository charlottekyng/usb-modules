include usb-modules/Makefile.inc
include usb-modules/config.inc
include usb-modules/variant_callers/somatic/somaticVariantCaller.inc

LOGDIR ?= log/tvc_somtic.$(NOW)

$(info TARGETS_FILE_INTERVALS $(TARGETS_FILE_INTERVALS))
PHONY += tvc_somatic tvc_somatic_vcfs tvc_somatic_tables tvc_somatic_vcfs_sets tvc_somatic_tables_sets
tvc_somatic : tvc_somatic_vcfs tvc_somatic_tables tvc_somatic_vcfs_sets tvc_somatic_tables_sets

VARIANT_TYPES ?= tvc_snps tvc_indels
tvc_somatic_vcfs : $(foreach type,$(VARIANT_TYPES),$(call SOMATIC_VCFS,$(type)) $(addsuffix .idx,$(call SOMATIC_VCFS,$(type))))
tvc_somatic_tables : $(foreach type,$(VARIANT_TYPES),$(call SOMATIC_TABLES,$(type)))

VARIANT_TYPES_SETS ?= tvc_snps_sufam tvc_indels_sufam
tvc_somatic_vcfs_sets : $(foreach type,$(VARIANT_TYPES_SETS),$(call SOMATIC_VCFS_SETS,$(type)) $(addsuffix .idx,$(call SOMATIC_VCFS_SETS,$(type))))
tvc_somatic_tables_sets : $(foreach type,$(VARIANT_TYPES_SETS),$(call SOMATIC_TABLES_SETS,$(type)))

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



define tvc-somatic-sufam-get-positions
tvc/sufam/$1/sufampos.vcf : $2
	$$(call LSCRIPT_MEM,22G,03:59:59,"$$(LOAD_JAVA8_MODULE); $$(call COMBINE_VARIANTS,21G) \
		$$(foreach vcf,$$^,--variant $$(vcf) ) -o $$@ --genotypemergeoption UNSORTED -R $$(REF_FASTA)")
endef
$(foreach set,$(SAMPLE_SETS),\
	$(eval $(call tvc-somatic-sufam-get-positions,$(set),\
		$(foreach tumor,$(tumor.$(set)),vcf/$(tumor)_$(normal.$(set)).$(SOMATIC_SUFFIX).vcf) )) )

define tvc-somatic-sufam
tvc/sufam/$1/$2.vcf : $2.bam tvc/sufam/$1/sufampos.vcf
	$$(call LSCRIPT_PARALLEL_MEM,4,10G,05:59:59,"$$(TVC) -s $$(<<) -i $$(<<) -r $$(REF_FASTA) -o $$(@D) -N 4 \
		-m $$(TVC_MOTIF) -t $$(TVC_ROOT_DIR) --primer-trim-bed $$(PRIMER_TRIM_BED)")
endef
$(foreach set,$(SAMPLE_SETS),\
	$(foreach sample,$(subst _,$( ),$(set)),\
		$(eval $(call tvc-somatic-sufam-merge,$(set),$(sample) )) ) )

$(info INFO $(foreach sample,$(subst _,$( ),pt4PL_pt4RABx_pt4U_pt4BCnormal),tvc/sufam/pt4PL_pt4RABx_pt4U_pt4BCnormal/$(sample).tvc_snps_sufam.vcf))
define tvc-somatic-sufam-merge
tvc/sufam/$1.tvc_snps_sufam.vcf : $(foreach sample,$(subst _,$( ),$1),tvc/sufam/$1/$(sample).tvc_snps_sufam.vcf)
	$$(call LSCRIPT_MEM,22G,03:59:59,"$$(LOAD_JAVA8_MODULE); $$(call COMBINE_VARIANTS,21G) \
		$$(foreach vcf,$$^,--variant $$(vcf) ) -o $$@ --genotypemergeoption UNSORTED -R $$(REF_FASTA))")

tvc/sufam/$1.tvc_indels_sufam.vcf : $(foreach sample,$(subst _,$( ),$1),tvc/sufam/$1/$(sample).tvc_indels_sufam.vcf)
	$$(call LSCRIPT_MEM,22G,03:59:59,"$$(LOAD_JAVA8_MODULE); $$(call COMBINE_VARIANTS,21G) \
		$$(foreach vcf,$$^,--variant $$(vcf) ) -o $$@ --genotypemergeoption UNSORTED -R $$(REF_FASTA))")
endef
$(foreach set,$(SAMPLE_SETS),\
	$(eval $(call tvc-somatic-sufam-merge,$(set) )) )


include usb-modules/vcf_tools/vcftools.mk
include usb-modules/variant_callers/TVC.mk

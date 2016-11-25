include usb-modules/Makefile.inc
include usb-modules/config.inc

ifeq ($(findstring true,$(SOMATIC_ANN)),true)
include usb-modules/variant_callers/somatic/somaticVariantCaller.inc
else
include usb-modules/variant_callers/variantCaller.inc
endif


LOGDIR ?= log/ann_ext_vcf.$(NOW)

EXT_NAMES ?=

SOMATIC_FILTERS := 


ann_ext_vcf : ext_vcfs 
#ext_tables

ifeq ($(findstring true,$(SOMATIC_ANN)),true)
ext_vcfs : $(foreach ext_name,$(EXT_NAMES),$(call SOMATIC_VCFS,$(ext_name)))
ext_tables : $(foreach ext_name,$(EXT_NAMES),$(call SOMATIC_TABLES,$(ext_name)))
else
ext_vcfs : $(foreach ext_name,$(EXT_NAMES),$(call VCFS,$(EXT_NAME)))
ext_tables : $(foreach ext_name,$(EXT_NAMES),$(call TABLES,$(EXT_NAME)))
endif

#$(info $(EXT_NAMES))

.DELETE_ON_ERROR:
.SECONDARY:
.PHONY : ann_ext_vcf

include usb-modules/vcf_tools/vcftools.mk

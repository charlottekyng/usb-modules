include usb-modules/Makefile.inc
include usb-modules/config.inc

ifeq ($(findstring true,$(SOMATIC_ANN)),true)
include usb-modules/variant_callers/somatic/somaticVariantCaller.inc
else
include usb-modules/variant_callers/variantCaller.inc
endif


LOGDIR ?= log/ann_ext_vcf.$(NOW)

EXT_NAME ?=

SOMATIC_FILTERS := 


ann_ext_vcf : ext_vcfs
#ext_tables

ifeq ($(findstring true,$(SOMATIC_ANN)),true)
ext_vcfs : $(call SOMATIC_VCFS,$(EXT_NAME))
ext_tables : $(call SOMATIC_TABLES,$(EXT_NAME))
else
ext_vcfs : $(call VCFS,$(EXT_NAME))
ext_tables : $(call TABLES,$(EXT_NAME))
endif

$(info $(EXT_NAME))

.DELETE_ON_ERROR:
.SECONDARY:
.PHONY : ann_ext_vcf

include usb-modules/vcf_tools/vcftools.mk

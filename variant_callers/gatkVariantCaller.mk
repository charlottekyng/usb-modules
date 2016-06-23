# GATK variant caller module
# Author: Raymond Lim <raylim@mm.st> & Fong Chun Chan <fongchunchan@gmail.com>
# 

LOGDIR ?= log/gatk.$(NOW)

include usb-modules/Makefile.inc
include usb-modules/config.inc
include usb-modules/variant_callers/variantCaller.inc

VPATH ?= bam

.DELETE_ON_ERROR:

.SECONDARY: 

VARIANT_TYPES = gatk_snps gatk_indels

PHONY += gatk # gatk_vcfs gatk_tables gatk_reports

gatk : gatk_vcfs gatk_tables # reports

gatk_vcfs : $(foreach type,$(VARIANT_TYPES),$(call VCFS,$(type)) $(addsuffix .idx,$(call VCFS,$(type))))

gatk_tables : $(foreach type,$(VARIANT_TYPES),$(call TABLES,$(type)))

gatk_reports : $(foreach type,gatk_indels gatk_snps,reports/$(type).dp_ft.grp)

include usb-modules/variant_callers/gatk.mk

.PHONY : $(PHONY)

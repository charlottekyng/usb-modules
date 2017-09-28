include usb-modules/Makefile.inc
include usb-modules/config.inc

LOGDIR ?= log/rsem.$(NOW)

.DELETE_ON_ERROR:
.SECONDARY:
.PHONY: rsem

rsem : $(foreach sample,$(SAMPLES),rsem/$(sample).genes.results)

define rsem-calc-expression
rsem/$1.genes.results : star/secondpass/$1.Aligned.toTranscriptome.out.bam 
	$$(call LSCRIPT_PARALLEL_MEM,8,$$(RESOURCE_REQ_LOWMEM),$$(RESOURCE_REQ_SHORT),"$$(LOAD_RSEM_MODULE); \
	$$(RSEM_CALC_EXPR) $$(RSEM_OPTIONS) $$(if $$(findstring true,$$(PAIRED_END)),--paired-end) \
	$$< $$(RSEM_INDEX) $$(@D)/$1")
endef
$(foreach sample,$(SAMPLES),\
	$(eval $(call rsem-calc-expression,$(sample))))

include usb-modules/bam_tools/processBam.mk
include usb-modules/aligners/align.mk
include usb-modules/aligners/starAligner.mk


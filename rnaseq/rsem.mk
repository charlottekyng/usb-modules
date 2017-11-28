include usb-modules/Makefile.inc
include usb-modules/config.inc

LOGDIR ?= log/rsem.$(NOW)

.DELETE_ON_ERROR:
.SECONDARY:
.PHONY: rsem

rsem : $(foreach sample,$(SAMPLES),rsem/$(sample).genes.results) rsem/all.genes.expected_count.results rsem/all.genes.TPM.results rsem/all.genes.FPKM.results rsem/all.isoforms.expected_count.results rsem/all.isoforms.TPM.results rsem/all.isoforms.FPKM.results

define rsem-calc-expression
rsem/$1.genes.results : star/secondpass/$1.Aligned.toTranscriptome.out.bam 
	$$(call LSCRIPT_PARALLEL_MEM,4,$$(RESOURCE_REQ_LOWMEM),$$(RESOURCE_REQ_SHORT),"$$(LOAD_RSEM_MODULE); \
	$$(RSEM_CALC_EXPR) $$(RSEM_OPTIONS) $$(if $$(findstring true,$$(PAIRED_END)),--paired-end) \
	$$< $$(RSEM_INDEX) $$(@D)/$1")
endef
$(foreach sample,$(SAMPLES),\
	$(eval $(call rsem-calc-expression,$(sample))))

rsem/%.isoforms.results : rsem/%.genes.results

rsem/all.genes.expected_count.results : $(foreach sample,$(SAMPLES),rsem/$(sample).genes.results)
	$(call LSCRIPT_MEM,$(RESOURCE_REQ_LOWMEM),$(RESOURCE_REQ_VSHORT),"$(RSEM_GEN_DATA_MATRIX) expected_count $^ > $@")

rsem/all.genes.TPM.results : $(foreach sample,$(SAMPLES),rsem/$(sample).genes.results)
	$(call LSCRIPT_MEM,$(RESOURCE_REQ_LOWMEM),$(RESOURCE_REQ_VSHORT),"$(RSEM_GEN_DATA_MATRIX) TPM $^ > $@")

rsem/all.genes.FPKM.results : $(foreach sample,$(SAMPLES),rsem/$(sample).genes.results)
	$(call LSCRIPT_MEM,$(RESOURCE_REQ_LOWMEM),$(RESOURCE_REQ_VSHORT),"$(RSEM_GEN_DATA_MATRIX) FPKM $^ > $@")

rsem/all.isoforms.expected_count.results : $(foreach sample,$(SAMPLES),rsem/$(sample).isoforms.results)
	$(call LSCRIPT_MEM,$(RESOURCE_REQ_MEDIUM_MEM),$(RESOURCE_REQ_VSHORT),"$(RSEM_GEN_DATA_MATRIX) expected_count $^ > $@")

rsem/all.isoforms.TPM.results : $(foreach sample,$(SAMPLES),rsem/$(sample).isoforms.results)
	$(call LSCRIPT_MEM,$(RESOURCE_REQ_MEDIUM_MEM),$(RESOURCE_REQ_VSHORT),"$(RSEM_GEN_DATA_MATRIX) TPM $^ > $@")

rsem/all.isoforms.FPKM.results : $(foreach sample,$(SAMPLES),rsem/$(sample).isoforms.results)
	$(call LSCRIPT_MEM,$(RESOURCE_REQ_MEDIUM_MEM),$(RESOURCE_REQ_VSHORT),"$(RSEM_GEN_DATA_MATRIX) FPKM $^ > $@")

include usb-modules/bam_tools/processBam.mk
include usb-modules/aligners/align.mk
include usb-modules/aligners/starAligner.mk


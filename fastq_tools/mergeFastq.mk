# add in new fastq to existing bams
include usb-modules/Makefile.inc

LOGDIR ?= log/merge_fastq.$(NOW)

.PHONY: merge_fastq
.DELETE_ON_ERROR:
.SECONDARY:

MERGE_SAMPLE_FILE ?= merge_samples.txt
ifneq ($(wildcard $(MERGE_SAMPLE_FILE)),)
  MERGE_SAMPLES ?= $(shell sed '/^\#/d' $(MERGE_SAMPLE_FILE))
endif

ALIGNER ?= bwamem

merge_fastq : $(foreach sample,$(MERGE_SAMPLES),$(if $(wildcard bam/$(sample).bam),merged_bam/$(sample).bam,bam/$(sample).bam))
	$(call LSCRIPT_MEM,$(RESOURCE_REQ_MEDIUM_MEM),$(RESOURCE_REQ_SHORT),"$(LOAD_SAMTOOLS_MODULE); \
		for bam in $(filter merged_bam/%.bam,$^); do \
		ln -f \$${bam} bam/\$$(basename \$${bam}) && \
		$(SAMTOOLS) index \$${bam}; \
		done")

include usb-modules/aligners/$(ALIGNER)Aligner.mk

merged_bam/%.1.bam merged_bam/%.2.bam : $(ALIGNER)/bam/%.$(ALIGNER).$(BAM_SUFFIX)
	$(INIT)  ln -f $(<M) merged_bam/$*.1.bam && \
	ln -f bam/$*.bam merged_bam/$*.2.bam

merged_bam/%.header.sam : merged_bam/%.1.bam merged_bam/%.2.bam
	$(INIT) { $(LOAD_SAMTOOLS_MODULE); $(SAMTOOLS) view -H $(<M) | grep -v '^@RG'; \
	for bam in $(^); do \
	$(SAMTOOLS) view -H $$bam | grep '^@RG'; \
	done | sort | uniq; } > $@

merged_bam/%.bam : merged_bam/%.header.sam merged_bam/%.1.bam merged_bam/%.2.bam
	$(call LSCRIPT_MEM,$(RESOURCE_REQ_HIGHMEM),$(RESOURCE_REQ_SHORT),"$(LOAD_SAMTOOLS_MODULE); \
	$(SAMTOOLS) merge -f -h $< $(@) $(filter %.bam,$(^))")


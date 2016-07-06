# This module extract fastq files from a bam file.  It will use either the Picard (SamToFastq.jar) or bam2fastq programs to extract the fastq.  You can specify which program to use with the EXTRACT_TOOL variable
# input: $(SAMPLES)
# Author: Fong Chun Chan <fongchunchan@gmail.com>

include usb-modules/Makefile.inc
include usb-modules/config.inc

LOGDIR ?= log/extract_fastq.$(NOW)

.DELETE_ON_ERROR:
.SECONDARY:
.PHONY: extract_fastq

extract_fastq : $(foreach sample,$(SAMPLES),fastq/$(sample).1.fastq.gz)

ifeq (${EXTRACT_TOOL},PICARD)
fastq/%.1.fastq.gz fastq/%.2.fastq.gz : bam/%.bam
	$(call LSCRIPT_MEM,10G,00:59:59,"$(LOAD_JAVA8_MODULE); $(call SAM_TO_FASTQ,9G) \
		I=$< FASTQ=>(gzip -c > fastq/$*.1.fastq.gz) SECOND_END_FASTQ=>(gzip -c > fastq/$*.2.fastq.gz)")
else
fastq/%.1.fastq.gz fastq/%.2.fastq.gz : bam/%.bam
	$(call LSCRIPT_PARALLEL_MEM,4,4G,5G,"$(BAM_TO_FASTQ) -i <($(SAMTOOLS2) sort -T bam/$* -O bam -n -@ 4 -m 4G $<) \
		-fq >(gzip -c > fastq/$*.1.fastq.gz) -fq2 >(gzip -c > fastq/$*.2.fastq.gz)")
endif

include usb-modules/Makefile.inc
include usb-modules/config.inc

LOGDIR ?= log/poolednorm_bam.$(NOW)

poolednorm_bam : bam/poolednorm.bam $(addsuffix .bai,bam/poolednorm.bam)


unprocessed_bam/poolednorm.bam : $(foreach normal,$(POOLED_NORM_SAMPLES),bam/$(normal).downsampled.bam) $(addsuffix .bai,$(foreach normal,$(POOLED_NORM_SAMPLES),bam/$(normal).downsampled.bam))
	$(MKDIR) unprocessed_bam; $(call LSCRIPT_MEM,$(RESOURCE_REQ_HIGHMEM),$(RESOURCE_REQ_MEDIUM),"$(LOAD_SAMTOOLS_MODULE); \
	$(SAMTOOLS) merge -f $@ $(filter %.bam,$^) && $(RM) $^")

include usb-modules/bam_tools/processBam.mk
include usb-modules/bam_tools/fixRG.mk




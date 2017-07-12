# various bam processing steps
##### MAKE INCLUDES #####
include usb-modules/Makefile.inc

LOGDIR ?= log/endtoendbam.$(NOW)

BAMS = $(foreach sample,$(SAMPLES),bam/$(sample).bam)
fixed_bams : $(BAMS) $(addsuffix .bai,$(BAMS))

bam/%.bam : unprocessed_bam/%.endtoend.bam
	$(INIT) ln -f $(<) $(@)


include usb-modules/bam_tools/processBam.mk

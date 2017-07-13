# This module fixes the mate information of a bam file to prepare it for proper extraction
# input: $(SAMPLES)
# Author: Fong Chun Chan <fongchunchan@gmail.com>
# CN: 2170713 absolutely does not work

include usb-modules/Makefile.inc
include usb-modules/config.inc

SAMPLE_FILE = samplesToFixMate.txt
SAMPLES = $(shell cat $(SAMPLE_FILE))

LOGDIR = gsc_bam/logs

.DELETE_ON_ERROR:

.SECONDARY:

all : $(foreach sample,$(SAMPLES),gsc_bam/$(sample).fixmate.bam)

gsc_bam/%.fixmate.bam : gsc_bam/%.bam
	SGE_RREQ="$(SGE_RREQ) $(call MEM_FREE,10G,40G)" $(MKDIR) $(LOGDIR);\
	$(FIXMATE) I=$< O=$@ &> ${LOGDIR}/$(@F).fixmate.log

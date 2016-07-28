include usb-modules/Makefile.inc
include usb-modules/config.inc

LOGDIR ?= log/pon_bam.$(NOW)

pon_bam : bam/pon.bam $(addsuffix .bai,bam/pon.bam)


bam/pon.header.sam : $(foreach normal,$(NORMAL_SAMPLES),bam/$(normal).downsampled.bam)
	$(INIT) $(LOAD_SAMTOOLS_MODULE); $(SAMTOOLS) view -H $< | grep -v '^@RG' > $@.tmp; \
	for bam in $^; do $(SAMTOOLS) view -H $$bam | grep '^@RG' >> $@.tmp; done; \
	uniq $@.tmp > $@ && $(RM) $@.tmp

bam/pon.bam : bam/pon.header.sam $(foreach normal,$(NORMAL_SAMPLES),bam/$(normal).downsampled.bam) $(addsuffix .bai,$(foreach normal,$(NORMAL_SAMPLES),bam/$(normal).downsampled.bam))
	$(call LSCRIPT_MEM,20G,02:29:29,"$(LOAD_SAMTOOLS_MODULE); $(SAMTOOLS) merge -f -h $< $@ $(filter %.bam,$^) && $(RM) $^")


include usb-modules/bam_tools/processBam.mk





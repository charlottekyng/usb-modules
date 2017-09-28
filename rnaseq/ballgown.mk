include usb-modules/Makefile.inc
include usb-modules/config.inc

LOGDIR ?= log/ballgown.$(NOW)

.PHONY: ballgown
.DELETE_ON_ERROR:

ballgown : $(foreach sample,$(SAMPLES),stringtie/$(sample).ballgown.gtf)

stringtie/%.stringtie.gtf : hisat2/bam/%.hisat2.bam
	$(call LSCRIPT_PARALLEL_MEM,8,$(RESOURCE_REQ_MEDIUM_MEM),$(RESOURCE_REQ_SHORT),"$(LOAD_STRINGTIE_MODULE); \
	$(STRINGTIE) $(STRINGTIE_OPTS) -G $(GENCODE_GTF) -o $@ -l $* $<")

stringtie/all.stringtie_merged.gtf : $(foreach sample,$(SAMPLES),stringtie/$(sample).stringtie.gtf)
	ls $^ > stringtie/all.stringtie_merge_list.txt; \	
	$(call LSCRIPT_PARALLEL_MEM,8,$(RESOURCE_REQ_MEDIUM_MEM),$(RESOURCE_REQ_SHORT),"$(LOAD_STRINGTIE_MODULE); \
	$(STRINGTIE) --merge $(STRINGTIE_OPTS) -G $(GENCODE_GTF) -o $@ stringtie/all.stringtie_merge_list.txt")

stringtie/%.ballgown.gtf : stringtie/all.stringtie_merged.gtf hisat2/bam/%.hisat2.bam
	$(call LSCRIPT_PARALLEL_MEM,8,$(RESOURCE_REQ_LOW_MEM),$(RESOURCE_REQ_SHORT),"$(LOAD_STRINGTIE_MODULE); \
	$(STRINGTIE) -e -B $(STRINGTIE_OPTS) -G $< -o $@ $(<<)")

include usb-modules/aligners/hisat2Aligner.mk
include usb-modules/bam_tools/processBam.mk

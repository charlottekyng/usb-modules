# vim: set ft=make :
# Run Fastqc on bam files

include usb-modules/Makefile.inc
include usb-modules/config.inc

LOGDIR ?= log/fastqc.$(NOW)

.PHONY: fastqc
.SECONDARY: 

#fastqc : $(foreach sample,$(SAMPLES),fastqc/$(sample)_fastqc/summary.txt) fastqc/all_summary.txt
fastqc : $(foreach sample,$(SAMPLES),fastqc/$(sample)_fastqc.zip)
fastqc/%_fastqc.zip : bam/%.bam
	$(call LSCRIPT_NAMED_MEM,$*_fastqc,$(RESOURCE_REQ_MEDIUM_MEM),$(RESOURCE_REQ_SHORT),"$(LOAD_FASTQC_MODULE); \
	$(FASTQC) -o fastqc $^")

fastqc/%_fastqc/summary.txt : fastqc/%_fastqc.zip
	$(INIT) $(UNZIP) -o -d fastqc $< &> $(LOG) && touch $@

fastqc/all_summary.txt : $(foreach sample,$(SAMPLES),fastqc/$(sample)_fastqc/summary.txt)
	$(INIT) $(LOAD_R_MODULE); $(FASTQC_SUMMARY_PLOT) --outPrefix fastqc/all_summary $^ &> $(LOG)

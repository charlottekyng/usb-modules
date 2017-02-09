# vim: set ft=make :
# BWA-mem alignment of short reads
# OPTIONS: NO_MARKDUP = true/false (default: false)
# 		   EXTRACT_FASTQ = true/false (default: false)
# 		   BAM_NO_RECAL = true/false (default: false)
ALIGNER := bwamem

include usb-modules/Makefile.inc
include usb-modules/config.inc
include usb-modules/aligners/align.inc

LOGDIR ?= log/bwamem.$(NOW)

VPATH ?= unprocessed_bam

..DUMMY := $(shell mkdir -p version; $(BWA) &> version/bwamem.txt; echo "options: $(BWA_ALN_OPTS)" >> version/bwamem.txt )
.SECONDARY:
.DELETE_ON_ERROR: 
.PHONY: bwamem

BWA_BAMS = $(foreach sample,$(SAMPLES),bam/$(sample).bam)

bwamem : $(BWA_BAMS) $(addsuffix .bai,$(BWA_BAMS))

$(info BAM_SUFFIX-BWA is $(BAM_SUFFIX))

bam/%.bam : bwamem/bam/%.bwamem.$(BAM_SUFFIX)
	$(call LSCRIPT," ln -f $(<) $(@)")

#$(call align-split-fastq,name,split-name,fastqs)
define align-split-fastq
bwamem/bam/$2.bwamem.bam : $3
	$$(call LSCRIPT_PARALLEL_MEM,$$(BWA_NUM_CORES),4G,01:29:59,"$$(LOAD_BWA_MODULE); $$(LOAD_SAMTOOLS_MODULE); \
		$$(BWA_MEM) -t $$(BWA_NUM_CORES) $$(BWA_ALN_OPTS) \
		-R \"@RG\tID:$2\tLB:$1\tPL:$${SEQ_PLATFORM}\tSM:$1\" $$(REF_FASTA) $$^ | $$(SAMTOOLS) view -bhS - > $$@")
endef
$(foreach ss,$(SPLIT_SAMPLES),$(if $(fq.$(ss)),$(eval $(call align-split-fastq,$(split.$(ss)),$(ss),$(fq.$(ss))))))

bwamem/bam/%.bwamem.bam : fastq/%.1.fastq.gz fastq/%.2.fastq.gz
	LBID=`echo "$*" | sed 's/_[A-Za-z0-9]\+//'`; \
	$(call LSCRIPT_PARALLEL_MEM,$(BWA_NUM_CORES),4G,01:29:59,"$(LOAD_BWA_MODULE); $(LOAD_SAMTOOLS_MODULE); \
		$(BWA_MEM) -t $(BWA_NUM_CORES) $(BWA_ALN_OPTS) \
		-R \"@RG\tID:$*\tLB:$${LBID}\tPL:${SEQ_PLATFORM}\tSM:$${LBID}\" $(REF_FASTA) $^ | \
		$(SAMTOOLS) view -bhS - > $@")

bwamem/bam/%.bwamem.bam : fastq/%.fastq.gz
	LBID=`echo "$*" | sed 's/_[A-Za-z0-9]\+//'`; \
	$(call LSCRIPT_PARALLEL_MEM,$(BWA_NUM_CORES),4G,01:29:59,"$(LOAD_BWA_MODULE); $(LOAD_SAMTOOLS_MODULE); \
		$(BWA_MEM) -t $(BWA_NUM_CORES) $(BWA_ALN_OPTS) \
		-R \"@RG\tID:$*\tLB:$${LBID}\tPL:${SEQ_PLATFORM}\tSM:$${LBID}\" $(REF_FASTA) $^ | \
		$(SAMTOOLS) view -bhS - > $@")

fastq/%.fastq.gz : fastq/%.fastq
	$(call LSCRIPT,"gzip -c $< > $(@) && $(RM) $<")

include usb-modules/fastq_tools/fastq.mk
include usb-modules/bam_tools/processBam.mk
include usb-modules/aligners/align.mk

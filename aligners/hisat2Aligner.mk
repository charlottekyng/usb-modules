# This module is for the hisat aligner
# input: $(SAMPLES) 
include usb-modules/Makefile.inc
include usb-modules/config.inc

BAM_NO_REALN = true
BAM_NO_RECAL = true
BAM_NO_FILTER = true
BAM_DUP_TYPE = none

ALIGNER := hisat2

LOGDIR = log/hisat2.$(NOW)

..DUMMY := $(shell mkdir -p version; $(HISAT2) --version &> version/hisat2.txt; echo "options: $(HISAT2_OPTS)" >> version/hisat2.txt)
.SECONDARY:
.DELETE_ON_ERROR:
.PHONY : hisat_bams
	
BAMS = $(foreach sample,$(SAMPLES),bam/$(sample).bam)
hisat2_bams : $(BAMS) $(addsuffix .bai,$(BAMS))

bam/%.bam : hisat2/bam/%.hisat2.$(BAM_SUFFIX)
	$(INIT) ln -f $(<) $(@) 

hisat2/bam/%.hisat2.bam hisat2/unmapped_bam/%.hisat2_unmapped.bam : fastq/%.1.fastq.gz fastq/%.2.fastq.gz
	$(call LSCRIPT_PARALLEL_MEM,8,$(RESOURCE_REQ_MEDIUM_MEM),$(RESOURCE_REQ_SHORT),"$(MKDIR) hisat2/bam hisat2/unmapped_bam; \
		$(LOAD_SAMTOOLS_MODULE); $(LOAD_HISAT2_MODULE); \
		LBID=`echo \"$*\" | sed 's/_.\+//'`; \
		$(HISAT2) $(HISAT2_OPTS) \
		--rg-id $* --rg \"SM:\$${LBID}\" \
		--rg PL:${SEQ_PLATFORM} --rg \"LB:\$${LBID}\" \
		-S hisat2/bam/$*.hisat.sam --un hisat2/unmapped_bam/$*.hisat_unmapped.sam \
		-1 $(<) -2 $(<<) && \
		$(SAMTOOLS) view -Sbh hisat2/bam/$*.hisat.sam > hisat2/bam/$*.hisat.bam && \
		$(SAMTOOLS) view -Sbh hisat2/unmapped_bam/$*.hisat_unmapped.sam > hisat2/unmapped_bam/$*.hisat_unmapped.bam && \
		$(RM) hisat2/bam/$*.hisat.sam hisat2/unmapped_bam/$*.hisat_unmapped.sam")

#define align-merged-split-fastq
#hisat/bam/$2.hisat.bam : $3
#	$$(call LSCRIPT_NAMED_PARALLEL_MEM,$2_hisat,8,1G,1.5G,"mkdir -p hisat/bam hisat/unmapped_bam; \
#		$$(HISAT) $$(HISAT_OPTS) \
#		--rg-id $2 --rg \"SM:$1\" \
#		--rg PL:$${SEQ_PLATFORM} --rg \"LB:\$1\" \
#		-S >($$(SAMTOOLS2) view -bh - > hisat/bam/$2.hisat.bam) \
#		--un >($$(SAMTOOLS2) view -bh - > hisat/unmapped_bam/$2.hisat_unmapped.bam) \
#		$$(if $$(<<),-1 $$(<) -2 $$(<<),-U $$(<))")
#hisat/unmapped_bam/$2.hisat_unmapped.bam : hisat/bam/$2.hisat.bam
#endef
#$(foreach ss,$(SPLIT_SAMPLES),\
#	$(if $(fq.$(ss)),\
#	$(eval $(call align-merged-split-fastq,$(split.$(ss)),$(ss),$(fq.$(ss))))))


include usb-modules/fastq_tools/fastq.mk
include usb-modules/bam_tools/processBam.mk
include usb-modules/aligners/align.mk

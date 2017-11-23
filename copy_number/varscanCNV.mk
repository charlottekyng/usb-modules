# Run VarScan to detect copynumber
# Detect copy number
##### DEFAULTS ######

LOGDIR = log/varscanCNV.$(NOW)

##### MAKE INCLUDES #####
include usb-modules/Makefile.inc
include usb-modules/config.inc

.DELETE_ON_ERROR:
.SECONDARY: 
.PHONY: all
#copycalls segments cnv

all : copynum copycalls segments geneCN 

copynum : $(foreach pair,$(SAMPLE_PAIRS),$(foreach pool,$(TARGETS_FILE_INTERVALS_POOLS),varscan/copynum/$(pair).$(notdir $(pool)).copynumber))
copycalls : $(foreach pair,$(SAMPLE_PAIRS),$(foreach pool,$(TARGETS_FILE_INTERVALS_POOLS),varscan/copycall/$(pair).$(notdir $(pool)).copycall))
segments : $(foreach pair,$(SAMPLE_PAIRS),varscan/segment/$(pair).segment.Rdata)
geneCN : varscan/segment/geneCN.txt


ifeq ($(findstring ILLUMINA,$(SEQ_PLATFORM)),ILLUMINA)
define varscan-copynum-tumor-normal
varscan/copynum/$1_$2.copynumber :  bam/$1.bam bam/$2.bam
	$$(call LSCRIPT_CHECK_MEM,$$(RESOURCE_REQ_HIGHMEM),$$(RESOURCE_REQ_SHORT),"$$(LOAD_SAMTOOLS_MODULE); $$(LOAD_JAVA7_MODULE); \
	$$(SAMTOOLS) mpileup $$(CBS_MPILEUP_OPTS) -l $$(TARGETS_FILE_INTERVALS) -f $$(REF_FASTA) $$(word 2,$$^) $$< | awk 'NF == 9 { print }' |  \
	$$(VARSCAN) copynumber - $$(basename $$@) --mpileup 1")
endef
$(foreach pair,$(SAMPLE_PAIRS),\
	$(eval $(call varscan-copynum-tumor-normal,$(tumor.$(pair)),$(normal.$(pair)))))
endif

ifeq ($(findstring IONTORRENT,$(SEQ_PLATFORM)),IONTORRENT)
define varscan-copynum-tumor-normal
varscan/copynum/$1_$2.$$(notdir $3).copynumber :  bam/$1.bam bam/$2.bam $3
	$$(call LSCRIPT_CHECK_MEM,$$(RESOURCE_REQ_MEDIUM_MEM),$$(RESOURCE_REQ_VSHORT),"$$(LOAD_SAMTOOLS_MODULE); $$(LOAD_JAVA7_MODULE); $$(LOAD_BEDTOOLS_MODULE); \
	TMP1=`mktemp`.bam && $$(BEDTOOLS) intersect -wa -F 1 -a $$< -b $$(word 3, $$^) > \$$$$TMP1 && $$(SAMTOOLS) index \$$$$TMP1 && \
	TMP2=`mktemp`.bam && $$(BEDTOOLS) intersect -wa -F 1 -a $$(word 2,$$^) -b $$(word 3, $$^) > \$$$$TMP2 && $$(SAMTOOLS) index \$$$$TMP2 && \
	$$(SAMTOOLS) mpileup $$(CBS_MPILEUP_OPTS) -l $$(word 3,$$^) -f $$(REF_FASTA) \$$$$TMP2 \$$$$TMP1 | awk 'NF == 9 { print }' |  \
	$$(VARSCAN) copynumber - $$(basename $$@) --mpileup 1 --max-segment-size $$(VARSCAN_CNV_MAX_SEG_SIZE) && $$(RM) \$$$$TMP1 \$$$$TMP2")
endef
$(foreach pair,$(SAMPLE_PAIRS),\
	$(foreach pool,$(TARGETS_FILE_INTERVALS_POOLS),\
		$(eval $(call varscan-copynum-tumor-normal,$(tumor.$(pair)),$(normal.$(pair)),$(pool)))))
endif

varscan/copycall/%.copycall : varscan/copynum/%.copynumber
	$(call LSCRIPT_CHECK_MEM,$(RESOURCE_REQ_MEDIUM_MEM),$(RESOURCE_REQ_SHORT),"n=`awk '{ total += \$$7 } END { printf \"\%.6f\"$(,) total / NR }' $<`; \
	if [ \$$(bc <<< \"\$$n > 0\") -eq 1 ]; then \
		recenter_opt=\"--recenter-up \$$n\"; \
	else \
		n=\$$(bc <<< \"\$$n*-1\"); \
		recenter_opt=\"--recenter-down \$$n\"; \
	fi; \
	$(VARSCAN) copyCaller $< --output-file $@ \$$recenter_opt")

ifeq ($(findstring ILLUMINA,$(SEQ_PLATFORM)),ILLUMINA)
varscan/segment/%.segment.Rdata : varscan/copycall/%.copycall
	$(call LSCRIPT_CHECK_MEM,$(RESOURCE_REQ_LOWMEM),$(RESOURCE_REQ_VSHORT),"$(LOAD_R_MODULE); \
	$(CBS_SEGMENTCNV) --alpha $(CBS_SEG_ALPHA) --smoothRegion $(CBS_SEG_SMOOTH) \
	--undoSD $(CBS_SEG_SD) $(if $(CENTROMERE_TABLE),--centromereFile=$(CENTROMERE_TABLE)) --prefix=$(@D)/$* $^")
endif

ifeq ($(findstring IONTORRENT,$(SEQ_PLATFORM)),IONTORRENT)
define varscan-segment
varscan/segment/$1_$2.segment.Rdata : $$(foreach pool,$$(TARGETS_FILE_INTERVALS_POOLS),varscan/copycall/$1_$2.$$(notdir $$(pool)).copycall)
	$$(call LSCRIPT_CHECK_MEM,$$(RESOURCE_REQ_LOWMEM),$$(RESOURCE_REQ_VSHORT),"$$(LOAD_R_MODULE); \
	$$(CBS_SEGMENTCNV) --alpha $$(CBS_SEG_ALPHA) --smoothRegion $$(CBS_SEG_SMOOTH) --trim $$(CBS_TRIM) --clen $$(CBS_CLEN) \
	--undoSD $$(CBS_SEG_SD) $$(if $$(CENTROMERE_TABLE),--centromereFile=$$(CENTROMERE_TABLE)) --prefix=$$(@D)/$1_$2 $$^")
endef
$(foreach pair,$(SAMPLE_PAIRS),\
	$(eval $(call varscan-segment,$(tumor.$(pair)),$(normal.$(pair)))))
endif

varscan/segment/%.collapsed_seg.txt : varscan/segment/%.segment.Rdata
	

varscan/segment/geneCN.txt : $(foreach pair,$(SAMPLE_PAIRS),varscan/segment/$(pair).collapsed_seg.txt)
	$(call LSCRIPT_CHECK_MEM,$(RESOURCE_REQ_LOWMEM),$(RESOURCE_REQ_VSHORT),"$(LOAD_R_MODULE); \
	$(VARSCAN_GENE_CN) $(VARSCAN_GENE_CN_OPTS) --outFile $@ $^")	

define varscan-segment-sd-alpha-smooth
varscan/segment_sd$1_alpha$2_smooth$3/%.segment.Rdata : varscan/copycall/%.copycall
	$$(call LSCRIPT_CHECK_NAMED_MEM,$$*_seg_$1_$2_$3,$(RESOURCE_REQ_LOWMEM),$(RESOURCE_REQ_VSHORT),"$$(LOAD_R_MODULE); \
	$$(CBS_SEGMENTCNV) --undoSD $1 --alpha $2 --smoothRegion $3 \
	$$(if $$(CENTROMERE_TABLE),--centromereFile=$$(CENTROMERE_TABLE)) --prefix=$$(@D)/$$* $$<")
endef
$(foreach sd,$(SEG_SDS),\
	$(foreach alpha,$(SEG_ALPHAS),\
		$(foreach smooth,$(SEG_SMOOTHS),\
			$(eval $(call varscan-segment-sd-alpha-smooth,$(sd),$(alpha),$(smooth))))))


# Run VarScan to detect copynumber
# Detect copy number
##### DEFAULTS ######

LOGDIR = log/varscan.$(NOW)

##### MAKE INCLUDES #####
include usb-modules/Makefile.inc

.DELETE_ON_ERROR:
.SECONDARY: 
.PHONY: all copycalls segments cnv

all : copycalls segments # geneCN 

CGHCALLS := $(foreach pair,$(SAMPLE_PAIRS),varscan/segment/$(pair).cgh_call.txt)
ifeq ($(MULTIPARAM_SEGMENT),true) 
CGHCALLS += $(foreach pair,$(SAMPLE_PAIRS),\
				$(foreach sd,$(SEG_SDS),\
					$(foreach alpha,$(SEG_ALPHAS),\
						$(foreach smooth,$(SEG_SMOOTHS),\
							varscan/segment_sd$(sd)_alpha$(alpha)_smooth$(smooth)/$(pair).cgh_call.txt))))
endif

segments : $(foreach pair,$(SAMPLE_PAIRS),varscan/segment/$(pair).segment.Rdata)
copycalls : $(foreach pair,$(SAMPLE_PAIRS),varscan/copycall/$(pair).copycall)
geneCN : varscan/segment/geneCN.txt

define varscan-copynum-tumor-normal
varscan/copynum/$1_$2.copynumber :  bam/$1.bam bam/$2.bam
	$$(call LSCRIPT_CHECK_MEM,9G,01:59:59,"$$(LOAD_SAMTOOLS_MODULE); $$(LOAD_JAVA7_MODULE); \
	$$(SAMTOOLS) mpileup $$(CBS_MPILEUP_OPTS) -f $$(REF_FASTA) $$(word 2,$$^) $$< | awk 'NF == 9 { print }' |  \
	$$(VARSCAN) copynumber - $$(basename $$@) --mpileup 1")
endef
$(foreach pair,$(SAMPLE_PAIRS),\
		$(eval $(call varscan-copynum-tumor-normal,$(tumor.$(pair)),$(normal.$(pair)))))

varscan/copycall/%.copycall : varscan/copynum/%.copynumber
	$(call LSCRIPT_CHECK_MEM,9G,01:59:59,"n=`awk '{ total += \$$7 } END { print total / NR }' $<`; \
	if [ \$$(bc <<< \"\$$n > 0\") -eq 1 ]; then \
		recenter_opt=\"--recenter-up \$$n\"; \
	else \
		n=\$$(bc <<< \"\$$n*-1\"); \
		recenter_opt=\"--recenter-down \$$n\"; \
	fi; \
	$(VARSCAN) copyCaller $< --output-file $@ \$$recenter_opt")

varscan/segment/%.segment.Rdata : varscan/copycall/%.copycall
	$(call LSCRIPT_CHECK_MEM,4G,00:29:29,"$(LOAD_R_MODULE); $(CBS_SEGMENTCNV) --alpha $(CBS_SEG_ALPHA) --smoothRegion $(CBS_SEG_SMOOTH) \
	--undoSD $(CBS_SEG_SD) --centromereFile=$(CENTROMERE_TABLE) --prefix=$(@D)/$* $<")

varscan/segment/geneCN.txt : $(foreach pair,$(SAMPLE_PAIRS),varscan/segment/$(pair).collapsed_seg.txt)
	$(call LSCRIPT_CHECK_MEM,4G,00:29:29,"$(LOAD_R_MODULE); $(VARSCAN_GENE_CN) $(VARSCAN_GENE_CN_OPTS) --outFile $@ $^")	

define varscan-segment-sd-alpha-smooth
varscan/segment_sd$1_alpha$2_smooth$3/%.segment.Rdata : varscan/copycall/%.copycall
	$$(call LSCRIPT_CHECK_NAMED_MEM,$$*_seg_$1_$2_$3,4G,00:29:29,"$$(LOAD_R_MODULE); $$(CBS_SEGMENTCNV) --undoSD $1 \
	--alpha $2 --smoothRegion $3 --centromereFile=$$(CENTROMERE_TABLE) --prefix=$$(@D)/$$* $$<")
endef
$(foreach sd,$(SEG_SDS),\
	$(foreach alpha,$(SEG_ALPHAS),\
		$(foreach smooth,$(SEG_SMOOTHS),\
			$(eval $(call varscan-segment-sd-alpha-smooth,$(sd),$(alpha),$(smooth))))))


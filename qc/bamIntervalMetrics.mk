# generate bam interval metrics per sample

#NO_RM := true

include usb-modules/Makefile.inc
include usb-modules/variant_callers/gatk.inc
# picard format intervals file, needs requires sam format header

VPATH ?= bam

LOGDIR ?= log/metrics.$(NOW)

EXOME ?= false

PLOT_HS_METRICS = $(RSCRIPT) usb-modules/qc/plotHsMetrics.R
NON_REF_FREQ = usb-modules/qc/nonRefFreqFromPileup.pl
NON_REF_FREQ_BIN_SIZE = 0.01

.DELETE_ON_ERROR:

.SECONDARY: 

.PHONY: bam_interval_metrics hs_metrics amplicon_metrics wgs_metrics interval_report non_ref_metrics

ifeq ($(CAPTURE_METHOD),NONE)
bam_interval_metrics : wgs_metrics non_ref_metrics
endif
ifeq ($(CAPTURE_METHOD),BAITS)
bam_interval_metrics : hs_metrics non_ref_metrics
endif
ifeq ($(CAPTURE_METHOD),PCR)
bam_interval_metrics : amplicon_metrics non_ref_metrics
endif

#bam_interval_metrics : hs_metrics interval_report non_ref_metrics
non_ref_metrics : $(foreach sample,$(SAMPLES),metrics/$(sample).interval_nonref_freq.txt)
hs_metrics : metrics/hs_metrics.txt metrics/interval_hs_metrics.txt
amplicon_metrics : $(foreach sample,$(SAMPLES),metrics/$(sample).amplicon_metrics.txt)
wgs_metrics : $(foreach sample,$(SAMPLES), metrics/$(sample).wgs_metrics.txt)
interval_report : metrics/interval_report/index.html

# interval metrics per sample
metrics/%.hs_metrics.txt metrics/%.interval_hs_metrics.txt : bam/%.bam bam/%.bam.bai
	$(call LSCRIPT_MEM,10G,00:59:59,"TMP=`mktemp`.intervals; \
	$(SAMTOOLS) view -H $< | grep '^@SQ' > \$$TMP &&  grep -P \"\t\" $(TARGETS_FILE_INTERVALS) | awk 'BEGIN {OFS = \"\t\"} { print \$$1$(,)\$$2+1$(,)\$$3$(,)\"+\"$(,)NR }' >> \$$TMP; \
	$(CALC_HS_METRICS) INPUT=$< OUTPUT=metrics/$*.hs_metrics.txt METRIC_ACCUMULATION_LEVEL=ALL_READS PER_TARGET_COVERAGE=metrics/$*.interval_hs_metrics.txt TARGET_INTERVALS=\$$TMP BAIT_SET_NAME=hs BAIT_INTERVALS=\$$TMP")

# not sure how this differs from above, see picard doc
metrics/%.amplicon_metrics.txt metrics/%.interval_amplicon_metrics.txt : bam/%.bam bam/%.bam.bai
	$(call LSCRIPT_MEM,10G,00:59:59,"TMP=`mktemp`.intervals; \
	$(SAMTOOLS) view -H $< | grep '^@SQ' > \$$TMP && grep -P \"\t\"  $(TARGETS_FILE_INTERVALS) | awk 'BEGIN {OFS = \"\t\"} { print \$$1$(,)\$$2+1$(,)\$$3$(,)\"+\"$(,)NR }' >> \$$TMP; \
	$(COLLECT_TARGETED_METRICS) INPUT=$< OUTPUT=$@ AMPLICON_INTERVALS=\$$TMP TARGET_INTERVALS=\$$TMP METRIC_ACCUMULATION_LEVEL=ALL_READS PER_TARGET_COVERAGE=metrics/$*.interval_amplicon_metrics.txt")

# summarize metrics into one file
metrics/hs_metrics.txt : $(foreach sample,$(SAMPLES),metrics/$(sample).hs_metrics.txt)
	$(INIT) \
	{ \
	sed '/^$$/d; /^#/d; s/SAMPLE.*//; s/BAIT_SET/SAMPLE/; s/\s$$//' $< | head -1; \
	for metrics in $^; do \
		samplename=$$(basename $${metrics%%.hs_metrics.txt}); \
	    sed "/^#/d; /^BAIT/d; /^\$$/d; s/^hs/$$samplename/; s/\t\+$$//" $$metrics; \
	done; \
	} > $@
# print mean, min, max for average bait coverage
# awk 'BEGIN {min = 10000000; max = 0; } NR > 1 {if ($21 > max) {max = $21; } if ($21 < min) {min = $21;} total += $21} END {print total / (NR-1), min, max }'

# summarize interval metrics into one file
metrics/interval_hs_metrics.txt : $(foreach sample,$(SAMPLES),metrics/$(sample).interval_hs_metrics.txt)
	$(INIT) \
	sed '/^#/d; /^$$/d' $< | cut -f 1-6 > $@.tmp; \
	for metrics in $^; do \
		samplename=$$(basename $${metrics%%.interval_hs_metrics.txt}); \
		sed '/^#/d; /^$$/d' $$metrics | cut -f 7,8 | sed "s/mean_coverage/$${samplename}_mean_coverage/; s/normalized_coverage/$${samplename}_normalized_coverage/" | paste $@.tmp - > $@; \
		cp $@ $@.tmp; \
	done; \
	rm -f $@.tmp

metrics/interval_report/index.html : metrics/hs_metrics.txt
	$(call LSCRIPT_MEM,3G,00:29:29,"$(PLOT_HS_METRICS) --outDir $(@D) $<")

metrics/%.interval_nonref_freq.txt : bam/%.bam
	$(call LSCRIPT_MEM,8G,00:29:29,"$(SAMTOOLS) mpileup -l $(TARGETS_FILE_INTERVALS) -f $(REF_FASTA) $< | $(NON_REF_FREQ) -b $(NON_REF_FREQ_BIN_SIZE) > $@")

include usb-modules/bam_tools/processBam.mk

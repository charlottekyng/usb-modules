# generate bam interval metrics per sample

include usb-modules/Makefile.inc
include usb-modules/config.inc

LOGDIR ?= log/metrics.$(NOW)

.DELETE_ON_ERROR:
.SECONDARY: 
.PHONY: bam_metrics #hs_metrics amplicon_metrics wgs_metrics rna_metrics #interval_report #non_ref_metrics

ifeq ($(CAPTURE_METHOD),NONE)
bam_metrics : wgs_metrics oxog_wgs flagstats alignment_summary_metrics dup
endif
ifeq ($(CAPTURE_METHOD),BAITS)
bam_metrics : hs_metrics oxog flagstats alignment_summary_metrics #dup
endif
ifeq ($(CAPTURE_METHOD),PCR)
bam_metrics : amplicon_metrics flagstats alignment_summary_metrics
endif
ifeq ($(CAPTURE_METHOD),RNA)
bam_metrics : rna_metrics flagstats alignment_summary_metrics
#oxog_wgs flagstats alignment_summary_metrics dup
endif

hs_metrics : $(shell rm -f metrics/all.hs_metrics.txt metrics/all.interval_hs_metrics.txt) metrics/all.hs_metrics.txt metrics/all.interval_hs_metrics.txt
amplicon_metrics : $(shell rm -f metrics/all.hs_metrics.txt metrics/all.interval_hs_metrics.txt) metrics/all.amplicon_metrics.txt metrics/all.interval_amplicon_metrics.txt
wgs_metrics : $(foreach sample,$(SAMPLES), metrics/$(sample).wgs_metrics.txt)
rna_metrics : $(shell rm -f metrics/all.hs_metrics.txt metrics/all.interval_hs_metrics.txt) metrics/all.rnaseq_metrics.txt metrics/all.normalized_coverage.rnaseq_metrics.txt
#metrics/all.rnaseq_report/index.html
flagstats : $(shell rm -f metrics/all.flagstats.txt) metrics/all.flagstats.txt
alignment_summary_metrics : $(shell rm -f metrics/all.alignment_summary_metrics.txt) metrics/all.alignment_summary_metrics.txt
#gc : $(foreach sample,$(SAMPLES),metrics/$(sample).gc_bias_metrics.txt)
artifacts : $(foreach sample,$(SAMPLES),metrics/$(sample).artifact_metrics.bait_bias_summary_metrics)
artifacts_wgs : $(foreach sample,$(SAMPLES),metrics/$(sample).wgs.artifact_metrics.bait_bias_summary_metrics)
oxog : metrics/all.oxog_metrics.txt
oxog_wgs : $(foreach sample,$(SAMPLES),metrics/$(sample).wgs.oxog_metrics.txt)
dup : $(shell rm -f metrics/all.dup_metrics.txt) metrics/all.dup_metrics.txt

#interval_report : metrics/interval_report/index.html
#non_ref_metrics : $(foreach sample,$(SAMPLES),metrics/$(sample).interval_nonref_freq.txt)

# interval metrics per sample
metrics/%.hs_metrics.txt metrics/%.interval_hs_metrics.txt : bam/%.bam bam/%.bam.bai
	$(call LSCRIPT_MEM,$(RESOURCE_REQ_MEDIUM_MEM),$(RESOURCE_REQ_SHORT),"$(LOAD_SAMTOOLS_MODULE); $(LOAD_JAVA8_MODULE); \
	TMP=`mktemp`.intervals; TMPCOVERED=`mktemp`.covered_intervals; \
	$(SAMTOOLS) view -H $< | grep '^@SQ' > \$$TMP &&  grep -P \"\t\" $(TARGETS_FILE_INTERVALS) | \
	awk 'BEGIN {OFS = \"\t\"} { print \$$1$(,)\$$2+1$(,)\$$3$(,)\"+\"$(,)NR }' >> \$$TMP; \
	$(SAMTOOLS) view -H $< | grep '^@SQ' > \$$TMPCOVERED &&  grep -P \"\t\" $(TARGETS_FILE_COVERED_INTERVALS) | \
	awk 'BEGIN {OFS = \"\t\"} { print \$$1$(,)\$$2+1$(,)\$$3$(,)\"+\"$(,)NR }' >> \$$TMPCOVERED; \
	$(call PICARD,CollectHsMetrics,$(RESOURCE_REQ_MEDIUM_MEM)) INPUT=$< OUTPUT=metrics/$*.hs_metrics.txt \
	PER_TARGET_COVERAGE=metrics/$*.interval_hs_metrics.txt TARGET_INTERVALS=\$$TMPCOVERED BAIT_SET_NAME=hs BAIT_INTERVALS=\$$TMP")

metrics/%.amplicon_metrics.txt metrics/%.interval_amplicon_metrics.txt : bam/%.bam bam/%.bam.bai
	$(call LSCRIPT_MEM,$(RESOURCE_REQ_MEDIUM_MEM),$(RESOURCE_REQ_SHORT),"$(LOAD_SAMTOOLS_MODULE); $(LOAD_JAVA8_MODULE); \
	TMP=`mktemp`.intervals; \
	$(SAMTOOLS) view -H $< | grep '^@SQ' > \$$TMP && grep -P \"\t\" $(TARGETS_FILE_INTERVALS) | \
	awk 'BEGIN {OFS = \"\t\"} { print \$$1$(,)\$$2+1$(,)\$$3$(,)\"+\"$(,)NR }' >> \$$TMP; \
	$(call PICARD,CollectTargetedPcrMetrics,$(RESOURCE_REQ_MEDIUM_MEM)) INPUT=$< OUTPUT=$@ AMPLICON_INTERVALS=\$$TMP TARGET_INTERVALS=\$$TMP \
	PER_TARGET_COVERAGE=metrics/$*.interval_amplicon_metrics.txt COVERAGE_CAP=50000")

define amplicon-metrics-pools
POOLNAME=$$(shell basename $2)
metrics/$1.amplicon_metrics_$$(POOLNAME).txt : bam/$1.bam bam/$1.bam.bai $2
	$$(call LSCRIPT_MEM,$(RESOURCE_REQ_MEDIUM_MEM),$(RESOURCE_REQ_SHORT),"$$(LOAD_SAMTOOLS_MODULE); $$(LOAD_JAVA8_MODULE); \
	TMP=`mktemp`.intervals; \
	$$(SAMTOOLS) view -H $$< | grep '^@SQ' > \$$$$TMP && grep -P \"\t\" $2 | \
	awk 'BEGIN {OFS = \"\t\"} { print \$$$$1$$(,)\$$$$2+1$$(,)\$$$$3$$(,)\"+\"$$(,)NR }' >> \$$$$TMP; \
	$$(call PICARD,CollectTargetedPcrMetrics,$$(RESOURCE_REQ_MEDIUM_MEM)) INPUT=$$< OUTPUT=$$@ AMPLICON_INTERVALS=\$$$$TMP TARGET_INTERVALS=\$$$$TMP \
	COVERAGE_CAP=50000")
endef
$(if $(TARGETS_FILE_INTERVALS_POOLS),\
	$(foreach pool,$(TARGETS_FILE_INTERVALS_POOLS),\
		$(foreach sample,$(SAMPLES),\
			$(eval $(call amplicon-metrics-pools,$(sample),$(pool))))))			

metrics/%.wgs_metrics.txt : bam/%.bam bam/%.bam.bai
	$(call LSCRIPT_MEM,$(RESOURCE_REQ_MEDIUM_MEM),$(RESOURCE_REQ_SHORT),"$(LOAD_JAVA8_MODULE); \
		$(call PICARD,CollectWgsMetrics,$(RESOURCE_REQ_MEDIUM_MEM)) \
		INPUT=$< OUTPUT=$@ COUNT_UNPAIRED=true")

metrics/%.rnaseq_metrics.txt : bam/%.bam bam/%.bam.bai
	$(call LSCRIPT_CHECK_MEM,$(RESOURCE_REQ_MEDIUM_MEM),$(RESOURCE_REQ_SHORT),"$(LOAD_R_MODULE); $(LOAD_JAVA8_MODULE); \
		$(call PICARD,CollectRnaSeqMetrics,$(RESOURCE_REQ_MEDIUM_MEM)) \
		REF_FLAT=$(GENE_REF_FLAT) RIBOSOMAL_INTERVALS=$(RIBOSOMAL_INTERVALS) \
		STRAND_SPECIFICITY=$(STRAND_SPECIFICITY) \
		INPUT=$< OUTPUT=$@ CHART_OUTPUT=$@.pdf VERBOSITY=ERROR")

metrics/%.alignment_summary_metrics.txt : bam/%.bam bam/%.bam.bai
	$(call LSCRIPT_MEM,$(RESOURCE_REQ_MEDIUM_MEM),$(RESOURCE_REQ_SHORT),"$(LOAD_JAVA8_MODULE); \
	$(call PICARD,CollectAlignmentSummaryMetrics,$(RESOURCE_REQ_MEDIUM_MEM)) INPUT=$< OUTPUT=$@")

# does not work, there's a conflict of java versions
#metrics/%.gc_bias_metrics.txt metrics/%.gc_bias_metrics_summary.txt : bam/%.bam bam/%.bam.bai
#	$(call LSCRIPT_MEM,12G,02:59:59,"$(LOAD_JAVA6_MODULE); $(LOAD_R_MODULE); \
#		$(call COLLECT_GCBIAS_METRICS,11G) \
#		INPUT=$< OUTPUT=$@ CHART_OUTPUT=$(addsuffix .pdf,$@) \
#		SUMMARY_OUTPUT=metrics/$*.gc_bias_metrics_summary.txt")

metrics/%.artifact_metrics.bait_bias_summary_metrics : bam/%.bam bam/%.bam.bai
	$(call LSCRIPT_MEM,$(RESOURCE_REQ_MEDIUM_MEM),$(RESOURCE_REQ_SHORT),"$(LOAD_SAMTOOLS_MODULE); $(LOAD_JAVA8_MODULE); \
	TMP=`mktemp`.intervals; \
	$(SAMTOOLS) view -H $< | grep '^@SQ' > \$$TMP &&  grep -P \"\t\" $(TARGETS_FILE_INTERVALS) | \
	awk 'BEGIN {OFS = \"\t\"} { print \$$1$(,)\$$2+1$(,)\$$3$(,)\"+\"$(,)NR }' >> \$$TMP; \
	$(call PICARD,CollectSequencingArtifactMetrics,$(RESOURCE_REQ_MEDIUM_MEM)) INPUT=$< OUTPUT=$@ \
	DB_SNP=$(DBSNP) INTERVALS=\$$TMP")

metrics/%.wgs.artifact_metrics.bait_bias_summary_metrics : bam/%.bam bam/%.bam.bai
	$(call LSCRIPT_MEM,$(RESOURCE_REQ_MEDIUM_MEM),$(RESOURCE_REQ_SHORT),"$(LOAD_JAVA8_MODULE); \
	$(call PICARD,CollectSequencingArtifactMetrics,$(RESOURCE_REQ_MEDIUM_MEM)) \	
	INPUT=$< OUTPUT=$@ DB_SNP=$(DBSNP)")

metrics/%.oxog_metrics.txt : bam/%.bam bam/%.bam.bai
	$(call LSCRIPT_MEM,$(RESOURCE_REQ_MEDIUM_MEM),$(RESOURCE_REQ_SHORT),"$(LOAD_SAMTOOLS_MODULE); $(LOAD_JAVA8_MODULE); \
	TMP=`mktemp`.intervals; \
	$(SAMTOOLS) view -H $< | grep '^@SQ' > \$$TMP &&  grep -P \"\t\" $(TARGETS_FILE_INTERVALS) | \
	awk 'BEGIN {OFS = \"\t\"} { print \$$1$(,)\$$2+1$(,)\$$3$(,)\"+\"$(,)NR }' >> \$$TMP; \
	$(call PICARD,CollectOxoGMetrics,$(RESOURCE_REQ_MEDIUM_MEM)) INPUT=$< OUTPUT=$@ \
	DB_SNP=$(DBSNP) INTERVALS=\$$TMP")

metrics/%.wgs.oxog_metrics.txt : bam/%.bam bam/%.bam.bai
	$(call LSCRIPT_MEM,$(RESOURCE_REQ_MEDIUM_MEM),$(RESOURCE_REQ_SHORT),"$(LOAD_JAVA8_MODULE); \
	$(call PICARD,CollectOxoGMetrics,$(RESOURCE_REQ_MEDIUM_MEM)) INPUT=$< OUTPUT=$@ DB_SNP=$(DBSNP)")

metrics/%.flagstats.txt : bam/%.bam bam/%.bam.bai
	$(call LSCRIPT_MEM,$(RESOURCE_REQ_LOWMEM),$(RESOURCE_REQ_SHORT),"$(LOAD_SAMTOOLS_MODULE); $(SAMTOOLS) flagstat $< > $@")

# summarize metrics into one file
metrics/all.hs_metrics.txt : $(foreach sample,$(SAMPLES),metrics/$(sample).hs_metrics.txt)
	$(INIT) \
	{ \
	sed '/^$$/d; /^#/d; s/SAMPLE.*//; s/BAIT_SET/SAMPLE/; s/\s$$//' $< | head -1; \
	for metrics in $^; do \
		samplename=$$(basename $${metrics%%.hs_metrics.txt}); \
	    sed "/^#/d; /^BAIT/d; /^\$$/d; s/^hs/$$samplename/; s/\t\+$$//" $$metrics | grep "^$$samplename"; \
	done; \
	} > $@

# summarize interval metrics into one file
metrics/all.interval_hs_metrics.txt : $(foreach sample,$(SAMPLES),metrics/$(sample).interval_hs_metrics.txt)
	$(INIT) \
	sed '/^#/d; /^$$/d' $< | cut -f 1-6 > $@.tmp; \
	for metrics in $^; do \
		samplename=$$(basename $${metrics%%.interval_hs_metrics.txt}); \
		sed '/^#/d; /^$$/d' $$metrics | cut -f 7,8 | \
		sed "s/mean_coverage/$${samplename}_mean_coverage/; s/normalized_coverage/$${samplename}_normalized_coverage/" | \
		paste $@.tmp - > $@; \
		cp $@ $@.tmp; \
	done; \
	rm -f $@.tmp

# summarize metrics into one file
metrics/all.amplicon_metrics.txt : $(foreach sample,$(SAMPLES),metrics/$(sample).amplicon_metrics.txt) $(if $(TARGETS_FILE_INTERVALS_POOLS),$(foreach pool,$(TARGETS_FILE_INTERVALS_POOLS),$(foreach sample,$(SAMPLES),metrics/$(sample).amplicon_metrics_$(shell basename $(pool)).txt)))
	$(INIT) \
	{ \
	sed '/^$$/d; /^#/d; s/CUSTOM_AMPLICON_SET/SAMPLE/; s/\s$$//' $< | head -1; \
	for metrics in $^; do \
		samplename=$$(basename $${metrics%%.amplicon_metrics.txt}); \
		head -8 $$metrics | sed "/^#/d; /^CUSTOM_AMPLICON_SET/d; /^\$$/d; s/^tmp/$$samplename/; s/\t\+$$//" ; \
	done; \
	} > $@

# summarize interval metrics into one file
metrics/all.interval_amplicon_metrics.txt : $(foreach sample,$(SAMPLES),metrics/$(sample).interval_amplicon_metrics.txt)
	$(INIT) \
	sed '/^#/d; /^$$/d' $< | cut -f 1-6 > $@.tmp; \
	for metrics in $^; do \
	samplename=$$(basename $${metrics%%.interval_amplicon_metrics.txt}); \
	sed '/^#/d; /^$$/d' $$metrics | cut -f 7,8 | \
	sed "s/mean_coverage/$${samplename}_mean_coverage/; s/normalized_coverage/$${samplename}_normalized_coverage/" | \
		paste $@.tmp - > $@; \
		cp $@ $@.tmp; \
	done; \
	rm -f $@.tmp

metrics/all.rnaseq_metrics.txt : $(foreach sample,$(SAMPLES),metrics/$(sample).rnaseq_metrics.txt)
	$(INIT) \
	grep '^PF_BASES' $< > $@ && \
	for i in $^; do sample=`echo $$i | sed 's:.*/::; s/\..*//'`; \
		grep -A1 '^PF_BASES' $$i | tail -1 | awk -v sample=$$sample 'BEGIN { OFS = "\t" } { $$23=sample; print $$0 }'  >> $@; \
	done; 

metrics/all.normalized_coverage.rnaseq_metrics.txt : $(foreach sample,$(SAMPLES),metrics/$(sample).rnaseq_metrics.txt)
	$(INIT) \
	grep -A101 '^normalized_position' $< | cut -f1 > $@ && \
	for i in $^; do sample=`echo $$i | sed 's:.*/::; s/\..*//'`; \
		grep -A101 '^normalized_position' $$i | cut -f2 | sed "s/All_Reads/$$sample/" | paste $@ - > $@.tmp && mv $@.tmp $@; \
	done;

metrics/all.rnaseq_report/index.html : metrics/all.rnaseq_metrics.txt metrics/all.normalized_coverage.rnaseq_metrics.txt
	$(INIT) $(LOAD_R_MODULE); $(PLOT_RNASEQ_METRICS) --outDir $(@D) $^
      
metrics/all.alignment_summary_metrics.txt : $(foreach sample,$(SAMPLES),metrics/$(sample).alignment_summary_metrics.txt)
	$(INIT) \
	{ \
	sed '/^$$/d; /^#/d; s/SAMPLE.*//; s/\s$$//; s/^/SAMPLE\t/;' $< | head -1; \
	for metrics in $^; do \
		samplename=$$(basename $${metrics%%.alignment_summary_metrics.txt}); \
		sed "/^#/d; /^CATEGORY/d; /^\$$/d; s/^/$$samplename\t/; s/\t\+$$//" $$metrics | grep "^$$samplename"; \
	done; \
	} >$@

metrics/all.dup_metrics.txt : $(foreach sample,$(SAMPLES),metrics/$(sample).dup_metrics.txt)
	$(INIT) \
	{ \
	sed '/^$$/d; /^#/d; s/SAMPLE.*//; s/\s$$//; s/^/SAMPLE\t/;' $< | head -1; \
	for metrics in $^; do \
		samplename=$$(basename $${metrics%%.dup_metrics.txt}); \
		sed "/^#/d; /^LIBRARY/d; /^\$$/d; s/^/$$samplename\t/; s/\t\+$$//" $$metrics | grep "^$$samplename"; \
	done; \
	} >$@

metrics/all.oxog_metrics.txt : $(foreach sample,$(SAMPLES),metrics/$(sample).oxog_metrics.txt)
	$(INIT) \
	{ \
	sed '/^$$/d; /^#/d; s/\s$$//;' $< | head -1; \
	for metrics in $^; do \
		samplename=$$(basename $${metrics%%.oxog_metrics.txt}); \
		sed "/^#/d; /^SAMPLE_ALIAS/d; /^\$$/d; s/^[^\t]\+\t/$$samplename\t/;" $$metrics | \
		grep "^$$samplename" | sort -r -k 11 | head -1; \
	done; \
	} >$@

metrics/all.flagstats.txt : $(foreach sample,$(SAMPLES),metrics/$(sample).flagstats.txt)
	$(INIT) \
	{ \
	echo -ne "category\t"; sed 's/^[0-9]\+ + [0-9]\+ //;' $< | tr '\n' '\t';  echo "";\
	for metrics in $^; do \
		samplename=$$(basename $${metrics%%.flagstats.txt}); echo -ne "$$samplename\t"; \
		cut -f1 -d' ' $$metrics | tr '\n' '\t'; echo ""; \
	done; \
	} >$@
	
#metrics/interval_report/index.html : metrics/hs_metrics.txt
#	$(call LSCRIPT_MEM,3G,00:29:29,"$(PLOT_HS_METRICS) --outDir $(@D) $<")

#metrics/%.interval_nonref_freq.txt : bam/%.bam
#	$(call LSCRIPT_MEM,2G,01:59:29,"$(LOAD_SAMTOOLS_MODULE); $(LOAD_JAVA8_MODULE); \
#	$(SAMTOOLS) mpileup -l $(TARGETS_FILE_INTERVALS) -f $(REF_FASTA) $< | \
#	$(NON_REF_FREQ) -b $(NON_REF_FREQ_BIN_SIZE) > $@")

include usb-modules/bam_tools/processBam.mk

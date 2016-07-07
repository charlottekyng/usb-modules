# generate bam interval metrics per sample
SEQ_PLATFORM = ILLUMINA
REF = b37
PANEL = AGILENT_CLINICAL_EXOME
CAPTURE_METHOD = BAITS

include usb-modules/Makefile.inc
include usb-modules/config.inc

LOGDIR ?= log/metrics.$(NOW)

.DELETE_ON_ERROR:
.SECONDARY: 
.PHONY: bam_metrics #hs_metrics amplicon_metrics wgs_metrics rna_metrics #interval_report #non_ref_metrics

ifeq ($(CAPTURE_METHOD),NONE)
bam_metrics : wgs_metrics artifacts_wgs oxog_wgs flagstats gc alignment_summary_metrics
endif
ifeq ($(CAPTURE_METHOD),BAITS)
bam_metrics : hs_metrics artifacts oxog flagstats gc alignment_summary_metrics
endif
ifeq ($(CAPTURE_METHOD),PCR)
bam_metrics : amplicon_metrics artifacts oxog flagstats gc alignment_summary_metrics
endif
ifeq ($(CAPTURE_METHOD),RNA)
bam_metrics : rna_metrics artifacts_wgs oxog_wgs flagstats gc alignment_summary_metrics
endif

hs_metrics : metrics/hs_metrics.txt metrics/interval_hs_metrics.txt
amplicon_metrics : metrics/amplicon_metrics.txt metrics/interval_amplicon_metrics.txt
wgs_metrics : $(foreach sample,$(SAMPLES), metrics/$(sample).wgs_metrics.txt)
rna_metrics : metrics/all.rnaseq_metrics.txt metrics/all.normalized_coverage.rnaseq_metrics.txt metrics/rnaseq_report/index.html
flagstats : $(foreach sample,$(SAMPLES),metrics/$(sample).flagstats.txt)
alignment_summary_metrics : $(foreach sample,$(SAMPLES),metrics/$(sample).alignment_summary_metrics.txt)
gc : $(foreach sample,$(SAMPLES),metrics/$(sample).gc_bias_metrics.txt)
artifacts : $(foreach sample,$(SAMPLES),metrics/$(sample).artifact_metrics.txt)
artifacts_wgs : $(foreach sample,$(SAMPLES),metrics/$(sample).wgs.artifact_metrics.txt)
oxog : $(foreach sample,$(SAMPLES),metrics/$(sample).oxog_metrics.txt)
oxog_wgs : $(foreach sample,$(SAMPLES),metrics/$(sample).wgs.oxog_metrics.txt)

#interval_report : metrics/interval_report/index.html
#non_ref_metrics : $(foreach sample,$(SAMPLES),metrics/$(sample).interval_nonref_freq.txt)

# interval metrics per sample
metrics/%.hs_metrics.txt metrics/%.interval_hs_metrics.txt : bam/%.bam bam/%.bam.bai
	$(call LSCRIPT_MEM,10G,00:59:59,"$(LOAD_SAMTOOLS_MODULE); $(LOAD_JAVA8_MODULE); \
	TMP=`mktemp`.intervals; TMPBAITS=`mktemp`.baits_intervals; \
	$(SAMTOOLS) view -H $< | grep '^@SQ' > \$$TMP &&  grep -P \"\t\" $(TARGETS_FILE_INTERVALS) | \
	awk 'BEGIN {OFS = \"\t\"} { print \$$1$(,)\$$2+1$(,)\$$3$(,)\"+\"$(,)NR }' >> \$$TMP; \
	$(SAMTOOLS) view -H $< | grep '^@SQ' > \$$TMP &&  grep -P \"\t\" $(TARGETS_FILE_BAITS_INTERVALS) | \
	awk 'BEGIN {OFS = \"\t\"} { print \$$1$(,)\$$2+1$(,)\$$3$(,)\"+\"$(,)NR }' >> \$$TMPBAITS; \
	$(call COLLECT_HS_METRICS,9G) INPUT=$< OUTPUT=metrics/$*.hs_metrics.txt \
	PER_TARGET_COVERAGE=metrics/$*.interval_hs_metrics.txt TARGET_INTERVALS=\$$TMP BAIT_SET_NAME=hs BAIT_INTERVALS=\$$TMPBAITS")

# not sure how this differs from above, see picard doc
metrics/%.amplicon_metrics.txt metrics/%.interval_amplicon_metrics.txt : bam/%.bam bam/%.bam.bai
	$(call LSCRIPT_MEM,10G,00:59:59,"$(LOAD_SAMTOOLS_MODULE); $(LOAD_JAVA8_MODULE); \
	TMP=`mktemp`.intervals; \
	$(SAMTOOLS) view -H $< | grep '^@SQ' > \$$TMP && grep -P \"\t\" $(TARGETS_FILE_INTERVALS) | \
	awk 'BEGIN {OFS = \"\t\"} { print \$$1$(,)\$$2+1$(,)\$$3$(,)\"+\"$(,)NR }' >> \$$TMP; \
	$(call COLLECT_TARGETEDPCR_METRICS,9G) INPUT=$< OUTPUT=$@ AMPLICON_INTERVALS=\$$TMP TARGET_INTERVALS=\$$TMP \
	PER_TARGET_COVERAGE=metrics/$*.interval_amplicon_metrics.txt COVERAGE_CAP=50000")

metrics/%.wgs_metrics.txt : bam/%.bam bam/%.bam.bai
	$(call LSCRIPT_MEM,12G,02:59:59,"$(LOAD_JAVA8_MODULE); $(call COLLECT_WGS_METRICS,11G) \
		INPUT=$< OOUTPUT=$@ REFERENCE_SEQUENCE=$(REF_FASTA) COUNT_UNPAIRED=true")

metrics/%.rnaseq_metrics.txt : bam/%.bam bam/%.bam.bai
	$(call LSCRIPT_CHECK_MEM,8G,01:59:59,"$(LOAD_JAVA8_MODULE); $(call COLLECT_RNASEQ_METRICS,8G) \
		REF_FLAT=$(GENE_REF_FLAT) RIBOSOMAL_INTERVALS=$(RIBOSOMAL_INTERVALS) \
		STRAND_SPECIFICITY=$(STRAND_SPECIFICITY) REFERENCE_SEQUENCE=$(REF_FASTA) \
		INPUT=$< OUTPUT=$@ CHART_OUTPUT=$@.pdf VERBOSITY=ERROR")

metrics/%.alignment_summary_metrics.txt : bam/%.bam bam/%.bam.bai
	$(call LSCRIPT_MEM,12G,02:59:59,"$(LOAD_JAVA8_MODULE); $(call COLLECT_ALIGNMENT_METRICS,11G) \
		INPUT=$< OUTPUT=metrics/$@ REFERENCE_SEQUENCE=$(REF_FASTA)")

metrics/%.gc_bias_metrics.txt metrics/%.gc_bias_metrics_summary.txt : bam/%.bam bam/%.bam.bai
	$(call LSCRIPT_MEM,12G,02:59:59,"$(LOAD_JAVA8_MODULE); $(call COLLECT_GCBIAS_METRICS,11G) \
		INPUT=$< OUTPUT=$@ CHART_OUTPUT=$(addsuffix .pdf,$@) \
		SUMMARY_OUTPUT=metrics/$*.gc_bias_metrics_summary.txt REFERENCE_SEQUENCE=$(REF_FASTA)")

metrics/%.artifact_metrics.txt : bam/%.bam bam/%.bam.bai
	$(call LSCRIPT_MEM,12G,02:59:59,"$(LOAD_SAMTOOLS_MODULE); $(LOAD_JAVA8_MODULE); \
	TMP=`mktemp`.intervals; \
	$(SAMTOOLS) view -H $< | grep '^@SQ' > \$$TMP &&  grep -P \"\t\" $(TARGETS_FILE_INTERVALS) | \
	awk 'BEGIN {OFS = \"\t\"} { print \$$1$(,)\$$2+1$(,)\$$3$(,)\"+\"$(,)NR }' >> \$$TMP; \
	$(call COLLECT_SEQ_ARTIFACT_METRICS,11G) INPUT=$< OUTPUT=$@ REFERENCE_SEQUENCE=$(REF_FASTA) \
	DB_SNP=$(DBSNP) TARGET_INTERVALS=\$$TMP")

metrics/%.wgs.artifact_metrics.txt : bam/%.bam bam/%.bam.bai
	$(call LSCRIPT_MEM,12G,02:59:59,"$(LOAD_JAVA8_MODULE); $(call COLLECT_SEQ_ARTIFACT_METRICS,11G) \
		INPUT=$< OUTPUT=$@ REFERENCE_SEQUENCE=$(REF_FASTA) DB_SNP=$(DBSNP)")

metrics/%.oxog_metrics.txt : bam/%.bam bam/%.bam.bai
	$(call LSCRIPT_MEM,12G,02:59:59,"$(LOAD_SAMTOOLS_MODULE); $(LOAD_JAVA8_MODULE); \
	TMP=`mktemp`.intervals; \
	$(SAMTOOLS) view -H $< | grep '^@SQ' > \$$TMP &&  grep -P \"\t\" $(TARGETS_FILE_INTERVALS) | \
	awk 'BEGIN {OFS = \"\t\"} { print \$$1$(,)\$$2+1$(,)\$$3$(,)\"+\"$(,)NR }' >> \$$TMP; \
	$(call COLLECT_OXOG_METRICS,11G) INPUT=$< OUTPUT=$@ REFERENCE_SEQUENCE=$(REF_FASTA) \
	DB_SNP=$(DBSNP) TARGET_INTERVALS=\$$TMP")

metrics/%.wgs.oxog_metrics.txt : bam/%.bam bam/%.bam.bai
	$(call LSCRIPT_MEM,12G,02:59:59,"$(LOAD_JAVA8_MODULE); $(call COLLECT_OXOG_METRICS,11G) \
		INPUT=$< OUTPUT=$@ REFERENCE_SEQUENCE=$(REF_FASTA) DB_SNP=$(DBSNP)")

metrics/%.flagstats.txt : bam/%.bam bam/%.bam.bai
	$(call LSCRIPT_MEM,2G,00:29:29,"$(LOAD_SAMTOOLS_MODULE); $(SAMTOOLS) flagstat $< > $@")

# summarize metrics into one file
metrics/hs_metrics.txt : $(foreach sample,$(SAMPLES),metrics/$(sample).hs_metrics.txt)
	$(INIT) \
	{ \
	sed '/^$$/d; /^#/d; s/SAMPLE.*//; s/BAIT_SET/SAMPLE/; s/\s$$//' $< | head -1; \
	for metrics in $^; do \
		samplename=$$(basename $${metrics%%.hs_metrics.txt}); \
	    sed "/^#/d; /^BAIT/d; /^\$$/d; s/^hs/$$samplename/; s/\t\+$$//" $$metrics | grep "^$$samplename"; \
	done; \
	} > $@

# summarize interval metrics into one file
metrics/interval_hs_metrics.txt : $(foreach sample,$(SAMPLES),metrics/$(sample).interval_hs_metrics.txt)
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
metrics/amplicon_metrics.txt : $(foreach sample,$(SAMPLES),metrics/$(sample).amplicon_metrics.txt)
	$(INIT) \
	{ \
	sed '/^$$/d; /^#/d; s/CUSTOM_AMPLICON_SET/SAMPLE/; s/\s$$//' $< | head -1; \
	for metrics in $^; do \
		samplename=$$(basename $${metrics%%.amplicon_metrics.txt}); \
		head -8 $$metrics | sed "/^#/d; /^CUSTOM_AMPLICON_SET/d; /^\$$/d; s/^tmp/$$samplename/; s/\t\+$$//" ; \
	done; \
	} > $@

# summarize interval metrics into one file
metrics/interval_amplicon_metrics.txt : $(foreach sample,$(SAMPLES),metrics/$(sample).interval_amplicon_metrics.txt)
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
	done'

metrics/rnaseq_report/index.html : metrics/all.rnaseq_metrics.txt metrics/all.normalized_coverage.rnaseq_metrics.txt
	$(INIT) $(LOAD_R_MODULE); $(PLOT_RNASEQ_METRICS) --outDir $(@D) $^

#metrics/interval_report/index.html : metrics/hs_metrics.txt
#	$(call LSCRIPT_MEM,3G,00:29:29,"$(PLOT_HS_METRICS) --outDir $(@D) $<")

#metrics/%.interval_nonref_freq.txt : bam/%.bam
#	$(call LSCRIPT_MEM,2G,01:59:29,"$(LOAD_SAMTOOLS_MODULE); $(LOAD_JAVA8_MODULE); \
#	$(SAMTOOLS) mpileup -l $(TARGETS_FILE_INTERVALS) -f $(REF_FASTA) $< | \
#	$(NON_REF_FREQ) -b $(NON_REF_FREQ_BIN_SIZE) > $@")

include usb-modules/bam_tools/processBam.mk

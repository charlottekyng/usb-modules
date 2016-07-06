## defaults
LOGDIR = metrics/log

## includes
include modules/Makefile.inc
include modules/variant_callers/gatk.inc


.DELETE_ON_ERROR:

.SECONDARY: 

PHONY += bam_metrics
bam_metrics : summary_metrics gc flagstats wgs_metrics

PHONY += flagstats
flagstats : $(foreach sample,$(SAMPLES),metrics/$(sample).flagstats)
PHONY += summary_metrics
summary_metrics : $(foreach sample,$(SAMPLES),metrics/$(sample).alignment_summary_metrics)
PHONY += wgs_metrics
wgs_metrics : $(foreach sample,$(SAMPLES),metrics/$(sample).wgs_metrics)
PHONY += dup
dup : $(foreach sample,$(SAMPLES),metrics/$(sample).dup_metrics)
PHONY += gc
gc : $(foreach sample,$(SAMPLES),metrics/$(sample).gc_bias_metrics)

metrics/%.alignment_summary_metrics : bam/%.bam
	$(call LSCRIPT_MEM,12G,02:59:59,"$(call COLLECT_ALIGNMENT_METRICS,11G) I=$< O=metrics/$* REFERENCE_SEQUENCE=$(REF_FASTA)")

metrics/%.wgs_metrics : bam/%.bam
	$(call LSCRIPT_MEM,12G,02:59:59,"$(call COLLECT_WGS_METRICS,11G) I=$< O=$@ REFERENCE_SEQUENCE=$(REF_FASTA)")

metrics/%.gc_bias_metrics : bam/%.bam
	$(call LSCRIPT_MEM,12G,02:59:59,"$(call COLLECT_GC_METRICS,11G) I=$< O=$@ CHART_OUTPUT=$(addsuffix .pdf,$@) REFERENCE_SEQUENCE=$(REF_FASTA)")

metrics/%.flagstats : bam/%.bam
	$(call LSCRIPT_MEM,2G,00:29:29,"$(SAMTOOLS) flagstat $< > $@")
	
bam/%.markdup.bam metrics/%.dup_metrics : bam/%.bam
	$(call LSCRIPT_MEM,18G,02:59:59,"$(call MARK_DUP,17G) I=$< O=bam/$*.markdup.bam METRICS_FILE=metrics/$*.dup_metrics")

metrics/dup_metrics.txt : $(foreach sample,$(SAMPLES),metrics/$(sample).dup_metrics.txt)
	$(INIT) grep '^LIBRARY' $< > $@ && \
	for metrics in $^; do \
	    grep -A1 '^LIBRARY' $$metrics | sed '1d' >> $@; \
	done

.PHONY: $(PHONY)

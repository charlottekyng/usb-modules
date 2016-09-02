include usb-modules/Makefile.inc
include usb-modules/config.inc
include usb-modules/variant_callers/variantCaller.inc 

LOGDIR ?= log/facets.$(NOW)

.SECONDARY:
.DELETE_ON_ERROR:
.PHONY : facets

facets : $(foreach pair,$(SAMPLE_PAIRS),facets/cncf/$(pair).cncf.txt) facets/geneCN.txt
#	facets/geneCN.txt facets/geneCN.fill.txt facets/geneCN.heatmap.pdf facets/geneCN.fill.heatmap.pdf

ifeq ($(FACETS_GATK_VARIANTS),true)
GET_BASE_COUNTS_POS2 = facets/base_pos/snploc_chr
else
GET_BASE_COUNTS_POS2 = $(GET_BASE_COUNTS_POS)
endif

#facets/base_pos/gatk.vcf : $(foreach normalsample,$(NORMAL_SAMPLES),vcf/$(normalsample).$(call VCF_SUFFIXES,gatk_snps).vcf)
facets/base_pos/gatk.vcf : $(foreach normalsample,$(NORMAL_SAMPLES),vcf/$(normalsample).gatk_snps.vcf)
	$(call LSCRIPT_MEM,22G,03:59:59,"$(LOAD_JAVA8_MODULE); $(call COMBINE_VARIANTS,21G) \
		$(foreach vcf,$^,--variant $(vcf) ) -o $@ --genotypemergeoption UNSORTED -R $(REF_FASTA)")

define base-count-pos
facets/base_pos/snploc_chr$1 : facets/base_pos/gatk.pass.vcf $$(GET_BASE_COUNTS_POS)$1
	$(INIT) \
	$(MKDIR) facets/base_pos/ && \
	grep -v "^\#" $$< | perl -p -e "/^$$$$1\t/;" | cut -f2 > $$@ && \
	cat $$(word 2,$$^) >> $$@ && \
	sort -n $$@ > $$@.tmp && mv $$@.tmp $$@;	
endef
$(foreach chr,$(CHROMOSOMES),$(eval $(call base-count-pos,$(chr))))

define snp-pileup-tumor-normal
facets/snp_pileup/$1_$2.bc.gz : bam/$1.bam bam/$2.bam $(foreach chr,$(CHROMOSOMES),$(GET_BASE_COUNTS_POS2)$(chr))
	$$(call LSCRIPT_CHECK_MEM,3G,00:59:59,"$$(LOAD_PERL_MODULE); $$(GET_BASE_COUNTS) bam/$2.bam bam/$1.bam $$@ \
	$$(REF_FASTA) $$(TARGETS_FILE_INTERVALS) $$(GET_BASE_COUNTS_MIN_DEPTH) $$(GET_BASE_COUNTS_MAX_DEPTH) $$(GET_BASE_COUNTS_POS2)")
endef
$(foreach pair,$(SAMPLE_PAIRS),$(eval $(call snp-pileup-tumor-normal,$(tumor.$(pair)),$(normal.$(pair)))))

#facets/params.txt : $(foreach pair,$(SAMPLE_PAIRS),facets/snp_pileup/$(pair).bc.gz)
	

facets/cncf/%.cncf.txt : facets/snp_pileup/%.bc.gz
	$(call LSCRIPT_CHECK_MEM,3G,00:29:59,"$(LOAD_R_MODULE); $(FACETS) --minNDepth $(GET_BASE_COUNTS_MIN_DEPTH) --snp_nbhd $(FACETS_WINDOW_SIZE) \
	--cval2 $(FACETS_CVAL2) --cval1 $(FACETS_CVAL1) --genome $(REF) --min_nhet $(FACETS_MIN_NHET) --pre_cval $(FACETS_PRE_CVAL)  \
	--outPrefix $(@D)/$* $<")

facets/geneCN.txt : $(foreach pair,$(SAMPLE_PAIRS),facets/cncf/$(pair).cncf.txt)
	$(call LSCRIPT_CHECK_MEM,8G,00:29:59,"$(LOAD_R_MODULE); $(FACETS_GENE_CN) $(FACETS_GENE_CN_OPTS) --outFile $@ $^")

facets/geneCN.fill.txt : facets/geneCN.txt $(foreach pair,$(SAMPLE_PAIRS),facets/cncf/$(pair).cncf.txt)
	$(call LSCRIPT_CHECK_MEM,4G,00:29:59,"$(LOAD_R_MODULE); $(FACETS_FILL_GENE_CN) --outFile $@ --geneCNFile $< \
		$(filter %.cncf.txt,$^)")

facets/geneCN%heatmap.pdf  : facets/geneCN%txt
	$(call LSCRIPT_MEM,4G,00:29:29,"$(LOAD_R_MODULE); $(FACETS_PLOT_GENE_CN) $(FACETS_PLOT_GENE_CN_OPTS) $< $@")

include usb-modules/variant_callers/gatk.mk
include usb-modules/vcf_tools/vcftools.mk
include usb-modules/variant_callers/gatkVariantCaller.mk
include usb-modules/bam_tools/processBam.mk

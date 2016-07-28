include usb-modules/Makefile.inc

LOGDIR ?= log/facets.$(NOW)

.SECONDARY:
.DELETE_ON_ERROR:
.PHONY : facets

facets : $(foreach pair,$(SAMPLE_PAIRS),facets/cncf/$(pair).cncf.txt) facets/geneCN.txt
#	facets/geneCN.txt facets/geneCN.fill.txt facets/geneCN.heatmap.pdf facets/geneCN.fill.heatmap.pdf

facets/vcf/dbsnp_het_gatk.snps.vcf : $(FACETS_DBSNP:.gz=) $(foreach sample,$(SAMPLES),gatk/vcf/$(sample).variants.snps.het.pass.vcf) 
	$(call LSCRIPT_CHECK_MEM,4G,6G,"$(call GATK_MEM,3G) -T CombineVariants \
		--minimalVCF $(foreach i,$^, --variant $i) -R $(REF_FASTA) -o $@")

# flag homozygous calls
%.het.vcf : %.vcf
	$(call LSCRIPT_CHECK_MEM,9G,12G,"$(call GATK_MEM,8G) -V $< -T VariantFiltration -R $(REF_FASTA) \
		--genotypeFilterName 'hom' --genotypeFilterExpression 'isHet == 0' -o $@")

%.vcf.gz : %.vcf
	$(INIT) cat $< | gzip -c > $@

%.vcf : %.vcf.gz
	$(INIT) zcat $< > $@

# no flag target definitions
facets/vcf/targets_dbsnp.vcf.gz : $(TARGETS_FILE)
	$(INIT) $(BEDTOOLS) intersect -header -u -a $(DBSNP) -b $< | gzip -c > $@

define snp-pileup-tumor-normal
facets/snp_pileup/$1_$2.bc.gz : bam/$1.bam bam/$2.bam
	$$(call LSCRIPT_CHECK_MEM,3G,00:59:59,"$$(LOAD_PERL_MODULE); $$(GET_BASE_COUNTS) bam/$2.bam bam/$1.bam $$@ $$(GET_BASE_COUNTS_PARAMS)")
endef
$(foreach pair,$(SAMPLE_PAIRS),$(eval $(call snp-pileup-tumor-normal,$(tumor.$(pair)),$(normal.$(pair)))))

#facets/params.txt : $(foreach pair,$(SAMPLE_PAIRS),facets/snp_pileup/$(pair).bc.gz)
	

facets/cncf/%.cncf.txt : facets/snp_pileup/%.bc.gz
	$(call LSCRIPT_CHECK_MEM,3G,00:29:59,"$(LOAD_R_MODULE); $(FACETS) $(FACETS_OPTS) --outPrefix $(@D)/$* $<")

facets/geneCN.txt : $(foreach pair,$(SAMPLE_PAIRS),facets/cncf/$(pair).cncf.txt)
	$(call LSCRIPT_CHECK_MEM,8G,00:29:59,"$(LOAD_R_MODULE); $(FACETS_GENE_CN) $(FACETS_GENE_CN_OPTS) --outFile $@ $^")

facets/geneCN.fill.txt : facets/geneCN.txt $(foreach pair,$(SAMPLE_PAIRS),facets/cncf/$(pair).cncf.txt)
	$(call LSCRIPT_CHECK_MEM,4G,00:29:59,"$(LOAD_R_MODULE); $(FACETS_FILL_GENE_CN) --outFile $@ --geneCNFile $< \
		$(filter %.cncf.txt,$^)")

facets/geneCN%heatmap.pdf  : facets/geneCN%txt
	$(call LSCRIPT_MEM,4G,00:29:29,"$(LOAD_R_MODULE); $(FACETS_PLOT_GENE_CN) $(FACETS_PLOT_GENE_CN_OPTS) $< $@")

include usb-modules/variant_callers/gatk.mk
include usb-modules/bam_tools/processBam.mk

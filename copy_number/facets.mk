include usb-modules/Makefile.inc
include usb-modules/config.inc
include usb-modules/variant_callers/variantCaller.inc 

LOGDIR ?= log/facets.$(NOW)

.SECONDARY:
.DELETE_ON_ERROR:
.PHONY : facets

facets : $(foreach pair,$(SAMPLE_PAIRS),facets/cncf/$(pair).cncf.txt) facets/geneCN.txt
#	facets/geneCN.txt facets/geneCN.fill.txt facets/geneCN.heatmap.pdf facets/geneCN.fill.heatmap.pdf

facets/base_pos/%.gatk.vcf : vcf/%.gatk_snps.vcf
	$(INIT) ln $< $@

facets/base_pos/%.gatk.dbsnp.vcf : facets/base_pos/%.gatk.vcf $(DBSNP)
	$(call LSCRIPT_MEM,22G,03:59:59,"$(LOAD_JAVA8_MODULE); $(call COMBINE_VARIANTS,21G) \
	--variant $(word 1,$^) --variant $(word 2,$^) -o $@ --genotypemergeoption UNSORTED -R $(REF_FASTA)")

define snp-pileup-tumor-normal
facets/snp_pileup/$1_$2.bc.gz : bam/$1.bam bam/$2.bam $$(if $$(findstring true,$$(FACETS_GATK_VARIANTS)),facets/base_pos/$2.gatk.dbsnp.vcf,$$(DBSNP))
	$$(call LSCRIPT_CHECK_MEM,3G,00:59:59,"$$(FACETS_SNP_PILEUP) -A -d $$(FACETS_SNP_PILEUP_MAX_DEPTH) -g \
	-q $$(FACETS_SNP_PILEUP_MINMAPQ) -Q $$(FACETS_SNP_PILEUP_MINBASEQ) -r $$(FACETS_SNP_PILEUP_MIN_DEPTH) \
	$$(word 3,$$^) $$@ $$(word 1,$$^) $$(word 2,$$^) && $$(RM) $$(word 3,$$^)")
endef
$(foreach pair,$(SAMPLE_PAIRS),$(eval $(call snp-pileup-tumor-normal,$(tumor.$(pair)),$(normal.$(pair)))))

facets/cncf/%.cncf.txt : facets/snp_pileup/%.bc.gz
	$(call LSCRIPT_CHECK_MEM,3G,00:29:59,"$(LOAD_R_MODULE); $(FACETS) --minNDepth $(FACETS_SNP_PILEUP_MIN_DEPTH) --snp_nbhd $(FACETS_WINDOW_SIZE) \
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

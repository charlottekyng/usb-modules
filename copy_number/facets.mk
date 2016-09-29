include usb-modules/Makefile.inc
include usb-modules/config.inc

LOGDIR ?= log/facets.$(NOW)

.SECONDARY:
.DELETE_ON_ERROR:
.PHONY : facets

facets : $(foreach pair,$(SAMPLE_PAIRS),facets/cncf/$(pair).cncf.txt)
#facets/geneCN.txt
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

facets/vcf/targets_dbsnp.vcf.gz : $(TARGETS_FILE_INTERVALS)
	$(INIT) $(BEDTOOLS) intersect -header -u -a $(DBSNP) -b $< | gzip -c > $@

ifeq ($(findstring ILLUMINA,$(SEQ_PLATFORM)),ILLUMINA)
define snp-pileup-tumor-normal
facets/snp_pileup/$1_$2.bc.gz : bam/$1.bam bam/$2.bam
	$$(call LSCRIPT_CHECK_MEM,3G,00:59:59,"$$(LOAD_PERL_MODULE); $$(GET_BASE_COUNTS) bam/$2.bam bam/$1.bam $$@ \
	$$(REF_FASTA) $$(TARGETS_FILE_INTERVALS) $$(GET_BASE_COUNTS_MIN_DEPTH) $$(GET_BASE_COUNTS_MAX_DEPTH)")
endef
$(foreach pair,$(SAMPLE_PAIRS),$(eval $(call snp-pileup-tumor-normal,$(tumor.$(pair)),$(normal.$(pair)))))
endif

ifeq ($(findstring IONTORRENT,$(SEQ_PLATFORM)),IONTORRENT)
facets/snp_pileup/%/TSVC_variants.vcf : bam/%.bam
	$(call LSCRIPT_CHECK_MEM,20G,00:59:59,"$(TVC) -s $(DBSNP) -i $< -r $$(REF_FASTA) -o $(@D) -N 8 -t $(TVC_ROOT_DIR) --primer-trim-bed $(PRIMER_TRIM_BED)

define snp-pileup-tumor-normal
facets/snp_pileup/$1_$2.bc.gz : facets/snp_pileup/$2/TSVC_variants.vcf facets/snp_pileup/$1/TSVC_variants.vcf
	$$(call LSCRIPT_CHECK_MEM,20G,03:59:59,"$$(LOAD_SNP_EFF_MODULE); $$(LOAD_JAVA8_MODULE); $$(call COMBINE_VARIANTS,19G) \
	$$(foreach vcf,$$^,--variant $$(vcf) ) --genotypemergeoption UNSORTED -R $$(REF_FASTA) | \
	$$(SNP_SIFT) extractFields - CHROM POS REF ALT GEN[0].FRO GEN[0].FAO GEN[0].FXX GEN[0].DP GEN[1].FRO GEN[1].FAO GEN[1].FXX GEN[1].DP | \
	perl -p -e \"s/$$(,)[\w]+//g;\" | perl -p -e \"s/\#CHROM.+$$$$/Chromosome$$(,)Position$$(,)Ref$$(,)Alt$$(,)File1R$$(,)File1A$$(,)File1E$$(,)File1D$$(,)File2R$$(,)File2A$$(,)File2E$$(,)File2D/g;\" | \
	sed 's/\t/$$(,)/g;' | gzip > $$@")
endef
$(foreach pair,$(SAMPLE_PAIRS),$(eval $(call snp-pileup-tumor-normal,$(tumor.$(pair)),$(normal.$(pair)))))
endif

#facets/params.txt : $(foreach pair,$(SAMPLE_PAIRS),facets/snp_pileup/$(pair).bc.gz)
	

facets/cncf/%.cncf.txt : facets/snp_pileup/%.bc.gz
	$(call LSCRIPT_CHECK_MEM,3G,00:29:59,"$(LOAD_R_MODULE); $(FACETS) \
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
include usb-modules/bam_tools/processBam.mk

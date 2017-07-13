include usb-modules/Makefile.inc
include usb-modules/config.inc
include usb-modules/variant_callers/variantCaller.inc

LOGDIR ?= log/facets.$(NOW)

.SECONDARY:
.DELETE_ON_ERROR:
.PHONY : facets

facets : $(foreach cval1,$(FACETS_CVAL1),$(foreach pair,$(SAMPLE_PAIRS),facets/cncfTN_$(cval1)/$(pair).out) facets/cncfTN_$(cval1)/summary.txt facets/cncfTN_$(cval1)/geneCN.filled.txt)
#	facets/geneCN.txt facets/geneCN.fill.txt facets/geneCN.heatmap.pdf facets/geneCN.fill.heatmap.pdf

ifeq ($(findstring ILLUMINA,$(SEQ_PLATFORM)),ILLUMINA)
facets/base_pos/%.gatk.dbsnp.vcf : gatk/dbsnp/%.gatk_snps.vcf gatk/vcf/%.variants.vcf
	$(call LSCRIPT_MEM,$(RESOURCE_REQ_HIGHMEM),$(RESOURCE_REQ_SHORT),"$(LOAD_JAVA8_MODULE); \
	$(call GATK,CombineVariants,$(RESOURCE_REQ_HIGHMEM)) \
	$(foreach vcf,$^,--variant $(vcf) ) -o $@ --genotypemergeoption UNSORTED -R $(REF_FASTA)")

define snp-pileup-tumor-normal
facets/snp_pileup/$1_$2.bc.gz : bam/$1.bam bam/$2.bam $$(if $$(findstring true,$$(FACETS_GATK_VARIANTS)),facets/base_pos/$1.gatk.dbsnp.vcf,$$(FACETS_TARGETS_INTERVALS))
	$$(call LSCRIPT_CHECK_MEM,$$(RESOURCE_REQ_LOWMEM),$$(RESOURCE_REQ_SHORT),"$$(FACETS_SNP_PILEUP) \
	-A -d $$(FACETS_SNP_PILEUP_MAX_DEPTH) -g -q $$(FACETS_SNP_PILEUP_MINMAPQ) \
	-Q $$(FACETS_SNP_PILEUP_MINBASEQ) -r $$(FACETS_SNP_PILEUP_MIN_DEPTH)$$(,)0 \
	$$(word 3,$$^) $$@ $$(word 2,$$^) $$(word 1,$$^)")
endef
$(foreach pair,$(SAMPLE_PAIRS),$(eval $(call snp-pileup-tumor-normal,$(tumor.$(pair)),$(normal.$(pair)))))
endif

ifeq ($(findstring IONTORRENT,$(SEQ_PLATFORM)),IONTORRENT)

define snp-pileup-tumor-normal
facets/snp_pileup/$1_$2.bc.gz : tvc/dbsnp/$2/TSVC_variants.vcf tvc/dbsnp/$1/TSVC_variants.vcf
	$$(call LSCRIPT_CHECK_MEM,$$(RESOURCE_REQ_HIGHMEM),$$(RESOURCE_REQ_SHORT),"$$(LOAD_SNP_EFF_MODULE); $$(LOAD_JAVA8_MODULE); \
	$$(call GATK,CombineVariants,$$(RESOURCE_REQ_HIGHMEM)) \
	$$(foreach vcf,$$^,--variant $$(vcf) ) --genotypemergeoption UNSORTED -R $$(REF_FASTA) | \
	$$(SNP_SIFT) extractFields - CHROM POS REF ALT GEN[0].FRO GEN[0].FAO GEN[0].FXX GEN[0].FXX GEN[1].FRO GEN[1].FAO GEN[1].FXX GEN[1].FXX | \
	perl -p -e \"s/$$(,)[\w]+//g;\" | sed 's/^chr//g;' | sed 's/\t\t/\t0\t/g; s/\t$$$$/\t0/g; s/\t\t/\t0\t/g; s/\t\t/\t0\t/g;' | \
	perl -p -e \"s/\#CHROM.+$$$$/Chromosome$$(,)Position$$(,)Ref$$(,)Alt$$(,)File1R$$(,)File1A$$(,)File1E$$(,)File1D$$(,)File2R$$(,)File2A$$(,)File2E$$(,)File2D/g;\" | \
	grep -v '\.' | \
	sed 's/\t/$$(,)/g;' | gzip > $$@")
endef
$(foreach pair,$(SAMPLE_PAIRS),$(eval $(call snp-pileup-tumor-normal,$(tumor.$(pair)),$(normal.$(pair)))))
endif

define facets-cval1-tumor-normal
facets/cncfTN_$1/$2_$3.out : facets/snp_pileup/$2_$3.bc.gz
	$$(call LSCRIPT_CHECK_MEM,$$(RESOURCE_REQ_LOWMEM),$$(RESOURCE_REQ_VSHORT),"$$(LOAD_R_MODULE); \
	$$(FACETS) --minNDepth $$(FACETS_SNP_PILEUP_MIN_DEPTH) \
	--maxNDepth $$(FACETS_SNP_PILEUP_MAX_DEPTH) --snp_nbhd $$(FACETS_WINDOW_SIZE) --minGC $$(FACETS_MINGC) --maxGC $$(FACETS_MAXGC) \
	--cval1 $1 --genome $$(REF) --min_nhet $$(FACETS_MIN_NHET) \
	--outPrefix $$* $$< ")

facets/cncfTN_$1/$2_$3.Rdata : facets/cncfTN_$1/$2_$3.out
facets/cncfTN_$1/$2_$3.cncf.txt : facets/cncfTN_$1/$2_$3.out

endef
$(foreach cval1,$(FACETS_CVAL1),$(foreach pair,$(SAMPLE_PAIRS),$(eval $(call facets-cval1-tumor-normal,$(cval1),$(tumor.$(pair)),$(normal.$(pair))))))

define facets-merge
facets/cncfTN_$1/summary.txt : $$(foreach pair,$$(SAMPLE_PAIRS),facets/cncfTN_$1/$$(pair).out)
	$$(INIT) paste $$^ > $$@;

facets/cncfTN_$1/geneCN.filled.txt : $$(foreach pair,$$(SAMPLE_PAIRS),facets/cncfTN_$1/$$(pair).cncf.txt)
	$$(call LSCRIPT_CHECK_MEM,$$(RESOURCE_REQ_MEDIUM_MEM),$$(RESOURCE_REQ_SHORT),"$$(LOAD_R_MODULE); \
	$$(FACETS_GENE_CN) $$(FACETS_GENE_CN_OPTS) --outFile $$@ $$^")
endef
$(foreach cval1,$(FACETS_CVAL1),$(eval $(call facets-merge,$(cval1))))

facets/cncf/geneCN%heatmap.pdf  : facets/geneCN%txt
	$(call LSCRIPT_MEM,$(RESOURCE_REQ_LOWMEM),$(RESOURCE_REQ_VSHORT),"$(LOAD_R_MODULE); \
	$(FACETS_PLOT_GENE_CN) $(FACETS_PLOT_GENE_CN_OPTS) $< $@")

include usb-modules/variant_callers/gatk.mk
include usb-modules/variant_callers/TVC.mk
include usb-modules/bam_tools/processBam.mk
include usb-modules/vcf_tools/vcftools.mk
include usb-modules/variant_callers/gatkVariantCaller.mk

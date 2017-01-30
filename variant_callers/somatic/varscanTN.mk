# Run VarScan on tumour-normal matched pairs
# Detect point mutations
##### DEFAULTS ######

LOGDIR ?= log/varscanTN.$(NOW)

##### MAKE INCLUDES #####
include usb-modules/Makefile.inc
include usb-modules/config.inc
include usb-modules/variant_callers/somatic/somaticVariantCaller.inc

VPATH ?= bam

VARIANT_TYPES = varscan_indels varscan_snps

PHONY += varscan varscan_vcfs varscan_tables
varscan : varscan_vcfs varscan_tables
varscan_vcfs : $(foreach type,$(VARIANT_TYPES),$(call SOMATIC_VCFS,$(type)))
varscan_tables : $(foreach type,$(VARIANT_TYPES),$(call SOMATIC_TABLES,$(type)))

define varscan-somatic-tumor-normal-chr
varscan/chr_vcf/$1_$2.$3.varscan_timestamp : bam/$1.bam bam/$2.bam bam/$1.bam.bai bam/$2.bam.bai
	$$(call LSCRIPT_MEM,4G,00:29:59,"$(LOAD_SAMTOOLS_MODULE); $(LOAD_JAVA8_MODULE); \
	$$(SAMTOOLS) mpileup -r $3 -q $$(MIN_MQ) -f $$(REF_FASTA) $$(word 2,$$^) $$< | \
	$$(VARSCAN) somatic - varscan/chr_vcf/$1_$2.$3 $$(VARSCAN_OPTS) && \
	$$(VARSCAN) processSomatic varscan/chr_vcf/$1_$2.$3.snp.vcf --min-tumor-freq $$(MIN_AF_SNP) && \
	$$(VARSCAN) processSomatic varscan/chr_vcf/$1_$2.$3.indel.vcf --min-tumor-freq $$(MIN_AF_INDEL) && \
	$$(RM) varscan/chr_vcf/$1_$2.$3.snp.vcf varscan/chr_vcf/$1_$2.$3.indel.vcf && \
	$$(RM) varscan/chr_vcf/$1_$2.$3.indel.Germline* varscan/chr_vcf/$1_$2.$3.indel.LOH* && \
	$$(RM) varscan/chr_vcf/$1_$2.$3.snp.Germline* varscan/chr_vcf/$1_$2.$3.snp.LOH* && touch $$@")
varscan/chr_vcf/$1_$2.$3.snp.Somatic.hc.vcf : varscan/chr_vcf/$1_$2.$3.varscan_timestamp
varscan/chr_vcf/$1_$2.$3.indel.Somatic.hc.vcf : varscan/chr_vcf/$1_$2.$3.varscan_timestamp

varscan/chr_vcf/$1_$2.$3.snp.Somatic.hc.fpft.vcf : varscan/chr_vcf/$1_$2.$3.snp.Somatic.hc.vcf bam/$1.bam bam/$1.bam.bai
	$$(call LSCRIPT_MEM,4G,00:29:59,"awk '! /\#/' $$< | \
	awk '{if(length($$$$4) > length($$$$5)) print $$$$1\"\t\"($$$$2-1)\"\t\"($$$$2+length($$$$4)-1); \
	else print $$$$1\"\t\"($$$$2-1)\"\t\"($$$$2+length($$$$5)-1)}' > varscan/chr_vcf/$1_$2.$3.snp.Somatic.hc.vcf.region && \
	$$(BAM_READCOUNT) -f $$(REF_FASTA) -l varscan/chr_vcf/$1_$2.$3.snp.Somatic.hc.vcf.region $$(word 2,$$^) > \
	varscan/chr_vcf/$1_$2.$3.snp.Somatic.hc.vcf.bamrc && \
	$$(VARSCAN) fpfilter $$< varscan/chr_vcf/$1_$2.$3.snp.Somatic.hc.vcf.bamrc \
	--output-file $$@ --filtered-file varscan/chr_vcf/$1_$2.$3.snp.Somatic.hc.fpfail.vcf \
	--min-var-freq $$(MIN_AF_SNP) --min-ref-readpos 0 --min-var-readpos 0 --min-ref-dist3 0 --min-var-dist3 0 && \
	rm varscan/chr_vcf/$1_$2.$3.snp.Somatic.hc.vcf.region varscan/chr_vcf/$1_$2.$3.snp.Somatic.hc.vcf.bamrc")

varscan/chr_vcf/$1_$2.$3.indel.Somatic.hc.fpft.vcf : varscan/chr_vcf/$1_$2.$3.indel.Somatic.hc.vcf bam/$1.bam bam/$1.bam.bai
	$$(call LSCRIPT_MEM,4G,00:29:59,"awk '! /\#/' $$< | \
	awk '{if(length($$$$4) > length($$$$5)) print $$$$1\"\t\"($$$$2-1)\"\t\"($$$$2+length($$$$4)-1); \
	else print $$$$1\"\t\"($$$$2-1)\"\t\"($$$$2+length($$$$5)-1)}' > varscan/chr_vcf/$1_$2.$3.indel.Somatic.hc.vcf.region && \
	$$(BAM_READCOUNT) -f $$(REF_FASTA) -l varscan/chr_vcf/$1_$2.$3.indel.Somatic.hc.vcf.region $$(word 2,$$^) > \
	varscan/chr_vcf/$1_$2.$3.indel.Somatic.hc.vcf.bamrc && \
	$$(VARSCAN) fpfilter $$< varscan/chr_vcf/$1_$2.$3.indel.Somatic.hc.vcf.bamrc \
	--output-file $$@ --filtered-file varscan/chr_vcf/$1_$2.$3.indel.Somatic.hc.fpfail.vcf \
	--min-var-freq $$(MIN_AF_INDEL) --min-ref-readpos 0 --min-var-readpos 0 --min-ref-dist3 0 --min-var-dist3 0 && \
	rm varscan/chr_vcf/$1_$2.$3.indel.Somatic.hc.vcf.region varscan/chr_vcf/$1_$2.$3.indel.Somatic.hc.vcf.bamrc")

varscan/chr_vcf/$1_$2.%.fixed.vcf : varscan/chr_vcf/$1_$2.%.vcf
	$$(INIT) perl $$(FIX_VARSCAN_VCF) -t $1 -n $2 $$< > $$@

endef
$(foreach chr,$(CHROMOSOMES), \
	$(foreach pair,$(SAMPLE_PAIRS), \
	$(eval $(call varscan-somatic-tumor-normal-chr,$(tumor.$(pair)),$(normal.$(pair)),$(chr)))))

define varscan-somatic-tumor-normal-merge
varscan/vcf/$1_$2.snp.Somatic.hc.fpft.fixed.vcf : $$(foreach chr,$$(CHROMOSOMES),varscan/chr_vcf/$1_$2.$$(chr).snp.Somatic.hc.fpft.vcf)
	$$(INIT) $$(LOAD_PERL_MODULE); grep '^#' $$< > $$@.tmp; cat $$^ | grep -v '^#' >>$@.tmp; $$(VCF_SORT) $$(REF_DICT) $$@.tmp > $$@ 2> $$(LOG) && $$(RM) $$@.tmp

varscan/vcf/$1_$2.indel.Somatic.hc.fpft.fixed.vcf : $$(foreach chr,$$(CHROMOSOMES),varscan/chr_vcf/$1_$2.$$(chr).indel.Somatic.hc.fpft.vcf)
	$$(INIT) $$(LOAD_PERL_MODULE); grep '^#' $$< > $$@.tmp; cat $$^ | grep -v '^#' >>$@.tmp; $$(VCF_SORT) $$(REF_DICT) $$@.tmp > $$@ 2> $$(LOG) && $$(RM) $$@.tmp
endef
$(foreach pair,$(SAMPLE_PAIRS),\
	$(eval $(call varscan-somatic-tumor-normal-merge,$(tumor.$(pair)),$(normal.$(pair)))))

vcf/%.varscan_indels.vcf : varscan/vcf/%.indel.Somatic.hc.fpft.fixed.vcf
	$(INIT) ln -f $< $@

vcf/%.varscan_snps.vcf : varscan/vcf/%.snp.Somatic.hc.fpft.fixed.vcf
	$(INIT) ln -f $< $@

include usb-modules/vcf_tools/vcftools.mk

.DELETE_ON_ERROR:
.SECONDARY: 
.PHONY: $(PHONY)


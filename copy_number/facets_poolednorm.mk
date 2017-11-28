include usb-modules/Makefile.inc
include usb-modules/config.inc
include usb-modules/variant_callers/variantCaller.inc

LOGDIR ?= log/facets_poolednorm.$(NOW)

.SECONDARY:
.DELETE_ON_ERROR:
.PHONY : facets_poolednorm

facets_poolednorm : $(foreach cval1,$(FACETS_CVAL1),$(foreach sample,$(SAMPLES),facets/cncf_poolednorm_$(cval1)/$(sample)_poolednorm.out) facets/cncf_poolednorm_$(cval1)/summary.txt facets/cncf_poolednorm_$(cval1)/geneCN.filled.txt)

ifeq ($(findstring ILLUMINA,$(SEQ_PLATFORM)),ILLUMINA)
define snp-pileup-tumor-poolednorm
facets/snp_pileup/$1_poolednorm.bc.gz : bam/$1.bam bam/poolednorm.bam $$(if $$(findstring true,$$(FACETS_GATK_VARIANTS)),facets/base_pos/$1.gatk.dbsnp.vcf,$$(FACETS_TARGETS_INTERVALS)) bam/$1.bam.bai bam/poolednorm.bam.bai
	$$(call LSCRIPT_CHECK_MEM,$$(RESOURCE_REQ_LOWMEM),$$(RESOURCE_REQ_SHORT),"$$(FACETS_SNP_PILEUP) \
	-A -d $$(FACETS_SNP_PILEUP_MAX_DEPTH) -g -q $$(FACETS_SNP_PILEUP_MINMAPQ) \
	-Q $$(FACETS_SNP_PILEUP_MINBASEQ) -r $$(FACETS_SNP_PILEUP_MIN_DEPTH)$$(,)0 \
	$$(word 3,$$^) $$@ $$(word 2,$$^) $$(word 1,$$^)")
endef
$(foreach sample,$(SAMPLES),$(eval $(call snp-pileup-tumor-poolednorm,$(sample))))
endif

ifeq ($(findstring IONTORRENT,$(SEQ_PLATFORM)),IONTORRENT)
define snp-pileup-tumor-poolednorm
facets/snp_pileup/$1_poolednorm.bc.gz : tvc/dbsnp/$1/TSVC_variants.vcf tvc/dbsnp/poolednorm/TSVC_variants.vcf
	$$(call LSCRIPT_CHECK_MEM,$$(RESOURCE_REQ_HIGHMEM),$$(RESOURCE_REQ_SHORT),"$$(LOAD_SNP_EFF_MODULE); $$(LOAD_JAVA8_MODULE); \
	$$(call GATK,CombineVariants,$$(RESOURCE_REQ_HIGHMEM)) \
	$$(foreach vcf,$$^,--variant $$(vcf) ) --genotypemergeoption UNSORTED -R $$(REF_FASTA) | \
	$$(SNP_SIFT) extractFields - CHROM POS REF ALT GEN[0].FRO GEN[0].FAO GEN[0].FXX GEN[0].FXX GEN[1].FRO GEN[1].FAO GEN[1].FXX GEN[1].FXX | \
	perl -p -e \"s/$$(,)[\w]+//g;\" | sed 's/^chr//g;' | sed 's/\t\t/\t0\t/g; s/\t$$$$/\t0/g; s/\t\t/\t0\t/g; s/\t\t/\t0\t/g;' | \
	perl -p -e \"s/\#CHROM.+$$$$/Chromosome$$(,)Position$$(,)Ref$$(,)Alt$$(,)File1R$$(,)File1A$$(,)File1E$$(,)File1D$$(,)File2R$$(,)File2A$$(,)File2E$$(,)File2D/g;\" | \
	grep -v '\.' | \
	sed 's/\t/$$(,)/g;' | gzip > $$@")
endef
$(foreach sample,$(SAMPLES),$(eval $(call snp-pileup-tumor-poolednorm,$(sample))))
endif

define facets-cval1-sample
facets/cncf_poolednorm_$1/$2_poolednorm.out : facets/snp_pileup/$2_poolednorm.bc.gz
	$$(call LSCRIPT_CHECK_MEM,$$(RESOURCE_REQ_LOWMEM),$$(RESOURCE_REQ_VSHORT),"$$(LOAD_R_MODULE); \
	$$(FACETS) --minNDepth $$(FACETS_SNP_PILEUP_MIN_DEPTH) \
	--maxNDepth $$(FACETS_SNP_PILEUP_MAX_DEPTH) --snp_nbhd $$(FACETS_WINDOW_SIZE) --minGC $$(FACETS_MINGC) --maxGC $$(FACETS_MAXGC) --unmatched TRUE \
	--cval1 $1 --genome $$(REF) --min_nhet $$(FACETS_MIN_NHET) \
	--outPrefix $$* $$<")

facets/cncf_poolednorm_$1/$2_poolednorm.Rdata : facets/cncf_poolednorm_$1/$2_poolednorm.out
facets/cncf_poolednorm_$1/$2_poolednorm.cncf.txt : facets/cncf_poolednorm_$1/$2_poolednorm.out

endef
$(foreach cval1,$(FACETS_CVAL1),$(foreach sample,$(SAMPLES),$(eval $(call facets-cval1-sample,$(cval1),$(sample)))))

define facets-merge-poolednorm
facets/cncf_poolednorm_$1/summary.txt : $$(foreach sample,$$(SAMPLES),facets/cncf_poolednorm_$1/$$(sample)_poolednorm.out)
	$$(INIT) paste $$^ > $$@;

facets/cncf_poolednorm_$1/geneCN.filled.txt : $$(foreach sample,$$(SAMPLES),facets/cncf_poolednorm_$1/$$(sample)_poolednorm.cncf.txt)
	$$(call LSCRIPT_CHECK_MEM,$$(RESOURCE_REQ_MEDIUM_MEM),$$(RESOURCE_REQ_SHORT),"$$(LOAD_R_MODULE); \
	$$(FACETS_GENE_CN) $$(FACETS_GENE_CN_OPTS) --outFile $$@ $$^")
endef
$(foreach cval1,$(FACETS_CVAL1),$(eval $(call facets-merge-poolednorm,$(cval1))))

#--cval2 $(FACETS_CVAL2) --pre_cval $(FACETS_PRE_CVAL)
include usb-modules/copy_number/facets.mk
include usb-modules/variant_callers/gatk.mk
include usb-modules/variant_callers/TVC.mk
include usb-modules/bam_tools/processBam.mk
include usb-modules/bam_tools/poolednormBam.mk
include usb-modules/vcf_tools/vcftools.mk
include usb-modules/variant_callers/gatkVariantCaller.mk


# common GATK steps for variant calling
# Author: Raymond Lim <raylim@mm.st> & Fong Chun Chan <fongchunchan@gmail.com>
#
# GATK needs sorted, properly ordered bams
# INPUTS: sample bams
# OUTPUTS: vcf file for SNPs and indels
# OPTIONS: GATK_HARD_FILTER_SNPS = true/false (default: true)
# 		   GATK_POOL_SNP_RECAL = true/false (default: false)
# 		   SPLIT_CHR = true/false (default: true)
#

ifndef GATK_MK

include usb-modules/Makefile.inc
include usb-modules/config.inc

###### RECIPES #######

%.intervals : %.vcf
	$(INIT) sed '/^#/d' $< | awk '{print $$1":"$$2 }' > $@

#### if SPLIT_CHR is set to true, we will split gatk processing by chromosome
ifeq ($(SPLIT_CHR),true)
## call sample sets
ifdef SAMPLE_SET_PAIRS
define hapcall-vcf-sets-chr
gatk/chr_vcf/$1.$2.variants.vcf : $$(foreach sample,$$(samples.$1),gatk/chr_vcf/$$(sample).$2.variants.intervals) $$(foreach sample,$$(samples.$1),bam/$$(sample).bam bam/$$(sample).bai)
	$$(call LSCRIPT_CHECK_MEM,$$(RESOURCE_REQ_MEDIUM_MEM),$$(RESOURCE_REQ_SHORT),"$$(LOAD_JAVA8_MODULE); \
	$$(call GATK,HaplotypeCaller,$$(RESOURCE_REQ_MEDIUM_MEM)) $$(HAPLOTYPE_CALLER_OPTS) \
	$$(foreach bam,$$(filter %.bam,$$^),-I $$(bam) ) $$(foreach intervals,$$(filter %.intervals,$$^),-L $$(intervals) ) -o $$@")
endef
$(foreach chr,$(CHROMOSOMES),$(foreach set,$(SAMPLE_SET_PAIRS),$(eval $(call hapcall-vcf-sets-chr,$(set),$(chr)))))

define merge-chr-variants-sets
gatk/vcf/$1.variants.vcf : $$(foreach chr,$$(CHROMOSOMES),gatk/chr_vcf/$1.$$(chr).variants.vcf)
	$$(call LSCRIPT_CHECK_MEM,$$(RESOURCE_REQ_MEDIUM_MEM),$$(RESOURCE_REQ_SHORT),"$$(LOAD_JAVA8_MODULE); \
	$$(call GATK,CombineVariants,$$(RESOURCE_REQ_MEDIUM_MEM)) \
	--assumeIdenticalSamples $$(foreach i,$$^, --variant $$i) -R $$(REF_FASTA) -o $$@")
endef
$(foreach set,$(SAMPLE_SET_PAIRS),$(eval $(call merge-chr-variants-sets,$(set))))
endif # def SAMPLE_SET_PAIRS

define chr-variants
gatk/chr_vcf/%.$1.variants.vcf : bam/%.bam bam/%.bai
	$$(call LSCRIPT_CHECK_MEM,$$(RESOURCE_REQ_MEDIUM_MEM),$$(RESOURCE_REQ_SHORT),"$$(LOAD_JAVA8_MODULE); \
	$$(call GATK37,HaplotypeCaller,$$(RESOURCE_REQ_MEDIUM_MEM)) $$(HAPLOTYPE_CALLER_OPTS) \
	-L $1 -I $$< -o $$@")
endef
$(foreach chr,$(CHROMOSOMES),$(eval $(call chr-variants,$(chr))))

define merge-chr-variants
gatk/vcf/$1.variants.vcf : $$(foreach chr,$$(CHROMOSOMES),gatk/chr_vcf/$1.$$(chr).variants.vcf)
	$$(call LSCRIPT_CHECK_MEM,$$(RESOURCE_REQ_MEDIUM_MEM),$$(RESOURCE_REQ_SHORT),"$$(LOAD_JAVA8_MODULE); \
		$$(call GATK,CombineVariants,$$(RESOURCE_REQ_MEDIUM_MEM)) \
		--assumeIdenticalSamples $$(foreach i,$$^, --variant $$i) -R $$(REF_FASTA) -o $$@")
endef
$(foreach sample,$(SAMPLES),$(eval $(call merge-chr-variants,$(sample))))

else #### no splitting by chr ####
## call sample sets
ifdef SAMPLE_SETS
define hapcall-vcf-sets
gatk/vcf/$1.variants.vcf : $$(foreach sample,$$(samples.$1),gatk/vcf/$$(sample).variants.vcf) $$(foreach sample,$$(samples.$1),bam/$$(sample).bam bam/$$(sample).bai)
	$$(call LSCRIPT_CHECK_MEM,$$(RESOURCE_REQ_MEDIUM_MEM),00:59:59,"$$(LOAD_JAVA8_MODULE); \
	$$(call GATK,HaplotypeCaller,$$(RESOURCE_REQ_MEDIUM_MEM)) $$(HAPLOTYPE_CALLER_OPTS) \
	$$(foreach bam,$$(filter %.bam,$$^),-I $$(bam) ) $$(foreach vcf,$$(filter %.vcf,$$^),-L $$(vcf) ) -o $$@")
endef
$(foreach set,$(SAMPLE_SET_PAIRS),$(eval $(call hapcall-vcf-sets,$(set))))
endif

define hapcall-vcf
gatk/vcf/$1.variants.vcf : bam/$1.bam bam/$1.bai
	$$(call LSCRIPT_CHECK_MEM,$$(RESOURCE_REQ_MEDIUM_MEM),$$(RESOURCE_REQ_MEDIUM),"$$(LOAD_JAVA8_MODULE); \
	$$(call GATK,HaplotypeCaller,$$(RESOURCE_REQ_MEDIUM_MEM)) $$(HAPLOTYPE_CALLER_OPTS) -I $$< -o $$@")
endef
$(foreach sample,$(SAMPLES),$(eval $(call hapcall-vcf,$(sample))))

endif # split by chr

gatk/vcf/%.variants.snps.vcf : gatk/vcf/%.variants.vcf gatk/vcf/%.variants.vcf.idx
	$(call LSCRIPT_CHECK_MEM,$(RESOURCE_REQ_MEDIUM_MEM),$(RESOURCE_REQ_VSHORT),"$(LOAD_JAVA8_MODULE); \
	$(call GATK,SelectVariants,$(RESOURCE_REQ_MEDIUM_MEM)) \
	-R $(REF_FASTA) --variant $<  -o $@ -selectType SNP")

gatk/vcf/%.variants.indels.vcf : gatk/vcf/%.variants.vcf gatk/vcf/%.variants.vcf.idx
	$(call LSCRIPT_CHECK_MEM,$(RESOURCE_REQ_MEDIUM_MEM),$(RESOURCE_REQ_VSHORT),"$(LOAD_JAVA8_MODULE); \
	$(call GATK,SelectVariants,$(RESOURCE_REQ_MEDIUM_MEM)) \
	-R $(REF_FASTA) --variant $<  -o $@ -selectType INDEL")

gatk/dbsnp/%.gatk_snps.vcf : bam/%.bam bam/%.bai
	$(call LSCRIPT_PARALLEL_MEM,8,$(RESOURCE_REQ_LOWMEM),$(RESOURCE_REQ_SHORT),"$(LOAD_JAVA8_MODULE); \
	$(call GATK,UnifiedGenotyper,$(RESOURCE_REQ_LOWMEM)) \
	-nt 8 -R $(REF_FASTA) --dbsnp $(DBSNP_TARGETS_INTERVALS) $(foreach bam,$(filter %.bam,$<),-I $(bam) ) \
	--genotyping_mode GENOTYPE_GIVEN_ALLELES -alleles $(DBSNP_TARGETS_INTERVALS) -o $@ --output_mode EMIT_ALL_SITES")

$(REF_FASTA).fai : $(REF_FASTA)
	$(call LSCRIPT_CHECK_MEM,$(RESOURCE_REQ_LOWMEM),$(RESOURCE_REQ_MEDIUM),"$(LOAD_SAMTOOLS_MODULE); $(SAMTOOLS) faidx $<")

$(REF_FASTA:.fasta=.dict) : $(REF_FASTA)
	$(call LSCRIPT_CHECK_MEM,$(RESOURCE_REQ_LOWMEM),$(RESOURCE_REQ_VSHORT)),"$(LOAD_JAVA8_MODULE); \
		$(call PICARD,CreateSequenceDictionary,$(RESOURCE_REQ_LOWMEM)) REFERENCE=$< OUTPUT=$@")

define VARIANT_RECAL
	$$(call LSCRIPT_CHECK_PARALLEL_MEM,6,$$(RESOURCE_REQ_LOWMEM),$$(RESOURCE_REQ_VSHORT),"$$(LOAD_JAVA8_MODULE); \
	$$(call GATK,VariantRecalibrator,$$(RESOURCE_REQ_LOWMEM)) \
	-R $$(REF_FASTA) -nt 6 \
	-resource:hapmap$(,)known=false$(,)training=true$(,)truth=true$(,)prior=15.0 $$(HAPMAP) \
	-resource:omni$(,)known=false$(,)training=true$(,)truth=false$(,)prior=12.0 $$(OMNI) \
	-resource:dbsnp$(,)known=true$(,)training=false$(,)truth=false$(,)prior=8.0 $$(DBSNP) \
	$$(foreach i,$$(VARIANT_RECAL_ANNOTATIONS), -an $$i) \
	$$(foreach i,$$(filter %.vcf,$2), -input $$i) \
	-recalFile $1 -tranchesFile $$(basename $1).tranches -rscriptFile $$(basename $1).snps.plots.R")
endef

define APPLY_VARIANT_RECAL
	$$(call LSCRIPT_CHECK_MEM,$$(RESOURCE_REQ_MEDIUM_MEM),$$(RESOURCE_REQ_SHORT),"$$(LOAD_JAVA8_MODULE); \
	$$(call GATK,ApplyRecalibration,$$(RESOURCE_REQ_MEDIUM_MEM)) \
	-R $$(REF_FASTA) -input $2 -recalFile $3 \
	--ts_filter_level $$(VARIANT_RECAL_TRUTH_SENSITIVITY_LEVEL) \
	-tranchesFile $$(basename $3).tranches -o $1")
endef

ifeq ($(GATK_HARD_FILTER_SNPS),true)
gatk/vcf/%.variants.snps.filtered.vcf : gatk/vcf/%.variants.snps.vcf gatk/vcf/%.variants.snps.vcf.idx
	$(call LSCRIPT_CHECK_MEM,$(RESOURCE_REQ_MEDIUM_MEM),$(RESOURCE_REQ_VSHORT),"$(LOAD_JAVA8_MODULE); \
	$(call GATK,VariantFiltration,$(RESOURCE_REQ_MEDIUM_MEM)) -R $(REF_FASTA) $(SNP_FILTERS) -o $@ \
	--variant $<")
else 

# pool sample vcfs for recalibration
ifeq ($(GATK_POOL_SNP_RECAL),true)
gatk/vcf/samples.snps.recal.vcf : $(foreach sample,$(SAMPLES),gatk/vcf/$(sample).variants.snps.vcf) $(foreach sample,$(SAMPLES),gatk/vcf/$(sample).variants.snps.vcf.idx)
	$(call VARIANT_RECAL,$@,$^)

define sample-apply-recal
gatk/vcf/$1.variants.snps.filtered.vcf : gatk/vcf/$1.variants.snps.vcf gatk/vcf/samples.snps.recal.vcf gatk/vcf/samples.snps.recal.vcf.idx gatk/vcf/$1.variants.snps.vcf.idx 
	$$(call APPLY_VARIANT_RECAL,$$@,$$<,$$(word 2,$$^))
endef
$(foreach sample,$(SAMPLES),$(eval $(call sample-apply-recal,$(sample))))

ifdef SAMPLE_SETS
gatk/vcf/sets.snps.recal.vcf : $(foreach set,$(SAMPLE_SET_PAIRS),gatk/vcf/$(set).variants.snps.vcf gatk/vcf/$(set).variants.snps.vcf.idx )
	$(call VARIANT_RECAL,$@,$^)
			
define sets-apply-recal
gatk/vcf/$1.variants.snps.filtered.vcf : gatk/vcf/$1.variants.snps.vcf gatk/vcf/sets.snps.recal.vcf gatk/vcf/sets.snps.recal.vcf.idx gatk/vcf/$1.variants.snps.vcf.idx 
	$$(call APPLY_VARIANT_RECAL,$$@,$$<,$$(word 2,$$^))
endef
$(foreach set,$(SAMPLE_SET_PAIRS),$(eval $(call sets-apply-recal,$(set))))
endif
else 
gatk/vcf/%.variants.snps.recal.vcf : gatk/vcf/%.variants.snps.vcf gatk/vcf/%.variants.snps.vcf.idx
	$(call VARIANT_RECAL,$@,$^)

gatk/vcf/%.variants.snps.filtered.vcf : gatk/vcf/%.variants.snps.vcf gatk/vcf/%.variants.snps.recal.vcf gatk/vcf/%.variants.snps.vcf.idx gatk/vcf/%.variants.snps.recal.vcf.idx
	$(call APPLY_VARIANT_RECAL,$@,$<,$(word 2,$^))
endif
endif

# hard filter indels %=sample
gatk/vcf/%.variants.indels.filtered.vcf : gatk/vcf/%.variants.indels.vcf gatk/vcf/%.variants.indels.vcf.idx
	$(call LSCRIPT_CHECK_MEM,$(RESOURCE_REQ_MEDIUM_MEM),$(RESOURCE_REQ_VSHORT),"$(LOAD_JAVA8_MODULE); \
	$(call GATK,VariantFiltration,$(RESOURCE_REQ_MEDIUM_MEM)) -R $(REF_FASTA) $(INDEL_FILTERS) -o $@ \
	--variant $<")

# filter for only novel snps/indels
%.novel.txt : %.txt
	$(INIT) /bin/awk 'NR == 1 || $$4 == "."' $< > $@


vcf/%.gatk_snps.vcf : gatk/vcf/%.variants.snps.filtered.vcf
	$(INIT) $(FIX_GATK_VCF) $< > $@

vcf/%.gatk_indels.vcf : gatk/vcf/%.variants.indels.filtered.vcf
	$(INIT) $(FIX_GATK_VCF) $< > $@


reports/%/index.html : reports/%.dp_ft.grp metrics/hs_metrics.txt
	$(call LSCRIPT,"$(VARIANT_EVAL_GATK_REPORT) --metrics $(word 2,$^) --outDir $(@D) $<")


# merge variants 
include usb-modules/bam_tools/processBam.mk
include usb-modules/vcf_tools/vcftools.mk

endif
GATK_MK = true


##### DEFAULTS ######
LOGDIR = log/screen_hotspot.$(NOW)

##### MAKE INCLUDES #####
include usb-modules/Makefile.inc
include usb-modules/config.inc

VPATH ?= bam

.DELETE_ON_ERROR:
.SECONDARY: 
.PHONY : all

all : $(foreach pair,$(SAMPLE_PAIRS),hotspots/$(pair).snps.screened.target_ft.dp_ft.hotspot.pass.eff.vcf)

define hotspot-samplepairs
hotspots/$1_$2.snps.vcf : hotspots/$1.snps.vcf hotspots/$2.snps.vcf
	$$(call LSCRIPT_MEM,22G,03:59:59,"$$(LOAD_JAVA8_MODULE); $$(call COMBINE_VARIANTS,21G) \
		$$(foreach vcf,$$^,--variant $$(vcf) ) -o $$@ --genotypemergeoption UNSORTED -R $$(REF_FASTA)")
endef
$(foreach pair,$(SAMPLE_PAIRS),\
	$(eval $(call hotspot-samplepairs,$(tumor.$(pair)),$(normal.$(pair)))))

define hotspot-samplepairs-screen
hotspots/$1_$2.snps.screened.vcf : hotspots/$1_$2.snps.vcf
	$$(call LSCRIPT_MEM,8G,00:29:59,"$$(LOAD_JAVA8_MODULE); $$(call VARIANT_FILTRATION,7G) -R $$(REF_FASTA) -V $$< -o $$@ \
		--filterExpression 'vc.getGenotype(\"$1\").getAD().1 > 0 || vc.getGenotype(\"$1\").getAD().2 > 0 || vc.getGenotype(\"$1\").getAD().3 > 0' \
		--filterName presentInTumor")
endef
$(foreach pair,$(SAMPLE_PAIRS),\
	$(eval $(call hotspot-samplepairs-screen,$(tumor.$(pair)),$(normal.$(pair)))))

hotspots/%.snps.vcf : bam/%.bam hotspots/sites.to.screen.vcf
	$(call LSCRIPT_PARALLEL_MEM,4,5G,03:59:59,"$(LOAD_JAVA8_MODULE); $(call UNIFIED_GENOTYPER,4G) \
		-nt 4 -R $(REF_FASTA) --dbsnp $(DBSNP) $(foreach bam,$(filter %.bam,$<),-I $(bam) ) \
		--genotyping_mode GENOTYPE_GIVEN_ALLELES -alleles $(word 2,$^) -o $@ --output_mode EMIT_ALL_SITES")

ifndef $(TARGETS_FILE_INTERVALS)
hotspots/sites.to.screen.vcf : $(TARGETS_FILE_INTERVALS) $(CANCER_HOTSPOT_VCF)
	$(INIT) grep "^\#" $(CANCER_HOTSPOT_VCF) > $@ && \
	$(LOAD_BEDTOOLS_MODULE); $(BEDTOOLS) intersect -b $(TARGETS_FILE_INTERVALS) -a $(CANCER_HOTSPOT_VCF) -header >>$@
else
hotspots/sites.to.screen.vcf : $(CANCER_HOTSPOT_VCF)
	$(INIT) ln $(CANCER_HOTSPOT_VCF) $@
endif


include usb-modules/vcf_tools/vcftools.mk

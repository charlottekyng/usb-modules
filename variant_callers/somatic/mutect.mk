#### MAKE INCLUDES #####
include usb-modules/Makefile.inc
include usb-modules/config.inc
include usb-modules/variant_callers/somatic/somaticVariantCaller.inc

LOGDIR ?= log/mutect.$(NOW)

PHONY += mutect mutect_vcfs mutect_tables ext_output #mut_report
..DUMMY := $(shell mkdir -p version; echo "$(MUTECT) &> version/mutect.txt")

mutect : mutect_vcfs mutect_tables ext_output
mutect_vcfs : $(call SOMATIC_VCFS,mutect) $(addsuffix .idx,$(call SOMATIC_VCFS,mutect))
mutect_tables : $(call SOMATIC_TABLES,mutect)
ext_output : $(foreach pair,$(SAMPLE_PAIRS),mutect/tables/$(pair).mutect.txt)
#mut_report : mutect/report/index.html mutect/lowAFreport/index.html mutect/highAFreport/index.html

.DELETE_ON_ERROR:
.SECONDARY:
.PHONY : $(PHONY)

# run mutect on each chromosome
#$(call mutect-tumor-normal-chr,tumor,normal,chr)
define mutect-tumor-normal-chr
mutect/chr_vcf/$1_$2.$3.mutect%vcf mutect/chr_tables/$1_$2.$3.mutect%txt mutect/coverage/$1_$2.$3.mutect_cov%txt: bam/$1%bam bam/$2%bam
	$$(MKDIR) mutect/chr_tables mutect/chr_vcf mutect/coverage; $$(call LSCRIPT_CHECK_MEM,12G,01:59:59,"$$(LOAD_JAVA6_MODULE); $$(MUTECT) \
		--enable_extended_output $$(MUTECT_OPTS) --cosmic $$(COSMIC) \
		--intervals $3 --reference_sequence $$(REF_FASTA) --dbsnp $$(DBSNP) \
		--input_file:tumor $$< --input_file:normal $$(word 2,$$^) \
		--vcf mutect/chr_vcf/$1_$2.$3.mutect.vcf --out mutect/chr_tables/$1_$2.$3.mutect.txt \
		--coverage_file mutect/coverage/$1_$2.$3.mutect_cov.txt")
endef
$(foreach chr,$(CHROMOSOMES), \
	$(foreach pair,$(SAMPLE_PAIRS), \
			$(eval $(call mutect-tumor-normal-chr,$(tumor.$(pair)),$(normal.$(pair)),$(chr)))))

# merge variant tables 
define ext-mutect-tumor-normal
mutect/tables/$1.mutect.txt : $$(foreach chr,$$(CHROMOSOMES),mutect/chr_tables/$1.$$(chr).mutect.txt)
	$$(INIT) head -2 $$< > $$@; for table in $$^; do sed '1,2d' $$$$table >> $$@; done
endef
$(foreach pair,$(SAMPLE_PAIRS),$(eval $(call ext-mutect-tumor-normal,$(pair))))

#mutect/report/index.html: $(foreach pair,$(SAMPLE_PAIRS),mutect/tables/$(pair).mutect.txt)
#	$(call LSCRIPT_NAMED_MEM,mutect_report,6G,35G,"$(KNIT) $(MUT_FREQ_REPORT) $(@D) $^")

#mutect/lowAFreport/index.html: $(foreach pair,$(SAMPLE_PAIRS),mutect/tables/$(pair).mutect.txt)
#	$(call LSCRIPT_NAMED_MEM,mutect_lowaf_report,6G,35G,"$(KNIT) $(MUT_FREQ_REPORT) $(@D) --lowAF $^")

#mutect/highAFreport/index.html: $(foreach pair,$(SAMPLE_PAIRS),mutect/tables/$(pair).mutect.txt)
#	$(call LSCRIPT_NAMED_MEM,mutect_highaf_report,6G,35G,"$(KNIT) $(MUT_FREQ_REPORT) $(@D) --highAF $^")

# merge variants 
define mutect-tumor-normal
vcf/$1_$2.mutect.vcf : $$(foreach chr,$$(CHROMOSOMES),mutect/chr_vcf/$1_$2.$$(chr).mutect.vcf)
	$$(INIT) $$(LOAD_PERL_MODULE); grep '^#' $$< > $$@; cat $$^ | grep -v '^#' | $$(VCF_SORT) $$(REF_DICT) - >> $$@ 2> $$(LOG)
endef
$(foreach pair,$(SAMPLE_PAIRS),\
		$(eval $(call mutect-tumor-normal,$(tumor.$(pair)),$(normal.$(pair)))))

include usb-modules/vcf_tools/vcftools.mk
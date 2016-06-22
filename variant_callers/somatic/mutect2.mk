#### MAKE INCLUDES #####
include usb-modules/Makefile.inc
include usb-modules/config.inc
include usb-modules/variant_callers/somatic/somaticVariantCaller.inc

LOGDIR ?= log/mutect2.$(NOW)

PHONY += mutect2 #mutect_vcfs mutect_tables ext_output #mut_report
..DUMMY := $(shell mkdir -p version; echo "$(MUTECT2) &> version/mutect2.txt")

mutect2 : pon mutect2_vcfs mutect_tables #ext_output
pon : $(foreach chr,$(CHROMOSOMES),mutect2/chr_vcf_pon/pon.$(chr).mutect2.vcf)
mutect2_vcfs : $(call SOMATIC_VCFS,mutect2) $(addsuffix .idx,$(call SOMATIC_VCFS,mutect2))
mutect_tables : $(call SOMATIC_TABLES,mutect2)
#ext_output : $(foreach pair,$(SAMPLE_PAIRS),mutect/tables/$(pair).mutect.txt)
#mut_report : mutect/report/index.html mutect/lowAFreport/index.html mutect/highAFreport/index.html

.DELETE_ON_ERROR:
.SECONDARY:
.PHONY : $(PHONY)

define mutect2-normal-chr
mutect2/chr_vcf_pon/$1.$2.mutect2.vcf : bam/$1.bam
	$$(MKDIR) mutect2/chr_vcf_pon; $$(call LSCRIPT_CHECK_MEM,12G,09:59:59,"$$(LOAD_JAVA8_MODULE); $$(MUTECT2) \
		--reference_sequence $$(REF_FASTA) --input_file:tumor $$< --artifact_detection_mode \
		--dbsnp $$(DBSNP) --cosmic $$(COSMIC) --intervals $2 \
		--annotation TandemRepeatAnnotator --annotation OxoGReadCounts \
		--out mutect2/chr_vcf_pon/$1.$2.mutect2.vcf")
endef
$(foreach chr,$(CHROMOSOMES), \
	$(foreach normal,$(NORMAL_SAMPLES), \
		$(eval $(call mutect2-normal-chr,$(normal),$(chr)))))

define mutect2-pon-chr
mutect2/chr_vcf_pon/pon.$1.mutect2.vcf : $$(foreach normal,$$(NORMAL_SAMPLES),mutect2/chr_vcf_pon/$$(normal).$1.mutect2.vcf)
	$$(call LSCRIPT_CHECK_MEM,8G,00:29:59,"$$(LOAD_JAVA8_MODULE); $$(COMBINE_VARIANTS) \
	--reference_sequence $$(REF_FASTA) -minN 2 --setKey \"null\" --filteredAreUncalled --filteredrecordsmergetype KEEP_IF_ANY_UNFILTERED \
	$$(foreach normal,$$(NORMAL_SAMPLES),--variant mutect2/chr_vcf_pon/$$(normal).$1.mutect2.vcf) --intervals $1 \
	--out mutect2/chr_vcf_pon/pon.$1.mutect2.vcf")
endef
$(foreach chr,$(CHROMOSOMES), \
	$(eval $(call mutect2-pon-chr,$(chr))))

define mutect2-tumor-normal-chr
mutect2/chr_vcf/$1_$2.$3.mutect2%vcf : bam/$1%bam bam/$2%bam mutect2/chr_vcf_pon/pon.$3.mutect2.vcf
	$$(MKDIR) mutect2/chr_vcf; $$(call LSCRIPT_CHECK_MEM,12G,09:59:59,"$$(LOAD_JAVA8_MODULE); $$(MUTECT2) \
		--reference_sequence $$(REF_FASTA) --input_file:tumor $$< --input_file:normal $$(word 2,$$^) \
		--dbsnp $$(DBSNP) --cosmic $$(COSMIC) --intervals $3 $$(MUTECT_OPTS) \
		--annotation TandemRepeatAnnotator --annotation OxoGReadCounts --normal_panel $$(word 3,$$^) \
		--out mutect2/chr_vcf/$1_$2.$3.mutect2.vcf")
endef
$(foreach chr,$(CHROMOSOMES), \
	$(foreach pair,$(SAMPLE_PAIRS), \
			$(eval $(call mutect2-tumor-normal-chr,$(tumor.$(pair)),$(normal.$(pair)),$(chr)))))

# merge variants 
define mutect2-tumor-normal
vcf/$1_$2.mutect2.vcf : $$(foreach chr,$$(CHROMOSOMES),mutect2/chr_vcf/$1_$2.$$(chr).mutect2.vcf)
	$$(INIT) $$(LOAD_PERL_MODULE); grep '^#' $$< > $$@; cat $$^ | grep -v '^#' | $$(VCF_SORT) $$(REF_DICT) - >> $$@ 2> $$(LOG)
endef
$(foreach pair,$(SAMPLE_PAIRS),\
		$(eval $(call mutect2-tumor-normal,$(tumor.$(pair)),$(normal.$(pair)))))

include usb-modules/vcf_tools/vcftools.mk

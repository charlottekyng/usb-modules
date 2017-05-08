include usb-modules/Makefile.inc
include usb-modules/config.inc
include usb-modules/variant_callers/somatic/somaticVariantCaller.inc

LOGDIR ?= log/pon.$(NOW)

#PHONY : pon

#pon : $(foreach chr,$(CHROMOSOMES),mutect2/chr_vcf_pon/pon.$(chr).mutect2.vcf)

define mutect2-pon-chr
mutect2/chr_vcf_pon/$1.$2.mutect2.vcf : bam/$1.bam
	$$(MKDIR) mutect2/chr_vcf_pon; $$(call LSCRIPT_CHECK_MEM,12G,09:59:59,"$$(LOAD_JAVA8_MODULE); $$(call MUTECT2,18G) \
		--reference_sequence $$(REF_FASTA) --input_file:tumor $$< --artifact_detection_mode \
		--dbsnp $$(DBSNP) $$(if $$(findstring GRCm38,$$(REF)),,--cosmic $$(COSMIC)) --intervals $2 \
		--annotation TandemRepeatAnnotator --annotation OxoGReadCounts \
		--out mutect2/chr_vcf_pon/$1.$2.mutect2.vcf")
endef
$(foreach chr,$(CHROMOSOMES), \
	$(foreach normal,$(NORMAL_SAMPLES), \
		$(eval $(call mutect2-pon-chr,$(normal),$(chr)))))

define mutect2-pon-chr
mutect2/chr_vcf_pon/pon.$1.mutect2.vcf : $$(foreach normal,$$(NORMAL_SAMPLES),mutect2/chr_vcf_pon/$$(normal).$1.mutect2.vcf)
	$$(call LSCRIPT_CHECK_MEM,20G,03:29:59,"$$(LOAD_JAVA8_MODULE); $$(call COMBINE_VARIANTS,19G) \
		--reference_sequence $$(REF_FASTA) -minN 2 --setKey \"null\" --filteredAreUncalled --filteredrecordsmergetype KEEP_IF_ANY_UNFILTERED \
		$$(foreach normal,$$(NORMAL_SAMPLES),--variant mutect2/chr_vcf_pon/$$(normal).$1.mutect2.vcf) --intervals $1 \
		--out $$@")
endef
$(foreach chr,$(CHROMOSOMES), \
	$(eval $(call mutect2-pon-chr,$(chr))))

mutect2/pon.mutect2.vcf : $(foreach chr,$(CHROMOSOMES),mutect2/chr_vcf_pon/pon.$(chr).mutect2.vcf)
	$(call LSCRIPT_CHECK_MEM,20G,03:29:59,"$(LOAD_JAVA8_MODULE); $(call COMBINE_VARIANTS,19G) \
	--reference_sequence $(REF_FASTA) $(foreach chr,$(CHROMOSOMES),--variant mutect2/chr_vcf_pon/pon.$(chr).mutect2.vcf) \
	--genotypemergeoption UNIQUIFY --out $@")

tvc/pon.tvc.vcf : $(foreach sample,$(NORMAL_SAMPLES),tvc/vcf_pon/$(sample)/TSVC_variants.vcf)
	$(call LSCRIPT_CHECK_MEM,20G,03:29:59,"$(LOAD_JAVA8_MODULE); $(call COMBINE_VARIANTS,8G) \
	--reference_sequence $(REF_FASTA) $(foreach sample,$(NORMAL_SAMPLES),--variant tvc/vcf_pon/$(sample)/TSVC_variants.vcf) \
	--genotypemergeoption UNIQUIFY --out $@")

tvc/vcf_pon/%/TSVC_variants.vcf : bam/%.bam bam/%.bam.bai
	$(call LSCRIPT_PARALLEL_MEM,4,10G,05:59:59,"$(TVC) -i $< -r $(REF_FASTA) -o $(@D) -N 4 \
	$(if $(TARGETS_FILE_INTERVALS),-b $(TARGETS_FILE_INTERVALS)) -m $(TVC_MOTIF) $(TVC_SENSITIVE_JSON) \
	-t $(TVC_ROOT_DIR) --primer-trim-bed $(PRIMER_TRIM_BED)")

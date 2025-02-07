include usb-modules/Makefile.inc
include usb-modules/config.inc
include usb-modules/variant_callers/somatic/somaticVariantCaller.inc

LOGDIR ?= log/pon.$(NOW)

#PHONY : pon

#pon : $(foreach chr,$(CHROMOSOMES),mutect2/chr_vcf_pon/pon.$(chr).mutect2.vcf)

define mutect2-pon-chr
mutect2/chr_vcf_pon/$1.$2.mutect2.vcf : bam/$1.bam
	$$(MKDIR) mutect2/chr_vcf_pon; $$(call LSCRIPT_CHECK_MEM,$$(RESOURCE_REQ_HIGHMEM),$$(RESOURCE_REQ_LONG),"$$(LOAD_JAVA8_MODULE); \
		$$(call GATK,MuTect2,$$(RESOURCE_REQ_HIGHMEM)) \
		--reference_sequence $$(REF_FASTA) --input_file:tumor $$< --artifact_detection_mode \
		--dbsnp $$(DBSNP) $$(if $$(findstring GRCm38,$$(REF)),,--cosmic $$(COSMIC)) --intervals $2 \
		--annotation TandemRepeatAnnotator --annotation OxoGReadCounts \
		--out mutect2/chr_vcf_pon/$1.$2.mutect2.vcf")
endef
$(foreach chr,$(CHROMOSOMES), \
	$(foreach normal,$(PANEL_OF_NORMAL_SAMPLES), \
		$(eval $(call mutect2-pon-chr,$(normal),$(chr)))))

define mutect2-pon-chr-merge
mutect2/chr_vcf_pon/pon.$1.mutect2.vcf : $$(foreach normal,$$(PANEL_OF_NORMAL_SAMPLES),mutect2/chr_vcf_pon/$$(normal).$1.mutect2.vcf)
	$$(call LSCRIPT_CHECK_MEM,$$(RESOURCE_REQ_HIGHMEM),$$(RESOURCE_REQ_SHORT),"$$(LOAD_JAVA8_MODULE); \
	$$(call GATK,CombineVariants,$$(RESOURCE_REQ_HIGHMEM)) --reference_sequence $(REF_FASTA) \
	-minN 2 --setKey \"null\" --filteredAreUncalled --filteredrecordsmergetype KEEP_IF_ANY_UNFILTERED \
	$$(foreach normal,$$(NORMAL_SAMPLES),--variant mutect2/chr_vcf_pon/$$(normal).$1.mutect2.vcf) --intervals $1 \
	--out $$@")
endef
$(foreach chr,$(CHROMOSOMES), \
	$(eval $(call mutect2-pon-chr-merge,$(chr))))

mutect2/pon.mutect2.vcf : $(foreach chr,$(CHROMOSOMES),mutect2/chr_vcf_pon/pon.$(chr).mutect2.vcf)
	$(call LSCRIPT_CHECK_MEM,$(RESOURCE_REQ_HIGHMEM),$(RESOURCE_REQ_SHORT),"$(LOAD_JAVA8_MODULE); \
	$(call GATK,CombineVariants,$(RESOURCE_REQ_HIGHMEM)) --reference_sequence $(REF_FASTA) \
	$(foreach chr,$(CHROMOSOMES),--variant mutect2/chr_vcf_pon/pon.$(chr).mutect2.vcf) \
	--genotypemergeoption UNIQUIFY --out $@")

tvc/pon.tvc.vcf : $(foreach sample,$(PANEL_OF_NORMAL_SAMPLES),tvc/vcf_pon/$(sample)/TSVC_variants.vcf)
	$(call LSCRIPT_CHECK_MEM,$(RESOURCE_REQ_HIGHMEM),$(RESOURCE_REQ_SHORT),"$(LOAD_JAVA8_MODULE); \
	$(call GATK,CombineVariants,$(RESOURCE_REQ_HIGHMEM)) --reference_sequence $(REF_FASTA) \
	$(foreach sample,$(NORMAL_SAMPLES),--variant tvc/vcf_pon/$(sample)/TSVC_variants.vcf) \
	--genotypemergeoption UNIQUIFY --out $@")

tvc/vcf_pon/%/TSVC_variants.vcf : bam/%.bam bam/%.bam.bai
	$(call LSCRIPT_PARALLEL_MEM,4,$(RESOURCE_REQ_MEDIUM_MEM),$(RESOURCE_REQ_SHORT),"$(TVC) \
	-i $< -r $(REF_FASTA) -o $(@D) -N 4 \
	$(if $(TARGETS_FILE_INTERVALS),-b $(TARGETS_FILE_INTERVALS)) -m $(TVC_MOTIF) $(TVC_SOMATIC_JSON) \
	-t $(TVC_ROOT_DIR) --primer-trim-bed $(PRIMER_TRIM_BED)")

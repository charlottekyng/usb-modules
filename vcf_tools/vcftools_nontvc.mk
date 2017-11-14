%.altad_ft.vcf : %.vcf
	$(call LSCRIPT_CHECK_MEM,$(RESOURCE_REQ_MEDIUM_MEM),$(RESOURCE_REQ_SHORT),"$(LOAD_JAVA8_MODULE); \
		$(call GATK,VariantFiltration,$(RESOURCE_REQ_MEDIUM_MEM)) \
		-R $(REF_FASTA) -V $< -o $@ --filterExpression 'vc.getGenotype(\"$1\").getAD().1 == 0' --filterName zeroAD && $(RM) $< $<.idx")

%.dp_ft.vcf : %.vcf
	$(call LSCRIPT_CHECK_MEM,$(RESOURCE_REQ_MEDIUM_MEM),$(RESOURCE_REQ_SHORT),"$(LOAD_JAVA8_MODULE); \
		$(call GATK,VariantFiltration,$(RESOURCE_REQ_MEDIUM_MEM)) \
		-R $(REF_FASTA) -V $< -o $@ --filterExpression 'DP < $(MIN_NORMAL_DEPTH) || DP == \".\"' --filterName DepthO")


ifdef SAMPLE_PAIRS
define som-ad-ft-tumor-normal
vcf/$1_$2.%.som_ad_ft.vcf : vcf/$1_$2.%.vcf
	$$(call LSCRIPT_CHECK_MEM,$$(RESOURCE_REQ_MEDIUM_MEM),$$(RESOURCE_REQ_SHORT),"$$(LOAD_JAVA8_MODULE); \
		$$(call GATK,VariantFiltration,$$(RESOURCE_REQ_MEDIUM_MEM)) -R $$(REF_FASTA) -V $$< -o $$@ \
		--filterExpression 'vc.getGenotype(\"$1\").getAD().1 < $(MIN_TUMOR_AD)' \
		--filterName tumorVarAlleleDepth \
		--filterExpression 'if (vc.getGenotype(\"$2\").getDP() > $(MIN_NORMAL_DEPTH)) { \
		( vc.getGenotype(\"$2\").getAD().1 * 1.0 / vc.getGenotype(\"$2\").getDP()) > ( vc.getGenotype(\"$1\").getAD().1 * 1.0 / vc.getGenotype(\"$1\").getDP()) / $(MIN_TN_AD_RATIO) \
			} else { ( vc.getGenotype(\"$2\").getAD().1 * 1.0 / vc.getGenotype(\"$2\").getDP()) > ( vc.getGenotype(\"$1\").getAD().1 * 1.0 / vc.getGenotype(\"$1\").getDP()) / $(MIN_TN_AD_RATIO) && \
			vc.getGenotype(\"$1\").getAD().1 * 1.0 < 1 && vc.getGenotype(\"$2\").getAD().1 > 1 }' \
		--filterName somaticAlleleDepth \
		--filterExpression 'vc.getGenotype(\"$1\").getDP() <= $$(MIN_TUMOR_DEPTH) || vc.getGenotype(\"$2\").getDP() <= $$(MIN_NORMAL_DEPTH)' \
		--filterName depthFilter && sed -i 's/getGenotype(\"\([^\"]*\)\")/getGenotype(\1)/g' $$@ && $$(RM) $$< $$<.idx")
endef
$(foreach pair,$(SAMPLE_PAIRS),$(eval $(call som-ad-ft-tumor-normal,$(tumor.$(pair)),$(normal.$(pair)))))

define sufam
ifeq ($3,$2)
vcf/$1_$2.%.sufam.vcf : vcf/$1_$2.%.vcf
	$$(INIT) ln -f $$< $$@
else
vcf/$3.%.sufam.tmp : $$(foreach tumor,$$(wordlist 1,$$(shell expr $$(words $$(subst _,$$( ),$3)) - 1),$$(subst _,$$( ),$3)),vcf/$$(tumor)_$$(lastword $$(subst _,$$( ),$3)).%.vcf)
	$$(call LSCRIPT_MEM,$$(RESOURCE_REQ_HIGHMEM),$$(RESOURCE_REQ_SHORT),"$$(LOAD_JAVA8_MODULE); \
		$$(call GATK,CombineVariants),$$(RESOURCE_REQ_HIGHMEM)) \
		$$(foreach vcf,$$^,--variant $$(vcf) ) -o $$@ --genotypemergeoption UNSORTED -R $$(REF_FASTA)")

ifeq ($$(findstring varscan,$$(MUT_CALLER)),varscan)
vcf/$1_$2.%.sufam.vcf : vcf/$1_$2.%.vcf vcf/$3.%.sufam.tmp bam/$1.bam bam/$2.bam
	$$(call CHECK_VCF_CMD,$$(word 2,$$^),cp $$(word 1,$$^) $$@,
		$$(call LSCRIPT_CHECK_MEM,$$(RESOURCE_REQ_HIGHMEM),$$(RESOURCE_REQ_SHORT),"$$(LOAD_SNP_EFF_MODULE); $$(LOAD_JAVA8_MODULE); $(LOAD_SAMTOOLS_MODULE); \
		$$(call GATK,SelectVariants,$$(RESOURCE_REQ_MEDIUM_MEM)) -R $$(REF_FASTA) \
		--variant $$(word 2,$$^) --discordance $$(word 1,$$^) | grep -v \"#\" | \
		awk 'BEGIN { OFS=\"\t\"} { print \$$$$1$$(,)\$$$$2-10$$(,)\$$$$2+10 }' > $$@.tmp1 && \
		$$(SAMTOOLS) mpileup -l $$@.tmp1 -d 10000 -q $$(MIN_MQ) -f $$(REF_FASTA) $$(word 4,$$^) $$(word 3,$$^) | \
		$$(VARSCAN) somatic - $$@.tmp2 --output-vcf 1 --mpileup 1 --validation && \
		$$(FIX_VARSCAN_VCF) -t $1 -n $2 $$@.tmp2.validation > $$@.tmp3 && \
		$$(call GATK,SelectVariants,$$(RESOURCE_REQ_HIGHMEM)) -R $$(REF_FASTA) \
		--variant $$@.tmp3 --concordance $$(word 2,$$^) > $$@.tmp4 && \
		$$(call GATK,VariantFiltration,$$(RESOURCE_REQ_MEDIUM_MEM)) -R $$(REF_FASTA) -V $$@.tmp4 -o $$@.tmp5 \
		--filterExpression 'vc.getGenotype(\"$1\").getAD().1 * 1.0 > 0' --filterName interrogation && \
		$$(SNP_SIFT) filter $$(SNP_SIFT_OPTS) -f $$@.tmp5 \"(FILTER has 'interrogation')\" > $$@.tmp6 && \
		$$(call GATK,CombineVariants,$$(RESOURCE_REQ_HIGHMEM)) --variant $$< --variant $$@.tmp6 -o $$@ \
		--genotypemergeoption UNSORTED -R $$(REF_FASTA) && \
		$$(RM) $$@.tmp* $$(word 2,$$^) $$(word 2,$$^).idx"))
else
vcf/$1_$2.%.sufam.vcf : vcf/$1_$2.%.vcf vcf/$3.%.sufam.tmp bam/$1.bam bam/$2.bam
	$$(call CHECK_VCF_CMD,$$(word 2,$$^),cp $$(word 1,$$^) $$@,\
		$$(call LSCRIPT_CHECK_MEM,$$(RESOURCE_REQ_HIGHMEM),$$(RESOURCE_REQ_SHORT),"$$(LOAD_SNP_EFF_MODULE); $$(LOAD_JAVA8_MODULE); \
		$$(call GATK,SelectVariants,$$(RESOURCE_REQ_MEDIUM_MEM)) -R $$(REF_FASTA) \
		--variant $$(word 2,$$^) --discordance $$(word 1,$$^) -o $$@.tmp1 && \
		$$(call GATK,HaplotypeCaller,$$(RESOURCE_REQ_MEDIUM_MEM)) -R $$(REF_FASTA) -I $$(word 3,$$^) -I $$(word 4,$$^) \
		--downsampling_type NONE --dbsnp $(DBSNP_TARGETS_INTERVALS) \
		--genotyping_mode GENOTYPE_GIVEN_ALLELES --output_mode EMIT_ALL_SITES -alleles $$@.tmp1 -o $$@.tmp2 && \
		$$(FIX_GATK_VCF) $$@.tmp2 > $$@.tmp3 && mv $$@.tmp3 $$@.tmp2 && \
		$$(call GATK,VariantFiltration,$$(RESOURCE_REQ_MEDIUM_MEM)) -R $$(REF_FASTA) -V $$@.tmp2 -o $$@.tmp3 \
		--filterExpression 'vc.getGenotype(\"$1\").getAD().1 * 1.0 > 0' --filterName interrogation && \
		$$(SNP_SIFT) filter $$(SNP_SIFT_OPTS) -f $$@.tmp3 \"(FILTER has 'interrogation')\" > $$@.tmp4 && \
		$$(call GATK,CombineVariants,$$(RESOURCE_REQ_HIGHMEM)) --variant $$< --variant $$@.tmp4 -o $$@ \
		--genotypemergeoption UNSORTED -R $$(REF_FASTA) && \
		$$(RM) $$@.tmp1 $$@.tmp2 $$@.tmp3 $$@.tmp4 $$(word 2,$$^) $$@.tmp1.idx $$@.tmp2.idx $$@.tmp3.idx $$@.tmp4.idx $$(word 2,$$^).idx"))
endif
endif #ifeq ($3,$2)
endef #define sufam

$(foreach set,$(SAMPLE_SETS),\
	$(foreach tumor,$(wordlist 1,$(shell expr $(words $(subst _,$( ),$(set))) - 1),$(subst _,$( ),$(set))),\
		$(eval $(call sufam,$(tumor),$(lastword $(subst _,$( ),$(set))),$(subst $(tumor)_,,$(set))))))
endif #ifdef SAMPLE_PAIRS

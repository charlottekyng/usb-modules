ifdef SAMPLE_PAIRS
define som-ad-ft-tumor-normal
vcf/$1_$2.%.som_ad_ft.vcf : vcf/$1_$2.%.vcf
	$$(call LSCRIPT_CHECK_MEM,8G,00:59:59,"$$(LOAD_JAVA8_MODULE); $$(call VARIANT_FILTRATION,7G) -R $$(REF_FASTA) -V $$< -o $$@ \
		--filterExpression 'vc.getGenotype(\"$1\").getAnyAttribute(\"FAO\") < $(MIN_TUMOR_AD)' \
		--filterName tumorVarAlleleDepth \
		--filterExpression 'if ((vc.getGenotype(\"$2\").getAnyAttribute(\"FAO\") * 1.0 + vc.getGenotype(\"$2\").getAnyAttribute(\"FRO\") * 1.0) > $(MIN_NORMAL_DEPTH)) { \
		( vc.getGenotype(\"$2\").getAnyAttribute(\"FAO\") * 1.0 / (vc.getGenotype(\"$2\").getAnyAttribute(\"FAO\") * 1.0 + vc.getGenotype(\"$2\").getAnyAttribute(\"FRO\") * 1.0) ) \
			> ( vc.getGenotype(\"$1\").getAnyAttribute(\"FAO\") * 1.0 / (vc.getGenotype(\"$1\").getAnyAttribute(\"FAO\") * 1.0 + vc.getGenotype(\"$1\").getAnyAttribute(\"FRO\") * 1.0)) / $(MIN_TN_AD_RATIO) && \
				( vc.getGenotype(\"$2\").getAnyAttribute(\"AO\") * 1.0 / (vc.getGenotype(\"$2\").getAnyAttribute(\"AO\") * 1.0 + vc.getGenotype(\"$2\").getAnyAttribute(\"RO\") * 1.0) ) \
				> ( vc.getGenotype(\"$1\").getAnyAttribute(\"AO\") * 1.0 / (vc.getGenotype(\"$1\").getAnyAttribute(\"AO\") * 1.0 + vc.getGenotype(\"$1\").getAnyAttribute(\"RO\") * 1.0)) / $(MIN_TN_AD_RATIO) \
			} else {  ( vc.getGenotype(\"$2\").getAnyAttribute(\"FAO\") * 1.0 / (vc.getGenotype(\"$2\").getAnyAttribute(\"FAO\") * 1.0 + vc.getGenotype(\"$2\").getAnyAttribute(\"FRO\") * 1.0) ) \
				> ( vc.getGenotype(\"$1\").getAnyAttribute(\"FAO\") * 1.0 / (vc.getGenotype(\"$1\").getAnyAttribute(\"FAO\") * 1.0 + vc.getGenotype(\"$1\").getAnyAttribute(\"FRO\") * 1.0)) / $(MIN_TN_AD_RATIO) && \
				( vc.getGenotype(\"$2\").getAnyAttribute(\"AO\") * 1.0 / (vc.getGenotype(\"$2\").getAnyAttribute(\"AO\") * 1.0 + vc.getGenotype(\"$2\").getAnyAttribute(\"RO\") * 1.0) ) \
				> ( vc.getGenotype(\"$1\").getAnyAttribute(\"AO\") * 1.0 / (vc.getGenotype(\"$1\").getAnyAttribute(\"AO\") * 1.0 + vc.getGenotype(\"$1\").getAnyAttribute(\"RO\") * 1.0)) / $(MIN_TN_AD_RATIO) && \
				vc.getGenotype(\"$1\").getAnyAttribute(\"FAO\") * 1.0 < 1 && vc.getGenotype(\"$2\").getAnyAttribute(\"FAO\") > 1 && \
				vc.getGenotype(\"$1\").getAnyAttribute(\"AO\") * 1.0 < 1 && vc.getGenotype(\"$2\").getAnyAttribute(\"AO\") > 1 }' \
		--filterName somaticAlleleDepth \
		--filterExpression '(vc.getGenotype(\"$1\").getAnyAttribute(\"FAO\") * 1.0 + vc.getGenotype(\"$1\").getAnyAttribute(\"FRO\") * 1.0) <= $$(MIN_TUMOR_DEPTH) || \
			(vc.getGenotype(\"$2\").getAnyAttribute(\"FAO\") * 1.0 + vc.getGenotype(\"$2\").getAnyAttribute(\"FRO\") * 1.0) <= $$(MIN_NORMAL_DEPTH)' \
		--filterName depthFilter \
		&& sed -i 's/getGenotype(\"\([^\"]*\)\")/getGenotype(\1)/g' $$@ && $$(RM) $$< $$<.idx")
endef
$(foreach pair,$(SAMPLE_PAIRS),$(eval $(call som-ad-ft-tumor-normal,$(tumor.$(pair)),$(normal.$(pair)))))



define sufam
ifeq ($3,$2)
vcf/$1_$2.%.sufam.vcf : vcf/$1_$2.%.vcf
	$$(INIT) ln -f $$< $$@
else
vcf/$3.%.sufam.tmp : $$(foreach tumor,$$(wordlist 1,$$(shell expr $$(words $$(subst _,$$( ),$3)) - 1),$$(subst _,$$( ),$3)),vcf/$$(tumor)_$$(lastword $$(subst _,$$( ),$3)).%.vcf)
	$$(call LSCRIPT_MEM,22G,03:59:59,"$$(LOAD_JAVA8_MODULE); $$(call COMBINE_VARIANTS,21G) \
		$$(foreach vcf,$$^,--variant $$(vcf) ) -o $$@ --genotypemergeoption UNSORTED -R $$(REF_FASTA)")


vcf/$1_$2.%.sufam.tmp1.vcf : vcf/$1_$2.%.vcf vcf/$3.%.sufam.tmp
	$$(call CHECK_VCF,$$(word 2,$$^),$$@,\
		$$(call LSCRIPT_CHECK_MEM,12G,00:29:59,"$$(LOAD_JAVA8_MODULE); \
		$$(call SELECT_VARIANTS,10G) -R $$(REF_FASTA) --variant $$(word 2,$$^) --discordance $$(word 1,$$^) -o $$@ && $$(RM) $$(word 2,$$^) $$(word 2,$$^).idx"))

vcf/$1_$2.%.sufam.tmp1.T.vcf : vcf/$1_$2.%.sufam.tmp1.vcf bam/$1.bam
	$$(call CHECK_VCF,$$(word 1,$$^),$$@,\
		$$(call LSCRIPT_CHECK_MEM,10G,03:59:59,"$$(LOAD_BCFTOOLS_MODULE); $$(LOAD_JAVA8_MODULE); \
			$$(TVC) -s $$(word 1,$$^) -i $$(word 2,$$^) -r $$(REF_FASTA) -N 1 -m $$(TVC_MOTIF) -o $$(@D)/$1_$2.$$*.sufam.tmp1.T \
			-t $$(TVC_ROOT_DIR) --primer-trim-bed $$(PRIMER_TRIM_BED) && \
			$$(BCFTOOLS) norm -f $$(REF_FASTA) -m-both $$(@D)/$1_$2.$$*.sufam.tmp1.T/TSVC_variants.vcf.gz | grep -v \"##contig\" | \
			$$(FIX_TVC_VCF) > $$(@D)/$1_$2.$$*.sufam.tmp1.T/TSVC_variants.norm.vcf && \
			$$(call SELECT_VARIANTS,6G) -R $$(REF_FASTA) --variant $$(@D)/$1_$2.$$*.sufam.tmp1.T/TSVC_variants.norm.vcf -o $$@ --concordance $$(word 1,$$^) && $$(RMR) $$(@D)/$1_$2.$$*.sufam.tmp1.T"))

vcf/$1_$2.%.sufam.tmp1.N.vcf : vcf/$1_$2.%.sufam.tmp1.vcf bam/$2.bam
	$$(call CHECK_VCF,$$(word 1,$$^),$$@,\
		$$(call LSCRIPT_CHECK_MEM,10G,03:59:59,"$$(LOAD_BCFTOOLS_MODULE); $$(LOAD_JAVA8_MODULE); \
			$$(TVC) -s $$(word 1,$$^) -i $$(word 2,$$^) -r $$(REF_FASTA) -N 1 -m $$(TVC_MOTIF) -o $$(@D)/$1_$2.$$*.sufam.tmp1.N \
			-t $$(TVC_ROOT_DIR) --primer-trim-bed $$(PRIMER_TRIM_BED) && \
			$$(BCFTOOLS) norm -f $$(REF_FASTA) -m-both $$(@D)/$1_$2.$$*.sufam.tmp1.N/TSVC_variants.vcf.gz | grep -v \"##contig\" | \
			$$(FIX_TVC_VCF) > $$(@D)/$1_$2.$$*.sufam.tmp1.N/TSVC_variants.norm.vcf && \
			$$(call SELECT_VARIANTS,6G) -R $$(REF_FASTA) --variant $$(@D)/$1_$2.$$*.sufam.tmp1.N/TSVC_variants.norm.vcf -o $$@ --concordance $$(word 1,$$^) && $$(RMR) $$(@D)/$1_$2.$$*.sufam.tmp1.N"))

vcf/$1_$2.%.sufam.tmp2.vcf : vcf/$1_$2.%.sufam.tmp1.T.vcf vcf/$1_$2.%.sufam.tmp1.N.vcf
	$$(call CHECK_VCF,$$(word 1,$$^),$$@,\
		$$(call LSCRIPT_CHECK_MEM,12G,00:29:59,"$$(LOAD_JAVA8_MODULE); \
			$$(call COMBINE_VARIANTS,10G) -R $$(REF_FASTA) \
			--variant $$(word 1,$$^) --variant $$(word 2,$$^) -o $$@"))

vcf/$1_$2.%.sufam.tmp3.vcf : vcf/$1_$2.%.sufam.tmp2.vcf
	$$(call CHECK_VCF,$$<,$$@,\
		$$(call LSCRIPT_CHECK_MEM,12G,00:29:59,"$$(LOAD_JAVA8_MODULE); \
			$$(call VARIANT_FILTRATION,7G) -R $$(REF_FASTA) -V $$<  -o $$@  --filterName interrogation \
			--filterExpression 'vc.getGenotype(\"$1\").getAnyAttribute(\"FSAF\") > 0 || vc.getGenotype(\"$1\").getAnyAttribute(\"FSAR\") > 0' && $$(RM) $$< $$<.idx"))

vcf/$1_$2.%.sufam.tmp4.vcf : vcf/$1_$2.%.sufam.tmp3.vcf
	$$(call CHECK_VCF,$$<,$$@,\
		$$(call LSCRIPT_CHECK_MEM,12G,00:29:59,"$$(LOAD_SNP_EFF_MODULE); $$(SNP_SIFT) filter $$(SNP_SIFT_OPTS) -f $$< \"(FILTER has 'interrogation')\"  > $$@ && $$(RM) $$< $$<.idx"))

vcf/$1_$2.%.sufam.vcf : vcf/$1_$2.%.vcf vcf/$1_$2.%.sufam.tmp4.vcf
	$$(call CHECK_VCF_CMD,$$(word 2,$$^),cp $$< $$@,\
		$$(call LSCRIPT_CHECK_MEM,12G,00:29:59,"$$(LOAD_JAVA8_MODULE); \
			$$(call COMBINE_VARIANTS,21G) --variant $$< --variant $$(word 2,$$^) -o $$@ --genotypemergeoption UNSORTED -R $$(REF_FASTA) && $$(RM) $$(word 2,$$^) $$(word 2,$$^).idx"))


endif
endef

$(foreach set,$(SAMPLE_SETS),\
	$(foreach tumor,$(wordlist 1,$(shell expr $(words $(subst _,$( ),$(set))) - 1),$(subst _,$( ),$(set))),\
		$(eval $(call sufam,$(tumor),$(lastword $(subst _,$( ),$(set))),$(subst $(tumor)_,,$(set))))))


endif

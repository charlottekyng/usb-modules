%.altad_ft.vcf : %.vcf
	$(call LSCRIPT_CHECK_MEM,$(RESOURCE_REQ_MEDIUM_MEM),$(RESOURCE_REQ_SHORT),"$(LOAD_JAVA8_MODULE); \
		$(call GATK,VariantFiltration,$(RESOURCE_REQ_MEDIUM_MEM)) \
		-R $(REF_FASTA) -V $< -o $@ \
		--filterExpression 'vc.getGenotype(\"\$$1\").getAnyAttribute(\"FAO\") == 0 && vc.getGenotype(\"\$$1\").getAnyAttribute(\"AO\") == 0' \
		--filterName zeroAD && $(RM) $< $<.idx")

%.dp_ft.vcf : %.vcf
	$(call LSCRIPT_CHECK_MEM,$(RESOURCE_REQ_MEDIUM_MEM),$(RESOURCE_REQ_SHORT),"$(LOAD_JAVA8_MODULE); \
		$(call GATK,VariantFiltration,$(RESOURCE_REQ_MEDIUM_MEM)) \
		-R $(REF_FASTA) -V $< -o $@ \
		--filterExpression 'vc.getGenotype(\"\$$1\").getAnyAttribute(\"DP\") * 1.0 <= $(MIN_NORMAL_DEPTH) || \
		vc.getGenotype(\"\$$1\").getAnyAttribute(\"FDP\") * 1.0 <= $(MIN_NORMAL_DEPTH)' --filterName Depth && \
		$(RM) $< $<.idx")

ifdef SAMPLE_PAIRS
#define som-ad-ft-tumor-normal
#vcf/$1_$2.%.som_ad_ft.vcf : vcf/$1_$2.%.vcf
#	$$(call LSCRIPT_CHECK_MEM,$$(RESOURCE_REQ_MEDIUM_MEM),$$(RESOURCE_REQ_SHORT),"$$(LOAD_JAVA8_MODULE); \
#		$$(call GATK,VariantFiltration,$$(RESOURCE_REQ_MEDIUM_MEM)) -R $$(REF_FASTA) -V $$< -o $$@ \
#		--filterExpression 'vc.getGenotype(\"$1\").getAnyAttribute(\"FAO\") < $(MIN_TUMOR_AD) && vc.getGenotype(\"$1\").getAnyAttribute(\"AO\") < $(MIN_TUMOR_AD)' \
#		--filterName tumorVarAlleleDepth \
#		--filterExpression 'vc.getGenotype(\"$1\").getAnyAttribute(\"DP\") * 1.0 <= $(MIN_TUMOR_DEPTH) || vc.getGenotype(\"$1\").getAnyAttribute(\"FDP\") * 1.0 <= $(MIN_TUMOR_DEPTH) || \
#			vc.getGenotype(\"$2\").getAnyAttribute(\"DP\") * 1.0 <= $(MIN_NORMAL_DEPTH) || vc.getGenotype(\"$2\").getAnyAttribute(\"FDP\") * 1.0 <= $(MIN_NORMAL_DEPTH)' \
#		--filterName depthFilter \
#		--filterExpression 'if (vc.getGenotype(\"$1\").getAnyAttribute(\"FDP\") * 1.0 > $(MIN_TUMOR_DEPTH) && vc.getGenotype(\"$2\").getAnyAttribute(\"FDP\") * 1.0 > $(MIN_NORMAL_DEPTH)) { \
#			 ( vc.getGenotype(\"$2\").getAnyAttribute(\"FAO\") * 1.0 / vc.getGenotype(\"$2\").getAnyAttribute(\"FDP\") * 1.0  ) \
#				> ( vc.getGenotype(\"$1\").getAnyAttribute(\"FAO\") * 1.0 / vc.getGenotype(\"$1\").getAnyAttribute(\"FDP\") * 1.0  ) / $(MIN_TN_AD_RATIO) \
#			} else { vc.getGenotype(\"$2\").getAnyAttribute(\"FAO\") > 0 }' \
#		--filterName FsomaticAlleleDepth \
#		--filterExpression 'if (vc.getGenotype(\"$1\").getAnyAttribute(\"DP\") * 1.0 > $(MIN_TUMOR_DEPTH) && vc.getGenotype(\"$2\").getAnyAttribute(\"DP\") * 1.0 > $(MIN_NORMAL_DEPTH)) { \
#			( vc.getGenotype(\"$2\").getAnyAttribute(\"AO\") * 1.0 / vc.getGenotype(\"$2\").getAnyAttribute(\"DP\") * 1.0  ) \
#				> ( vc.getGenotype(\"$1\").getAnyAttribute(\"AO\") * 1.0 / vc.getGenotype(\"$1\").getAnyAttribute(\"DP\") * 1.0  ) / $(MIN_TN_AD_RATIO) \
#			} else { vc.getGenotype(\"$2\").getAnyAttribute(\"AO\") > 0 }'	\
#		--filterName somaticAlleleDepth \
#		&& sed -i 's/getGenotype(\"\([^\"]*\)\")/getGenotype(\1)/g' $$@ && $$(RM) $$< $$<.idx")
#endef

define som-ad-ft-tumor-normal
vcf/$1_$2.%.som_ad_ft.vcf : vcf/$1_$2.%.vcf
	$$(call LSCRIPT_CHECK_MEM,$$(RESOURCE_REQ_MEDIUM_MEM),$$(RESOURCE_REQ_SHORT),"$$(LOAD_JAVA8_MODULE); \
		$$(call GATK,VariantFiltration,$$(RESOURCE_REQ_MEDIUM_MEM)) -R $$(REF_FASTA) -V $$< -o $$@ \
		--filterExpression 'vc.getGenotype(\"$1\").getAnyAttribute(\"FAO\") < $(MIN_TUMOR_AD)' \
		--filterName tumorVarAD_raw \
		--filterExpression 'vc.getGenotype(\"$1\").getAnyAttribute(\"AO\") < $(MIN_TUMOR_AD)' \
		--filterName tumorVarAD_flow \
		--filterExpression 'vc.getGenotype(\"$1\").getAnyAttribute(\"DP\") * 1.0 <= $(MIN_TUMOR_DEPTH)' \
		--filterName tumorDepthFilter_raw \
		--filterExpression 'vc.getGenotype(\"$1\").getAnyAttribute(\"FDP\") * 1.0 <= $(MIN_TUMOR_DEPTH)' \
		--filterName tumorDepthFilter_flow \
		--filterExpression 'vc.getGenotype(\"$2\").getAnyAttribute(\"DP\") * 1.0 <= $(MIN_NORMAL_DEPTH)' \
		--filterName normalDepthFilter_raw \
		--filterExpression 'vc.getGenotype(\"$2\").getAnyAttribute(\"FDP\") * 1.0 <= $(MIN_NORMAL_DEPTH)' \
		--filterName normalDepthFilter_flow \
		--filterExpression '( vc.getGenotype(\"$2\").getAnyAttribute(\"AO\") * 1.0 / vc.getGenotype(\"$2\").getAnyAttribute(\"DP\") * 1.0  ) \
			> ( vc.getGenotype(\"$1\").getAnyAttribute(\"AO\") * 1.0 / vc.getGenotype(\"$1\").getAnyAttribute(\"DP\") * 1.0  ) / $(MIN_TN_AD_RATIO)' \
		--filterName somAD_raw \
		--filterExpression '( vc.getGenotype(\"$2\").getAnyAttribute(\"AO\") * 1.0 / vc.getGenotype(\"$2\").getAnyAttribute(\"DP\") * 1.0  ) \
			> ( vc.getGenotype(\"$1\").getAnyAttribute(\"FAO\") * 1.0 / vc.getGenotype(\"$1\").getAnyAttribute(\"FDP\") * 1.0  ) / $(MIN_TN_AD_RATIO)' \
		--filterName somAD_flow \
		&& sed -i 's/getGenotype(\"\([^\"]*\)\")/getGenotype(\1)/g' $$@ && $$(RM) $$< $$<.idx")
endef
$(foreach pair,$(SAMPLE_PAIRS),$(eval $(call som-ad-ft-tumor-normal,$(tumor.$(pair)),$(normal.$(pair)))))

#define ad-ft-tumor-normal
#%.altad_ft_tn.vcf : %.vcf
#	$$(call LSCRIPT_CHECK_MEM,$$(RESOURCE_REQ_MEDIUM_MEM),$$(RESOURCE_REQ_SHORT),"$$(LOAD_JAVA8_MODULE); \
#		$$(call GATK,VariantFiltration,$$(RESOURCE_REQ_MEDIUM_MEM)) \
#		-R $(REF_FASTA) -V $$< -o $$@ \
#		--filterExpression 'vc.getGenotype(\"\$$1\").getAnyAttribute(\"FAO\") == 0 && vc.getGenotype(\"\$$1\").getAnyAttribute(\"AO\") == 0' \
#		--filterName zeroAD && $(RM) $< $<.idx")
#endef
#$(foreach pair,$(SAMPLE_PAIRS),$(eval $(call ad-ft-tumor-normal,$(tumor.$(pair)),$(normal.$(pair)))))

define sufam
ifeq ($3,$2)
vcf/$1_$2.%.sufam.vcf : vcf/$1_$2.%.vcf
	$$(INIT) ln -f $$< $$@
else
#vcf/$3.%.sufam.tmp : $$(foreach tumor,$$(wordlist 1,$$(shell expr $$(words $$(subst _,$$( ),$3)) - 1),$$(subst _,$$( ),$3)),vcf/$$(tumor)_$$(lastword $$(subst _,$$( ),$3)).%.vcf)
#	$$(call LSCRIPT_MEM,$$(RESOURCE_REQ_HIGHMEM),$$(RESOURCE_REQ_VSHORT),"$$(LOAD_JAVA8_MODULE); \
#		$$(call GATK,CombineVariants,$$(RESOURCE_REQ_HIGHMEM)) \
#		$$(foreach vcf,$$^,--variant $$(vcf) ) -o $$@ --genotypemergeoption UNSORTED -R $$(REF_FASTA)")

#vcf/$1_$2.%.sufam.tmp1.vcf : vcf/$1_$2.%.vcf vcf/$3.%.sufam.tmp
#	$$(call CHECK_VCF,$$(word 2,$$^),$$@,\
#		$$(call LSCRIPT_CHECK_MEM,$$(RESOURCE_REQ_HIGHMEM),$$(RESOURCE_REQ_VSHORT),"$$(LOAD_JAVA8_MODULE); \
#		$$(call GATK,SelectVariants,$$(RESOURCE_REQ_HIGHMEM)) -R $$(REF_FASTA) --variant $$(word 2,$$^) \
#			--discordance $$(word 1,$$^) -o $$@ && $$(RM) $$(word 2,$$^) $$(word 2,$$^).idx"))

#vcf/$1_$2.%.sufam.tmp1.T.vcf : vcf/$1_$2.%.sufam.tmp1.vcf bam/$1.bam vcf/$1_$2.%.sufam.tmp1.vcf.gz vcf/$1_$2.%.sufam.tmp1.vcf.gz.tbi
#	$$(call CHECK_VCF,$$(word 1,$$^),$$@,\
#		$$(call LSCRIPT_CHECK_MEM,$$(RESOURCE_REQ_HIGHMEM),$$(RESOURCE_REQ_SHORT),"$$(LOAD_BCFTOOLS_MODULE); $$(LOAD_JAVA8_MODULE); $$(LOAD_TABIX_MODULE); \
#			$$(TVC) -s $$(word 1,$$^) -i $$(word 2,$$^) -r $$(REF_FASTA) -N 1 -m $$(TVC_MOTIF) -o $$(@D)/$1_$2.$$*.sufam.tmp1.T \
#			-t $$(TVC_ROOT_DIR) --primer-trim-bed $$(PRIMER_TRIM_BED) && \
#			$$(BCFTOOLS) norm -f $$(REF_FASTA) -m-both $$(@D)/$1_$2.$$*.sufam.tmp1.T/TSVC_variants.vcf.gz | \
#			grep -v \"##contig\" > $$(@D)/$1_$2.$$*.sufam.tmp1.T/TSVC_variants.vcf.tmp && \
#			$$(call GATK,LeftAlignAndTrimVariants,$$(RESOURCE_REQ_MEDIUM_MEM)) -R $$(REF_FASTA) \
#			--variant $$(@D)/$1_$2.$$*.sufam.tmp1.T/TSVC_variants.vcf.tmp -o $$(@D)/$1_$2.$$*.sufam.tmp1.T/TSVC_variants.vcf.tmp2 && \
#			$$(FIX_TVC_VCF) $$(@D)/$1_$2.$$*.sufam.tmp1.T/TSVC_variants.vcf.tmp2 > $$(@D)/$1_$2.$$*.sufam.tmp1.T/TSVC_variants.vcf.tmp3 && \
#			$$(BGZIP) -c $$(@D)/$1_$2.$$*.sufam.tmp1.T/TSVC_variants.vcf.tmp3 > $$(@D)/$1_$2.$$*.sufam.tmp1.T/TSVC_variants.vcf.tmp3.gz && \
#			$$(TABIX) -p vcf $$(@D)/$1_$2.$$*.sufam.tmp1.T/TSVC_variants.vcf.tmp3.gz && \
#			$$(BCFTOOLS) isec -O v -p $$(@D)/$1_$2.$$*.sufam.tmp1.T/isec $$(@D)/$1_$2.$$*.sufam.tmp1.T/TSVC_variants.vcf.tmp3.gz $$(<<<) && \
#			mv $$(@D)/$1_$2.$$*.sufam.tmp1.T/isec/0002.vcf $$@ && $$(RMR) $$(@D)/$1_$2.$$*.sufam.tmp1.T"))

#vcf/$1_$2.%.sufam.tmp1.N.vcf : vcf/$1_$2.%.sufam.tmp1.vcf bam/$2.bam vcf/$1_$2.%.sufam.tmp1.vcf.gz vcf/$1_$2.%.sufam.tmp1.vcf.gz.tbi
#	$$(call CHECK_VCF,$$(word 1,$$^),$$@,\
#		$$(call LSCRIPT_CHECK_MEM,$$(RESOURCE_REQ_HIGHMEM),$$(RESOURCE_REQ_SHORT),"$$(LOAD_BCFTOOLS_MODULE); $$(LOAD_JAVA8_MODULE); $$(LOAD_TABIX_MODULE);\
#			$$(TVC) -s $$(word 1,$$^) -i $$(word 2,$$^) -r $$(REF_FASTA) -N 1 -m $$(TVC_MOTIF) -o $$(@D)/$1_$2.$$*.sufam.tmp1.N \
#			-t $$(TVC_ROOT_DIR) --primer-trim-bed $$(PRIMER_TRIM_BED) && \
#			$$(BCFTOOLS) norm -f $$(REF_FASTA) -m-both $$(@D)/$1_$2.$$*.sufam.tmp1.N/TSVC_variants.vcf.gz | \
#			grep -v \"##contig\" > $$(@D)/$1_$2.$$*.sufam.tmp1.N/TSVC_variants.vcf.tmp && \
#			$$(call GATK,LeftAlignAndTrimVariants,$$(RESOURCE_REQ_MEDIUM_MEM)) -R $$(REF_FASTA) \
#			--variant $$(@D)/$1_$2.$$*.sufam.tmp1.N/TSVC_variants.vcf.tmp -o $$(@D)/$1_$2.$$*.sufam.tmp1.N/TSVC_variants.vcf.tmp2 && \
#			$$(FIX_TVC_VCF) $$(@D)/$1_$2.$$*.sufam.tmp1.N/TSVC_variants.vcf.tmp2 > $$(@D)/$1_$2.$$*.sufam.tmp1.N/TSVC_variants.vcf.tmp3 && \
#			$$(BGZIP) -c $$(@D)/$1_$2.$$*.sufam.tmp1.N/TSVC_variants.vcf.tmp3 > $$(@D)/$1_$2.$$*.sufam.tmp1.N/TSVC_variants.vcf.tmp3.gz && \
#			$$(TABIX) -p vcf $$(@D)/$1_$2.$$*.sufam.tmp1.N/TSVC_variants.vcf.tmp3.gz && \
#			$$(BCFTOOLS) isec -O v -p $$(@D)/$1_$2.$$*.sufam.tmp1.N/isec $$(@D)/$1_$2.$$*.sufam.tmp1.N/TSVC_variants.vcf.tmp3.gz $$(<<<) && \
#			mv $$(@D)/$1_$2.$$*.sufam.tmp1.N/isec/0002.vcf $$@ && $$(RMR) $$(@D)/$1_$2.$$*.sufam.tmp1.N"))

#vcf/$1_$2.%.sufam.tmp2.vcf : vcf/$1_$2.%.sufam.tmp1.T.vcf vcf/$1_$2.%.sufam.tmp1.N.vcf vcf/$1_$2.%.sufam.tmp1.vcf
#	$$(call CHECK_VCF,$$(word 1,$$^),$$@,\
#		$$(call LSCRIPT_CHECK_MEM,$$(RESOURCE_REQ_MEDIUM_MEM),$$(RESOURCE_REQ_VSHORT),"$$(LOAD_JAVA8_MODULE); \
#			$$(call GATK,CombineVariants,$$(RESOURCE_REQ_MEDIUM_MEM)) -R $$(REF_FASTA) \
#			--variant $$(word 1,$$^) --variant $$(word 2,$$^) -o $$@ && \
#			$$(RM) $$< $$<.idx $$(word 2,$$^) $$(word 2,$$^).idx $$(word 3,$$^)*"))

#vcf/$1_$2.%.sufam.tmp3.vcf : vcf/$1_$2.%.sufam.tmp2.vcf
#	$$(call CHECK_VCF,$$<,$$@,\
#		$$(call LSCRIPT_CHECK_MEM,$$(RESOURCE_REQ_MEDIUM_MEM),$$(RESOURCE_REQ_VSHORT),"$$(LOAD_JAVA8_MODULE); \
#			$$(call GATK,VariantFiltration,$$(RESOURCE_REQ_MEDIUM_MEM)) -R $$(REF_FASTA) -V $$<  -o $$@  \
#			--filterName interrogation \
#			--filterExpression 'vc.getGenotype(\"$1\").getAnyAttribute(\"FAO\") > 0 || vc.getGenotype(\"$1\").getAnyAttribute(\"AO\") > 0' \
#			--filterName interrogation_Absent \
#			--filterExpression 'vc.getGenotype(\"$1\").getAnyAttribute(\"FAO\") == 0 && vc.getGenotype(\"$1\").getAnyAttribute(\"AO\") == 0' \
#			--maskName HOTSPOT --mask $(CANCER_HOTSPOT_VCF) \
#			&& $$(RM) $$< $$<.idx"))

#vcf/$1_$2.%.sufam.tmp4.vcf : vcf/$1_$2.%.sufam.tmp3.vcf
#	$$(call CHECK_VCF,$$<,$$@,\
#		$$(call LSCRIPT_CHECK_MEM,12G,00:29:59,"$$(LOAD_SNP_EFF_MODULE); $$(SNP_SIFT) filter $$(SNP_SIFT_OPTS) \
#			-f $$< \"(FILTER has 'interrogation' || FILTER has 'interrogation_Absent')\"  > $$@ && $$(RM) $$< $$<.idx"))

#vcf/$1_$2.%.sufam.vcf : vcf/$1_$2.%.vcf vcf/$1_$2.%.sufam.tmp4.vcf
#	$$(call CHECK_VCF_CMD,$$(word 2,$$^),cp $$< $$@,\
#		$$(call LSCRIPT_CHECK_MEM,$$(RESOURCE_REQ_HIGHMEM),$$(RESOURCE_REQ_VSHORT),"$$(LOAD_JAVA8_MODULE); $$(LOAD_BCFTOOLS_MODULE); \
#			$$(call GATK,CombineVariants,$$(RESOURCE_REQ_HIGHMEM)) --variant $$< --variant $$(word 2,$$^) -o $$(word 2,$$^).tmp \
#			--genotypemergeoption UNSORTED -R $$(REF_FASTA) -assumeIdenticalSamples && \
#			$$(BCFTOOLS) norm -f $$(REF_FASTA) -m-both $$(word 2,$$^).tmp | grep -v \"##contig\" > $$(word 2,$$^).tmp2 && \
#			$$(call GATK,LeftAlignAndTrimVariants,$$(RESOURCE_REQ_MEDIUM_MEM)) -R $$(REF_FASTA) --variant $$(word 2,$$^).tmp2 -o $$@ && \
#			$$(RM) $$(word 2,$$^) $$(word 2,$$^).idx $$(word 2,$$^).tmp*"))

vcf/$3.%.sufam.tmp1.vcf : $$(foreach tumor,$$(wordlist 1,$$(shell expr $$(words $$(subst _,$$( ),$3)) - 1),$$(subst _,$$( ),$3)),vcf/$$(tumor)_$$(lastword $$(subst _,$$( ),$3)).%.vcf)
	$$(call LSCRIPT_MEM,$$(RESOURCE_REQ_MEDIUM_MEM),$$(RESOURCE_REQ_VSHORT),"$$(LOAD_JAVA8_MODULE); \
	$$(call GATK,CombineVariants,$$(RESOURCE_REQ_HIGHMEM)) \
	$$(foreach vcf,$$^,--variant $$(vcf) ) -o $$@ --genotypemergeoption UNSORTED -R $$(REF_FASTA)")

vcf/$1_$2.%/isec/0001.vcf : vcf/$1_$2.%.vcf vcf/$3.%.sufam.tmp1.vcf vcf/$1_$2.%.vcf.gz vcf/$3.%.sufam.tmp1.vcf.gz vcf/$1_$2.%.vcf.gz.tbi vcf/$3.%.sufam.tmp1.vcf.gz.tbi
	$(MKDIR) $$(dir $$@); $$(call CHECK_VCF,$$(word 2,$$^),$$@,\
	$$(call LSCRIPT_MEM,$$(RESOURCE_REQ_LOWMEM),$$(RESOURCE_REQ_VSHORT),"$$(LOAD_BCFTOOLS_MODULE); \
	$$(BCFTOOLS) isec -O v -p $$(dir $$@) $$(word 3,$$^) $$(word 4,$$^)"))

vcf/$1_$2.%/TSVC_variants.vcf.gz : vcf/$1_$2.%/isec/0001.vcf bam/$1.bam bam/$1.bam.bai
	$$(call LSCRIPT_PARALLEL_MEM,4,$$(RESOURCE_REQ_MEDIUM_MEM),$$(RESOURCE_REQ_SHORT),"$$(TVC) \
	-s $$< -i $$(word 2,$$^) -r $$(REF_FASTA) -o $$(@D) -N 4 \
	$$(if $$(TARGETS_FILE_INTERVALS),-b $$(TARGETS_FILE_INTERVALS)) -m $$(TVC_MOTIF) \
	-t $$(TVC_ROOT_DIR) --primer-trim-bed $$(PRIMER_TRIM_BED)")

vcf/$1_$2.%/isec2/0002.vcf : vcf/$1_$2.%/TSVC_variants.biallelic_ft.vcf vcf/$1_$2.%/isec/0001.vcf vcf/$1_$2.%/TSVC_variants.biallelic_ft.vcf.gz vcf/$1_$2.%/isec/0001.vcf.gz vcf/$1_$2.%/TSVC_variants.biallelic_ft.vcf.gz.tbi vcf/$1_$2.%/isec/0001.vcf.gz.tbi
	$$(call LSCRIPT_MEM,$$(RESOURCE_REQ_LOWMEM),$$(RESOURCE_REQ_VSHORT),"$$(LOAD_BCFTOOLS_MODULE); \
	$$(BCFTOOLS) isec -O v -p $$(dir $$@) $$(word 3,$$^) $$(word 4,$$^)")

vcf/$1_$2.%/isec2/0001.vcf : vcf/$1_$2.%/isec2/0002.vcf
	

vcf/$1_$2.%/isec3/0002.vcf : vcf/$1_$2.%/TSVC_variants.multiallelic_ft.norm.left_align.vcf vcf/$1_$2.%/isec2/0001.vcf vcf/$1_$2.%/TSVC_variants.multiallelic_ft.norm.left_align.vcf.gz vcf/$1_$2.%/isec2/0001.vcf.gz vcf/$1_$2.%/TSVC_variants.multiallelic_ft.norm.left_align.vcf.gz.tbi vcf/$1_$2.%/isec2/0001.vcf.gz.tbi
	$$(call LSCRIPT_MEM,$$(RESOURCE_REQ_LOWMEM),$$(RESOURCE_REQ_VSHORT),"$$(LOAD_BCFTOOLS_MODULE); \
	$$(BCFTOOLS) isec -O v -p $$(dir $$@) $$(word 3,$$^) $$(word 4,$$^)")

vcf/$1_$2.%.sufam.tmp2.vcf : vcf/$1_$2.%/isec2/0002.post_bcftools.vcf vcf/$1_$2.%/isec3/0002.post_bcftools.vcf
	$$(call LSCRIPT_MEM,$$(RESOURCE_REQ_MEDIUM_MEM),$$(RESOURCE_REQ_VSHORT),"$$(LOAD_JAVA8_MODULE); \
	$$(call GATK,CombineVariants,$$(RESOURCE_REQ_MEDIUM_MEM)) -R $$(REF_FASTA) \
	$$(foreach vcf,$$^,-V $$(vcf)) -o $$@ --assumeIdenticalSamples && $$(RMR) $$(dir $$<)")

vcf/$1_$2.%.sufam.tmp3.vcf : vcf/$1_$2.%.sufam.tmp2.vcf
	$$(call CHECK_VCF,$$<,$$@,\
	$$(call LSCRIPT_CHECK_MEM,$$(RESOURCE_REQ_MEDIUM_MEM),$$(RESOURCE_REQ_VSHORT),"$$(LOAD_JAVA8_MODULE); \
	$$(call GATK,VariantFiltration,$$(RESOURCE_REQ_MEDIUM_MEM)) -R $$(REF_FASTA) -V $$<  -o $$@  \
	--filterName interrogation \
	--filterExpression 'vc.getGenotype(\"$1\").getAnyAttribute(\"FAO\") > 0 || vc.getGenotype(\"$1\").getAnyAttribute(\"AO\") > 0' \
	--filterName interrogation_Absent \
	--filterExpression 'vc.getGenotype(\"$1\").getAnyAttribute(\"FAO\") == 0 && vc.getGenotype(\"$1\").getAnyAttribute(\"AO\") == 0' \
	--maskName HOTSPOT --mask $(CANCER_HOTSPOT_VCF) \
	&& $$(RM) $$< $$<.idx"))

vcf/$1_$2.%.sufam.tmp4.vcf : vcf/$3.%.sufam.tmp1.vcf
	$$(call LSCRIPT_CHECK_MEM,$$(RESOURCE_REQ_MEDIUM_MEM),$$(RESOURCE_REQ_VSHORT),"$$(LOAD_JAVA8_MODULE); \
	$$(call GATK,SelectVariants,$$(RESOURCE_REQ_MEDIUM_MEM)) -R $$(REF_FASTA) -V $$< -o $$@ -sn $2")

vcf/$1_$2.%.sufam.vcf : vcf/$1_$2.%.sufam.tmp3.vcf vcf/$1_$2.%.sufam.tmp4.vcf vcf/$1_$2.%.vcf
	$$(call LSCRIPT_CHECK_MEM,$$(RESOURCE_REQ_MEDIUM_MEM),$$(RESOURCE_REQ_VSHORT),"$$(LOAD_JAVA8_MODULE); \
	$$(call GATK,CombineVariants,$$(RESOURCE_REQ_MEDIUM_MEM)) -R $$(REF_FASTA) $$(foreach vcf,$$^,-V $$(vcf)) -o $$@ \
	--genotypemergeoption UNSORTED -filteredRecordsMergeType KEEP_IF_ALL_UNFILTERED && $(RM) vcf/$1_$2.%.sufam.tmp*")

endif
endef

$(foreach set,$(SAMPLE_SETS),\
	$(foreach tumor,$(wordlist 1,$(shell expr $(words $(subst _,$( ),$(set))) - 1),$(subst _,$( ),$(set))),\
		$(eval $(call sufam,$(tumor),$(lastword $(subst _,$( ),$(set))),$(subst $(tumor)_,,$(set))))))


endif

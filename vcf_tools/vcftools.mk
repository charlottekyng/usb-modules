# vim: set ft=make :
# sub module containing vcf related tools

ifndef VCFTOOLS_MK

include usb-modules/Makefile.inc

LOGDIR ?= log/vcf.$(NOW)

.SECONDEXPANSION:
..DUMMY := $(shell mkdir -p version; echo "$(SNP_EFF) &> version/snp_eff.txt")


%.vcf.idx : %.vcf
	$(call LSCRIPT_CHECK_MEM,4G,00:29:29,"$(LOAD_IGVTOOLS_MODULE); $(IGVTOOLS) index $< && sleep 10")

%.vcf.gz : %.vcf
	$(call LSCRIPT_MEM,3G,00:29:29,"$(LOAD_TABIX_MODULE); $(BGZIP) -c -f $< >$@")

%.vcf.gz.tbi : %.vcf.gz
	$(call LSCRIPT_MEM,3G,00:29:29,"$(LOAD_TABIX_MODULE); $(TABIX) $<")

#%.vcf : %.vcf.gz
#	$(call LSCRIPT_MEM,3G,00:29:29,"$(LOAD_TABIX_MODULE); $(BGZIP) -c -d $< > $@")

############ FILTERS #########


ifeq ($(findstring ILLUMINA,$(SEQ_PLATFORM)),ILLUMINA)
%.nft.vcf : %.vcf mutect2/pon.mutect2.vcf
	$(call LSCRIPT_CHECK_MEM,8G,01:59:59,"$(LOAD_JAVA8_MODULE); $(call VARIANT_FILTRATION,7G) \
		-R $(REF_FASTA) -V $< -o $@ --maskName 'PoN' --mask $(word 2,$^) && $(RM) $< $<.idx")
endif
ifeq ($(findstring IONTORRENT,$(SEQ_PLATFORM)),IONTORRENT)
%.nft.vcf : %.vcf tvc/pon.tvc.vcf
	$(call LSCRIPT_CHECK_MEM,8G,01:59:59,"$(LOAD_JAVA8_MODULE); $(call VARIANT_FILTRATION,7G) \
		-R $(REF_FASTA) -V $< -o $@ --maskName 'PoN' --mask $(word 2,$^) && $(RM) $< $<.idx")
endif

%.target_ft.vcf : %.vcf
	$(call LSCRIPT_CHECK_MEM,8G,00:59:59,"$(LOAD_JAVA8_MODULE); $(call VARIANT_FILTRATION,7G) \
		-R $(REF_FASTA) -V $< -o $@ --mask $(TARGETS_FILE_INTERVALS) --maskName targetInterval --filterNotInMask && $(RM) $< $<.idx")

# process snp eff output with gatk %=sample.indels/snps
%.annotated.vcf : %.vcf %.gatk_eff.vcf %.gatk_eff.vcf.idx %.vcf.idx 
	$(call LSCRIPT_PARALLEL_MEM,5,3G,00:29:29,"$(LOAD_JAVA8_MODULE); $(call VARIANT_FILTRATION,2G) \
	-R $(REF_FASTA) -nt 5 -A SnpEff --variant $< --snpEffFile $(word 2,$^) -o $@ &> $(LOGDIR)/$@.log")

%.dp_ft.vcf : %.vcf
	$(call LSCRIPT_CHECK_MEM,8G,00:59:59,"$(LOAD_JAVA8_MODULE); $(call VARIANT_FILTRATION,7G) \
		 -R $(REF_FASTA) -V $< -o $@ --filterExpression 'DP < $(MIN_NORMAL_DEPTH)' --filterName Depth && $(RM) $< $<.idx")

%.het_ft.vcf : %.vcf
	$(call LSCRIPT_CHECK_MEM,8G,00:59:59,"$(LOAD_JAVA8_MODULE); $(call VARIANT_FILTRATION,7G) \
		-R $(REF_FASTA) -V $< -o $@ --genotypeFilterExpression 'isHet == 1' --genotypeFilterName 'Heterozygous positions'")
%.altad_ft.vcf : %.vcf
	$(call LSCRIPT_CHECK_MEM,8G,00:59:59,"$(LOAD_JAVA8_MODULE); $(call VARIANT_FILTRATION,7G) \
		-R $(REF_FASTA) -V $< -o $@ --filterExpression 'vc.getGenotype(\"$1\").getAD().1 < 0' --filterName nonZeroAD && $(RM) $< $<.idx")

#%.mapq_ft.vcf : %.vcf.gz %.vcf.gz.tbi $(foreach sample,$(SAMPLES),bam/$(sample).bam)
#	$(call 

%.encode_ft.vcf : %.vcf
	$(call LSCRIPT_CHECK_MEM,8G,01:59:59,"$(LOAD_JAVA8_MODULE); $(call VARIANT_FILTRATION,7G) \
		-R $(REF_FASTA) -V $< -o $@ --maskName 'encode' --mask $(ENCODE_BED) && $(RM) $< $<.idx")

%.hotspot.vcf : %.vcf
	$(call LSCRIPT_CHECK_MEM,8G,01:59:59,"$(LOAD_JAVA8_MODULE); $(call VARIANT_FILTRATION,7G) \
		-R $(REF_FASTA) -V $< -o $@ --maskName HOTSPOT --mask $(CANCER_HOTSPOT_VCF) && $(RM) $< $<.idx")

%.pass.vcf : %.vcf
	$(call CHECK_VCF,$<,$@,$(call LSCRIPT_CHECK_MEM,2G,00:29:29,"$(LOAD_SNP_EFF_MODULE); $(SNP_SIFT) filter \
		$(SNP_SIFT_OPTS) -f $< \"( na FILTER ) | (FILTER = 'PASS') | (FILTER has 'HOTSPOT')\" > $@"))

############ ANNOTATION #########

%.eff.vcf : %.vcf %.vcf.idx
	$(call CHECK_VCF,$<,$@,$(call LSCRIPT_CHECK_MEM,9G,01:59:59,"$(LOAD_SNP_EFF_MODULE); $(SNP_EFF) ann \
		$(SNP_EFF_OPTS) $(SNP_EFF_GENOME) $< > $@"))

%.nsfp.vcf : %.vcf %.vcf.idx
	$(call CHECK_VCF,$<,$@,$(call LSCRIPT_CHECK_MEM,9G,01:59:59,"$(LOAD_SNP_EFF_MODULE); $(SNP_SIFT) dbnsfp \
		$(SNP_SIFT_OPTS) -db $(DB_NSFP) $< | sed '/^##INFO=<ID=dbNSFP/ s/Character/String/' > $@ && $(RM) $^"))
#                $(SNP_SIFT_OPTS) -f $(subst $( ),$(,),$(NSFP_FIELDS)) -db $(DB_NSFP) $< | sed '/^##INFO=<ID=dbNSFP/ s/Character/String/' > $@ && $(RM) $^"))

%.gatk_eff.vcf : %.vcf %.vcf.idx
	$(call CHECK_VCF,$<,$@,$(call LSCRIPT_CHECK_MEM,5G,00:59:59,"$(LOAD_SNP_EFF_MODULE); $(SNP_EFF) eff \
		-i vcf -o gatk $(SNP_EFF_GENOME) $< > $@"))

define annotate-sample
vcf/$1.%.ann.vcf : vcf/$1.%.vcf bam/$1.bam bam/$1.bai
	$$(call LSCRIPT_PARALLEL_MEM,4,2G,00:29:29,"$$(LOAD_JAVA8_MODULE); $(call VARIANT_ANNOTATOR,1.5G) -nt 4 -R $$(REF_FASTA) \
		$$(foreach ann,$$(VCF_ANNOTATIONS),-A $$(ann) ) --dbsnp $$(DBSNP) $$(foreach bam,$$(filter %.bam,$$^),-I $$(bam) ) -L $$< -V $$< -o $$@ && $$(RM) $$< $$<.idx")
endef
$(foreach sample,$(SAMPLES),$(eval $(call annotate-sample,$(sample))))

define hrun-sample
vcf/$1.%.hrun.vcf : vcf/$1.%.vcf bam/$1.bam bam/$1.bai
	$$(call LSCRIPT_CHECK_PARALLEL_MEM,4,2G,00:29:29,"$$(LOAD_JAVA8_MODULE); $(call VARIANT_ANNOTATOR,1.5G) -nt 4 -R $$(REF_FASTA) \
		-L $$< -A HomopolymerRun --dbsnp $$(DBSNP) $$(foreach bam,$$(filter %.bam,$$^),-I $$(bam) ) -V $$< -o $$@ && $$(RM) $$< $$<.idx")
endef
$(foreach sample,$(SAMPLES),$(eval $(call hrun-sample,$(sample))))

# Here first select the uncalled variants
# then genotype them
# then select the ones supported by at least 1 read
# then set filter to "interrogation"
# then merge

define sufam
ifeq ($3,$2)
vcf/$1_$2.%.sufam.vcf : vcf/$1_$2.%.vcf
	$$(INIT) ln -f $$< $$@
else
vcf/$3.%.sufam.tmp : $$(foreach tumor,$$(wordlist 1,$$(shell expr $$(words $$(subst _,$$( ),$3)) - 1),$$(subst _,$$( ),$3)),vcf/$$(tumor)_$$(lastword $$(subst _,$$( ),$3)).%.vcf)
	$$(call LSCRIPT_MEM,22G,03:59:59,"$$(LOAD_JAVA8_MODULE); $$(call COMBINE_VARIANTS,21G) \
		$$(foreach vcf,$$^,--variant $$(vcf) ) -o $$@ --genotypemergeoption UNSORTED -R $$(REF_FASTA)")
#ifeq ($(findstring ILLUMINA,$(SEQ_PLATFORM)),ILLUMINA)
#$(info SUFAM_METHOD,$(SUFAM_METHOD))
ifeq ($(findstring gatk,$(SUFAM_METHOD)),gatk)
vcf/$1_$2.%.sufam.vcf : vcf/$1_$2.%.vcf vcf/$3.%.sufam.tmp bam/$1.bam bam/$2.bam
	$$(call CHECK_VCF_CMD,$$(word 2,$$^),cp $$(word 1,$$^) $$@,$$(call LSCRIPT_CHECK_MEM,10G,02:59:59,"$$(LOAD_SNP_EFF_MODULE); $$(LOAD_JAVA8_MODULE); \
		$$(call SELECT_VARIANTS,21G) -R $$(REF_FASTA) --variant $$(word 2,$$^) --discordance $$(word 1,$$^) -o $$@.tmp1 && \
		$$(call UNIFIED_GENOTYPER,4G) -R $$(REF_FASTA) -I $$(word 3,$$^) -I $$(word 4,$$^) --downsampling_type NONE --dbsnp $(DBSNP_TARGETS_INTERVALS) \
		--genotyping_mode GENOTYPE_GIVEN_ALLELES --output_mode EMIT_ALL_SITES -alleles $$@.tmp1 -o $$@.tmp2 && \
		$$(call VARIANT_FILTRATION,7G) -R $$(REF_FASTA) -V $$@.tmp2 -o $$@.tmp3 \
		--filterExpression 'vc.getGenotype(\"$1\").getAnyAttribute(\"FAO\").1 > 0' --filterName interrogation && \
		$$(SNP_SIFT) filter $$(SNP_SIFT_OPTS) -f $$@.tmp3 \"(FILTER has 'interrogation')\" | $$(FIX_GATK_VCF) > $$@.tmp4 && \
		$$(call COMBINE_VARIANTS,21G) --variant $$< --variant $$@.tmp4 -o $$@ \
		--genotypemergeoption UNSORTED -R $$(REF_FASTA) && \
		$$(RM) $$@.tmp1 $$@.tmp2 $$@.tmp3 $$@.tmp4 $$(word 2,$$^) $$@.tmp1.idx $$@.tmp2.idx $$@.tmp3.idx $$@.tmp4.idx $$(word 2,$$^).idx"))
endif
#ifeq ($(findstring IONTORRENT,$(SEQ_PLATFORM)),IONTORRENT)
ifeq ($(findstring tvc,$(SUFAM_METHOD)),tvc)
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
endif
endef
$(foreach set,$(SAMPLE_SETS),\
	$(foreach tumor,$(wordlist 1,$(shell expr $(words $(subst _,$( ),$(set))) - 1),$(subst _,$( ),$(set))),\
		$(eval $(call sufam,$(tumor),$(lastword $(subst _,$( ),$(set))),$(subst $(tumor)_,,$(set))))))


%.dbsnp.vcf : %.vcf %.vcf.idx 
	$(call CHECK_VCF,$<,$@,$(call LSCRIPT_CHECK_MEM,9G,00:59:59,"$(LOAD_SNP_EFF_MODULE); $(SNP_SIFT) annotate \
		$(SNP_SIFT_OPTS) $(DBSNP) $< > $@ && $(RM) $^"))

%.cosmic.vcf : %.vcf %.vcf.idx 
	$(call CHECK_VCF,$<,$@,$(call LSCRIPT_CHECK_MEM,9G,00:29:29,"$(LOAD_SNP_EFF_MODULE); $(SNP_SIFT) annotate \
		$(SNP_SIFT_OPTS) $(COSMIC) $< > $@ && $(RM) $^"))

%.clinvar.vcf : %.vcf %.vcf.idx 
	$(call CHECK_VCF,$<,$@,$(call LSCRIPT_CHECK_MEM,9G,00:29:29,"$(LOAD_SNP_EFF_MODULE); $(SNP_SIFT) annotate \
		$(SNP_SIFT_OPTS) $(CLINVAR) $< > $@ && $(RM) $^"))

%.exac_nontcga.vcf : %.vcf %.vcf.idx 
	$(call CHECK_VCF,$<,$@,$(call LSCRIPT_CHECK_MEM,9G,00:29:29,"$(LOAD_SNP_EFF_MODULE); $(SNP_SIFT) annotate \
		$(SNP_SIFT_OPTS) -info ExAC_AF $(EXAC_NONTCGA) $< > $@ && $(RM) $^"))

%.mut_taste.vcf : %.vcf
	$(INIT) $(call CHECK_VCF,$<,$@,$(MUTATION_TASTER) $< > $@ 2> $(LOG))

%.chasm.vcf : %.vcf
	$(call CHECK_VCF,$<,$@,$(call LSCRIPT_CHECK_MEM,8G,00:59:59,"unset PYTHONPATH && source $(CHASM_PYTHON_ENV)/bin/activate $(CHASM_PYTHON_ENV) && \
		$(CHASM) --genome $(REF) --classifier $(subst $( ),$(,),$(CHASM_CLASSIFIER)) --chasmDir $(CHASM_DIR) --python $(shell which python) --outFile $@ $< && $(RM) $< $<.idx"))

%.fathmm.vcf : %.vcf
	$(call CHECK_VCF,$<,$@,$(call LSCRIPT_CHECK_MEM,8G,00:59:59,"PYTHONPATH=$(FATHMM_PYTHONPATH) $(FATHMM) $(FATHMM_OPTS) --outFile $@ $< && $(RM) $< $<.idx"))

%.mutass.vcf : %.vcf
	$(call LSCRIPT_MEM,12G,00:59:59,$(MUT_ASS) --outFile $@ --maData $(MUT_ASS_RDATA) $<)

%.transfic.vcf : %.vcf
	$(call CHECK_VCF,$<,$@,$(call LSCRIPT_MEM,9G,00:59:59,"$(TRANSFIC) --genome $(REF) --transfic $(TRANSFIC_PERL_SCRIPT) --outFile $@ $< && $(RM) $< $<.idx"))

%.provean.vcf : %.vcf
	$(call LSCRIPT_MEM,8G,00:29:29,"$(PROVEAN) $(PROVEAN_OPTS) --outFile $@ $<")

%.pathogenic.vcf : %.vcf
	$(INIT) $(call CHECK_VCF,$<,$@,$(CLASSIFY_PATHOGENICITY) $< > $@ 2> $(LOG))

%.gene_ann.vcf : %.vcf
	$(call LSCRIPT_CHECK_MEM,4G,00:29:29,"$(LOAD_R_MODULE); $(ADD_GENE_LIST_ANNOTATION) --genome $(REF) \
		--geneBed $(HAPLOTYPE_INSUF_BED)$(,)$(CANCER_GENE_CENSUS_BED)$(,)$(KANDOTH_BED)$(,)$(LAWRENCE_BED)$(,)$(DGD_BED) \
		--name hap_insuf$(,)cancer_gene_census$(,)kandoth$(,)lawrence$(,)duplicatedGenesDB --outFile $@ $< && $(RM) $< $<.idx")

%.$(ANNOVAR_REF)_multianno.vcf : %.vcf
	$(call LSCRIPT_CHECK_MEM,7G,9G,"$(ANNOVAR) -out $* $(ANNOVAR_OPTS) $< $(ANNOVAR_DB) && $(RM) $< $<.idx")

ifdef SAMPLE_PAIRS
define annotate-facets-pair
vcf/$1_$2.%.facets.vcf : vcf/$1_$2.%.vcf facets/cncf/$1_$2.Rdata
	$$(call LSCRIPT_MEM,4G,00:29:59,"$$(LOAD_R_MODULE); $$(ANNOTATE_FACETS_VCF) --genome \"$$(REF)\" --tumor \"$1\" --facetsRdata $$(<<) --outFile $$@ $$<")
endef
$(foreach pair,$(SAMPLE_PAIRS),$(eval $(call annotate-facets-pair,$(tumor.$(pair)),$(normal.$(pair)))))
endif


######### WHAT THE HELL IS THIS???###########
%.norm.vcf.gz : %.vcf
	$(call LSCRIPT_MEM,9G,12G,"sed '/^##GATKCommandLine/d;/^##MuTect/d;' $< | \
		$(VT) view -h -f PASS - | \
		$(VT) decompose -s - | \
		$(VT) normalize -r $(REF_FASTA) - | \
		$(call SNP_EFF_MEM,8G) ann -c $(SNP_EFF_CONFIG) $(SNP_EFF_GENOME) -formatEff -classic | \
		bgzip -c > $@")






######### WHAT THE HELL IS THIS???###########
ifdef SAMPLE_SET_PAIRS
define somatic-filter-vcf-set
vcf/$1.%.som_ft.vcf : vcf/$1.%.vcf
	$$(INIT) $$(SOMATIC_FILTER_VCF) -n $(normal.$1) -f 0.03 $$< > $$@ 2> $$(LOG) && $$(RM) $$< $$<.idx
endef
$(foreach set,$(SAMPLE_SET_PAIRS),$(eval $(call somatic-filter-vcf-set,$(set))))
endif

      
      

      
      
ifdef SAMPLE_PAIRS
ifeq ($(findstring ILLUMINA,$(SEQ_PLATFORM)),ILLUMINA)
define som-ad-ft-tumor-normal
vcf/$1_$2.%.som_ad_ft.vcf : vcf/$1_$2.%.vcf
	$$(call LSCRIPT_CHECK_MEM,8G,00:59:59,"$$(LOAD_JAVA8_MODULE); $$(call VARIANT_FILTRATION,7G) -R $$(REF_FASTA) -V $$< -o $$@ \
		--filterExpression 'vc.getGenotype(\"$1\").getAD().1 < $(MIN_TUMOR_AD)' \
  		--filterName tumorVarAlleleDepth \
		--filterExpression 'if (vc.getGenotype(\"$2\").getDP() > $(MIN_NORMAL_DEPTH)) { \
			( vc.getGenotype(\"$2\").getAD().1 * 1.0 / vc.getGenotype(\"$2\").getDP()) > ( vc.getGenotype(\"$1\").getAD().1 * 1.0 / vc.getGenotype(\"$1\").getDP()) / $(MIN_TN_AD_RATIO) \
			} else { ( vc.getGenotype(\"$2\").getAD().1 * 1.0 / vc.getGenotype(\"$2\").getDP()) > ( vc.getGenotype(\"$1\").getAD().1 * 1.0 / vc.getGenotype(\"$1\").getDP()) / $(MIN_TN_AD_RATIO) && \
			vc.getGenotype(\"$1\").getAD().1 * 1.0 < 1 && vc.getGenotype(\"$2\").getAD().1 > 1 }' \
		--filterName somaticAlleleDepth \
		--filterExpression 'vc.getGenotype(\"$1\").getDP() <= $$(MIN_TUMOR_DEPTH) || vc.getGenotype(\"$2\").getDP() <= $$(MIN_NORMAL_DEPTH)' \
		--filterName depthFilter && sed -i 's/getGenotype(\"\([^\"]*\)\")/getGenotype(\1)/g' $$@ && $$(RM) $$< $$<.idx")

# somatic filter for structural variants
vcf/$1_$2.%.sv_som_ft.vcf : vcf/$1_$2.%.vcf
	$$(call LSCRIPT_CHECK_MEM,8G,00:59:59,"$$(LOAD_JAVA8_MODULE); $$(call VARIANT_FILTRATION,7G) -R $$(REF_FASTA) -V $$< -o $$@ \
		--filterExpression 'vc.getGenotype(\"$1\").getAnyAttribute(\"SU\") <= $$(DEPTH_FILTER)' \
		--filterName svSupport \
		--filterExpression 'vc.getGenotype(\"$1\").getAnyAttribute(\"SU\") < 5 * vc.getGenotype(\"$2\").getAnyAttribute(\"SU\")' \
		--filterName somaticSvSupport \
		&& sed -i 's/getGenotype(\"\([^\"]*\)\")/getGenotype(\1)/g' $$@ && $$(RM) $$< $$<.idx")
endef
endif

ifeq ($(findstring IONTORRENT,$(SEQ_PLATFORM)),IONTORRENT)
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
endif
$(foreach pair,$(SAMPLE_PAIRS),$(eval $(call som-ad-ft-tumor-normal,$(tumor.$(pair)),$(normal.$(pair)))))

define rename-samples-tumor-normal
vcf/$1_$2.%.rn.vcf : vcf/$1_$2.%.vcf
	$$(INIT) $(LOAD_PERL_MODULE); perl -ne 'if (/^#CHROM/) { s/NORMAL/$2/; s/TUMOR/$1/; } print;' $$< > $$@ && $$(RM) $$< $$<.idx
endef
$(foreach pair,$(SAMPLE_PAIRS),$(eval $(call rename-samples-tumor-normal,$(tumor.$(pair)),$(normal.$(pair)))))


define ad-tumor-normal
vcf/$1_$2.%.ad.vcf : vcf/$1_$2.%.vcf bam/$1.bam bam/$2.bam bam/$1.bai bam/$2.bai
	$$(call LSCRIPT_CHECK_PARALLEL_MEM,4,2G,00:29:29,"$$(LOAD_JAVA8_MODULE); $$(call VARIANT_ANNOTATOR,1.5G) -nt 4 -R $$(REF_FASTA) \
		-A DepthPerAlleleBySample --dbsnp $$(DBSNP) $$(foreach bam,$$(filter %.bam,$$^),-I $$(bam) ) -V $$< -o $$@ -L $$<")
endef
$(foreach pair,$(SAMPLE_PAIRS),\
	$(eval $(call ad-tumor-normal,$(tumor.$(pair)),$(normal.$(pair)))))

define annotate-tumor-normal
vcf/$1_$2.%.ann.vcf : vcf/$1_$2.%.vcf bam/$1.bam bam/$2.bam bam/$1.bai bam/$2.bai
	$$(call LSCRIPT_CHECK_PARALLEL_MEM,4,2G,00:29:29,"$$(LOAD_JAVA8_MODULE); $$(call VARIANT_ANNOTATOR,1.5G) -nt 4 -R $$(REF_FASTA) \
		$$(foreach ann,$$(VCF_ANNOTATIONS),-A $$(ann) ) --dbsnp $$(DBSNP) $$(foreach bam,$$(filter %.bam,$$^),-I $$(bam) ) -V $$< -o $$@ -L $$< && $$(RM) $$< $$<.idx")
endef
$(foreach pair,$(SAMPLE_PAIRS),\
		$(eval $(call annotate-tumor-normal,$(tumor.$(pair)),$(normal.$(pair)))))

define hrun-tumor-normal
vcf/$1_$2.%.hrun.vcf : vcf/$1_$2.%.vcf bam/$1.bam bam/$2.bam bam/$1.bai bam/$2.bai
	$$(call LSCRIPT_CHECK_PARALLEL_MEM,4,2G,00:29:29,"$$(LOAD_JAVA8_MODULE); $$(call VARIANT_ANNOTATOR,1.5G) -nt 4 -R $$(REF_FASTA) \
		-A HomopolymerRun --dbsnp $$(DBSNP) $$(foreach bam,$$(filter %.bam,$$^),-I $$(bam) ) -V $$< -L $$< -o $$@ && $$(RM) $$< $$<.idx")
endef
$(foreach pair,$(SAMPLE_PAIRS),$(eval $(call hrun-tumor-normal,$(tumor.$(pair)),$(normal.$(pair)))))
endif

# extract vcf to table
tables/%.opl_tab.txt : vcf/%.vcf
	$(call LSCRIPT_CHECK_MEM,2G,00:29:29,"$(LOAD_SNP_EFF_MODULE); \
	format_fields=\$$(grep '^##FORMAT=<ID=' $< | sed 's/dbNSFP_GERP++/dbNSFP_GERP/g' | sed 's/.*ID=//; s/,.*//;' | tr '\n' ' '); \
	N=\$$(expr \$$(grep '^#CHROM' $< | wc -w) - 10); \
	fields='$(VCF_FIELDS)'; \
	for f in \$$format_fields; do \
		for i in \$$(seq 0 \$$N); do \
			fields+=' 'GEN[\$$i].\$$f; \
		done; \
	done; \
	fields+=' '\$$(grep '^##INFO=<ID=' $< | grep -v '=REF,' | sed 's/dbNSFP_GERP++/dbNSFP_GERP/g' | \
		sed 's/.*ID=//; s/,.*//; s/\bANN\b/$(ANN_FIELDS)/; ' | tr '\n' ' '); \
	$(LOAD_PERL_MODULE); $(VCF_EFF_ONE_PER_LINE) < $< | sed 's/dbNSFP_GERP++/dbNSFP_GERP/g' | \
		$(SNP_SIFT) extractFields - \$$fields > $@; \
	for i in \`seq 0 \$$N\`; do \
	S=\$$(grep '^#CHROM' $< | cut -f \$$((\$$i + 10))); \
	sed -i \"1s/GEN\[\$$i\]/\$$S/g;\" $@; \
	done")

%.tab.txt : %.opl_tab.txt
	$(call LSCRIPT_MEM,8G,00:59:59,"$(LOAD_PERL_MODULE); $(VCF_JOIN_EFF) < $< > $@")
	
%.pass.txt : %.txt
	$(INIT) head -1 $< > $@ && awk '$$6 == "PASS" { print }' $< >> $@ || true


# merge tables
alltables/all.%.txt : $(foreach sample,$(SAMPLES),tables/$(sample).%.txt)
	$(call LSCRIPT_MEM,5G,00:29:29,"$(LOAD_R_MODULE); $(RSCRIPT) $(RBIND) --sampleName $< $^ > $@")
ifdef SAMPLE_SETS
alltables/allSS.%.txt : $(foreach set,$(SAMPLE_SETS),tables/$(set).%.txt)
	$(call LSCRIPT_MEM,5G,00:29:29,"$(LOAD_R_MODULE); $(RSCRIPT) $(RBIND) --normalLast $^ > $@")
endif
ifdef SAMPLE_PAIRS
alltables/allTN.%.txt : $(foreach pair,$(SAMPLE_PAIRS),tables/$(pair).%.txt)
	$(call LSCRIPT_MEM,5G,00:29:29,"$(LOAD_R_MODULE); $(RSCRIPT) $(RBIND) --tumorNormal $^ > $@")
endif

%.high_moderate.txt : %.txt
	col=$$(head -1 $< | tr '\t' '\n' | grep -n "IMPACT" | sed 's/:.*//'); \
	$(INIT) head -1 $< > $@ && awk -v col=$$col 'match($$col, /MODERATE/) || match($$col, /HIGH/)' $< >> $@

%.low_modifier.txt : %.txt
	col=$$(head -1 $< | tr '\t' '\n' | grep -n "IMPACT" | sed 's/:.*//'); \
	$(INIT) head -1 $< > $@ && awk -v col=$$col '! (match($$col, /MODERATE/) || match($$col, /HIGH/)) && (match($$col, /LOW/) || match($$col,/MODIFIER/))' $< >> $@

%.synonymous.txt : %.txt
	col_imp=$$(head -1 $< | tr '\t' '\n' | grep -n "IMPACT" | sed 's/:.*//'); \
	col_eff=$$(head -1 $< | tr '\t' '\n' | grep -n "EFFECT" | sed 's/:.*//'); \
	$(INIT) head -1 $< > $@ && awk -v col_imp=$$col_imp -v col_eff=$$col_eff '! (match($$col_imp, /MODERATE/) || match($$col_imp, /HIGH/)) && (match($$col_imp, /LOW/) && (match($$col_eff, /synonymous_variant/)))' $< >> $@

%.nonsynonymous.txt : %.txt
	col=$$(head -1 $< | tr '\t' '\n' | grep -n "IMPACT" | sed 's/:.*//'); \
	$(INIT) head -1 $< > $@ && awk -v col=$$col 'match($$col, /MODERATE/) || match($$col, /HIGH/)' $< >> $@

%.nonsynonymous_synonymous_hotspot.txt : %.txt
	col_imp=$$(head -1 $< | tr '\t' '\n' | grep -n "IMPACT" | sed 's/:.*//'); \
	col_eff=$$(head -1 $< | tr '\t' '\n' | grep -n "EFFECT" | sed 's/:.*//'); \
	col_filter=$$(head -1 $< | tr '\t' '\n' | grep -n "FILTER"); \
	$(INIT) head -1 $< > $@ && awk -v col_imp=$$col_imp -v col_eff=$$col_eff -v col_filter=$$col_filter \
	'(match($$col_imp, /MODERATE/) || match($$col_imp, /HIGH/) || (match($$col_imp, /LOW/) && match($$col_eff, /synonymous_variant/)) || match($$col_filter, /HOTSPOT/))' $< >> $@

%.nonsynonymous_synonymous_hotspot_lincRNA.txt : %.txt
	col_imp=$$(head -1 $< | tr '\t' '\n' | grep -n "IMPACT" | sed 's/:.*//'); \
	col_eff=$$(head -1 $< | tr '\t' '\n' | grep -n "EFFECT" | sed 's/:.*//'); \
	col_biotype=$$(head -1 $< | tr '\t' '\n' | grep -n "BIOTYPE" | sed 's/:.*//'); \
	col_filter=$$(head -1 $< | tr '\t' '\n' | grep -n "FILTER"); \
	$(INIT) head -1 $< > $@ && awk -v col_imp=$$col_imp -v col_eff=$$col_eff -v col_biotype=$$col_biotype -v col_filter=$$col_filter \
	'(match($$col_imp, /MODERATE/) || match($$col_imp, /HIGH/) || (match($$col_imp, /LOW/) && match($$col_eff, /synonymous_variant/)) || (match($$col_biotype, /lincRNA/) && match($$col_eff, /non_coding_exon/)) || match($$col_filter, /HOTSPOT/))' $< >> $@







#################### MAF ###################
ifdef SAMPLE_PAIRS
define vcf2maf-tumor-normal
maf/$1_$2.%.maf : vcf/$1_$2.%.vcf
	$$(call LSCRIPT_MEM,9G,12G,"$$(VCF2MAF) --input-vcf $$< --tumor-id $1 --normal-id $2 --ref-fasta $$(REF_FASTA) --vep-path $$(VEP_PATH) --vep-data $$(VEP_DATA) --output-maf $$@")
endef
$(foreach pair,$(SAMPLE_PAIRS),$(eval $(call vcf2maf-tumor-normal,$(tumor.$(pair)),$(normal.$(pair)))))

allmaf/allTN.%.maf : $(foreach pair,$(SAMPLE_PAIRS),maf/$(pair).%.maf)
	$(INIT) \
	{ \
	grep -v '^#' $< | sed -n 1p; \
	for x in $^; do grep -v '^#' $$x | sed 1d; done \
	} > $@
endif

define vcf2maf-sample
maf/$1.%.maf : vcf/$1.%.vcf
	$$(call LSCRIPT_MEM,9G,12G,"$$(VCF2MAF) --input-vcf $$< --tumor-id $1 --ref-fasta $$(REF_FASTA) --vep-path $$(VEP_PATH) --vep-data $$(VEP_DATA) --output-maf $$@")
endef
$(foreach sample,$(SAMPLES),$(eval $(call vcf2maf-sample,$(sample))))

allmaf/all.%.maf : $(foreach sample,$(SAMPLES),maf/$(sample).%.maf)
	$(INIT) \
	{ \
	sed -n 2p $<; \
	sed 1,2d $^; \
	} > $@
endif

# mouse genome project dbsnp
#%.mgp_dbsnp.vcf : %.vcf %.vcf.idx 
#	$(call CHECK_VCF,$<,$@,$(call LSCRIPT_CHECK_MEM,33G,03:59:59,"$(LOAD_SNP_EFF_MODULE); $(SNP_SIFT) annotate \
#		-tabix $(SNP_SIFT_OPTS) $(MGP_SNP_DBSNP) $< | $(call SNP_SIFT_MEM,10G) annotate \
#		-tabix $(SNP_SIFT_OPTS) $(MGP_INDEL_DBSNP) > $@ && $(RM) $^"))

# post-annotation filter
#VCF_POST_ANN_FILTER_EXPRESSION ?= ExAC_AF > 0.1
#%.cft.vcf : %.vcf
#	$(call LSCRIPT_CHECK_MEM,8G,00:59:59,"$(LOAD_JAVA8_MODULE); $(VARIANT_FILTRATION) \
#		-R $(REF_FASTA) -V $< -o $@ --filterExpression '$(VCF_POST_ANN_FILTER_EXPRESSION)' --filterName customFilter && $(RM) $<")

%.sdp_ft.vcf : %.vcf
	$(call CHECK_VCF,$<,$@,$(call LSCRIPT_CHECK_MEM,4G,00:59:59,"$(LOAD_SNP_EFF_MODULE); $(SNP_SIFT) filter \
		$(SNP_SIFT_OPTS) -f $< '(exists GEN[*].DP) & (GEN[*].DP > 20)' > $@"))

# Copy number regulated genes annotated per subtype
# FYI Endometrioid_MSI-L has no copy number regulated genes
#CN_ENDOMETRIAL_SUBTYPES = CN_high CN_low Endometrioid_MSI_H Endometrioid_MSS Endometrioid MSI POLE Serous
#CN_BREAST_SUBTYPES = ER_negative ER_positive HER2_postitive Pam50_Basal Pam50_Her2 Pam50_LumA Pam50_LumB Pam50_Normal Triple_negative
#CN_ENDOMETRIAL_BED = $(foreach set,$(CN_ENDOMETRIAL_SUBTYPES), $(HOME)/share/reference/annotation_gene_lists/cn_reg/endometrial/copy_number_regulated_genes_subtype_$(set)_spearmanrsquare0.4_fdrbh_adjp_lt0.05.HUGO.bed)
#CN_BREAST_BED = $(foreach set,$(CN_BREAST_SUBTYPES), $(HOME)/share/reference/annotation_gene_lists/cn_reg/breast/metabric_subtype_$(set)_copy_number_regulated_genes_std0.5_spearmanrsquare0.4_fdrbh_adjp_lt0.05.HUGO.bed)
#%.cn_reg.vcf : %.vcf
#	$(call LSCRIPT_MEM,8G,12G,"$(ADD_GENE_LIST_ANNOTATION) --genome $(REF) --geneBed \
#        $(subst $(space),$(,),$(strip $(CN_ENDOMETRIAL_BED)) $(strip $(CN_BREAST_BED))) --name $(subst $(space),$(,),$(foreach set,$(strip $(CN_ENDOMETRIAL_SUBTYPES)),endometrial_#$(set)) $(foreach set,$(strip $(CN_BREAST_SUBTYPES)),breast_$(set))) --outFile $@ $< && $(RM) $< $<.idx")

#-cancer does nothing
#%.som_eff.vcf : %.vcf %.vcf.pair
#	$(call CHECK_VCF,$<,$@,$(call LSCRIPT_CHECK_MEM,9G,14G,"$(call SNP_EFF_MEM,8G) ann -cancer -cancerSamples $(<<) $(SNP_EFF_OPTS) $(SNP_EFF_GENOME) $< > $@ && $(RM) $^"))

# VariantEval: generate vcf report
#reports/%.grp : $(foreach sample,$(SAMPLES),vcf/$(sample).%.vcf) $(foreach sample,$(SAMPLES),vcf/$(sample).%.vcf.idx)
#	$(call LSCRIPT_MEM,2G,5G,"$(call GATK_MEM,2G) -T VariantEval $(foreach sm,$(REPORT_STRATIFICATION), --stratificationModule $(sm)) -R $(REF_FASTA) --dbsnp $(DBSNP) $(foreach eval,$(filter %.vcf,$^), --eval:$(call strip-suffix,$(notdir $(eval))) $(eval)) -o $@")
#ifdef SAMPLE_PAIRS
#reports/%.grp : $(foreach pair,$(SAMPLE_PAIRS),vcf/$(pair).%.vcf vcf/$(pair).%.vcf.idx)
#	$(call LSCRIPT_MEM,2G,5G,"$(call GATK_MEM,2G) -T VariantEval $(foreach sm,$(REPORT_STRATIFICATION), --stratificationModule $(sm)) -R $(REF_FASTA) --dbsnp $(DBSNP) $(foreach eval,$(filter %.vcf,$^), --eval:$(call strip-suffix,$(notdir $(eval))) $(eval)) -o $@")
#endif

# apply dp filter for somatic sniper
#%.ss_dp_ft.vcf : %.vcf
#	$(call LSCRIPT_CHECK_MEM,8G,12G,"$(call GATK_MEM,8G) -T VariantFiltration -R $(REF_FASTA) -V $< -o $@ --filterExpression 'vc.getGenotype(\"TUMOR\").getDP() < $(DEPTH_FILTER) || vc.getGenotype(\"NORMAL\").getDP() < $(DEPTH_FILTER)' --filterName depthFilter && $(RM) $< $<.idx")

# varscan TN variant allele frequency: min tumor freq > 5% ; max normal freq < 5%
#%.freq_ft.vcf : %.vcf
#	$(call LSCRIPT_CHECK_MEM,2G,5G,"sed '/##FORMAT=<ID=FREQ/ s/String/Float/; /^#/! s/%//g' $< | $(call SNP_SIFT_MEM,2G) filter $(SNP_SIFT_OPTS) '(exists GEN[*].FREQ) & (GEN[0].FREQ < 5) & (GEN[1].FREQ[0] > 5)' > $@ && $(RM) $< $<.idx")

#%.vaf_ft.vcf : %.vcf
#	$(call LSCRIPT_CHECK_MEM,2G,5G,"$(call SNP_SIFT_MEM,2G) filter $(SNP_SIFT_OPTS) '(exists GEN[*].VAF) & (GEN[0].VAF > 0.05) & (GEN[1].VAF < 0.05)' < $< > $@ && $(RM) $< $<.idx")


# varscan depth filter (b/c varscan is dumb and only gives variant depth)
#%.vdp_ft.vcf : %.vcf
#	$(call LSCRIPT_CHECK_MEM,2G,5G,"cat $< | $(call SNP_SIFT_MEM,2G) filter $(SNP_SIFT_OPTS) '(exists GEN[*].AD) & (GEN[*].AD > $(DEPTH_FILTER))' > $@ && $(RM) $< $<.idx")

# add exon distance
#%.exondist.vcf : %.vcf
#	$(call LSCRIPT_CHECK_MEM,2G,3G,"$(INTRON_POSN_LOOKUP) $< > $@")

#%.common_ft.vcf : %.vcf
#	$(call LSCRIPT_MEM,4G,5G,"$(COMMON_FILTER_VCF) $< > $@")


#%.fp_ft.vcf : %.vcf
#	$(call LSCRIPT_MEM,8G,12G,"$(call GATK_MEM,8G) -T VariantFiltration \
#	-R $(REF_FASTA) -V $< -o $@ --maskName 'FuentesFalsePositive' --mask $(FALSE_POSITIVE_BED) && $(RM) $< $<.idx")


#%.hrun_ft.vcf : %.vcf
#	$(call LSCRIPT_CHECK_MEM,8G,01:59:59,"$(LOAD_JAVA8_MODULE); $(VARIANT_FILTRATION) \
#		-R $(REF_FASTA) -V $< -o $@ --filterExpression 'HRun > $(HRUN_FILTER)' --filterName HRun && $(RM) $< $<.idx")

# workaround to do a double pass
#%.pass2.vcf : %.vcf
#	$(call CHECK_VCF,$<,$@,$(call LSCRIPT_CHECK_MEM,2G,00:29:29,"$(LOAD_SNP_EFF_MODULE); $(SNP_SIFT) filter \
#		$(SNP_SIFT_OPTS) -f $< \"( na FILTER ) | (FILTER = 'PASS')\" > $@"))


VCFTOOLS_MK = true

#include usb-modules/copy_number/facets.mk

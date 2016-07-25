ifndef SOMATIC_VARIANT_CALLER_INC

### Filters
SOMATIC_FILTERS = $(if $(TARGETS_FILE_INTERVALS),target_ft) som_ad_ft
#SOMATIC_FILTERS += $(if $(findstring strelka,$1),,som_ad_ft)
SOMATIC_FILTERS += hotspot
#SOMATIC_FILTERS += ifeq ($(shell test $(words $(NORMAL_SAMPLES)) -gt 4),nft)

SOMATIC_FILTER_SUFFIX = $(subst $( ),.,$(strip $(call SOMATIC_FILTERS,$1)))

### Annotation
SOMATIC_COMMON_ANNS := $(if $(findstring GRCm38,$(REF)),mgp_dbsnp, \
        eff dbsnp cosmic exac_nontcga clinvar gene_ann nsfp)

ifeq ($(shell test $(words $(NORMAL_SAMPLES)) -gt 4; echo $$?),0)
	SOMATIC_COMMON_ANNS += nft
endif

ifeq ($(findstring true,$(ANN_FACETS)),true)
	SOMATIC_COMMON_ANNS += facets
endif

SOMATIC_SNV_ANNS = $(SOMATIC_COMMON_ANNS)
SOMATIC_INDEL_ANNS = $(SOMATIC_COMMON_ANNS)

#hrun $(if $(findstring b37,$(REF)),provean mut_taste)

SOMATIC_ANNS = $(if $(findstring indel,$1),$(SOMATIC_INDEL_ANNS),$(SOMATIC_SNV_ANNS))

#pathogenic
SOMATIC_ANN_SUFFIX = $(subst $( ),.,$(strip pass $(call SOMATIC_ANNS,$1)))

# filters that are run after first round of annotations
#SOMATIC_POST_FILTERS += cft common_ft
#SOMATIC_POST_FILTER_SUFFIX = $(subst $( ),.,$(strip $(call SOMATIC_POST_FILTERS,$1) pass2))
#SOMATIC_POST_ANNS += $(if $(findstring b37,$(REF)),$(if $(findstring indel,$1),mut_taste,nsfp chasm fathmm))
#SOMATIC_POST_ANNS += $(if $(findstring b37,$(REF)),gene_ann clinvar hotspot)
#SOMATIC_POST_ANN_SUFFIX = $(subst $( ),.,$(strip $(call SOMATIC_POST_ANNS,$1)))

SOMATIC_SUFFIX = $(subst $( ),.,$(strip $1 $(SOMATIC_FILTER_SUFFIX) $(SOMATIC_ANN_SUFFIX)))
SOMATICGT_SUFFIX = $(subst $( ),.,$(strip $1 $(SOMATIC_FILTER_SUFFIX) pass eff))

SOMATIC_EFF_TYPES := high_moderate low_modifier synonymous nonsynonymous
SOMATIC_VCF_SUFFIXES = $(foreach type,$1,$(call SOMATIC_SUFFIX,$(type)))
SOMATICGT_VCF_SUFFIXES = $(foreach type,$1,$(call SOMATICGT_SUFFIX,$(type)))

SOMATIC_TABLE_EXTENSIONS = tab $(foreach eff,$(SOMATIC_EFF_TYPES),tab.$(eff))
SOMATIC_TABLE_SUFFIXES = $(foreach suff,$(call SOMATIC_VCF_SUFFIXES,$1),\
				 $(foreach ext,$(SOMATIC_TABLE_EXTENSIONS),$(suff).$(ext)))
SOMATICGT_TABLE_SUFFIXES = $(foreach suff,$(call SOMATICGT_VCF_SUFFIXES,$1),\
				 $(foreach ext,$(SOMATIC_TABLE_EXTENSIONS),$(suff).$(ext)))

SOMATIC_VCFS = $(foreach pair,$(SAMPLE_PAIRS),\
    $(foreach suff,$(call SOMATIC_VCF_SUFFIXES,$1), \
        vcf/$(pair).$(suff).vcf))
SOMATIC_TABLES = $(foreach pair,$(SAMPLE_PAIRS),\
    $(foreach suff,$(call SOMATIC_TABLE_SUFFIXES,$1),\
		tables/$(pair).$(suff).txt))
SOMATIC_TABLES += $(foreach suff,$(call SOMATIC_TABLE_SUFFIXES,$1),alltables/allTN.$(suff).txt)

SOMATIC_VCFS_HOTSPOTGT = $(foreach pair,$(SAMPLE_PAIRS),\
	$(foreach suff,$(call SOMATICGT_VCF_SUFFIXES,$1),\
		vcf/$(pair).$(suff).vcf))

SOMATIC_TABLES_HOTSPOTGT = $(foreach pair,$(SAMPLE_PAIRS),\
	$(foreach suff,$(call SOMATICGT_TABLE_SUFFIXES,$1),\
		tables/$(pair).$(suff).txt))
SOMATIC_TABLES_HOTSPOTGT += $(foreach suff,$(call SOMATICGT_TABLE_SUFFIXES,$1),alltables/allTN.$(suff).txt)

endif
SOMATIC_VARIANT_CALLER_INC = true

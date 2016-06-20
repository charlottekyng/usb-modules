ifndef VARIANT_CALLER_INC


### Filters
FILTERS = pass
FILTERS += $(if $(TARGETS_FILE_INTERVALS),target_ft)
FILTERS += dp_ft
FILTERS += $(if $(NORMAL_VCF),nft)
FILTERS += pass

### Annotation
COMMON_ANNS := $(if $(findstring GRCm38,$(REF)),mgp_dbsnp, \
        eff dbsnp cosmic exac_nontcga clinvar gene_ann nsfp)

SNV_ANNS = $(COMMON_ANNS) $(if $(findstring b37,$(REF)),chasm fathmm) 

INDEL_ANNS = $(COMMON_ANNS) hrun $(if $(findstring b37,$(REF)),provean mut_taste)

ANNOTATIONS = $(if $(findstring indel,$1),$(INDEL_ANNS),$(SNV_ANNS)) pathogenic

#POSTANN_FILTERS += cft common

FILTER_SUFFIX = $1.$(subst $( ),.,$(strip $(call FILTERS,$1)))
ANN_SUFFIX = pass.$(subst $( ),.,$(strip $(call ANNOTATIONS,$1)))
#POSTANN_SUFFIX = $(subst $( ),.,$(strip $(call POSTANN_FILTERS,$1))).pass2
SUFFIX = $(FILTER_SUFFIX).$(ANN_SUFFIX) #.$(POSTANN_SUFFIX)

EFF_TYPES := high_moderate low_modifier synonymous nonsynonymous

VCF_SUFFIXES = $(foreach type,$1,$(call SUFFIX,$(type)))

TABLE_EXTENSIONS = tab $(foreach eff,$(EFF_TYPES),tab.$(eff))

TABLE_SUFFIXES = $(foreach suff,$(call VCF_SUFFIXES,$1),\
				 $(foreach ext,$(TABLE_EXTENSIONS),$(suff).$(ext)))

VCFS = $(foreach sample,$(SAMPLES),$(foreach suff,$(call VCF_SUFFIXES,$1),vcf/$(sample).$(suff).vcf))
TABLES = $(foreach sample,$(SAMPLES),$(foreach suff,$(call TABLE_SUFFIXES,$1),tables/$(sample).$(suff).txt))
TABLES += $(foreach suff,$(call TABLE_SUFFIXES,$1),alltables/all.$(suff).txt)

endif
VARIANT_CALLER_INC = true

ifndef VARIANT_CALLER_INC

### Filters
FILTERS += $(if $(TARGETS_FILE_INTERVALS),target_ft) dp_ft
FILTERS += $(if $(findstring GRCm38,$(REF)),,hotspot)

FILTER_SUFFIX = $(subst $( ),.,$(strip $(call FILTERS,$1)))

### Annotation

COMMON_ANNS := $(if $(findstring GRCm38,$(REF)),eff dbsnp, \
	eff dbsnp cosmic exac_nontcga clinvar gene_ann nsfp)

SNV_ANNS = $(COMMON_ANNS)
INDEL_ANNS = $(COMMON_ANNS)

ANNS = $(if $(findstring indel,$1),$(INDEL_ANNS),$(SNV_ANNS))
ANN_SUFFIX = $(subst $( ),.,$(strip pass $(call ANNS,$1)))

#hrun $(if $(findstring b37,$(REF)),provean mut_taste)

SUFFIX = $(subst $( ),.,$(strip $1 $(FILTER_SUFFIX) $(ANN_SUFFIX)))

EFF_TYPES := high_moderate low_modifier synonymous nonsynonymous
VCF_SUFFIXES = $(foreach type,$1,$(call SUFFIX,$(type)))

TABLE_EXTENSIONS = tab $(foreach eff,$(EFF_TYPES),tab.$(eff))
TABLE_SUFFIXES = $(foreach suff,$(call VCF_SUFFIXES,$1),\
	$(foreach ext,$(TABLE_EXTENSIONS),$(suff).$(ext)))

VCFS = $(foreach sample,$(SAMPLES),\
	$(foreach suff,$(call VCF_SUFFIXES,$1), \
	vcf/$(sample).$(suff).vcf))
TABLES = $(foreach sample,$(SAMPLES),\
	$(foreach suff,$(call TABLE_SUFFIXES,$1),\
	tables/$(sample).$(suff).txt))
TABLES += $(foreach suff,$(call TABLE_SUFFIXES,$1),alltables/all.$(suff).txt)

endif
VARIANT_CALLER_INC = true

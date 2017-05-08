include usb-modules/Makefile.inc
include usb-modules/config.inc

LOGDIR = log/star_fusion.$(NOW)

PHONY += star_fusion
.PHONY : $(PHONY)
.DELETE_ON_ERROR:

star_fusion : $(foreach sample,$(SAMPLES),star_fusion/$(sample)/star-fusion.STAR-Fusion.filter.ok)

define star-fusion
star_fusion/$1/star-fusion.STAR-Fusion.filter.ok : star/secondpass/$1.Chimeric.out.junction
	$$(call LSCRIPT_MEM,8G,00:29:59,"$$(LOAD_STAR_FUSION_MODULE); $$(STAR_FUSION) \
		--genome_lib_dir $$(STAR_CTAT_DIR) \
		-J $$< --output_dir star_fusion/$1 --max_promiscuity $$(STAR_FUSION_MAX_PROMISCUITY) \
		--min_novel_junction_support $$(STAR_FUSION_MIN_JUNC_SUPP) --min_alt_pct_junction $$(STAR_FUSION_MIN_ALT_PCT_JUNC)")
endef
$(foreach sample,$(SAMPLES),$(eval $(call star-fusion,$(sample))))

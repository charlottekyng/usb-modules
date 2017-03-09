include usb-modules/Makefile.inc
include usb-modules/config.inc

LOGDIR = log/star_fusion.$(NOW)


PHONY += star_fusion
star_fusion : $(foreach sample,$(SAMPLES),star_fusion/$(sample).star_fusion_timestamp)

star_fusion/%.star_fusion_timestamp : star/secondpass/%.Chimeric.out.junction fastq/%.1.fastq.gz $(if $(findstring true,$(PAIRED_END)),fastq/%.2.fastq.gz)
	$(call LSCRIPT_MEM,8G,05:59:59,"$(LOAD_STAR_FUSION_MODULE); $(STAR_FUSION) \
		--genome_lib_dir $(STAR_CTAT_DIR) \
		-J $< --verbose_level 2 \
		--left_fq $(<<) $(if $(findstring true,$(PAIRED_END)),--right_fq $(<<<)) && touch $@")

.PHONY: $(PHONY)
.SECONDARY: 
.DELETE_ON_ERROR:

include usb-modules/aligners/starAligner.mk
include usb-modules/fastq_tools/mergeSplitFastq.mk

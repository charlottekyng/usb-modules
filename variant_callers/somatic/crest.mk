# Run somatic sniper on tumour-normal matched pairs
# Detect point mutations
##### DEFAULTS ######

##### MAKE INCLUDES #####
include modules/Makefile.inc

LOGDIR = log/crest.$(NOW)
CREST_DIR = $(HOME)/share/usr/crest
CREST = PERL5LIB=$(PERL5LIB):$(CREST_DIR) $(PERL) $(CREST_DIR)/CREST.pl

EXTRACT_SCLIP = PERL5LIB=$(PERL5LIB):$(CREST_DIR) $(PERL) $(CREST_DIR)/extractSClip.pl

.SECONDARY:
.DELETE_ON_ERROR:
.PHONY: all

ifdef SAMPLE_PAIRS
all : $(foreach pair,$(SAMPLE_PAIRS),crest/sv/$(pair).predSV.txt)
else
all : $(foreach sample,$(SAMPLES),crest/sv/$(sample).predSV.txt)
endif

define sclip-chr
crest/sclip/%.$1.sclip.txt crest/sclip/%.$1.cover : bam/%.bam
	$$(call LSCRIPT_CHECK_MEM,4G,6G,"$$(EXTRACT_SCLIP) -i $$< --ref_genome $$(REF_FASTA) -r $1 -o $$(@D) -p $$*")
endef
$(foreach chr,$(CHROMOSOMES),$(eval $(call sclip-chr,$(chr))))

crest/sclip/%.cover : $(foreach chr,$(CHROMOSOMES),crest/sclip/%.$(chr).cover)
	$(INIT) cat $^ > $@

crest/sclip/%.sclip.txt : $(foreach chr,$(CHROMOSOMES),crest/sclip/%.$(chr).sclip.txt)
	$(INIT) cat $^ > $@

define sv-tumor-normal-chr
crest/sv/$1_$2.$3.predSV.txt : bam/$1.bam bam/$2.bam crest/sclip/$1.cover 
	$$(call LSCRIPT_CHECK_MEM,4G,6G,"$$(CREST) -p $1_$2.$3 -f $$(<<<) -d $$< -g $$(<<) --ref_genome $$(REF_FASTA) -t $$(REF_2BIT) -r $3")
endef
$(foreach pair,$(SAMPLE_PAIRS),\
	$(foreach chr,$(CHROMOSOMES),\
		$(eval $(call sv-tumor-normal-chr,$(tumor.$(pair)),$(normal.$(pair)),$(chr)))))

define sv-chr
crest/sv/%.$1.predSV.txt : bam/%.bam crest/sclip/%.cover 
	$$(call LSCRIPT_CHECK_MEM,4G,6G,"$$(CREST) -p $$*.$1 -f $$(<<) -d $$< --ref_genome $$(REF_FASTA) -t $$(REF_2BIT) -r $1")
endef
$(foreach chr,$(CHROMOSOMES),$(eval $(call sv-chr,$(chr))))

crest/sv/%.predSV.txt : $(foreach chr,$(CHROMOSOMES),crest/sv/%.$(chr).predSV.txt)
	$(INIT) cat $^ > $@


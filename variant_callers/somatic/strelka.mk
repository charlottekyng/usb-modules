

# Run strelka on tumour-normal matched pairs

include usb-modules/Makefile.inc
include usb-modules/variant_callers/somatic/somaticVariantCaller.inc

##### DEFAULTS ######


LOGDIR ?= log/strelka.$(NOW)
PHONY += strelka strelka_vcfs strelka_tables

strelka : strelka_vcfs strelka_tables

VARIANT_TYPES := strelka_indels
#strelka_snps strelka_indels
strelka_vcfs : $(foreach type,$(VARIANT_TYPES),$(call SOMATIC_VCFS,$(type)))
strelka_tables : $(foreach type,$(VARIANT_TYPES),$(call SOMATIC_TABLES,$(type)))

define strelka-tumor-normal
strelka/$1_$2/Makefile : bam/$1.bam bam/$2.bam
	$$(call LSCRIPT_NAMED,strelka_$1_$2,"rm -rf $$(@D) && $$(LOAD_PERL_MODULE) && $$(CONFIGURE_STRELKA) --tumor=$$< --normal=$$(<<) \
		--ref=$$(REF_FASTA) --config=$$(STRELKA_CONFIG) --output-dir=$$(@D)")

#$$(INIT) qmake -inherit -q jrf.q -- -j 20 -C $$< > $$(LOG) && touch $$@
strelka/$1_$2/task.complete : strelka/$1_$2/Makefile
	$$(call LSCRIPT_NAMED_PARALLEL_MEM,$1_$2.strelka,8,4G,02:29:29,"make -j 8 -C $$(<D)")

vcf/$1_$2.%.vcf : strelka/vcf/$1_$2.%.vcf
	$$(INIT) perl -ne 'if (/^#CHROM/) { s/NORMAL/$2/; s/TUMOR/$1/; } print;' $$< > $$@ && $$(RM) $$<

strelka/vcf/$1_$2.strelka_snps.vcf : strelka/$1_$2/task.complete
	$$(INIT) cp -f strelka/$1_$2/results/all.somatic.snvs.vcf $$@

strelka/vcf/$1_$2.strelka_indels.vcf : strelka/$1_$2/task.complete
	$$(INIT) $$(FIX_STRELKA_VCF) strelka/$1_$2/results/all.somatic.indels.vcf > $$@

endef
$(foreach pair,$(SAMPLE_PAIRS),$(eval $(call strelka-tumor-normal,$(tumor.$(pair)),$(normal.$(pair)))))

include usb-modules/vcf_tools/vcftools.mk

.DELETE_ON_ERROR:
.SECONDARY:
.PHONY: $(PHONY)


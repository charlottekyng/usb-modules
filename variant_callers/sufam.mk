##### DEFAULTS ######
LOGDIR = log/sufamscreen.$(NOW)

##### MAKE INCLUDES #####
include usb-modules/Makefile.inc
include usb-modules/config.inc
#include usb-modules/variant_callers/somatic/somaticVariantCaller.inc

VPATH ?= bam

.DELETE_ON_ERROR:
.SECONDARY: 
.PHONY : all

all : $(foreach sample,$(SAMPLES),sufamscreen/$(sample).sufamscreen.target_ft.dp_ft.altad_ft.pass.eff.vcf)

ifeq ($(findstring ILLUMINA,$(SEQ_PLATFORM)),ILLUMINA)
sufamscreen/%.sufamscreen.vcf : bam/%.bam sufamscreen/sites.to.screen.vcf
	$(call LSCRIPT_PARALLEL_MEM,8,5G,03:59:59,"$(LOAD_JAVA8_MODULE); $(call UNIFIED_GENOTYPER,4G) \
		-nt 8 -R $(REF_FASTA) --dbsnp $(DBSNP) $(foreach bam,$(filter %.bam,$<),-I $(bam) ) --downsampling_type NONE \
		--genotyping_mode GENOTYPE_GIVEN_ALLELES -alleles $(word 2,$^) -o $@ --output_mode EMIT_ALL_SITES")
endif

ifeq ($(findstring IONTORRENT,$(SEQ_PLATFORM)),IONTORRENT)
sufamscreen/%/sufamscreen.vcf : bam/%.bam sufamscreen/sites.to.screen.vcf.gz sufamscreen/sites.to.screen.vcf.gz.tbi
	$(call LSCRIPT_PARALLEL_MEM,8,32G,11:59:59,"$(LOAD_BCFTOOLS_MODULE); $(LOAD_JAVA8_MODULE); $(LOAD_TABIX_MODULE); \
	$(TVC) -s $(word 2,$^) -i $< -r $(REF_FASTA) -o $(@D) -N 8 \
	$(if $(TARGETS_FILE_INTERVALS),-b $(TARGETS_FILE_INTERVALS)) -p $(TVC_SENSITIVE_JSON) -m $(TVC_MOTIF) \
	-t $(TVC_ROOT_DIR) --primer-trim-bed $(PRIMER_TRIM_BED) && \
	$(BCFTOOLS) norm -m -both $(@D)/TSVC_variants.vcf.gz | grep -v \"##contig\" > $(@D)/TSVC_variants.vcf.tmp && \
	$(call LEFT_ALIGN_VCF,6G) -R $(REF_FASTA) --variant $(@D)/TSVC_variants.vcf.tmp -o $(@D)/TSVC_variants.vcf.tmp2 && \
	$(FIX_TVC_VCF) $(@D)/TSVC_variants.vcf.tmp2 > $(@D)/TSVC_variants.vcf.tmp3 && \
	$(BGZIP) -c $(@D)/TSVC_variants.vcf.tmp3 > $(@D)/TSVC_variants.vcf.tmp3.gz && \
	$(TABIX) -p vcf $(@D)/TSVC_variants.vcf.tmp3.gz && \
	$(BCFTOOLS) isec -O v -p $(@D)/isec $(@D)/TSVC_variants.vcf.tmp3.gz $(word 2,$^) && \
	mv $(@D)/isec/0002.vcf $@ && $(RMR) $(@D)/isec && $(RM) $(@D)/*tmp*")

sufamscreen/%.sufamscreen.vcf : sufamscreen/%/sufamscreen.vcf
	ln -f $< $@
endif

ifndef $(TARGETS_FILE_INTERVALS)
sufamscreen/sites.to.screen.vcf : $(TARGETS_FILE_INTERVALS) $(SUFAM_SCREEN_VCF)
	$(INIT) grep "^\#" $(CANCER_HOTSPOT_VCF) > $@ && \
	$(LOAD_BEDTOOLS_MODULE); $(BEDTOOLS) intersect -b $(TARGETS_FILE_INTERVALS) -a $(SUFAM_SCREEN_VCF) -header >>$@
else
sufamscreen/sites.to.screen.vcf : $(SUFAM_SCREEN_VCF)
	$(INIT) ln $(SUFAM_SCREEN_VCF) $@
endif

include usb-modules/vcf_tools/vcftools.mk

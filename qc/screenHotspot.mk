
##### DEFAULTS ######
LOGDIR = log/screen_hotspot.$(NOW)

##### MAKE INCLUDES #####
include usb-modules/Makefile.inc
include usb-modules/config.inc
#include usb-modules/variant_callers/somatic/somaticVariantCaller.inc

VPATH ?= bam

.DELETE_ON_ERROR:
.SECONDARY: 
.PHONY : all

all : $(foreach sample,$(SAMPLES),hotspots/$(sample).hotspotscreen.target_ft.dp_ft.altad_ft.pass.eff.vcf) $(foreach sample,$(SAMPLES),hotspots/$(sample).hotspotscreen.target_ft.dp_ft.altad_ft.pass.eff.tab.txt)

ifeq ($(findstring ILLUMINA,$(SEQ_PLATFORM)),ILLUMINA)
hotspots/%.hotspotscreen.vcf : bam/%.bam hotspots/sites.to.screen.vcf
	$(call LSCRIPT_PARALLEL_MEM,8,$(RESOURCE_REQ_LOWMEM),$(RESOURCE_REQ_SHORT),"$(LOAD_JAVA8_MODULE); \
		$(call GATK,UnifiedGenotyper,$(RESOURCE_REQ_LOWMEM)) \
		-nt 8 -R $(REF_FASTA) --dbsnp $(DBSNP) $(foreach bam,$(filter %.bam,$<),-I $(bam) ) --downsampling_type NONE \
		--genotyping_mode GENOTYPE_GIVEN_ALLELES -alleles $(word 2,$^) -o $@ --output_mode EMIT_ALL_SITES")
endif

ifeq ($(findstring IONTORRENT,$(SEQ_PLATFORM)),IONTORRENT)
hotspots/%/hotspotscreen.vcf : bam/%.bam hotspots/sites.to.screen.vcf.gz hotspots/sites.to.screen.vcf.gz.tbi
	$(call LSCRIPT_PARALLEL_MEM,8,$(RESOURCE_REQ_HIGHMEM),$(RESOURCE_REQ_LONG),"$(LOAD_BCFTOOLS_MODULE); $(LOAD_JAVA8_MODULE); $(LOAD_TABIX_MODULE); \
	$(TVC) -s $(word 2,$^) -i $< -r $(REF_FASTA) -o $(@D) -N 8 \
	$(if $(TARGETS_FILE_INTERVALS),-b $(TARGETS_FILE_INTERVALS)) -p $(TVC_SENSITIVE_JSON) -m $(TVC_MOTIF) \
	-t $(TVC_ROOT_DIR) --primer-trim-bed $(PRIMER_TRIM_BED) && \
	$(BCFTOOLS) norm -m -both $(@D)/TSVC_variants.vcf.gz | grep -v \"##contig\" > $(@D)/TSVC_variants.vcf.tmp && \
	$(call GATK,LeftAlignAndTrimVariants,$(RESOURCE_REQ_MEDIUM_MEM)) -R $(REF_FASTA) \
	--variant $(@D)/TSVC_variants.vcf.tmp -o $(@D)/TSVC_variants.vcf.tmp2 && \
	$(FIX_TVC_VCF) $(@D)/TSVC_variants.vcf.tmp2 > $(@D)/TSVC_variants.vcf.tmp3 && \
	$(BGZIP) -c $(@D)/TSVC_variants.vcf.tmp3 > $(@D)/TSVC_variants.vcf.tmp3.gz && \
	$(TABIX) -p vcf $(@D)/TSVC_variants.vcf.tmp3.gz && \
	$(BCFTOOLS) isec -O v -p $(@D)/isec $(@D)/TSVC_variants.vcf.tmp3.gz $(word 2,$^) && \
	mv $(@D)/isec/0002.vcf $@ && $(RMR) $(@D)/isec && $(RM) $(@D)/*tmp*")

hotspots/%.hotspotscreen.vcf : hotspots/%/hotspotscreen.vcf
	ln -f $< $@
endif

ifndef $(TARGETS_FILE_INTERVALS)
hotspots/sites.to.screen.vcf : $(TARGETS_FILE_INTERVALS) $(CANCER_HOTSPOT_VCF)
	$(INIT) grep "^\#" $(CANCER_HOTSPOT_VCF) > $@ && \
	$(LOAD_BEDTOOLS_MODULE); $(BEDTOOLS) intersect -b $(TARGETS_FILE_INTERVALS) -a $(CANCER_HOTSPOT_VCF) -header >>$@
else
hotspots/sites.to.screen.vcf : $(CANCER_HOTSPOT_VCF)
	$(INIT) ln $(CANCER_HOTSPOT_VCF) $@
endif

include usb-modules/vcf_tools/vcftools.mk

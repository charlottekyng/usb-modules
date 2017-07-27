
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

all : $(shell rm -rf hotspots/all.hotspotscreen.target_ft.dp_ft.altad_ft.pass.eff.tab.txt) hotspots/all.hotspotscreen.target_ft.dp_ft.altad_ft.pass.eff.tab.txt 

#$(foreach sample,$(SAMPLES),hotspots/$(sample).hotspotscreen.target_ft.dp_ft.altad_ft.hotspot.pass.eff.vcf) $(foreach sample,$(SAMPLES),hotspots/$(sample).hotspotscreen.target_ft.dp_ft.altad_ft.hotspot.pass.eff.tab.txt)

hotspots/all.hotspotscreen.target_ft.dp_ft.altad_ft.hotspot.eff.tab.txt : $(foreach sample,$(SAMPLES),hotspots/$(sample).hotspotscreen.target_ft.dp_ft.altad_ft.pass.eff.tab.txt)
	$(call LSCRIPT_MEM,$(RESOURCE_REQ_MEDIUM_MEM),$(RESOURCE_REQ_VSHORT),"$(LOAD_R_MODULE); $(RSCRIPT) $(RBIND) --sampleName $< $^ > $@")

ifeq ($(findstring ILLUMINA,$(SEQ_PLATFORM)),ILLUMINA)
hotspots/%.hotspotscreen.vcf : bam/%.bam hotspots/sites.to.screen.vcf
	$(call LSCRIPT_PARALLEL_MEM,8,$(RESOURCE_REQ_LOWMEM),$(RESOURCE_REQ_SHORT),"$(LOAD_JAVA8_MODULE); \
		$(call GATK,UnifiedGenotyper,$(RESOURCE_REQ_LOWMEM)) \
		-nt 8 -R $(REF_FASTA) --dbsnp $(DBSNP) $(foreach bam,$(filter %.bam,$<),-I $(bam) ) --downsampling_type NONE \
		--genotyping_mode GENOTYPE_GIVEN_ALLELES -alleles $(word 2,$^) -o $@ --output_mode EMIT_ALL_SITES")
endif

ifeq ($(findstring IONTORRENT,$(SEQ_PLATFORM)),IONTORRENT)
MUT_CALLER = tvc
hotspots/%/TSVC_variants.vcf.gz : bam/%.bam hotspots/sites.to.screen.vcf
	$(call LSCRIPT_PARALLEL_MEM,4,$(RESOURCE_REQ_MEDIUM_MEM),$(RESOURCE_REQ_SHORT),"$(LOAD_BCFTOOLS_MODULE); $(LOAD_JAVA8_MODULE); $(LOAD_TABIX_MODULE); \
	$(TVC) -s $(word 2,$^) -i $< -r $(REF_FASTA) -o $(@D) -N 4 \
	$(if $(TARGETS_FILE_INTERVALS),-b $(TARGETS_FILE_INTERVALS)) -p $(TVC_SENSITIVE_JSON) -m $(TVC_MOTIF) \
	-t $(TVC_ROOT_DIR) --primer-trim-bed $(PRIMER_TRIM_BED)")

hotspots/%/hotspotscreen.vcf : hotspots/%/TSVC_variants.norm.left_align.vcf.gz hotspots/sites.to.screen.vcf.gz hotspots/%/TSVC_variants.norm.left_align.vcf.gz.tbi hotspots/sites.to.screen.vcf.gz.tbi
	$(call LSCRIPT_MEM,$(RESOURCE_REQ_LOWMEM),$(RESOURCE_REQ_VSHORT),"$(LOAD_BCFTOOLS_MODULE); \
	$(BCFTOOLS) isec -O v -p $(dir $@)/isec $(word 1,$^) $(word 2,$^) && mv $(dir $@)/isec/0002.vcf $@ \
	&& $(RMR) $(@D)/isec && $(RM) $(@D)/*tmp*")

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

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

all : $(shell rm -rf sufamscreen/all.sufamscreen.target_ft.dp_ft.altad_ft.pass.eff.tab.txt) sufamscreen/all.sufamscreen.target_ft.dp_ft.altad_ft.pass.eff.tab.txt

#all : $(foreach sample,$(SAMPLES),sufamscreen/$(sample).sufamscreen.target_ft.dp_ft.altad_ft.pass.eff.vcf)

sufamscreen/all.sufamscreen.target_ft.dp_ft.altad_ft.pass.eff.tab.txt : $(foreach sample,$(SAMPLES),sufamscreen/$(sample).sufamscreen.target_ft.dp_ft.altad_ft.pass.eff.tab.txt)
	$(call LSCRIPT_MEM,$(RESOURCE_REQ_MEDIUM_MEM),$(RESOURCE_REQ_VSHORT),"$(LOAD_R_MODULE); $(RSCRIPT) $(RBIND) --sampleName $< $^ > $@")

ifeq ($(findstring ILLUMINA,$(SEQ_PLATFORM)),ILLUMINA)
sufamscreen/%.sufamscreen.vcf : bam/%.bam sufamscreen/sites.to.screen.vcf
	$(call LSCRIPT_PARALLEL_MEM,8,$(RESOURCE_REQ_MEDIUM_MEM),$(RESOURCE_REQ_SHORT),"$(LOAD_JAVA8_MODULE); \
		$(call GATK,UnifiedGenotyper,$(RESOURCE_REQ_MEDIUM_MEM)) \
		-nt 8 -R $(REF_FASTA) --dbsnp $(DBSNP) $(foreach bam,$(filter %.bam,$<),-I $(bam) ) --downsampling_type NONE \
		--genotyping_mode GENOTYPE_GIVEN_ALLELES -alleles $(word 2,$^) -o $@ --output_mode EMIT_ALL_SITES")
endif

ifeq ($(findstring IONTORRENT,$(SEQ_PLATFORM)),IONTORRENT)
MUT_CALLER = tvc
sufamscreen/%/TSVC_variants.vcf.gz : bam/%.bam sufamscreen/sites.to.screen.vcf
	$(call LSCRIPT_PARALLEL_MEM,4,$(RESOURCE_REQ_MEDIUM_MEM),$(RESOURCE_REQ_SHORT),"$(LOAD_BCFTOOLS_MODULE); $(LOAD_JAVA8_MODULE); $(LOAD_TABIX_MODULE); \
	$(TVC) -s $(word 2,$^) -i $< -r $(REF_FASTA) -o $(@D) -N 4 \
	$(if $(TARGETS_FILE_INTERVALS),-b $(TARGETS_FILE_INTERVALS)) -p $(TVC_SENSITIVE_JSON) -m $(TVC_MOTIF) \
	-t $(TVC_ROOT_DIR) --primer-trim-bed $(PRIMER_TRIM_BED)")

sufamscreen/%/sufamscreen.vcf : sufamscreen/%/TSVC_variants.norm.left_align.vcf.gz sufamscreen/sites.to.screen.vcf.gz
	$(call LSCRIPT_MEM,$(RESOURCE_REQ_LOWMEM),$(RESOURCE_REQ_VSHORT),"$(LOAD_BCFTOOLS_MODULE); \
	$(BCFTOOLS) isec -O v -p $(dir $@)/isec $(word 1,$^) $(word 2,$^) && mv $(dir $@)/isec/0002.vcf $@ \
	&& $(RMR) $(@D)/isec && $(RM) $(@D)/*tmp*")

sufamscreen/%.sufamscreen.vcf : sufamscreen/%/sufamscreen.vcf
	perl -p -e "s/NOCALL/\./g" < $< > $@
endif

ifndef $(TARGETS_FILE_INTERVALS)
sufamscreen/sites.to.screen.vcf : $(TARGETS_FILE_INTERVALS) $(SUFAM_SCREEN_VCF)
	$(INIT) grep "^\#" $(SUFAM_SCREEN_VCF) > $@ && \
	$(LOAD_BEDTOOLS_MODULE); $(BEDTOOLS) intersect -b $(TARGETS_FILE_INTERVALS) -a $(SUFAM_SCREEN_VCF) -header >>$@
else
sufamscreen/sites.to.screen.vcf : $(SUFAM_SCREEN_VCF)
	$(INIT) ln $(SUFAM_SCREEN_VCF) $@
endif

include usb-modules/vcf_tools/vcftools.mk

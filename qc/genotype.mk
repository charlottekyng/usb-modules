
##### DEFAULTS ######
LOGDIR = log/genotype.$(NOW)

##### MAKE INCLUDES #####
include usb-modules/Makefile.inc
include usb-modules/config.inc

VPATH ?= bam

.DELETE_ON_ERROR:
.SECONDARY: 
.PHONY : all

CLUSTER_VCF = $(RSCRIPT) usb-modules/qc/clusterVcf.R

all : genotype/snps_filtered.sdp_ft.clust.png

ifeq ($(findstring ILLUMINA,$(SEQ_PLATFORM)),ILLUMINA)
define snps-chr
genotype/snps.$1.vcf : $$(foreach sample,$$(SAMPLES),gatk/dbsnp/$$(sample).gatk_snps.vcf)
	$$(call LSCRIPT_MEM,$$(RESOURCE_REQ_VVHIGHMEM),$$(RESOURCE_REQ_SHORT),"$$(LOAD_JAVA8_MODULE); \
	$$(call GATK,CombineVariants,$$(RESOURCE_REQ_VVHIGHMEM)) \
	$$(foreach vcf,$$^,--variant $$(vcf) ) -o $$@ --genotypemergeoption UNSORTED \
	-L $1 -R $$(REF_FASTA)")
endef
$(foreach chr,$(CHROMOSOMES),$(eval $(call snps-chr,$(chr))))

ifeq ($(findstring true,$(GENOTYPE_WITH_CHR2)),true)
genotype/snps.vcf : genotype/snps.2.vcf
else
genotype/snps.vcf : $(foreach chr,$(CHROMOSOMES),genotype/snps.$(chr).vcf)
	$(call LSCRIPT_MEM,$(RESOURCE_REQ_VVHIGHMEM),$(RESOURCE_REQ_SHORT),"$(LOAD_JAVA8_MODULE); \
		$(call GATK,CombineVariants,$(RESOURCE_REQ_VVHIGHMEM)) \
		$(foreach vcf,$^,--variant $(vcf) ) -o $@ --genotypemergeoption UNSORTED -R $(REF_FASTA)")
endif
endif #end illumina

ifeq ($(findstring IONTORRENT,$(SEQ_PLATFORM)),IONTORRENT)
genotype/snps.vcf : $(foreach sample,$(SAMPLES),tvc/dbsnp/$(sample)/TSVC_variants.vcf)
	$(call LSCRIPT_MEM,$(RESOURCE_REQ_HIGHMEM),$(RESOURCE_REQ_MEDIUM),"$(LOAD_JAVA8_MODULE); \
	$(call GATK,CombineVariants,$(RESOURCE_REQ_HIGHMEM)) \
	$(foreach vcf,$^,--variant $(vcf) ) -o $@ --genotypemergeoption UNSORTED -R $(REF_FASTA)")
endif

genotype/snps_filtered.vcf : genotype/snps.vcf
	$(INIT) grep '^#' $< > $@ && grep -e '0/1' -e '1/1' $< >> $@

genotype/%.clust.png : genotype/%.vcf
	$(call LSCRIPT_MEM,22G,03:59:59,"$(LOAD_R_MODULE); $(CLUSTER_VCF) --outPrefix genotype/$* \
		$(if $(findstring RNA,$(CAPTURE_METHOD)),--clusterType hetSameHom) $<")

include usb-modules/vcf_tools/vcftools.mk
include usb-modules/variant_callers/TVC.mk
include usb-modules/variant_callers/gatk.mk

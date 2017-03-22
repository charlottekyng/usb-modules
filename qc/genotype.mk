
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
#genotype/snps.vcf : $(foreach sample,$(SAMPLES),gatk/dbsnp/$(sample).gatk_snps.vcf)
#	$(call LSCRIPT_MEM,95G,03:59:59,"$(LOAD_JAVA8_MODULE); $(call COMBINE_VARIANTS,96G) \
#		$(foreach vcf,$^,--variant $(vcf) ) -o $@ --genotypemergeoption UNSORTED \
#		$(if $(findstring AGILENT_CLINICAL_EXOME,$(PANEL)),-L 2,) -R $(REF_FASTA)")
genotype/snps_A.vcf : $(foreach sample,$(filter A%,$(SAMPLES)),gatk/dbsnp/$(sample).gatk_snps.vcf)
	$(call LSCRIPT_MEM,95G,03:59:59,"$(LOAD_JAVA8_MODULE); $(call COMBINE_VARIANTS,96G) \
		$(foreach vcf,$^,--variant $(vcf) ) -o $@ --genotypemergeoption UNSORTED \
		$(if $(findstring AGILENT_CLINICAL_EXOME,$(PANEL)),-L 2,) -R $(REF_FASTA)")

genotype/snps_B.vcf : $(foreach sample,$(filter B%,$(SAMPLES)),gatk/dbsnp/$(sample).gatk_snps.vcf)
	$(call LSCRIPT_MEM,95G,03:59:59,"$(LOAD_JAVA8_MODULE); $(call COMBINE_VARIANTS,96G) \
		$(foreach vcf,$^,--variant $(vcf) ) -o $@ --genotypemergeoption UNSORTED \
		$(if $(findstring AGILENT_CLINICAL_EXOME,$(PANEL)),-L 2,) -R $(REF_FASTA)")

genotype/snps_C.vcf : $(foreach sample,$(filter C%,$(SAMPLES)),gatk/dbsnp/$(sample).gatk_snps.vcf)
	$(call LSCRIPT_MEM,95G,03:59:59,"$(LOAD_JAVA8_MODULE); $(call COMBINE_VARIANTS,96G) \
		$(foreach vcf,$^,--variant $(vcf) ) -o $@ --genotypemergeoption UNSORTED \
		$(if $(findstring AGILENT_CLINICAL_EXOME,$(PANEL)),-L 2,) -R $(REF_FASTA)")

genotype/snps_D.vcf : $(foreach sample,$(filter D%,$(SAMPLES)),gatk/dbsnp/$(sample).gatk_snps.vcf)
	$(call LSCRIPT_MEM,95G,03:59:59,"$(LOAD_JAVA8_MODULE); $(call COMBINE_VARIANTS,96G) \
		$(foreach vcf,$^,--variant $(vcf) ) -o $@ --genotypemergeoption UNSORTED \
		$(if $(findstring AGILENT_CLINICAL_EXOME,$(PANEL)),-L 2,) -R $(REF_FASTA)")

genotype/snps_S.vcf : $(foreach sample,$(filter S%,$(SAMPLES)),gatk/dbsnp/$(sample).gatk_snps.vcf)
	$(call LSCRIPT_MEM,95G,03:59:59,"$(LOAD_JAVA8_MODULE); $(call COMBINE_VARIANTS,96G) \
		$(foreach vcf,$^,--variant $(vcf) ) -o $@ --genotypemergeoption UNSORTED \
		$(if $(findstring AGILENT_CLINICAL_EXOME,$(PANEL)),-L 2,) -R $(REF_FASTA)")

genotype/snps.vcf : genotype/snps_A.vcf genotype/snps_B.vcf genotype/snps_C.vcf genotype/snps_D.vcf genotype/snps_S.vcf
	$(call LSCRIPT_MEM,95G,03:59:59,"$(LOAD_JAVA8_MODULE); $(call COMBINE_VARIANTS,96G) \
		$(foreach vcf,$^,--variant $(vcf) ) -o $@ --genotypemergeoption UNSORTED \
		$(if $(findstring AGILENT_CLINICAL_EXOME,$(PANEL)),-L 2,) -R $(REF_FASTA)")


endif
ifeq ($(findstring IONTORRENT,$(SEQ_PLATFORM)),IONTORRENT)
genotype/snps.vcf : $(foreach sample,$(SAMPLES),tvc/dbsnp/$(sample)/TSVC_variants.vcf)
	$(call LSCRIPT_MEM,22G,03:59:59,"$(LOAD_JAVA8_MODULE); $(call COMBINE_VARIANTS,21G) \
		$(foreach vcf,$^,--variant $(vcf) ) -o $@ --genotypemergeoption UNSORTED -R $(REF_FASTA)")
endif

#genotype/snps.vcf : $(foreach sample,$(SAMPLES),genotype/$(sample).snps.vcf)
#	$(call LSCRIPT_MEM,22G,03:59:59,"$(LOAD_JAVA8_MODULE); $(call COMBINE_VARIANTS,21G) \
#		$(foreach vcf,$^,--variant $(vcf) ) -o $@ --genotypemergeoption UNSORTED -R $(REF_FASTA)")

genotype/snps_filtered.vcf : genotype/snps.vcf
	$(INIT) grep '^#' $< > $@ && grep -e '0/1' -e '1/1' $< >> $@

genotype/%.clust.png : genotype/%.vcf
	$(call LSCRIPT_MEM,22G,03:59:59,"$(LOAD_R_MODULE); $(CLUSTER_VCF) --outPrefix genotype/$* $<")

include usb-modules/vcf_tools/vcftools.mk
include usb-modules/variant_callers/TVC.mk
include usb-modules/variant_callers/gatk.mk

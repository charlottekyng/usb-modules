SEQ_PLATFORM = ILLUMINA
REF = hg19
PANEL = AGILENT_CLINICAL_EXOME
CAPTURE_METHOD = BAITS
INCLUDE_CHR_Y = true

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

genotype/snps.vcf : $(foreach sample,$(SAMPLES),genotype/$(sample).snps.vcf)
	$(call LSCRIPT_MEM,16G,20G,"$(LOAD_JAVA8_MODULE); $(COMBINE_VARIANTS) \
		$(foreach vcf,$^,--variant $(vcf) ) -o $@ --genotypemergeoption UNSORTED -R $(REF_FASTA)")

genotype/snps_filtered.vcf : genotype/snps.vcf
	$(INIT) grep '^#' $< > $@ && grep -e '0/1' -e '1/1' $< >> $@

genotype/%.snps.vcf : bam/%.bam genotype/sites.to.genotype.vcf
	$(call LSCRIPT_PARALLEL_MEM,4,2.5G,00:59:59,"$(LOAD_JAVA8_MODULE); $(UNIFIED_GENOTYPER) \
		-nt 4 -R $(REF_FASTA) --dbsnp $(DBSNP) $(foreach bam,$(filter %.bam,$<),-I $(bam) ) \
		--genotyping_mode GENOTYPE_GIVEN_ALLELES -alleles $(word 2,$^) -o $@ --output_mode EMIT_ALL_SITES")

ifndef $(TARGETS_FILE_INTERVALS)
genotype/sites.to.genotype.vcf : $(TARGETS_FILE_INTERVALS) $(DBSNP)
	$(INIT) $(LOAD_BEDTOOLS_MODULE); $(BEDTOOLS) intersect -b $(TARGETS_FILE_INTERVALS) -a $(DBSNP) -header | \
	awk '$$5!~/N/'> $@
else
genotype/sites.to.genotype.vcf : $(DBSNP)
	$(INIT) awk '$$5!~/N/' $(DBSNP) > genotype/sites.to.genotype.vcf
endif

genotype/%.clust.png : genotype/%.vcf
	$(INIT) $(LOAD_R_MODULE); $(CLUSTER_VCF) --outPrefix genotype/$* $<

include usb-modules/vcf_tools/vcftools.mk

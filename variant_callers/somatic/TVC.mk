include usb-modules/Makefile.inc
include usb-modules/config.inc
include usb-modules/variant_callers/somatic/somaticVariantCaller.inc

LOGDIR ?= log/tvc_somatic.$(NOW)

PHONY += tvc_somatic tvc_somatic_vcfs tvc_somatic_tables 
#tvc_somatic_vcfs_sets tvc_somatic_tables_sets
tvc_somatic : tvc_somatic_vcfs tvc_somatic_tables 
#tvc_somatic_vcfs_sets tvc_somatic_tables_sets

VARIANT_TYPES ?= tvc_snps tvc_indels
tvc_somatic_vcfs : $(foreach type,$(VARIANT_TYPES),$(call SOMATIC_VCFS,$(type)) $(addsuffix .idx,$(call SOMATIC_VCFS,$(type))))
tvc_somatic_tables : $(foreach type,$(VARIANT_TYPES),$(call SOMATIC_TABLES,$(type)))

VARIANT_TYPES_SETS ?= tvc_snps_sufam tvc_indels_sufam
tvc_somatic_vcfs_sets : $(foreach type,$(VARIANT_TYPES_SETS),$(call SOMATIC_VCFS_SETS,$(type)) $(addsuffix .idx,$(call SOMATIC_VCFS_SETS,$(type))))
tvc_somatic_tables_sets : $(foreach type,$(VARIANT_TYPES_SETS),$(call SOMATIC_TABLES_SETS,$(type)))

.DELETE_ON_ERROR:
.SECONDARY:
.PHONY : $(PHONY)

MUT_CALLER = tvc

define tvc-somatic-vcf
tvc/vcf/$1_$2/TSVC_variants_preliminary.vcf : bam/$1.bam bam/$1.bam.bai bam/$2.bam bam/$2.bam.bai
	$$(call LSCRIPT_PARALLEL_MEM,8,$$(RESOURCE_REQ_HIGHMEM),$$(RESOURCE_REQ_LONG),"$$(LOAD_BCFTOOLS_MODULE); $$(LOAD_JAVA8_MODULE); \
	$$(TVC) -i $$< -n $$(word 3,$$^) -r $$(REF_FASTA) -o $$(@D) -N 8 \
	$$(if $$(TARGETS_FILE_INTERVALS),-b $$(TARGETS_FILE_INTERVALS)) -p $$(TVC_SOMATIC_JSON) -m $$(TVC_MOTIF) \
	-t $$(TVC_ROOT_DIR) --primer-trim-bed $$(PRIMER_TRIM_BED) -g $$(basename $$(notdir $$<)) && \
	$$(BCFTOOLS) norm -m -both $$(@D)/TSVC_variants.vcf.gz | grep -v \"##contig\" > $$(@D)/TSVC_variants.vcf.tmp && \
	$$(call GATK,LeftAlignAndTrimVariants,$$(RESOURCE_REQ_MEDIUM_MEM)) -R $$(REF_FASTA) \
	--variant $$(@D)/TSVC_variants.vcf.tmp -o $$(@D)/TSVC_variants.vcf.tmp2 && \
	$$(FIX_TVC_VCF) $$(@D)/TSVC_variants.vcf.tmp2 > $$@ && $$(RM) $$@.gz $$(@D)/TSVC_variants.vcf.tmp $$(@D)/TSVC_variants.vcf.tmp2")

tvc/vcf/$1_$2/TSVC_variants_preliminary.fpft.vcf : tvc/vcf/$1_$2/TSVC_variants_preliminary.vcf bam/$1.bam bam/$1.bam.bai
	$$(call LSCRIPT_MEM,$$(RESOURCE_REQ_LOWMEM),$$(RESOURCE_REQ_VSHORT),"$$(LOAD_JAVA8_MODULE); $$(LOAD_TABIX_MODULE); awk '! /\#/' $$< | \
	awk '{if(length(\$$$$4) > length(\$$$$5)) print \$$$$1\"\t\"(\$$$$2-1)\"\t\"(\$$$$2+length(\$$$$4)-1); \
	else print \$$$$1\"\t\"(\$$$$2-1)\"\t\"(\$$$$2+length(\$$$$5)-1)}' > $$<.region && \
	$$(BAM_READCOUNT) -f $$(REF_FASTA) -l $$<.region $$(word 2,$$^) > $$<.bamrc && \
	$$(VARSCAN) fpfilter $$< $$<.bamrc --output-file $$@ --filtered-file $$<.fail \
	--min-var-freq $$(MIN_AF_SNP) --min-ref-readpos 0 --min-var-readpos 0 --min-ref-dist3 0 --min-var-dist3 0 && \
	rm $$<.region $$<.bamrc")

tvc/vcf/$1_$2/tumor/TSVC_variants_final.vcf : bam/$1.bam bam/$1.bam.bai tvc/vcf/$1_$2/TSVC_variants_preliminary.$$(if $$(findstring true,$$(USE_FPFILTER_FOR_TVC)),fpft.)vcf tvc/vcf/$1_$2/TSVC_variants_preliminary.$$(if $$(findstring true,$$(USE_FPFILTER_FOR_TVC)),fpft.)vcf.gz tvc/vcf/$1_$2/TSVC_variants_preliminary.$$(if $$(findstring true,$$(USE_FPFILTER_FOR_TVC)),fpft.)vcf.gz.tbi
	$$(call LSCRIPT_PARALLEL_MEM,8,$$(RESOURCE_REQ_HIGHMEM),$$(RESOURCE_REQ_LONG),"$$(LOAD_BCFTOOLS_MODULE); $$(LOAD_JAVA8_MODULE); $$(LOAD_TABIX_MODULE) && \
		$$(TVC) -s $$(<<<) -i $$(<) -r $$(REF_FASTA) -o $$(@D) -N 8 -p $$(TVC_SENSITIVE_JSON) \
		-m $$(TVC_MOTIF) -t $$(TVC_ROOT_DIR) --primer-trim-bed $$(PRIMER_TRIM_BED) && \
		$$(BCFTOOLS) norm -m -both $$(@D)/TSVC_variants.vcf.gz | grep -v \"##contig\" > $$(@D)/TSVC_variants.vcf.tmp && \
		$$(call GATK,LeftAlignAndTrimVariants,$$(RESOURCE_REQ_MEDIUM_MEM)) -R $$(REF_FASTA) --variant $$(@D)/TSVC_variants.vcf.tmp -o $$(@D)/TSVC_variants.vcf.tmp2 && \
		$$(FIX_TVC_VCF) $$(@D)/TSVC_variants.vcf.tmp2 > $$(@D)/TSVC_variants.vcf.tmp3 && \
		$$(BGZIP) -c $$(@D)/TSVC_variants.vcf.tmp3 > $$(@D)/TSVC_variants.vcf.tmp3.gz && \
		$$(TABIX) -p vcf $$(@D)/TSVC_variants.vcf.tmp3.gz && \
		$$(BCFTOOLS) isec -O v -p $$(@D)/isec $$(@D)/TSVC_variants.vcf.tmp3.gz $$(<<<<) && \
		mv $$(@D)/isec/0002.vcf $$@ && $$(RMR) $$(@D)/isec && $$(RM) $$(@D)/*tmp*")

tvc/vcf/$1_$2/normal/TSVC_variants_final.vcf : bam/$2.bam bam/$2.bam.bai tvc/vcf/$1_$2/TSVC_variants_preliminary.$$(if $$(findstring true,$$(USE_FPFILTER_FOR_TVC)),fpft.)vcf tvc/vcf/$1_$2/TSVC_variants_preliminary.$$(if $$(findstring true,$$(USE_FPFILTER_FOR_TVC)),fpft.)vcf.gz tvc/vcf/$1_$2/TSVC_variants_preliminary.$$(if $$(findstring true,$$(USE_FPFILTER_FOR_TVC)),fpft.)vcf.gz.tbi
	$$(call LSCRIPT_PARALLEL_MEM,8,$$(RESOURCE_REQ_HIGHMEM),$$(RESOURCE_REQ_LONG),"$$(LOAD_BCFTOOLS_MODULE); $$(LOAD_JAVA8_MODULE) ; $$(LOAD_TABIX_MODULE) ; \
		$$(TVC) -s $$(<<<) -i $$(<) -r $$(REF_FASTA) -o $$(@D) -N 8 $$(TVC_SENSITIVE_JSON) \
		-m $$(TVC_MOTIF) -t $$(TVC_ROOT_DIR) --primer-trim-bed $$(PRIMER_TRIM_BED) && \
		$$(BCFTOOLS) norm -m -both $$(@D)/TSVC_variants.vcf.gz | grep -v \"##contig\" > $$(@D)/TSVC_variants.vcf.tmp && \
		$$(call GATK,LeftAlignAndTrimVariants,$$(RESOURCE_REQ_MEDIUM_MEM)) -R $$(REF_FASTA) --variant $$(@D)/TSVC_variants.vcf.tmp -o $$(@D)/TSVC_variants.vcf.tmp2 && \
		$$(FIX_TVC_VCF) $$(@D)/TSVC_variants.vcf.tmp2 > $$(@D)/TSVC_variants.vcf.tmp3 && \
		$$(BGZIP) -c $$(@D)/TSVC_variants.vcf.tmp3 > $$(@D)/TSVC_variants.vcf.tmp3.gz && \
		$$(TABIX) -p vcf $$(@D)/TSVC_variants.vcf.tmp3.gz && \
		$$(BCFTOOLS) isec -O v -p $$(@D)/isec $$(@D)/TSVC_variants.vcf.tmp3.gz $$(<<<<) && \
		mv $$(@D)/isec/0002.vcf $$@ && $$(RMR) $$(@D)/isec && $$(RM) $$(@D)/*tmp*")

tvc/vcf/$1_$2/TSVC_variants_final.vcf : tvc/vcf/$1_$2/tumor/TSVC_variants_final.vcf.gz tvc/vcf/$1_$2/normal/TSVC_variants_final.vcf.gz tvc/vcf/$1_$2/tumor/TSVC_variants_final.vcf.gz.tbi tvc/vcf/$1_$2/normal/TSVC_variants_final.vcf.gz.tbi
	$$(call LSCRIPT_MEM,$$(RESOURCE_REQ_LOWMEM),$$(RESOURCE_REQ_VSHORT),"$$(LOAD_VCFTOOLS_MODULE); $$(LOAD_TABIX_MODULE); \
		$$(VCFTOOLS_MERGE) -c none $$< $$(word 2,$$^) | perl -p -e \"s/NOCALL/PASS/g;\" > $$@")
endef
$(foreach pair,$(SAMPLE_PAIRS), \
	$(eval $(call tvc-somatic-vcf,$(tumor.$(pair)),$(normal.$(pair)))))

include usb-modules/vcf_tools/vcftools.mk
include usb-modules/variant_callers/TVC.mk
include usb-modules/variant_callers/somatic/pon.mk

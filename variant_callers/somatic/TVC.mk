include usb-modules/Makefile.inc
include usb-modules/config.inc
include usb-modules/variant_callers/somatic/somaticVariantCaller.inc

LOGDIR ?= log/tvc_somtic.$(NOW)

PHONY += tvc_somatic tvc_somatic_vcfs tvc_somatic_tables tvc_somatic_vcfs_sets tvc_somatic_tables_sets
tvc_somatic : tvc_somatic_vcfs tvc_somatic_tables tvc_somatic_vcfs_sets tvc_somatic_tables_sets

VARIANT_TYPES ?= tvc_snps tvc_indels
tvc_somatic_vcfs : $(foreach type,$(VARIANT_TYPES),$(call SOMATIC_VCFS,$(type)) $(addsuffix .idx,$(call SOMATIC_VCFS,$(type))))
tvc_somatic_tables : $(foreach type,$(VARIANT_TYPES),$(call SOMATIC_TABLES,$(type)))

VARIANT_TYPES_SETS ?= tvc_snps_sufam tvc_indels_sufam
tvc_somatic_vcfs_sets : $(foreach type,$(VARIANT_TYPES_SETS),$(call SOMATIC_VCFS_SETS,$(type)) $(addsuffix .idx,$(call SOMATIC_VCFS_SETS,$(type))))
tvc_somatic_tables_sets : $(foreach type,$(VARIANT_TYPES_SETS),$(call SOMATIC_TABLES_SETS,$(type)))

.DELETE_ON_ERROR:
.SECONDARY:
.PHONY : $(PHONY)

define tvc-somatic-vcf
tvc/vcf/$1_$2/TSVC_variants_preliminary.vcf : bam/$1.bam bam/$1.bam.bai bam/$2.bam bam/$2.bam.bai
	$$(call LSCRIPT_PARALLEL_MEM,4,10G,05:59:59,"$$(TVC) -i $$< -n $$(word 3,$$^) -r $$(REF_FASTA) -o $$(@D) -N 4 \
	$$(if $$(TARGETS_FILE_INTERVALS),-b $$(TARGETS_FILE_INTERVALS)) -p $$(TVC_SOMATIC_JSON) -m $$(TVC_MOTIF) \
	-t $$(TVC_ROOT_DIR) --primer-trim-bed $$(PRIMER_TRIM_BED) -g $$(basename $$(notdir $$<)) \
	&& mv $$(@D)/TSVC_variants.vcf $$@")

tvc/vcf/$1_$2/TSVC_variants_prelimiary.fpft.vcf : tvc/vcf/$1_$2/TSVC_variants_preliminary.vcf bam/$1.bam bam/$1.bam.bai
	$$(call LSCRIPT_MEM,4G,00:29:59,"awk '! /\#/' $$< | \
	awk '{if(length($$$$4) > length($$$$5)) print $$$$1\"\t\"($$$$2-1)\"\t\"($$$$2+length($$$$4)-1); \
	else print $$$$1\"\t\"($$$$2-1)\"\t\"($$$$2+length($$$$5)-1)}' > $$<.region && \
	$$(BAM_READCOUNT) -f $$(REF_FASTA) -l $$<.region $$(word 2,$$^) > $$<.bamrc && \
	$$(VARSCAN) fpfilter $$< $$<.bamrc --output-file $$@ --filtered-file $$@.fail \
	--min-var-freq $$(MIN_AF_SNP) --min-ref-readpos 0 --min-var-readpos 0 --min-ref-dist3 0 --min-var-dist3 0 && \
	rm $$<.region $$<.bamrc")

tvc/vcf/$1_$2/tumor/TSVC_variants.vcf : bam/$1.bam bam/$1.bam.bai tvc/vcf/$1_$2/TSVC_variants_preliminary.$$(if $$(findstring true,$$(USE_FPFILTER_FOR_TVC)),fpft.)vcf
	$$(call LSCRIPT_PARALLEL_MEM,4,10G,05:59:59,"$$(TVC) -s $$(<<<) -i $$(<) -r $$(REF_FASTA) -o $$(@D) -N 4 \
		-m $$(TVC_MOTIF) -t $$(TVC_ROOT_DIR) --primer-trim-bed $$(PRIMER_TRIM_BED) && \
	$$(LOAD_JAVA8_MODULE) && \
	$$(call SELECT_VARIANTS,6G) -R $$(REF_FASTA) --variant $$@ -o $$@.tmp --concordance $$(<<<) && \
	mv $$@.tmp $$@")

tvc/vcf/$1_$2/normal/TSVC_variants.vcf : bam/$2.bam bam/$2.bam.bai tvc/vcf/$1_$2/TSVC_variants_preliminary.$$(if $$(findstring true,$$(USE_FPFILTER_FOR_TVC)),fpft.)vcf
	$$(call LSCRIPT_PARALLEL_MEM,4,10G,05:59:59,"$$(TVC) -s $$(<<<) -i $$(<) -r $$(REF_FASTA) -o $$(@D) -N 4 \
		-m $$(TVC_MOTIF) -t $$(TVC_ROOT_DIR) --primer-trim-bed $$(PRIMER_TRIM_BED) && \
	$$(LOAD_JAVA8_MODULE) && \
	$$(call SELECT_VARIANTS,6G) -R $$(REF_FASTA) --variant $$@ -o $$@.tmp --concordance $$(<<<) && \
	mv $$@.tmp $$@")

tvc/vcf/$1_$2/TSVC_variants.vcf : tvc/vcf/$1_$2/tumor/TSVC_variants.vcf.gz tvc/vcf/$1_$2/normal/TSVC_variants.vcf.gz tvc/vcf/$1_$2/tumor/TSVC_variants.vcf.gz.tbi tvc/vcf/$1_$2/normal/TSVC_variants.vcf.gz.tbi
	$$(call LSCRIPT_MEM,5G,00:29:29,"$$(LOAD_VCFTOOLS_MODULE); $$(LOAD_TABIX_MODULE); $$(VCFTOOLS_MERGE) $$< $$(word 2,$$^) | perl -p -e \"s/NOCALL/PASS/g;\" > $$@")
endef
$(foreach pair,$(SAMPLE_PAIRS), \
	$(eval $(call tvc-somatic-vcf,$(tumor.$(pair)),$(normal.$(pair)))))


################## SUFAM from here ####################
define tvc-somatic-sufam-get-positions
tvc/sufam/$1/sufampos.vcf : $(foreach tumor,$(wordlist 1,$(shell expr $(words $(subst _,$( ),$1)) - 1),$(subst _,$( ),$1)),tvc/vcf/$(tumor)_$(lastword $(subst _,$( ),$1))/TSVC_variants.vcf)
	$$(call LSCRIPT_MEM,22G,03:59:59,"$$(LOAD_JAVA8_MODULE); $$(call COMBINE_VARIANTS,21G) \
		$$(foreach vcf,$$^,--variant $$(vcf) ) -o $$@ --genotypemergeoption UNSORTED -R $$(REF_FASTA)")
endef
$(foreach set,$(SAMPLE_SETS),\
	$(eval $(call tvc-somatic-sufam-get-positions,$(set))))
	
define tvc-somatic-sufam
tvc/sufam/$1/$2.tvc_sufam.vcf : bam/$2.bam tvc/sufam/$1/sufampos.vcf
	$$(call LSCRIPT_PARALLEL_MEM,4,10G,05:59:59,"$$(TVC) -s $$(<<) -i $$(<) -r $$(REF_FASTA) -o $$(@D)/$2 -N 4 \
		-m $$(TVC_MOTIF) -t $$(TVC_ROOT_DIR) --primer-trim-bed $$(PRIMER_TRIM_BED) && \
		$$(LOAD_JAVA8_MODULE) && \
		$$(call SELECT_VARIANTS,6G) -R $$(REF_FASTA) --variant $$(@D)/$2/TSVC_variants.vcf -o $$@ --concordance $$(<<)")
endef
$(foreach set,$(SAMPLE_SETS),\
	$(foreach sample,$(subst _,$( ),$(set)),\
		$(eval $(call tvc-somatic-sufam,$(set),$(sample)))))

define tvc-somatic-sufam-merge
tvc/sufam/$1.tvc_sufam.vcf : $(foreach sample,$(subst _,$( ),$1),tvc/sufam/$1/$(sample).tvc_sufam.vcf.gz) $(foreach sample,$(subst _,$( ),$1),tvc/sufam/$1/$(sample).tvc_sufam.vcf.gz.tbi)
	$$(call LSCRIPT_MEM,5G,00:29:29,"$$(LOAD_VCFTOOLS_MODULE); $$(LOAD_TABIX_MODULE); $$(VCFTOOLS_MERGE) $$(filter-out %.tbi,$$^) > $$@")
endef
$(foreach set,$(SAMPLE_SETS),\
	$(eval $(call tvc-somatic-sufam-merge,$(set))))


include usb-modules/vcf_tools/vcftools.mk
include usb-modules/variant_callers/TVC.mk

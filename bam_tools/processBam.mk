# various bam processing steps
# can be used to reprocess bam files, merge them, or merge and reprocess bam files
# possible post-processing steps are defined in modules/aligners/align.inc
##### MAKE INCLUDES #####

ifndef PROCESS_BAM_MK

include usb-modules/Makefile.inc
include usb-modules/config.inc

LOGDIR ?= log/process_bam.$(NOW)

# not primary alignment
# read fails platform/vendor quality checks
BAM_FILTER_FLAGS ?= 768

.DELETE_ON_ERROR:
.SECONDARY: 

BAMS = $(foreach sample,$(SAMPLES),bam/$(sample).bam)
ifeq ($(BAM_REPROCESS),true)
processed_bams : $(BAMS) $(addsuffix .bai,$(BAMS))
bam/%.bam : unprocessed_bam/%.$(BAM_SUFFIX)
	$(INIT) ln -f $< $@
else
ifeq ($(MERGE_SPLIT_BAMS),true)
merged_bams : $(BAMS) $(addsuffix .bai,$(BAMS))
bam/%.bam : unprocessed_bam/%$(if $(findstring true,$(BAM_FIX_RG)),.rg).bam
	$(INIT) ln -f $< $@
endif
endif

ifeq ($(MERGE_SPLIT_BAMS),true)
define bam-header
unprocessed_bam/$1.header.sam : $$(foreach split,$2,unprocessed_bam/$$(split).bam)
	$$(INIT) $$(LOAD_SAMTOOLS_MODULE); $$(SAMTOOLS) view -H $$< | grep -v '^@RG' > $$@.tmp; \
	for bam in $$(^M); do $$(SAMTOOLS) view -H $$$$bam | grep '^@RG' >> $$@.tmp; done; \
	uniq $$@.tmp > $$@ && $(RM) $$@.tmp
endef
$(foreach sample,$(SPLIT_SAMPLES),$(eval $(call bam-header,$(sample),$(split.$(sample)))))

define merged-bam
unprocessed_bam/$1.bam : unprocessed_bam/$1.header.sam $$(foreach split,$2,unprocessed_bam/$$(split).bam)
	$$(call LSCRIPT_MEM,12G,01:59:59,"$$(LOAD_SAMTOOLS_MODULE); $$(SAMTOOLS) merge -f -h $$< $$@ $$(filter %.bam,$$^)")
endef
$(foreach sample,$(SPLIT_SAMPLES),$(eval $(call merged-bam,$(sample),$(split.$(sample)))))
endif


# indices
# if bam file is a symlink, need to create a symlink to index
BAMS = $(foreach sample,$(SAMPLES),bam/$(sample).bam)
index : $(BAMS) $(addsuffix .bai,$(BAMS))

%.bam.bai : %.bam
	$(call LSCRIPT_CHECK_MEM,3G,00:29:59,"$(LOAD_SAMTOOLS_MODULE); $(SAMTOOLS) index $<")

%.bai : %.bam.bai
	$(INIT) cp $< $@

# limit coverage
#%.dcov.bam : %.bam
#	$(call LSCRIPT_MEM,18G,00:59:59,"$(call GATK_MEM,18G) -T PrintReads -R $(REF_FASTA) -I $< -dcov 50 -o $@")

# filter
%.filtered.bam : %.bam
	$(call LSCRIPT_MEM,6G,00:59:59,"$(LOAD_SAMTOOLS_MODULE); $(SAMTOOLS) view -bF $(BAM_FILTER_FLAGS) $< > $@ && $(RM) $<")

%.fixmate.bam : %.bam
	$(call LSCRIPT_MEM,9G,01:59:59,"$(LOAD_JAVA8_MODULE); $(FIX_MATE) I=$< O=$@ && $(RM) $<")

# recalibrate base quality
%.recal_report.grp : %.bam %.bai
	$(call LSCRIPT_MEM,11G,02:59:59,"$(LOAD_JAVA8_MODULE); $(BASE_RECALIBRATOR) -R $(REF_FASTA) $(BAM_BASE_RECAL_OPTS) -I $< -o $@")

%.sorted.bam : %.bam
	$(call LSCRIPT_MEM,20G,03:59:59,"$(LOAD_JAVA8_MODULE); $(SORT_SAM) I=$< O=$@ SO=coordinate VERBOSITY=ERROR && $(RM) $<")

%.markdup.bam : %.bam
	$(call LSCRIPT_MEM,14G,03:59:59,"$(MKDIR) metrics; $(LOAD_JAVA8_MODULE); $(MARK_DUP) I=$< O=$@ \
		METRICS_FILE=metrics/$(call strip-suffix,$(@F)).dup_metrics.txt && $(RM) $<")

%.rmdup.bam : %.bam
	$(call LSCRIPT_MEM,4G,02:59:59,"$(LOAD_SAMTOOLS_MODULE); $(SAMTOOLS) rmdup $< $@ && $(RM) $<")

# clean sam files
%.clean.bam : %.bam
	$(call LSCRIPT_MEM,6G,01:59:59,"$(LOAD_JAVA8_MODULE); $(CLEANBAM) I=$< O=$@")

# add rg
%.rg.bam : %.bam
	$(call LSCRIPT_MEM,12G,00:59:59,"$(LOAD_JAVA8_MODULE); $(ADD_RG) I=$< O=$@ RGLB=$(call strip-suffix,$(@F)) \
		RGPL=$(SEQ_PLATFORM) RGPU=00000000 RGSM=$(call strip-suffix,$(@F)) RGID=$(call strip-suffix,$(@F)) \
		VERBOSITY=ERROR && $(RM) $<")

# if SPLIT_CHR is set to true, we will split realn processing by chromosome
ifeq ($(SPLIT_CHR),true)
# indel realignment intervals (i.e. where to do MSA)
# split by samples and chromosomes
# %=sample
# $(eval $(call chr-target-aln,chromosome))
define chr-target-realn
%.$1.chr_split.intervals : %.bam %.bam.bai
	$$(call LSCRIPT_PARALLEL_MEM,4,3G,00:29:59,"$$(LOAD_JAVA8_MODULE); $$(REALIGN_TARGET_CREATOR) \
		-I $$(<) -L $1 -nt 4 -R $$(REF_FASTA) -o $$@ $$(BAM_REALN_TARGET_OPTS)")
endef
$(foreach chr,$(CHROMOSOMES),$(eval $(call chr-target-realn,$(chr))))

# indel realignment per chromosome
# only realign if intervals is non-empty
# %=sample
# $(eval $(call chr-aln,chromosome))
define chr-realn
%.$(1).chr_realn.bam : %.bam %.$(1).chr_split.intervals %.bam.bai
	$$(call LSCRIPT_MEM,9G,02:59:59,"$$(LOAD_JAVA8_MODULE); \
		if [[ -s $$(word 2,$$^) ]]; then $$(INDEL_REALIGN) \
			-I $$(<) -R $$(REF_FASTA) -L $1 -targetIntervals $$(word 2,$$^) \
			-o $$(@) $$(BAM_REALN_OPTS); \
		else $$(PRINT_READS) -R $$(REF_FASTA) -I $$< -L $1 -o $$@ ; fi")
endef
$(foreach chr,$(CHROMOSOMES),$(eval $(call chr-realn,$(chr))))

# merge sample realn chromosome bams
%.realn.bam : $(foreach chr,$(CHROMOSOMES),%.$(chr).chr_realn.bam) $(foreach chr,$(CHROMOSOMES),%.$(chr).chr_realn.bai)
	$(call LSCRIPT_PARALLEL_MEM,2,10G,02:59:59,"$(LOAD_JAVA8_MODULE); $(MERGE_SAMS) $(foreach i,$(filter %.bam,$^), I=$(i)) \
		SORT_ORDER=coordinate O=$@ USE_THREADING=true && $(RM) $^ $(@:.realn.bam=.bam)")

# merge sample recal chromosome bams
%.recal.bam : $(foreach chr,$(CHROMOSOMES),%.$(chr).chr_recal.bam) $(foreach chr,$(CHROMOSOMES),%.$(chr).chr_recal.bai)
	$(call LSCRIPT_PARALLEL_MEM,2,10G,02:59:59,"$(LOAD_JAVA8_MODULE); $(MERGE_SAMS) \
		$(foreach i,$(filter %.bam,$^), I=$(i)) SORT_ORDER=coordinate O=$@ USE_THREADING=true && $(RM) $^ $(@:.recal.bam=.bam)")

define chr-recal
%.$1.chr_recal.bam : %.bam %.recal_report.grp
	$$(call LSCRIPT_MEM,11G,02:59:59,"$$(LOAD_JAVA8_MODULE); $$(PRINT_READS) -L $1 -R $$(REF_FASTA) -I $$< -BQSR $$(<<) -o $$@")
endef
$(foreach chr,$(CHROMOSOMES),$(eval $(call chr-recal,$(chr))))

else # no splitting by chr

# recalibration
%.recal.bam : %.bam %.recal_report.grp
	$(call LSCRIPT_MEM,11G,02:59:29,"$(LOAD_JAVA8_MODULE); $(PRINT_READS) -R $(REF_FASTA) -I $< -BQSR $(word 2,$^) -o $@ && $(RM) $<")

%.realn.bam : %.bam %.intervals %.bam.bai
	if [[ -s $(word 2,$^) ]]; then $(call LSCRIPT_MEM,9G,02:59:59,"$(LOAD_JAVA8_MODULE); $(INDEL_REALIGN) \
	-I $< -R $(REF_FASTA) -targetIntervals $(<<) -o $@ $(BAM_REALN_OPTS) && $(RM) $<") ; \
	else mv $< $@ ; fi

%.intervals : %.bam %.bam.bai
	$(call LSCRIPT_PARALLEL_MEM,4,3G,00:29:59,"$(LOAD_JAVA8_MODULE); $(REALIGN_TARGET_CREATOR) \
	-I $< -nt 4 -R $(REF_FASTA) -o $@ $(BAM_REALN_TARGET_OPTS)")
endif

endif
PROCESS_BAM_MK = true



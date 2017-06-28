# various bam processing steps
# can be used to reprocess bam files, merge them, or merge and reprocess bam files
# possible post-processing steps are defined in modules/aligners/align.inc
##### MAKE INCLUDES #####

ifndef PROCESS_BAM_MK

include usb-modules/Makefile.inc
include usb-modules/config.inc

LOGDIR ?= log/process_bam.$(NOW)

.DELETE_ON_ERROR:
.SECONDARY: 

BAMS = $(foreach sample,$(SAMPLES),bam/$(sample).bam)
$(info BAM_REPROCESS is $(BAM_REPROCESS))
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
index : $(addsuffix .bai,$(BAMS))

%.bam.bai : %.bam
	$(call LSCRIPT_CHECK_MEM,3G,00:29:59,"$(LOAD_SAMTOOLS_MODULE); $(SAMTOOLS) index $< && ln -f $@ $*.bai")

%.bai : %.bam.bai
	$(INIT) ln -f $@ $<

# limit coverage
#%.dcov.bam : %.bam
#	$(call LSCRIPT_MEM,18G,00:59:59,"$(call GATK_MEM,18G) -T PrintReads -R $(REF_FASTA) -I $< -dcov 50 -o $@")

%.downsampled.bam : %.bam
	$(call LSCRIPT_CHECK_MEM,20G,05:29:59,"$(LOAD_SAMTOOLS_MODULE); $(SAMTOOLS) view -bh -s $(SAMTOOLS_DOWNSAMPLE_FACTOR) $< > $@")

# filter
%.filtered.bam : %.bam
	$(call LSCRIPT_MEM,6G,00:59:59,"$(LOAD_SAMTOOLS_MODULE); $(SAMTOOLS) view -bF $(BAM_FILTER_FLAGS) $< > $@ && $(RM) $<")

%.fixmate.bam : %.bam
	$(call LSCRIPT_MEM,9G,05:59:59,"$(LOAD_JAVA8_MODULE); $(call FIX_MATE,8G) I=$< O=$@ && $(RM) $<")

# recalibrate base quality
%.recal_report.grp : %.bam %.bai
	$(call LSCRIPT_MEM,20G,23:59:59,"$(LOAD_JAVA8_MODULE); $(call BASE_RECALIBRATOR,19G) \
		-R $(REF_FASTA) $(BAM_BASE_RECAL_OPTS) -I $< -o $@")

%.sorted.bam : %.bam
	$(call LSCRIPT_MEM,20G,05:59:59,"$(LOAD_JAVA8_MODULE); $(call SORT_SAM,19G) I=$< O=$@ SO=coordinate VERBOSITY=ERROR && $(RM) $<")

%.nsorted.bam : %.bam
	$(call LSCRIPT_MEM,20G,05:59:59,"$(LOAD_JAVA8_MODULE); $(call SORT_SAM,19G) I=$< O=$@ SO=queryname")

%.markdup.bam : %.bam
	$(call LSCRIPT_MEM,20G,23:59:59,"$(MKDIR) metrics; $(LOAD_JAVA8_MODULE); $(call MARK_DUP,19G) I=$< O=$@ \
		METRICS_FILE=metrics/$(call strip-suffix,$(@F)).dup_metrics.txt && $(RM) $<")

%.rmdup.bam : %.bam
	$(call LSCRIPT_MEM,4G,05:59:59,"$(LOAD_SAMTOOLS_MODULE); $(SAMTOOLS) rmdup $< $@ && $(RM) $<")

%.splitntrim.bam : %.bam
	$(call LSCRIPT_MEM,20G,05:59:59,"$(LOAD_JAVA8_MODULE); $(call SPLIT_N_TRIM,19G) -I $< -o $@ \
		-rf ReassignOneMappingQuality -RMQF 255 -RMQT 60 -U ALLOW_N_CIGAR_READS -R $(REF_FASTA) && $(RM) $<")

# clean sam files
%.clean.bam : %.bam
	$(call LSCRIPT_MEM,6G,01:59:59,"$(LOAD_JAVA8_MODULE); $(call CLEANBAM,5G) I=$< O=$@")

# add rg
ifeq ($(findstring ILLUMINA,$(SEQ_PLATFORM)),ILLUMINA)
%.rg.bam : %.bam
	$(call LSCRIPT_MEM,12G,00:59:59,"$(LOAD_JAVA8_MODULE); $(call ADD_RG,11G) I=$< O=$@ RGLB=$(call strip-suffix,$(@F)) \
		RGPL=$(SEQ_PLATFORM) RGPU=00000000 RGSM=$(call strip-suffix,$(@F)) RGID=$(call strip-suffix,$(@F)) \
		VERBOSITY=ERROR && $(RM) $<")
endif

ifeq ($(findstring IONTORRENT,$(SEQ_PLATFORM)),IONTORRENT)
%.rg.bam : %.bam
	$(INIT) $(LOAD_SAMTOOLS_MODULE); \
		samplename=`basename $< .bam` && \
		$(SAMTOOLS) view -H $< | sed "s/SM:[a-zA-Z0-9 _-\.]*/SM:$${samplename}/" > $<.header
	$(call LSCRIPT_MEM,4G,00:29:29,"$(LOAD_SAMTOOLS_MODULE); $(SAMTOOLS) reheader $<.header $< > $@")
endif



%.read_len : %.bam
	$(call LSCRIPT_MEM,4G,00:29:59,"$(LOAD_SAMTOOLS_MODULE); \
	$(SAMTOOLS) view $< | awk '{ print length($$10) }' | sort -n | uniq -c | sort -rn | sed 's/^ \+//' | awk ' > $@")

# if SPLIT_CHR is set to true, we will split realn processing by chromosome
ifeq ($(SPLIT_CHR),true)
$(info CHR is $(CHROMOSOMES))
# indel realignment intervals (i.e. where to do MSA)
# split by samples and chromosomes
# %=sample
# $(eval $(call chr-target-aln,chromosome))
define chr-target-realn
%.$1.chr_split.intervals : %.bam %.bam.bai
	$$(call LSCRIPT_PARALLEL_MEM,4,4G,00:29:59,"$$(LOAD_JAVA8_MODULE); $$(call REALIGN_TARGET_CREATOR,3G) \
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
	if [[ -s $$(word 2,$$^) ]]; then $$(call INDEL_REALIGN,8G) \
	-I $$(<) -R $$(REF_FASTA) -L $1 -targetIntervals $$(word 2,$$^) \
	-o $$(@) $$(BAM_REALN_OPTS); \
	else $$(call PRINT_READS,8G) -R $$(REF_FASTA) -I $$< -L $1 -o $$@ ; fi")
endef
$(foreach chr,$(CHROMOSOMES),$(eval $(call chr-realn,$(chr))))

# merge sample realn chromosome bams
%.realn.bam : $(foreach chr,$(CHROMOSOMES),%.$(chr).chr_realn.bam) $(foreach chr,$(CHROMOSOMES),%.$(chr).chr_realn.bai)
	$(call LSCRIPT_PARALLEL_MEM,2,16G,02:59:59,"$(LOAD_JAVA8_MODULE); $(call MERGE_SAMS,14G) $(foreach i,$(filter %.bam,$^), I=$(i)) \
	SORT_ORDER=coordinate O=$@ USE_THREADING=true && $(RM) $^ $(@:.realn.bam=.bam)")

# merge sample recal chromosome bams
%.recal.bam : $(foreach chr,$(CHROMOSOMES),%.$(chr).chr_recal.bam) $(foreach chr,$(CHROMOSOMES),%.$(chr).chr_recal.bai)
	$(call LSCRIPT_PARALLEL_MEM,4,10G,23:59:59,"$(LOAD_JAVA8_MODULE); $(call MERGE_SAMS,9G) \
	$(foreach i,$(filter %.bam,$^), I=$(i)) SORT_ORDER=coordinate O=$@ USE_THREADING=true && $(RM) $^ $(@:.recal.bam=.bam)")

define chr-recal
%.$1.chr_recal.bam : %.bam %.recal_report.grp
	$$(call LSCRIPT_MEM,11G,23:59:59,"$$(LOAD_JAVA8_MODULE); $$(call PRINT_READS,10G) -L $1 -R $$(REF_FASTA) -I $$< -BQSR $$(<<) -o $$@")
endef
$(foreach chr,$(CHROMOSOMES),$(eval $(call chr-recal,$(chr))))

else # no splitting by chr

# recalibration
%.recal.bam : %.bam %.recal_report.grp
	$(call LSCRIPT_MEM,11G,02:59:29,"$(LOAD_JAVA8_MODULE); $(call PRINT_READS,10G) -R $(REF_FASTA) -I $< -BQSR $(word 2,$^) -o $@ && $(RM) $<")

%.realn.bam : %.bam %.intervals %.bam.bai
	if [[ -s $(word 2,$^) ]]; then $(call LSCRIPT_MEM,9G,02:59:59,"$(LOAD_JAVA8_MODULE); $(call INDEL_REALIGN,8G) \
	-I $< -R $(REF_FASTA) -targetIntervals $(<<) -o $@ $(BAM_REALN_OPTS) && $(RM) $<") ; \
	else mv $< $@ ; fi

%.intervals : %.bam %.bam.bai
	$(call LSCRIPT_PARALLEL_MEM,4,3G,00:29:59,"$(LOAD_JAVA8_MODULE); $(call REALIGN_TARGET_CREATOR,2G) \
	-I $< -nt 4 -R $(REF_FASTA) -o $@ $(BAM_REALN_TARGET_OPTS)")
endif

endif
PROCESS_BAM_MK = true



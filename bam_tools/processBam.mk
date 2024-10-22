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
	$$(call LSCRIPT_MEM,$(RESOURCE_REQ_MEDIUM_MEM),($RESOURCE_REQ_SHORT),"$$(LOAD_SAMTOOLS_MODULE); $$(SAMTOOLS) merge -f -h $$< $$@ $$(filter %.bam,$$^)")
endef
#$$(call LSCRIPT_MEM,12G,01:59:59,"$$(LOAD_SAMTOOLS_MODULE); $$(SAMTOOLS) merge -f -h $$< $$@ $$(filter %.bam,$$^)")
$(foreach sample,$(SPLIT_SAMPLES),$(eval $(call merged-bam,$(sample),$(split.$(sample)))))
endif

# indices
# if bam file is a symlink, need to create a symlink to index
index : $(addsuffix .bai,$(BAMS))

%.bam.bai : %.bam
	$(call LSCRIPT_CHECK_MEM,$(RESOURCE_REQ_LOWMEM),$(RESOURCE_REQ_SHORT),"$(LOAD_SAMTOOLS_MODULE); sleep 5; $(SAMTOOLS) index $< && ln -f $@ $*.bai")

%.bai : %.bam.bai
	$(INIT) sleep 5; ln -f $@ $<

# limit coverage
#%.dcov.bam : %.bam
#	$(call LSCRIPT_MEM,18G,00:59:59,"$(call GATK_MEM,18G) -T PrintReads -R $(REF_FASTA) -I $< -dcov 50 -o $@")

%.downsampled.bam : %.bam
	$(call LSCRIPT_CHECK_MEM,$(RESOURCE_REQ_MEDIUM_MEM),$(RESOURCE_REQ_SHORT),"$(LOAD_SAMTOOLS_MODULE); $(SAMTOOLS) view -bh -s $(SAMTOOLS_DOWNSAMPLE_FACTOR) $< > $@")

# filter
%.filtered.bam : %.bam
	$(call LSCRIPT_MEM,$(RESOURCE_REQ_LOWMEM),$(RESOURCE_REQ_MEDIUM),"$(LOAD_SAMTOOLS_MODULE); \
		$(SAMTOOLS) view -bF $(BAM_FILTER_FLAGS) $< > $@ && $(RM) $<")

%.fixmate.bam : %.bam
	$(call LSCRIPT_MEM,$(RESOURCE_REQ_MEDIUM_MEM),$(RESOURCE_REQ_SHORT),"$(LOAD_JAVA8_MODULE); \
		$(call PICARD,FixMateInformation,$(RESOURCE_REQ_MEDIUM_MEM)) I=$< O=$@ && $(RM) $<")

# recalibrate base quality
%.recal_report.grp : %.bam %.bai
	$(call LSCRIPT_MEM,$(RESOURCE_REQ_MEDIUM_MEM),$(RESOURCE_REQ_LONG),"$(LOAD_JAVA8_MODULE); \
		$(call GATK,BaseRecalibrator,$(RESOURCE_REQ_HIGHMEM)) \
		-R $(REF_FASTA) $(BAM_BASE_RECAL_OPTS) -I $< -o $@")

%.reordered.bam : %.bam $(REF_DICT)
	$(call LSCRIPT_MEM,$(RESOURCE_REQ_HIGHMEM),$(RESOURCE_REQ_SHORT),"$(LOAD_JAVA8_MODULE); \
		$(call PICARD,ReorderSam,$(RESOURCE_REQ_HIGHMEM)) I=$< O=$@ REFERENCE=$(REF_FASTA) && $(RM) $<")

%.sorted.bam : %.bam
	$(call LSCRIPT_MEM,$(RESOURCE_REQ_HIGHMEM),$(RESOURCE_REQ_LONG),"$(LOAD_JAVA8_MODULE); \
		$(call PICARD,SortSam,$(RESOURCE_REQ_HIGHMEM)) I=$< O=$@ SO=coordinate VERBOSITY=ERROR && $(RM) $<")

%.nsorted.bam : %.bam
	$(call LSCRIPT_MEM,$(RESOURCE_REQ_HIGHMEM),$(RESOURCE_REQ_LONG),"$(LOAD_JAVA8_MODULE); \
		$(call PICARD,SortSam,$(RESOURCE_REQ_HIGHMEM)) I=$< O=$@ SO=queryname")

%.markdup.bam : %.bam
	$(call LSCRIPT_MEM,$(RESOURCE_REQ_VVHIGHMEM),$(RESOURCE_REQ_LONG),"$(MKDIR) metrics; $(LOAD_JAVA8_MODULE); \
		$(call PICARD,MarkDuplicates,$(RESOURCE_REQ_VVHIGHMEM)) I=$< O=$@ TMP_DIR=$(TMPDIR) \
		METRICS_FILE=metrics/$(call strip-suffix,$(@F)).dup_metrics.txt && $(RM) $<")

%.rmdup.bam : %.bam
	$(call LSCRIPT_MEM,$(RESOURCE_REQ_LOWMEM),$(RESOURCE_REQ_SHORT),"$(LOAD_SAMTOOLS_MODULE); $(SAMTOOLS) rmdup $< $@ && $(RM) $<")

%.splitntrim.bam : %.bam
	$(call LSCRIPT_MEM,$(RESOURCE_REQ_MEDIUM_MEM),$(RESOURCE_REQ_MEDIUM),"$(LOAD_JAVA8_MODULE); \
	$(call GATK,SplitNCigarReads,$(RESOURCE_REQ_MEDIUM_MEM)) -I $< -o $@ \
	-rf ReassignOneMappingQuality -RMQF 255 -RMQT 60 -U ALLOW_N_CIGAR_READS -R $(REF_FASTA) && $(RM) $<")

# clean sam files
%.clean.bam : %.bam
	$(call LSCRIPT_MEM,$(RESOURCE_REQ_LOWMEM),$(RESOURCE_REQ_SHORT),"$(LOAD_JAVA8_MODULE); $(call PICARD,CleanSam,$(RESOURCE_REQ_LOWMEM)) I=$< O=$@")

%.endtoend.bam : %.bam
	$(call LSCRIPT_MEM,$(RESOURCE_REQ_LOWMEM),$(RESOURCE_REQ_VSHORT),"$(LOAD_SAMTOOLS_MODULE); \
		TMPSTART=`mktemp` && TMPEND=`mktemp` && \
		cut -f1$(,)2 $(TARGETS_FILE_INTERVALS) | awk 'BEGIN { OFS = \"\t\" } {print \$$1\"\t\"\$$2\"\t\"\$$2+1}' > \$$TMPSTART && \
		cut -f1$(,)3 $(TARGETS_FILE_INTERVALS) | awk 'BEGIN { OFS = \"\t\" } {print \$$1\"\t\"\$$2-1\"\t\"\$$2}' > \$$TMPEND && \
		$(SAMTOOLS) view -bh -L \$$TMPSTART $< | $(SAMTOOLS) view -bh -L \$$TMPEND - > $@")

# add rg
ifeq ($(findstring ILLUMINA,$(SEQ_PLATFORM)),ILLUMINA)
%.rg.bam : %.bam
	$(call LSCRIPT_MEM,$(RESOURCE_REQ_MEDIUM_MEM),$(RESOURCE_REQ_SHORT),"$(LOAD_JAVA8_MODULE); \
		$(call PICARD,AddOrReplaceReadGroups,$(RESOURCE_REQ_MEDIUM_MEM)) I=$< O=$@ RGLB=$(call strip-suffix,$(@F)) \
		RGPL=$(SEQ_PLATFORM) RGPU=00000000 RGSM=$(call strip-suffix,$(@F)) RGID=$(call strip-suffix,$(@F)) \
		VERBOSITY=ERROR && $(RM) $<")
endif

ifeq ($(findstring IONTORRENT,$(SEQ_PLATFORM)),IONTORRENT)
%.rg.bam : %.bam
	$(INIT) $(LOAD_SAMTOOLS_MODULE); \
		samplename=`basename $< .bam` && \
		$(SAMTOOLS) view -H $< | sed "s/SM:[a-zA-Z0-9 _-\.]*/SM:$${samplename}/" > $<.header
	$(call LSCRIPT_MEM,$(RESOURCE_REQ_LOWMEM),$(RESOURCE_REQ_VSHORT),"$(LOAD_SAMTOOLS_MODULE); $(SAMTOOLS) reheader $<.header $< > $@")
endif

%.read_len : %.bam
	$(call LSCRIPT_MEM,$(RESOURCE_REQ_LOWMEM),$(RESOURCE_REQ_VSHORT),"$(LOAD_SAMTOOLS_MODULE); \
	$(SAMTOOLS) view $< | awk '{ print length($$10) }' | sort -n | uniq -c | sort -rn | sed 's/^ \+//' | awk ' > $@")

# if SPLIT_CHR is set to true, we will split realn processing by chromosome
ifeq ($(SPLIT_CHR),true)
$(info CHR is $(CHROMOSOMES))

define chr-target-realn
%.$1.chr_split.intervals : %.bam %.bam.bai
	$$(call LSCRIPT_PARALLEL_MEM,4,$$(RESOURCE_REQ_LOWMEM),$$(RESOURCE_REQ_VSHORT),"$$(LOAD_JAVA8_MODULE); \
		$$(call GATK,RealignerTargetCreator,$$(RESOURCE_REQ_LOWMEM)) \
		-I $$(<) -L $1 -nt 4 -R $$(REF_FASTA) -o $$@ $$(BAM_REALN_TARGET_OPTS)")
endef
$(foreach chr,$(CHROMOSOMES),$(eval $(call chr-target-realn,$(chr))))

# only realign if intervals is non-empty
define chr-realn
%.$(1).chr_realn.bam : %.bam %.$(1).chr_split.intervals %.bam.bai
	$$(call LSCRIPT_MEM,$$(RESOURCE_REQ_MEDIUM_MEM),$$(RESOURCE_REQ_SHORT),"$$(LOAD_JAVA8_MODULE); \
	if [[ -s $$(word 2,$$^) ]]; then $$(call GATK,IndelRealigner,$$(RESOURCE_REQ_MEDIUM_MEM)) \
	-I $$(<) -R $$(REF_FASTA) -L $1 -targetIntervals $$(word 2,$$^) \
	-o $$(@) $$(BAM_REALN_OPTS); \
	else $$(call GATK,PrintReads,$$(RESOURCE_REQ_MEDIUM_MEM)) -R $$(REF_FASTA) -I $$< -L $1 -o $$@ ; fi")
endef
$(foreach chr,$(CHROMOSOMES),$(eval $(call chr-realn,$(chr))))

# merge sample realn chromosome bams
%.realn.bam : $(foreach chr,$(CHROMOSOMES),%.$(chr).chr_realn.bam) $(foreach chr,$(CHROMOSOMES),%.$(chr).chr_realn.bai)
	$(call LSCRIPT_MEM,$(RESOURCE_REQ_HIGHMEM),$(RESOURCE_REQ_LONG),"$(LOAD_JAVA8_MODULE); \
	$(call PICARD,MergeSamFiles,$(RESOURCE_REQ_HIGHMEM)) \
	$(foreach i,$(filter %.bam,$^), I=$(i)) SORT_ORDER=coordinate O=$@ USE_THREADING=true && $(RM) $^ $(@:.realn.bam=.bam)")

# merge sample recal chromosome bams
%.recal.bam : $(foreach chr,$(CHROMOSOMES),%.$(chr).chr_recal.bam) $(foreach chr,$(CHROMOSOMES),%.$(chr).chr_recal.bai)
	$(call LSCRIPT_MEM,$(RESOURCE_REQ_VHIGHMEM),$(RESOURCE_REQ_VLONG),"$(LOAD_JAVA8_MODULE); \
	$(call PICARD,MergeSamFiles,$(RESOURCE_REQ_VHIGHMEM)) \
	$(foreach i,$(filter %.bam,$^), I=$(i)) SORT_ORDER=coordinate O=$@ USE_THREADING=true && $(RM) $^ $(@:.recal.bam=.bam)")

define chr-recal
%.$1.chr_recal.bam : %.bam %.recal_report.grp
	$$(call LSCRIPT_MEM,$$(RESOURCE_REQ_MEDIUM_MEM),$$(RESOURCE_REQ_LONG),"$$(LOAD_JAVA8_MODULE); \
		$$(call GATK,PrintReads,$$(RESOURCE_REQ_MEDIUM_MEM)) -L $1 -R $$(REF_FASTA) -I $$< -BQSR $$(<<) -o $$@")
endef
$(foreach chr,$(CHROMOSOMES),$(eval $(call chr-recal,$(chr))))


else # no splitting by chr

# recalibration
%.recal.bam : %.bam %.recal_report.grp
	$(call LSCRIPT_MEM,$(RESOURCE_REQ_MEDIUM_MEM),$(RESOURCE_REQ_SHORT),"$(LOAD_JAVA8_MODULE); \
	$(call GATK,PrintReads,$(RESOURCE_REQ_MEDIUM_MEM)) -R $(REF_FASTA) -I $< -BQSR $(word 2,$^) -o $@ && $(RM) $<")

%.realn.bam : %.bam %.intervals %.bam.bai
	if [[ -s $(word 2,$^) ]]; then $(call LSCRIPT_MEM,$(RESOURCE_REQ_MEDIUM_MEM),$(RESOURCE_REQ_SHORT),"$(LOAD_JAVA8_MODULE); \
	$(call GATK,IndelRealigner,$(RESOURCE_REQ_MEDIUM_MEM)) \
	-I $< -R $(REF_FASTA) -targetIntervals $(<<) -o $@ $(BAM_REALN_OPTS) && $(RM) $<") ; \
	else mv $< $@ ; fi

%.intervals : %.bam %.bam.bai
	$(call LSCRIPT_PARALLEL_MEM,4,$(RESOURCE_REQ_LOWMEM),$(RESOURCE_REQ_VSHORT),"$(LOAD_JAVA8_MODULE); \
	$(call GATK,RealignerTargetCreator,$(RESOURCE_REQ_LOWMEM)) \
	-I $< -nt 4 -R $(REF_FASTA) -o $@ $(BAM_REALN_TARGET_OPTS)")
endif

endif
PROCESS_BAM_MK = true



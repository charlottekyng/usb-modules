include usb-modules/Makefile.inc
include usb-modules/config.inc

LOGDIR ?= log/star.$(NOW)

.PHONY: star
.DELETE_ON_ERROR:

ALIGNER := star
include usb-modules/aligners/align.inc

star : $(foreach sample,$(SAMPLES),bam/$(sample).bam)

star/firstpass/%.SJ.out.tab : fastq/%.1.fastq.gz
	$(call LSCRIPT_PARALLEL_MEM,8,4G,00:29:59,"$(MKDIR) star star/firstpass/; \
	$(LOAD_STAR_MODULE); STAR --runMode alignReads \
	--runThreadN 8 --genomeDir $(STAR_GENOME_DIR) --readFilesIn $< --readFilesCommand gunzip -c \
	--alignSJoverhangMin 8 --alignSJDBoverhangMin 1 --alignIntronMax 1000000 --alignMatesGapMax 1000000\
	--outFilterType BySJout --outFilterMultimapNmax 20 --outFilterMismatchNmax 999 \
	--outFileNamePrefix $(@D)/$*. \
	--outSAMprimaryFlag AllBestScore --outSAMtype BAM SortedByCoordinate \
	--outReadsUnmapped Fastx --outMultimapperOrder Random --outSAMattrIHstart 0")


bam/%.bam : fastq/%.1.fastq.gz $(foreach sample,$(SAMPLES),star/firstpass/$(sample).SJ.out.tab)
	$(call LSCRIPT_PARALLEL_MEM,8,4G,00:29:59,"$(MKDIR) star star/secondpass/; \
	$(LOAD_STAR_MODULE); STAR --runMode alignReads \
	--runThreadN 8 --genomeDir $(STAR_GENOME_DIR) --readFilesIn $< --readFilesCommand gunzip -c \
	--alignSJoverhangMin 8 --alignSJDBoverhangMin 10 --alignIntronMax 200000 --alignMatesGapMax 200000 \
	--alignSJstitchMismatchNmax 5 -1 5 5 \
	--outFilterType BySJout --outFilterMultimapNmax 20 --outFilterMismatchNmax 999 \
	--outFileNamePrefix $(@D)/$*_ \
	--sjdbFileChrStartEnd $(filter %.SJ.out.tab,$^) \
	--outSAMprimaryFlag AllBestScore --outSAMtype BAM SortedByCoordinate \
	--outReadsUnmapped None --outMultimapperOrder Random --outSAMattrIHstart 0 \
	--chimSegmentMin 12 --chimJunctionOverhangMin 12 --chimSegmentReadGapMax parameter 3 \
	--quantMode TranscriptomeSAM && \
	ln star/secondpass/$*_Aligned.sortedByCoord.out.bam $@")

include usb-modules/fastq_tools/fastq.mk
include usb-modules/bam_tools/processBam.mk
include usb-modules/aligners/align.mk

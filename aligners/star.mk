include usb-modules/Makefile.inc
include usb-modules/config.inc

LOGDIR ?= log/star.$(NOW)

.PHONY: star
.DELETE_ON_ERROR:

star : $(foreach sample,$(SAMPLES),star/$(sample).bam)

trimmed_fastq/%_trimmed.fq.gz : fastq/%.fastq.gz
	$(call LSCRIPT_MEM,2G,00:29:59,"$(LOAD_TRIM_GALORE_MODULE); $(LOAD_FASTQC_MODULE);
	TRIM_GALORE -q 20 --output trimmed_fastq --clip_R1 $(CLIP5P) --clip_R2 $(CLIP3P)"

define merged-fastq
trimmed_fastq/$1.fastq.gz : $$(foreach split,$2,trimmed_fastq/$$(split).1_trimmed.fq.gz)
	$$(INIT) cat $$^ > $@
endef
$(foreach sample,$(SPLIT_SAMPLES),$(eval $(call merged-fastq,$(sample),$(split.$(sample)))))


star/%.firstpass.SJ.out.tab : trimmed_fastq/%.fastq.gz
	$(call LSCRIPT_PARALLEL_MEM,8,4G,00:29:59,"$(LOAD_STAR_MODULE); STAR --runMode alignReads \
	--runThreadN 8 --genomeDir $(STAR_GENOME_DIR) --readFilesIn $< --readFilesCommand gunzip -c \
	--alignSJoverhangMin 8 --alignSJDBoverhangMin 1 --alignIntronMax 1000000 --alignMatesGapMax 1000000\
	--outFilterType BySJout --outFilterMultimapNmax 20 --outFilterMismatchNmax 999 \
	--outFileNamePrefix star/`basename $< .fastq.gz`.firstpass/ \
	--outSAMprimaryFlag AllBestScore --outSAMtype BAM SortedByCoordinate \
	--outReadsUnmapped Fastx --outMultimapperOrder Random --outSAMattrIHstart 0"


bam/%.bam : trimmed_fastq/%.fastq.gz $(foreach sample,$(SAMPLES),star/$(sample).firstpass.SJ.out.tab)
	$(call LSCRIPT_PARALLEL_MEM,8,4G,00:29:59,"$(LOAD_STAR_MODULE); STAR --runMode alignReads \
	--runThreadN 8 --genomeDir $(STAR_GENOME_DIR) --readFilesIn $< --readFilesCommand gunzip -c \
	--alignSJoverhangMin 8 --alignSJDBoverhangMin 1 --alignIntronMax 1000000 --alignMatesGapMax 1000000\
	--outFilterType BySJout --outFilterMultimapNmax 20 --outFilterMismatchNmax 999 \
	--outFileNamePrefix star/`basename $< .fastq.gz`.secondpass/ \
	--sjdbFileChrStartEnd $(foreach sample,$(SAMPLES),star/$(sample).firstpass.SJ.out.tab) \
	--outSAMprimaryFlag AllBestScore --outSAMtype BAM SortedByCoordinate \
	--outReadsUnmapped Fastx --outMultimapperOrder Random --outSAMattrIHstart 0 \
	--chimSegmentMin 25 --quantMode TranscriptomeSAM && \
	ln star/`basename $< .fastq.gz`.secondpass/`basename $< .fastq.gz`Aligned.sortedByCoord.out.bam $@"

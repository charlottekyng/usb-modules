include usb-modules/Makefile.inc
include usb-modules/config.inc

LOGDIR ?= log/star.$(NOW)

.PHONY: star
.DELETE_ON_ERROR:

ALIGNER := star
include usb-modules/aligners/align.inc

BAM_NO_SORT = true
BAM_FIX_RG = true
BAM_DUP_TYPE = markdup

BAMS = $(foreach sample,$(SAMPLES),bam/$(sample).originalstar.bam) $(foreach sample,$(SAMPLES),bam/$(sample).bam)
star : $(BAMS) $(addsuffix .bai,$(BAMS)) star/all.ReadsPerGene.out.tab


star/firstpass/%.SJ.out.tab : fastq/%.1.fastq.gz $(if $(findstring true,$(PAIRED_END)),fastq/%.2.fastq.gz)
	$(call LSCRIPT_PARALLEL_MEM,8,4G,00:59:59,"$(MKDIR) star star/firstpass/; \
	$(LOAD_STAR_MODULE); STAR --runMode alignReads \
	--runThreadN 8 --genomeDir $(STAR_GENOME_DIR) --readFilesIn $< $(if $(findstring true,$(PAIRED_END)),$(word 2,$^)) \
	--readFilesCommand gunzip -c \
	--alignSJoverhangMin 8 --alignSJDBoverhangMin 1 --alignIntronMax 1000000 --alignMatesGapMax 1000000\
	--outFilterType BySJout --outFilterMultimapNmax 20 --outFilterMismatchNmax 999 \
	--outFileNamePrefix $(@D)/$*. \
	--outSAMprimaryFlag AllBestScore --outSAMtype BAM SortedByCoordinate \
	--outReadsUnmapped Fastx --outMultimapperOrder Random --outSAMattrIHstart 0")

define align-star-secondpass
star/secondpass/$1.star.bam : fastq/$1.1.fastq.gz $(if $(findstring true,$(PAIRED_END)),fastq/%.2.fastq.gz) $(foreach sample,$(SAMPLES),star/firstpass/$(sample).SJ.out.tab)
	$$(call LSCRIPT_PARALLEL_MEM,8,8G,00:59:59,"$$(MKDIR) star star/secondpass/; \
	$$(LOAD_STAR_MODULE); STAR --runMode alignReads \
	--runThreadN 8 --genomeDir $$(STAR_GENOME_DIR) --readFilesIn $$< $$(if $$(findstring true,$(PAIRED_END)),$$(word 2,$$^)) \
	--readFilesCommand gunzip -c \
	--alignSJoverhangMin 8 --alignSJDBoverhangMin 10 --alignIntronMax 200000 --alignMatesGapMax 200000 \
	--alignSJstitchMismatchNmax 5 -1 5 5 --limitSjdbInsertNsj 5000000 \
	--outFilterType BySJout --outFilterMultimapNmax 20 --outFilterMismatchNmax 999 \
	--outFileNamePrefix star/secondpass/$$*. \
	--sjdbFileChrStartEnd $$(filter %.SJ.out.tab,$$^) \
	--outSAMprimaryFlag AllBestScore --outSAMtype BAM SortedByCoordinate \
	--outReadsUnmapped None --outMultimapperOrder Random --outSAMattrIHstart 0 \
	--chimSegmentMin 12 --chimJunctionOverhangMin 12 --chimSegmentReadGapMax parameter 3 \
	--quantMode GeneCounts && mv star/secondpass/$1.Aligned.sortedByCoord.out.bam star/secondpass/$1.star.bam")
endef
$(foreach sample,$(SAMPLES),\
	$(eval $(call align-star-secondpass,$(sample))))

bam/%.originalstar.bam : star/secondpass/%.star.bam
	$(INIT) ln $< $@

star/secondpass/%.star-original.bam : star/secondpass/%.star.bam
	$(INIT) cp $< $@

bam/%.bam : star/secondpass/%.star-original.$(BAM_SUFFIX)
	$(INIT) ln $< $@ && rename "-original" "" $<

star/all.ReadsPerGene.out.tab : $(foreach sample,$(SAMPLES),bam/$(sample).bam)
	perl -p -e "s/N_unmapped/GENE\t\t\t\nN_unmapped/;" `ls star/secondpass/*.ReadsPerGene.out.tab|head -1` | cut -f 1 > $@; \
	for x in $^; do sample=`echo $$x | sed 's/.*\///; s/\..*//'`; perl -p -e "s/N_unmapped/\t$$sample\t\t\nN_unmapped/;" \
	star/secondpass/$$sample.ReadsPerGene.out.tab | cut -f 2 | paste $@ - > $@.tmp; mv $@.tmp $@; done

include usb-modules/fastq_tools/fastq.mk
include usb-modules/bam_tools/processBam.mk
include usb-modules/aligners/align.mk

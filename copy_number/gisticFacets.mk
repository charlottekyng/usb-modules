include usb-modules/Makefile.inc
include usb-modules/config.inc

LOGDIR = log/gisticFacets.$(NOW)

SHELL = usb-modules/scripts/Rshell
.SHELLFLAGS = -m $(MEM) -p $(PE) -n $(@F) -l $(LOGDIR) -e 

.ONESHELL:
.DELETE_ON_ERROR:
.SECONDARY:
.PHONY: all gistic_inputs gistic_heatmaps

MEM := 2G
PE := 1

CNV_SIZES = 100000 300000

all : gistic_inputs $(foreach size,$(CNV_SIZES),gistic/gistic_cnv$(size).timestamp)
gistic_inputs : gistic/markersfile.txt gistic/segmentationfile.txt $(foreach size,$(CNV_SIZES),gistic/cnv.$(size).txt)
gistic_heatmaps : $(foreach size,$(CNV_SIZES),gistic/gistic_cnv$(size)/gistic_cnv_heatmap.pdf)

gistic/markersfile.txt : gistic/segmentationfile.txt
	suppressPackageStartupMessages(library("rtracklayer"));
	seg <- read.table('$<', sep = '\t', stringsAsFactors = F, col.names = c('samplePair', 'chr', 'start', 'end', 'numMarkers', 'logRatio'))
	targets <- import('$(TARGETS_FILE_INTERVALS)')
	markers <- data.frame(chr = seqnames(targets), pos = start(targets))
	markers <- markers[markers$$chr %in% seg$$chr, ]
	dir.create('$(@D)', showWarnings = F)
	write.table(markers, col.names = F, file = "$@", sep = "\t", quote = F, na = "")
	
gistic/segmentationfile.txt : PE := 8
gistic/segmentationfile.txt : MEM := $(RESOURCE_REQ_LOWMEM)
gistic/segmentationfile.txt : $(foreach pair,$(SAMPLE_PAIRS),facets/cncf/$(pair).cncf.txt)
	suppressPackageStartupMessages(library("rtracklayer"));
	suppressPackageStartupMessages(library("foreach"));
	suppressPackageStartupMessages(library("doMC"));
	segFiles <- unlist(strsplit("$^", " "));
	segNames <- sub(".*/", "", sub("\\..*", "", segFiles))
	targets <- import('$(TARGETS_FILE_INTERVALS)')
	width(targets) <- 1
	registerDoMC(8)
	seg <- foreach (i = 1:length(segFiles), .combine = 'rbind') %dopar% {
		segFile <- segFiles[i]
		segName <- segNames[i]
		s <- read.delim(segFile, header = T, as.is = T, col.names=c("Chromosome","Seg","numMarkers","nhet","log2_ratio_seg","mafR","segclust","cnlr.median.clust","mafR.clust","Start","End","cf.em","tcn.em","lcn.em","clonal.cluster"))
		s <- cbind(segName, s)
		colnames(s)[1] <- "Sample"
		s[['Chromosome']][s[['Chromosome']] == 23] <- "X"
		s[['Chromosome']][s[['Chromosome']] == 24] <- "Y"
#		gr <- with(s, GRanges(seqnames = Chromosome, range = IRanges(start = Start, end = End), segmented = as.numeric(log2_ratio_seg)))

		out <- read.delim(gsub("cncf.txt", "out", segFile), as.is=T, sep="=", header=F)
		purity <- as.numeric(out[grep("Purity",out[,1]),2])
		ploidy <- as.numeric(out[grep("Ploidy",out[,1]),2])
		if(!is.na(purity) & !is.na(ploidy)) {
			s$isar_corrected <- s$log2_ratio_seg/(purity-(2*(1-purity)/(purity*ploidy)))
		} else { s$isar_corrected <- s$log2_ratio_seg }
		gr <- with(s, GRanges(seqnames = Chromosome, range = IRanges(start = Start, end = End), segmented = as.numeric(isar_corrected)))

		redGr <- reduce(gr)
		x <- findOverlaps(redGr, gr, select = 'first')
		redGr$$segmented <- gr[x]$$segmented
		# reduced the genomic range, need to intersect with targets
		numMarkers <- countOverlaps(redGr, targets)
		Start <- start(targets)[findOverlaps(redGr, targets, select = 'first')]
		End <- start(targets)[findOverlaps(redGr, targets, select = 'last')]
		seg <- data.frame(segName, chrom = seqnames(redGr), start = Start, end = End, numMarkers, segmented = redGr$$segmented)
		seg <- subset(seg, numMarkers > 0)
		seg[!duplicated(seg), ]
	}
#	splitSeg <- split(seg, list(as.factor(seg$$segName), as.factor(seg$$chrom)))
#	seg <- do.call('rbind', lapply(splitSeg, function(x) {
#		rx <- Rle(x$$segmented)
#		nx <- x[start(rx), ]
#		nx$$end <- x[end(rx), "end"]
#		nx$$numMarkers <- aggregate(x$$numMarkers, rx, sum)
#		nx
#	}))
	dir.create('$(@D)', showWarnings = F)
	write.table(seg, file = "$@", sep = "\t", row.names = F, col.names = F, quote = F)

gistic/cnv.%.txt : gistic/markersfile.txt
	suppressPackageStartupMessages(library("GenomicRanges"));
	dgv <- read.delim("$(DGV_FILE)", as.is=T)
	dgv <- dgv[which(dgv$$varianttype=="CNV"), ]
	dgv <- dgv[, c("variantaccession", "chr", "start", "end")]
	dgv$$size = dgv$$end-dgv$$start+1
	dgv <- dgv[which(dgv$$size <= $*), ]
	markers <- read.delim("$<", as.is=T, header=F)
	dgvGR <- GRanges(seqnames = dgv$$chr, ranges = IRanges(start = dgv$$start, end = dgv$$end))
	markersGR <- GRanges(seqnames = markers[,2], ranges = IRanges(start = markers[,3], end = markers[,3]))
	markers <- cbind(markers, countOverlaps(markersGR, dgvGR))
	#markers <- cbind(markers, apply(markers, 1, function(x, dgv) {
	#	length(which(dgv$$chr==x[2] & dgv$$start <= x[3] & dgv$$end >= x[3])) }, dgv))
	cnv <- markers[which(markers[,4] > 0),]
	cnv <- cbind(cnv[,1], cnv[,1])
	dir.create('$(@D)', showWarnings = F)
	write.table(cnv, file = "$@", sep = "\t", row.names = F, col.names = F, quote = F, na = "")

gistic/gistic_cnv%.timestamp : MEM := 12G
gistic/gistic_cnv%.timestamp : gistic/segmentationfile.txt gistic/markersfile.txt gistic/cnv.%.txt
	Sys.setenv(LD_LIBRARY_PATH = "/scicore/home/terracci/GROUP/usr_nobackup/local/MCR_R2014a/v83/runtime/glnxa64:/scicore/home/terracci/GROUP/usr_nobackup/local/MCR_R2014a/v83/bin/glnxa64:/scicore/home/terracci/GROUP/usr_nobackup/local/MCR_R2014a/v83/sys/os/glnxa64")
	Sys.setenv(MCR_DIR = "/scicore/home/terracci/GROUP/usr_nobackup/local/MCR_R2014a/")
	dir.create('$(@D)/gistic_cnv$*', showWarnings = F, recursive = T)
	system("umask 002; $(GISTIC) -b $(@D)/gistic_cnv$* -seg $< -mk $(<<) -refgene $(GISTIC_REF) -cnv $(<<<) $(GISTIC_OPTS) 2>&1 && touch $@")

#gistic/gistic_cnv%/gistic_cnv_heatmap.pdf : gistic/gistic_cnv%.timestamp
#	system("$(PLOT_GISTIC_HEATMAP) --out $@ gistic/gistic_cnv$*/all_thresholded.by_genes.txt")




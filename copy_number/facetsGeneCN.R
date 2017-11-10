#!/usr/bin/env Rscript
#### turn segmented copy number data to gene-based copy number with findOverlaps
## define HomDel as TCN=0, loss as TCN<ploidy, gain as TCN>ploidy, amp as TCN>=ploidy+4
## where ploidy= mode of TCN
### some variant of the below, also need one for the breast panel, IMPACT310 and exome

#---------------
# initialization
#---------------

# load base libraries
suppressMessages(pacman::p_load(optparse,RColorBrewer,GenomicRanges,plyr,dplyr,tibble,readr,stringr,tidyr,purrr,magrittr,rlist,crayon,foreach,Cairo,RMySQL,rtracklayer,colorspace,ggplot2,grid,gridExtra,RColorBrewer))
suppressPackageStartupMessages(library("facets"));

#--------------
# parse options
#--------------

optList <- list(
				make_option("--outFile", default = NULL, help = "output file"),
#				make_option("--mysqlHost", default = '10.0.200.48', help = "MySQL server hostname"),
#				make_option("--mysqlPort", default = 38493, help = "MySQL server port"),
#				make_option("--mysqlUser", default = 'embl', help = "MySQL server username"),
#				make_option("--mysqlPassword", default = 'embl', help = "MySQL server password"),
#				make_option("--mysqlDb", default = 'homo_sapiens_core_75_37', help = "MySQL server database"),
				make_option("--genesFile", default = NULL, help = "list of genes to include (hgnc symbols)"))
parser <- OptionParser(usage = "%prog [options] [facets files]", option_list = optList);

arguments <- parse_args(parser, positional_arguments = T);
opt <- arguments$options;

if (length(arguments$args) < 1) {
	cat("Need facets output files\n")
	print_help(parser);
	stop();
} else if (is.null(opt$outFile)) {
	cat("Need output prefix\n")
	print_help(parser);
	stop();
} else if (is.null(opt$genesFile)) {
        cat("Need genes files\n")
        print_help(parser);
        stop();
} else {
	facetsFiles <- arguments$args
}

#connect <- function() dbConnect(MySQL(), host = opt$mysqlHost, port = opt$mysqlPort, user = opt$mysqlUser, password = opt$mysqlPassword, dbname = opt$mysqlDb)
#cat('Connecting to ensembl ... ')
#mydb <- connect()
#on.exit(dbDisconnect(mydb))

#query <- "select r.name as chrom,
#g.seq_region_start as start,
#g.seq_region_end as end,
#x.display_label as hgnc,
#k.band as band
#from gene as g
#join seq_region as r on g.seq_region_id = r.seq_region_id
#join xref as x on g.display_xref_id = x.xref_id
#left join karyotype k on g.seq_region_id = k.seq_region_id
#and ((g.seq_region_start >= k.seq_region_start and g.seq_region_start <= k.seq_region_end)
#or (g.seq_region_end >= k.seq_region_start and g.seq_region_end <= k.seq_region_end))
#where x.external_db_id = 1100;"
#repeat {
#	rs <- try(dbSendQuery(mydb, query), silent = T)
#	if (is(rs, "try-error")) {
#		cat("Lost connection to mysql db ... ")
#		mydb <- connect()
#		cat("reconnected\n")
#	} else {
#		break
#	}
#}
#genes <- dbFetch(rs, -1)
#cat(paste("Found", nrow(genes), "records\n"))

#genes %<>% filter(chrom %in% as.character(c(1:22, "X", "Y"))) %>%
#	filter(!duplicated(hgnc)) %>% 
#	arrange(as.integer(chrom), start, end)

#if (!is.null(opt$genesFile)) {
#	g <- scan(opt$genesFile, what = 'character')
#	genes %<>% filter(hgnc %in% g)
#	absentGenes <- g[!g %in% genes$hgnc]
#	if (length(absentGenes) > 0) {
#		print("Unable to find", length(absentGenes), "in database\n");
#		cat(absentGenes, sep = '\n');
#	}
#}

#cat(paste("Filtering to", nrow(genes), "records\n"))

genes <- read.delim(opt$genesFile, as.is=T, check.names=F)
genes$chrom <- gsub("chr", "", genes$chrom)

genesGR <- genes %$% GRanges(seqnames = chrom, ranges = IRanges(start, end), band = band, hgnc = hgnc)
			
mm <- lapply(facetsFiles, function(f) {
	tab <- read.delim(f, as.is=T)
	tab$chrom[which(tab$chrom==23)] <- "X"

	tabGR <- tab %$% GRanges(seqnames = chrom, ranges = IRanges(start, end))
	mcols(tabGR) <- tab %>% select(num.mark,cnlr.median:mafR.clust,cf.em:clonal.cluster)

	fo <- findOverlaps(tabGR, genesGR)

	df <- as.data.frame(cbind(mcols(genesGR)[subjectHits(fo),], mcols(tabGR)[queryHits(fo),]))
	df %<>% group_by(hgnc) %>% top_n(1, abs(cnlr.median))

	ploidy <- median(unlist(apply(cbind(df$tcn.em, df$num.mark),1,function(x){rep(x[1], x[2])})))
#	ploidy <- as.numeric(names(ploidy)[which.max(ploidy)])

	df$GL <- 0
	df$GL[df$tcn.em < ploidy] <- -1
	df$GL[df$tcn.em == 0] <- -2
	df$GL[df$tcn.em > ploidy] <- 1
	df$GL[df$tcn.em >= ploidy + 4] <- 2

	load(gsub("cncf.txt", "Rdata", f, fixed=T))
	noise <- median(abs(out2$jointseg$cnlr-  unlist(apply(out2$out[,c("cnlr.median", "num.mark")], 1, function(x) {rep(x[1], each=x[2])}))))

	lrr <- sort(out2$jointseg$cnlr)
#	if (noise <= 0.2) { lrr <- lrr[round(0.25*length(lrr)):round(0.75*length(lrr))]
#	} else if ( noise <= 0.3 ) { lrr <- lrr[round(0.275*length(lrr)):round(0.725*length(lrr))]
#	} else { lrr <- lrr[round(0.3*length(lrr)):round(0.7*length(lrr))]}

#	if (noise <= 0.2) { lrr <- lrr[round(0.3*length(lrr)):round(0.7*length(lrr))]
#	} else if ( noise <= 0.3 ) { lrr <- lrr[round(0.275*length(lrr)):round(0.725*length(lrr))]
#	} else { 
lrr <- lrr[round(0.25*length(lrr)):round(0.75*length(lrr))] 
	df$GL2 <- 0
	df$GL2[df$cnlr.median < median(lrr)-(2.5*sd(lrr))] <- -1
	df$GL2[df$cnlr.median < median(lrr)-(7*sd(lrr))] <- -2
	df$GL2[df$cnlr.median > median(lrr)+(2*sd(lrr))] <- 1
	df$GL2[df$cnlr.median > median(lrr)+(6*sd(lrr))] <- 2

	df %>% select(hgnc, GL, GL2, tcn.em, lcn.em, cnlr.median) %>% ungroup
})
names(mm) <- facetsFiles
for (f in facetsFiles) {
	n <- sub('\\..*', '', sub('.*/', '', f))
	colnames(mm[[f]])[2:ncol(mm[[f]])] <- paste(n,  colnames(mm[[f]])[2:ncol(mm[[f]])], sep="_")
}

mm <- left_join(genes, join_all(mm, type = 'full', by="hgnc")) %>% arrange(as.integer(chrom), start, end)
#write.table(mm, file=opt$outFile, sep="\t", row.names=F, na="", quote=F)

seg_sample <- seg_chr <- seg_band <- seg_start <- seg_end <- seg_cnlr <- seg_genes <- seg_type <- seg_GLtype <- NA
for (i in grep("GL", colnames(mm))) {
	for(chr in intersect(c(1:22,"X"), unique(mm$chrom))) {
		tt <- mm[which(mm$chrom==chr),c(1:5,i), drop=F]
		tt[which(is.na(tt[,6])),6] <- -1000
		rr <- rle(tt[,6]); 
		if (rr$values[1]== -1000) {
			rr$values[1] <- rr$values[2]
		}
		if (rr$values[length(rr$values)]== -1000) {
			rr$values[length(rr$values)] <- rr$values[length(rr$values)-1]
		}
		for ( idx in which(rr$values== -1000)) {
			if (rr$values[idx-1]== rr$values[idx+1]) { rr$values[idx] <- rr$values[idx-1]}
			else {rr$values[idx] <- 0}
		}
		mm[which(mm$chrom==chr),i] <- as.vector(unlist(apply(cbind(rr$value,rr$length), 1, function(x){rep(x[1],x[2])})))

		tt <- mm[which(mm$chrom==chr),c(1:5,i), drop=F]
		rr <- rle(tt[,6]); 
		if (length(rr$length)>1) {
			cs <- cumsum(rr$lengths)
			start <- c(1,cs[1:(length(cs)-1)]+1)
			end <- cs
		} else {start <- 1; end <- rr$lengths[1] }

		for (idx in which(rr$values %in% c(-2,2))) {
			if (rr$values[idx] %in% c(-2,2)) {
				seg_sample <- c(seg_sample, colnames(mm)[i])
				seg_chr <- c(seg_chr, chr)
				seg_band <- c(seg_band, paste(tt[start[idx],"band"], tt[end[idx],"band"], sep="-"))
				seg_start <- c(seg_start, tt[start[idx],"start"])
				seg_end <- c(seg_end, tt[end[idx],"end"])
				seg_genes <- c(seg_genes, toString(mm[start[idx]:end[idx],"hgnc"]))
				seg_type <- c(seg_type, rr$values[idx])
				seg_GLtype <- c(seg_GLtype, colnames(mm)[i])
			}
		}		

	}
}
seg_type[which(seg_type==2)] <- "amp"
seg_type[which(seg_type== -2)] <- "del"
write.table(cbind(seg_sample, seg_chr, seg_band, seg_start, seg_end, seg_genes, seg_type, seg_GLtype), file=gsub("txt", "ampdel.txt", opt$outFile), sep="\t", row.names=F, na="", quote=F)
write.table(mm, file=opt$outFile, sep="\t", row.names=F, na="", quote=F)

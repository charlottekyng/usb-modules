# run the facets library

suppressPackageStartupMessages(library("optparse"));
suppressPackageStartupMessages(library("RColorBrewer"));
suppressPackageStartupMessages(library("plyr"));
suppressPackageStartupMessages(library("dplyr"));
suppressPackageStartupMessages(library("tidyr"));
suppressPackageStartupMessages(library("stringr"));
suppressPackageStartupMessages(library("magrittr"));
suppressPackageStartupMessages(library("facets"));
suppressPackageStartupMessages(library("foreach"));
#suppressPackageStartupMessages(library("Cairo"));

if (!interactive()) {
    options(warn = -1, error = quote({ traceback(); q('no', status = 1) }))
}

optList <- list(
	make_option("--seed", default = 1234),
	make_option("--snp_nbhd", default = 250, type = 'integer', help = "window size"),
	make_option("--minNDepth", default = 25, type = 'integer', help = "minimum depth in normal to keep the position"),
	make_option("--maxNDepth", default= 1000, type= 'integer', help = "maximum depth in normal to keep the position"),
	make_option("--pre_cval", default = 50, type = 'integer', help = "pre-processing critical value"),
	make_option("--cval1", default = 150, type = 'integer', help = "critical value for estimating diploid log Ratio"),
	make_option("--cval2", default = 50, type = 'integer', help = "starting critical value for segmentation (increases by 10 until success)"),
	make_option("--max_cval", default = 5000, type = 'integer', help = "maximum critical value for segmentation (increases by 10 until success)"),
	make_option("--min_nhet", default = 25, type = 'integer', help = "minimum number of heterozygote snps in a segment used for bivariate t-statistic during clustering of segment"),
	make_option("--gene_loc_file", default = '~/share/reference/IMPACT410_genes_for_copynumber.txt', type = 'character', help = "file containing gene locations"),
	make_option("--genome", default = 'b37', type = 'character', help = "genome of counts file"),
#	make_option("--chroms", default=NULL, type='character', help="chromosomes"),
	make_option("--outPrefix", default = NULL, help = "output prefix"))

parser <- OptionParser(usage = "%prog [options] [tumor-normal base counts file]", option_list = optList);

arguments <- parse_args(parser, positional_arguments = T);
opt <- arguments$options;

if (length(arguments$args) < 1) {
    cat("Need base counts file\n")
    print_help(parser);
    stop();
} else if (is.null(opt$outPrefix)) {
    cat("Need output prefix\n")
    print_help(parser);
    stop();
} else {
    baseCountFile <- arguments$args[1];
}

tumorName <- baseCountFile %>% sub('.*/', '', .) %>% sub('_.*', '', .)
normalName <- baseCountFile %>% sub('.*/', '', .) %>% sub('.*_', '', .) %>% sub('\\..*', '', .)

switch(opt$genome,
	b37={gbuild="hg19"},
	GRCh37={gbuild="hg19"},
	hg19={gbuild="hg19"},
	mm9={gbuild="mm9"},
	mm10={gbuild="mm10"},
	GRCm38={gbuild="mm10"},
       { stop(paste("Invalid Genome",opt$genome)) })

buildData=installed.packages()["facets",]
cat("#Module Info\n")
for(fi in c("Package","LibPath","Version","Built")){
    cat("#",paste(fi,":",sep=""),buildData[fi],"\n")
}
version=buildData["Version"]
cat("\n")

rcmat <- readSnpMatrix(gzfile(baseCountFile))
chromLevels=unique(rcmat[,1])
print(chromLevels)
if (gbuild %in% c("hg19", "hg18")) { chromLevels=intersect(chromLevels, c(1:22,"X"))
} else { chromLevels=intersect(chromLevels, c(1:19,"X"))}
print(chromLevels)
preOut=preProcSample(rcmat, snp.nbhd = opt$snp_nbhd, ndepth = opt$minNDepth, cval = opt$pre_cval, gbuild=gbuild, ndepthmax=opt$maxNDepth)
### Used this instead of preProc for wes_hall_pe
#    if (gbuild %in% c("hg19", "hg18"))
#        nX <- 23
#    if (gbuild %in% c("mm9", "mm10"))
#        nX <- 20
#pmat <- facets:::procSnps(rcmat, ndepth=opt$minNDepth, het.thresh = 0.25, snp.nbhd = opt$snp_nbhd, gbuild=gbuild, unmatched=F, ndepthmax=opt$maxNDepth)
#pmat$keep[which(pmat$chrom==2 & pmat$maploc>=28641233 & pmat$maploc<=28691172)] <- 1
#pmat$keep[which(pmat$chrom==19 & pmat$maploc>=32757577 & pmat$maploc<=32826160)] <- 1
#dmat <- facets:::counts2logROR(pmat[pmat$rCountT > 0, ], gbuild, unmatched=F)
#tmp <- facets:::segsnps(dmat, opt$pre_cval, hetscale=F)
#out <- list(pmat = pmat, gbuild=gbuild, n=nX)
#preOut <- c(out,tmp)


out1 <- preOut %>% procSample(cval = opt$cval1, min.nhet = opt$min_nhet)

print ("Completed preProc and proc")
cval <- opt$cval2
success <- F
while (!success && cval < opt$max_cval) {
    out2 <- preOut %>% procSample(cval = cval, min.nhet = opt$min_nhet, dipLogR = out1$dipLogR)
    print(str_c("attempting to run emncf() with cval2 = ", cval))
    fit <- tryCatch({
        out2 %>% emcncf
    }, error = function(e) {
        print(paste("Error:", e))
        return(NULL)
    })
    if (!is.null(fit)) {
        success <- T
	fit2 <- out2 %>% emcncf2
    } else {
        cval <- cval + 100
    }
}
if (!success) {
    stop("Failed to segment data\n")
} else { print ("Completed segmentation")}

formatSegmentOutput <- function(out,sampID) {
	seg=list()
	seg$ID=rep(sampID,nrow(out$out))
	seg$chrom=out$out$chr
	seg$loc.start=rep(NA,length(seg$ID))
	seg$loc.end=seg$loc.start
	seg$num.mark=out$out$num.mark
	seg$seg.mean=out$out$cnlr.median
	for(i in 1:nrow(out$out)) {
		lims=range(out$jointseg$maploc[(out$jointseg$chrom==out$out$chr[i] & out$jointseg$seg==out$out$seg[i])],na.rm=T)
		seg$loc.start[i]=lims[1]
		seg$loc.end[i]=lims[2]
	}	
	as.data.frame(seg)
}
id <- paste(tumorName, normalName, sep = '_')
out2$IGV = formatSegmentOutput(out2, id)
save(preOut, out1, out2, fit, fit2, file = str_c(opt$outPrefix, ".Rdata"), compress=T)

ff = str_c(opt$outPrefix, ".out")
cat("# Version =", version, "\n", file = ff, append = T)
cat("# Input =", basename(baseCountFile), "\n", file = ff, append = T)
cat("# tumor =", tumorName, "\n", file = ff, append = T)
cat("# normal =", normalName, "\n", file = ff, append = T)
cat("# snp.nbhd =", opt$snp_nbhd, "\n", file = ff, append = T)
cat("# cval1 =", opt$cval1, "\n", file = ff, append = T)
cat("# cval2 =", cval, "\n", file = ff, append = T)
cat("# min.nhet =", opt$min_nhet, "\n", file = ff, append = T)
cat("# genome =", opt$genome, "\n", file = ff, append = T)
cat("# Purity =", fit$purity, "\n", file = ff, append = T)
cat("# Ploidy =", fit$ploidy, "\n", file = ff, append = T)
cat("# dipLogR =", fit$dipLogR, "\n", file = ff, append = T)
cat("# dipt =", fit$dipt, "\n", file = ff, append = T)
cat("# loglik =", fit$loglik, "\n", file = ff, append = T)

#CairoPNG(file = str_c(opt$outPrefix, ".cncf.png"), height = 1100, width = 850)
#plotSample(out2, fit)
#dev.off()

if(sum(out2$out$num.mark)<=10000) { height=4; width=7} else { height=6; width=9}
pdf(file = str_c(opt$outPrefix, ".cncf.pdf"), height = height, width = width)
plotSample(out2, fit, chromlevels = chromLevels)
dev.off()

#tab <- cbind(out2$IGV[, 1:4], fit$cncf[, 2:ncol(fit$cncf)])
#write.table(tab, str_c(opt$outPrefix, ".cncf.txt"), row.names = F, quote = F, sep = '\t')

if("clonal.cluster" %in% colnames(fit2$cncf)) { clonal.cluster = fit2$cncf$clonal.cluster
} else { clonal.cluster <- rep(NA, nrow(fit2$cncf))}

writetable <- cbind(fit2$cncf[,c("chrom", "seg", "num.mark", "nhet", "cnlr.median", "mafR", "segclust", "cnlr.median.clust", "mafR.clust", "start", "end")],
	 fit2$cncf[,c("cf.em", "tcn.em", "lcn.em")], clonal.cluster=clonal.cluster)

write.table(writetable, str_c(opt$outPrefix, ".cncf.txt"), row.names = F, quote = F, sep = '\t')
warnings()


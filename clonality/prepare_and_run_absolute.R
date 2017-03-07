cat("List of functions: mutations_to_absolute, facets_to_absolute, make_param_template, read_params, run_absolute_step1/2/3\n")


mutations_to_absolute <- function(muts, sample_col="Tumor", gene_col="Gene", depth_col="Tumor.Depth", maf_col="MAF", ref_count_col=NULL, alt_count_col=NULL, chrom_col="Chromosome", pos_col="Position", ref_col="REF", alt_col="ALT", out.dir="inputmafs") {
	if (!file.exists(out.dir)) { dir.create(out.dir) }
	muts[which(muts[,chrom_col]=="X"), chrom_col] <- 23
	muts[which(muts[,chrom_col]=="Y"), chrom_col] <- 24

	if (is.null(ref_count_col) & is.null(alt_count_col)){
		muts[,maf_col] <- gsub("%", "", muts[,maf_col])
		if (max(as.numeric(muts[,maf_col]))>1) { muts[,maf_col] <- muts[,maf_col]/100 }

		muts$t_alt_count <- round(as.numeric(muts[,depth_col])*as.numeric(muts[,maf_col]))
		muts$t_ref_count <- as.numeric(muts[,depth_col])-muts$t_alt_count
	} else {
		muts$t_alt_count <- muts[,alt_count_col]
		muts$t_ref_count <- muts[,ref_count_col]
	}

	lapply(unique(muts[,sample_col]), function(x) {
		tab <- muts[which(muts[,sample_col]==x),]
		tab <- cbind(tab[,c(sample_col, gene_col, chrom_col, pos_col, ref_col, alt_col)], tab[,c( "t_ref_count", "t_alt_count")], "validated")
		colnames(tab) <- c("Tumor_Sample_Barcode", "Hugo_Symbol", "Chromosome", "Start_position", "Reference_Allele", "Tumor_Seq_Allele2", "t_ref_count", "t_alt_count", "dbSNP_Val_Status")

		write.table(tab, file=paste(out.dir, "/", x, ".txt", sep=""), row.names=F, na="", quote=F, sep="\t")
	})
}


facets_to_absolute <- function(facets_dir, absolute_segments_dir, chrom_col="chrom", start_col="start", end_col="end", nummark_col="num.mark", segmean_col="cnlr.median") {

	if (!file.exists(absolute_segments_dir)) { dir.create(absolute_segments_dir) }
	files <- dir(facets_dir, pattern="cncf.txt", full=T)
	lapply(files, function(facets_file) {
		print(facets_file)
		tab <- read.delim(facets_file, as.is=T)
		tab$ID <- gsub(".cncf.txt", "", basename(facets_file), fixed=T)
		tab <- tab[,c("ID", chrom_col, start_col, end_col, nummark_col, segmean_col)]
		colnames(tab) <- c("ID", "Chromosome", "Start", "End", "Num_Probes", "Segment_Mean")
		write.table(tab, file=paste(absolute_segments_dir, basename(facets_file), sep="/"), sep="\t", row.names=F, quote=F, na="")
	})
}

make_param_template <- function(samples, maf_dir="inputmafs", facets_dir="segments", params_file="absolute_params.txt") {
	library(gdata)

	facets_file <- dir(facets_dir, "cncf.txt", full=T)
	oo <- match(samples, unlist(lapply(facets_file, function(x) { strsplit(basename(x), split="_")[[1]][1]})))
	facets_file=facets_file[oo]

	params <- data.frame(include="Y", patient=samples, sample=samples, seg.dat.fn= facets_file, sigma.p=0.05, max.sigma.h=0.5, 
		min.ploidy=1.2, max.ploidy=6, primary.disease="BRCA", platform="Illumina_WES", sample.name=samples, 
		results.dir=paste("absolute/", samples, sep=""), max.as.seg.count=800, copy_num_type="total", max.neg.genome=0.05, max.non.clonal=0.5, 
		maf.fn=paste("inputmafs/", samples, ".txt", sep=""), min.mut.af=0, output.fn.base=samples)

	write.table(params, params_file, row.names=F, sep="\t", quote=F, na="")
}


read_params <- function(file) {
	if (length(grep("txt", file))==1) { 
		params <- read.delim(file, sep="\t", as.is=T)
	} else if (length(grep("xls", file))==1) {
		library(gdata)
		params <- read.xls(file, 1, stringsAsFactors=F)
	}
	params[which(params$include=="Y"), , drop=F]
}

run_absolute_step1 <- function(params, useSnow=T, numCores=3) {

	if (useSnow) {
		library(snow)
		cl <- makeCluster(numCores, "SOCK")
		parApply(cl, params, 1, function(x) {
			library(ABSOLUTE)
			RunAbsolute(seg.dat.fn=x[4],
				sigma.p=as.numeric(x[5]), max.sigma.h=as.numeric(x[6]),
				min.ploidy=as.numeric(x[7]), max.ploidy=as.numeric(x[8]), primary.disease=x[9],
				platform=x[10], sample.name=x[11],
				results.dir=x[12], max.as.seg.count=as.numeric(x[13]), copy_num_type=x[14],
				max.neg.genome=as.numeric(x[15]), max.non.clonal=as.numeric(x[16]), 
				maf.fn=x[17], min.mut.af=as.numeric(x[18]), output.fn.base=x[19], 
				verbose=TRUE)
		})
		stopCluster(cl)
	} else {
		library(ABSOLUTE)
		apply(params, 1, function(x) {
			RunAbsolute(seg.dat.fn=x[4],
			sigma.p=as.numeric(x[5]), max.sigma.h=as.numeric(x[6]),
			min.ploidy=as.numeric(x[7]), max.ploidy=as.numeric(x[8]), primary.disease=x[9],
			platform=x[10], sample.name=x[11],
			results.dir=x[12], max.as.seg.count=as.numeric(x[13]), copy_num_type=x[14],
			max.neg.genome=as.numeric(x[15]), max.non.clonal=as.numeric(x[16]), 
			maf.fn=x[17], min.mut.af=as.numeric(x[18]), output.fn.base=x[19], 
			verbose=TRUE)
		})
	}

}

run_absolute_step1_CNAonly <- function(params, useSnow=T) {

	if (useSnow) {
		library(snow)
		cl <- makeCluster(3, "SOCK")
		parApply(cl, params, 1, function(x) {
			library(ABSOLUTE)
			RunAbsolute(seg.dat.fn=x[4],
				sigma.p=as.numeric(x[5]), max.sigma.h=as.numeric(x[6]),
				min.ploidy=as.numeric(x[7]), max.ploidy=as.numeric(x[8]), primary.disease=x[9],
				platform=x[10], sample.name=x[11],
				results.dir=x[12], max.as.seg.count=as.numeric(x[13]), copy_num_type=x[14],
				max.neg.genome=as.numeric(x[15]), max.non.clonal=as.numeric(x[16]), 
				output.fn.base=x[17], 
				verbose=TRUE)
		})
		stopCluster(cl)
	} else {
		library(ABSOLUTE)
		apply(params, 1, function(x) {
			RunAbsolute(seg.dat.fn=x[4],
			sigma.p=as.numeric(x[5]), max.sigma.h=as.numeric(x[6]),
			min.ploidy=as.numeric(x[7]), max.ploidy=as.numeric(x[8]), primary.disease=x[9],
			platform=x[10], sample.name=x[11],
			results.dir=x[12], max.as.seg.count=as.numeric(x[13]), copy_num_type=x[14],
			max.neg.genome=as.numeric(x[15]), max.non.clonal=as.numeric(x[16]), 
			output.fn.base=x[17], 
			verbose=TRUE)
		})
	}

}

run_absolute_step2 <- function(obj.name, step1_results_dir="absolute", results.dir=file.path("absolute_summary", "summary")) {
	if (file.exists(results.dir)) { file.remove(results.dir) }
	absolute.files <- dir(step1_results_dir, pattern=".ABSOLUTE.RData", full=T, recursive=T)
	print (absolute.files)
	print (results.dir)
	

	if (file.exists(paste(results.dir, "/", obj.name, ".PP-calls_tab.txt", sep=""))) {
		i=1
		while (file.exists(paste(results.dir, "/", obj.name, ".PP-calls_tab", i, ".txt", sep=""))) {
			i=i+1
		}
		file.rename(paste(results.dir, "/", obj.name, ".PP-calls_tab.txt", sep=""), paste(results.dir, "/", obj.name, ".PP-calls_tab", i, ".txt", sep=""))
	}
			
	library(ABSOLUTE)
	CreateReviewObject(obj.name, absolute.files, results.dir, "total", verbose=TRUE)
}

run_absolute_step3 <- function(obj.name, results.dir=file.path("absolute_summary", "summary")) {
	library(ABSOLUTE)
	ExtractReviewedResults(reviewed.pp.calls.fn=paste(results.dir, "/", obj.name, ".PP-calls_tab.txt", sep=""), 
	analyst.id="CN", modes.fn=paste(results.dir, "/", obj.name, ".PP-modes.data.RData", sep=""), out.dir.base="reviewed", 
	obj.name=obj.name, copy_num_type="total")
}


merge_absolute_to_mutations <- function(muttab, absolute_dir="absolute/reviewed/SEG_MAF", outFile=NULL, sample_col="TUMOR_SAMPLE", chrom_col="CHROM", pos_col="POS", ref_col="REF", alt_col="ALT") {
	abs_files <- dir(absolute_dir, pattern="ABS_MAF.txt", full=T)
	cat ("Found ", length(abs_files), " files\n")
	abs <- do.call("rbind", lapply(abs_files, read.delim, as.is=T))
	abs$id <- paste(abs$sample, abs$Chromosome, abs$Start_position, abs$Reference_Allele, abs$Tumor_Seq_Allele2, sep="_")

	muttab$id <- paste(muttab[,sample_col], muttab[,chrom_col], muttab[,pos_col], muttab[,ref_col], muttab[,alt_col], sep="_")

	abs <- abs[match(muttab$id, abs$id),]
	muttab <- cbind(muttab, abs[,c("Pr_somatic_clonal", "Pr_subclonal", "clonal.ix", "subclonal.ix", "cell_mult", "cancer_cell_frac", "ccf_CI95_low", "ccf_CI95_high")])
	write.table(muttab, file=outFile, sep="\t", row.names=F, na="", quote=F)
}
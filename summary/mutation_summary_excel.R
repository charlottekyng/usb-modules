options(java.parameters = "-Xmx8000m")
suppressPackageStartupMessages(library("optparse"));
suppressPackageStartupMessages(library("xlsx"));

optList <- list(
        make_option("--outFile", default = NULL, help = "output excel file")
        )
parser <- OptionParser(usage = "%prog vcf.file", option_list = optList);
arguments <- parse_args(parser, positional_arguments = T);
opt <- arguments$options;

if (is.null(opt$outFile)) {
    cat("Need output file\n");
    print_help(parser);
    stop();
} else if (length(arguments$args) <= 1) {
    cat("Need vcf files\n");
    print_help(parser);
    stop();
}

files <- arguments$args;

output_fields = c("TUMOR_SAMPLE", "NORMAL_SAMPLE", "ANN[*].GENE", "ANN[*].HGVS_P", "ANN[*].HGVS_C", "ANN[*].EFFECT", 
	"ANN[*].IMPACT", "ANN[*].BIOTYPE", "ANN[*].FEATURE", "ANN[*].FEATUREID", "TUMOR.FA", "NORMAL.FA", "TUMOR.AF", "NORMAL.AF",
	"TUMOR.DP", "NORMAL.DP", "TUMOR.AD", "NORMAL.AD", "TUMOR.FAO", "NORMAL.FAO", "TUMOR.FDP", "NORMAL.FDP",
	"cancer_gene_census", "kandoth", "lawrence", "hap_insuf", 
	"facetsCF_EM", "facetsTCN_EM", "facetsLCN_EM",
	"duplicatedGenesDB", "dbNSFP_MutationTaster_pred", "dbNSFP_Polyphen2_HVAR_pred", "dbNSFP_Interpro_domain", 
	"CHROM", "POS", "ID", "REF", "ALT", "FILTER", "dbNSFP_ExAC_AC", "dbNSFP_ExAC_AF", "dbNSFP_ExAC_Adj_AF", "dbNSFP_Uniprot_acc")

for (file in files) {
	tab <- read.delim(file, as.is=T, check.names=F)
	print("These fields are not present in the input file - removing!")
	print(output_fields[which(!output_fields %in% colnames(tab))])

	output_fields2 <- output_fields[which(output_fields %in% colnames(tab))]

	tab <- tab[,output_fields2]

	impact <- tab[,"ANN[*].IMPACT"]
	impact_index <- lapply(impact, function(x) {
		xx <- unlist(strsplit(x, "|", fixed=T))
		if ("HIGH" %in% xx) { which(xx=="HIGH") 
		} else if ("MODERATE" %in% xx) { which(xx=="MODERATE")
		} else if ("LOW" %in% xx) { which(xx=="LOW") 
		} else if ("MODIFIER" %in% xx) {which (xx=="MODIFIER") }
	})

	aa <- tab[,"ANN[*].HGVS_P"]
	aa_index <- lapply(aa, function(x) {
		xx <- unlist(strsplit(x, "|", fixed=T))
		grep ("p\\.", xx)
	})

	for (i in grep("ANN[*]", colnames(tab), fixed=T)) {
		for (j in 1:nrow(tab)) {
			val <- tab[j,i]
			val <- unlist(strsplit(val, "|", fixed=T))
			if (length(intersect(impact_index[[j]], aa_index[[j]]))>0) {
				select <- intersect(impact_index[[j]], aa_index[[j]])
			} else {
				select <- impact_index[[j]]
			}
			tab[j,i] <- toString(val[select])
		}
	}	

	name <- strsplit (file, ".", fixed=T)[[1]]
	name <- toString(name[c(2,length(name)-1)])
	if (file.exists(opt$outFile)) {
		write.xlsx2(tab, opt$outFile, sheetName=name, append=TRUE, showNA=FALSE, row.names=F)
	} else { write.xlsx2(tab, opt$outFile, sheetName=name, append=FALSE, showNA=FALSE, row.names=F)}
}

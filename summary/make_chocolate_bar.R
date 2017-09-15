#### Usage:

### Currently this figure can handle 3 types of data, to be represented by different shades of blue, boxes around each cell, and diagonal lines
# e.g. color = CCF/VAF, boxes = clonal mutation, diagonals = LOH

# step 1: turn a table with mutations into input for the next step
# dat=make_CCF_matrix(mut)  
# where mut is a mutation table to input for figure
# use the various *col_name arguments to specify the appropriate columns.
# Z_col_name: a column with numeric data (eg. CCF, VAF)
# diags_col_name and boxes_col_name: columns with TRUE/FALSE, 1/0
# Z is mandatory, diags and boxes are optional
# Use sample_order if the samples are to be shown in specific order. Otherwise samples will be in the order of the input file
# gene_order not yet implemented
# output: dat is a list with names "Z", "diags" and "boxes" and this goes into step 2

# step 2: make the skeleton chocolate bar figure
# make_chocolate_bar(dat, filename)
# (optional) anno: list of annotation for each of the mutations (in the order of dat$Z) with two elements "pathogenic" and "membership", these should be numeric
# sort_genes will sort the genes in decreasing order of CCF/MAF or whatever the colours represent in decreasing order. Shared mutations go first.
# the rest of the options are mostly colours etc

# step 3 (optional): make a color key for the chocolate bar figure
# make_chocolate_bar_legend()
# make the legend, duh!



make_chocolate_bar <- function (dat, filename, height=NULL, width=NULL, 
	sort_genes=T, write_gene_names="all", z_breaks=c(-1e-10, c(0, 0.01, 0.05, 0.2, 0.4, 0.6, 0.8, 1)+1e-10), 
	z_cols=c("grey90", brewer.pal(9, "Blues")[3:9]), zlim=c(0,1), box.col="orange", diags_col="grey90", 
	patho_col=c(NA, "red", "darkorchid4", "grey"), mem_col=c(NA, "orange"), return_dat=F, rotate_fig=F, main="",
	mar = if (!rotate_fig) c(18,10,2,3) else c(3,18,10,2)) {

	library(RColorBrewer)
	if(sort_genes) {

		if (!is.null(dat$diags)){
			oo <- do.call(order, transform(dat$diags))
			dat <- lapply(dat, function(x) { x[oo,,drop=F] })
		}

		if (!is.null(dat$boxes)){
			oo <- do.call(order, transform(dat$boxes))
			dat <- lapply(dat, function(x) { x[oo,,drop=F] })
		}

		oo <- do.call(order, transform(dat$Z))
		dat <- lapply(dat, function(x) { x[oo,,drop=F] })

		oo <- rev(do.call(order, transform(dat$Z > 0)))
		dat <- lapply(dat, function(x) { x[oo,,drop=F] })

		oo <- order(apply(dat$Z, 1, function(y){length(which(y>0))}), decreasing=T)
		dat <- lapply(dat, function(x) { x[oo,,drop=F] })
	}
	if (!is.null(dat$anno)) { 
		anno <- data.frame(dat$anno)
		colnames(anno) <- colnames(dat$anno)
		dat <- dat[-which(names(dat)=="anno")]
	} else {anno <- NULL}

	if (!rotate_fig) {
		dat <- lapply(dat, function(x) x[,ncol(x):1, drop=F])
		# define the coordinates for grids/diagonal lines/boxes etc.
		if (nrow(dat$Z)>1) { mut_midpoint <- seq(0, 1, 1/(nrow(dat$Z)-1))} else { mut_midpoint <- 0.5 }
		if (ncol(dat$Z)>1) { sample_midpoint <- seq(0, 1, length=ncol(dat$Z)) } else { sample_midpoint <- 0.5 }
		if (length(mut_midpoint )>1) {
			vert_grid <- unlist(lapply(1:(length(mut_midpoint)-1), function(x, mut_midpoint) { 
				mean(mut_midpoint[x:(x+1)])}, mut_midpoint))
			vert_grid <- c(0-vert_grid[1], vert_grid, 1+ (1-vert_grid[length(vert_grid)]))
		} else { vert_grid <- c(-1,1) }

		if (length(sample_midpoint)>1) {
			horiz_grid <- unlist(lapply(1:(length(sample_midpoint)-1), function(x, sample_midpoint) { 
				mean(sample_midpoint[x:(x+1)])}, sample_midpoint))
			horiz_grid <- c(0-horiz_grid[1], horiz_grid, 1+(1-horiz_grid[length(horiz_grid)]))
		} else { horiz_grid <- c(-1,1) }

		if(is.null(height)) { height=(ncol(dat$Z)*0.6)+4 }
		if(is.null(width)) { width=(nrow(dat$Z)/5)+2 }
	} else if(rotate_fig) {
		dat <- lapply(dat, t)
		dat <- lapply(dat, function(x) x[,ncol(x):1, drop=F])
		if(!is.null(anno)){ anno <- anno[nrow(anno):1,,drop=F] }
		# define the coordinates for grids/diagonal lines/boxes etc.
		if (ncol(dat$Z)>1) { mut_midpoint <- seq(0, 1, 1/(ncol(dat$Z)-1))} else { mut_midpoint <- 0.5 }
		if (nrow(dat$Z)>1) { sample_midpoint <- seq(0, 1, length=nrow(dat$Z)) } else { sample_midpoint <- 0.5 }
		if (length(mut_midpoint )>1) {
			horiz_grid <- unlist(lapply(1:(length(mut_midpoint)-1), function(x, mut_midpoint) { 
				mean(mut_midpoint[x:(x+1)])}, mut_midpoint))
			horiz_grid <- c(0-horiz_grid[1], horiz_grid, 1+ (1-horiz_grid[length(horiz_grid)]))
		} else { horiz_grid <- c(-1,1) }

		if (length(sample_midpoint)>1) {
			vert_grid <- unlist(lapply(1:(length(sample_midpoint)-1), function(x, sample_midpoint) { 
				mean(sample_midpoint[x:(x+1)])}, sample_midpoint))
			vert_grid <- c(0-vert_grid[1], vert_grid, 1+(1-vert_grid[length(vert_grid)]))
		} else { vert_grid <- c(-1,1) }

		if(is.null(height)) { height=max((ncol(dat$Z)/5)+2,4) }
		if(is.null(width)) { width=max((nrow(dat$Z)*0.6)+4,4) }
	}
	if(!is.null(filename)){ 
		pdf(filename, height=height, width=width)
	}
#	if(!rotate_fig) {
#		par(mar=c(18,10,2,3), xpd=T)
#	} else {
#		par(mar=c(3,18,10,2), xpd=T)
#	}
	par(mar=mar, xpd=T)
	if (!rotate_fig) {
		# draw the heatmap
		image(as.matrix(dat$Z), ylab="", xaxt='n', yaxt='n', xlab="", bty='n', col= z_cols, zlim=zlim, breaks= z_breaks, main=main)
		axis(2, at=sample_midpoint, labels=colnames(dat$Z), las=2, tick=F, cex.axis = 1.2)
		if (write_gene_names=="all") { 
			axis(1, at= mut_midpoint, labels=rownames(dat$Z), las=2, tick=F, cex.axis= 1, line=1.8)
		} else if (write_gene_names=="annotated" & !is.null(anno)) {
			labels <- rownames(dat$Z)
			labels[which(apply(anno,1,function(x){all(x==1)}))] <- ""
			axis(1, at= mut_midpoint, labels=labels, las=2, tick=F, cex.axis= 1, line=1.8)
		}	


	} else if(rotate_fig) {
		# draw the heatmap
		image(as.matrix(dat$Z), ylab="", xaxt='n', yaxt='n', xlab="", bty='n', col= z_cols, zlim=zlim, breaks= z_breaks, main=main)
		axis(3, at=sample_midpoint, labels=rownames(dat$Z), las=2, tick=F, cex.axis = 1.2)
		if (write_gene_names=="all") { 
			axis(2, at= mut_midpoint, labels=colnames(dat$Z), las=2, tick=F, cex.axis= 1, line=1.8)
		} else if (write_gene_names=="annotated" & !is.null(anno)) {
			labels <- colnames(dat$Z)
			labels[which(apply(anno,1,function(x){all(x==1)}))] <- ""
			axis(2, at= mut_midpoint, labels=labels, las=2, tick=F, cex.axis= 1, line=1.8)
		}
	}
	# draw the grids
	if (!rotate_fig) {
		for (i in horiz_grid) { lines(c(min(vert_grid), max(vert_grid)), c(i,i), col="white", lwd=4)}
		for (i in vert_grid) { lines(c(i,i), c(min(horiz_grid), max(horiz_grid)), col="white", lwd=2)}
	}  else if(rotate_fig) {
		for (i in horiz_grid) { lines(c(min(vert_grid), max(vert_grid)), c(i,i), col="white", lwd=2)}
		for (i in vert_grid) { lines(c(i,i), c(min(horiz_grid), max(horiz_grid)), col="white", lwd=4)}
	}
	# draw the diagonal lines, if defined
	if (!is.null(dat$diags)) {	
		for (i in 1:ncol(dat$diags)) {
			for (j in 1:nrow(dat$diags)) {
				if (dat$diags[j,i]>0) { 
					lines(vert_grid[j:(j+1)], horiz_grid[i:(i+1)], col= diags_col, lwd=3)
				}
			}
		}
	}

	# draw the colored boxes around, if defined
	roffset <- 0.06*(1/nrow(dat$Z))
	coffset <- 0.06*(1/ncol(dat$Z))

	if (!is.null(dat$boxes)){
		for (i in 1:ncol(dat$boxes)) {
			for (j in 1:nrow(dat$boxes)) {
				if (dat$boxes[j,i]>0) { 
					polygon(c(vert_grid[j]+roffset, vert_grid[j+1]-roffset, vert_grid[j+1]-roffset, vert_grid[j]+roffset),
						c(horiz_grid[i]+coffset, horiz_grid[i]+coffset, horiz_grid[i+1]-coffset, horiz_grid[i+1]-coffset),
						border=box.col, lwd=2.2) 

				}
			}
		}
	}

	# draw the dots for the cancer genes, if defined
	if (!is.null(anno)) {
		if (!is.null(anno$pathogenicity)) { 
			if (!rotate_fig) {
				points(mut_midpoint, rep(horiz_grid[1]-(0.35*1/ncol(dat$Z)), length(mut_midpoint)), 
					col= patho_col[as.numeric(anno$pathogenicity)], pch=16, cex=2)
			} else if(rotate_fig) {
				points(rep(vert_grid[1]-(0.3*1/nrow(dat$Z)), length(mut_midpoint)), mut_midpoint, 
					col= patho_col[as.numeric(anno$pathogenicity)], pch=16, cex=2)
			}
		}
		if (!is.null(anno$membership)) {
			if (!rotate_fig) {
				points(mut_midpoint, rep(horiz_grid[1]- (0.8*1/ncol(dat$Z)), length(mut_midpoint)), 
					col= mem_col[as.numeric(anno$membership)], pch=16, cex=2)
			} else if(rotate_fig) {
				points(rep(vert_grid[1]- (0.9*1/nrow(dat$Z)), length(mut_midpoint)), mut_midpoint, 
					col= mem_col[as.numeric(anno$membership)], pch=16, cex=2)
			}
		}
	}

	if(!is.null(filename)){
		dev.off()
	}
	if(return_dat){return(dat)}

}

make_chocolate_bar_legend <- function(legend=c("0%", ">0%-1%", ">1%-5%", ">5%-20%", ">20%-40%", ">40%-60%", ">60%-80%", ">80%-100%"), 
	cols= c("grey90", brewer.pal(9, "Blues")[3:9]), file="CCF_colourkey.pdf") {
	print(length(legend))
	print(length(cols))
	pdf(file, height=10, width=4)
	plot(0,0, type='n', bty='n', xaxt='n', yaxt='n', xlab="", ylab="", xlim=c(-1,1), ylim=c(-2.2,1))
	leg <- cbind(legend,cols)
	legend(-1,1, legend=leg[,1], fill=leg[,2], bty='n', cex=0.8)
	dev.off()
}


make_CCF_matrix <- function(muts, sample_col_name="TUMOR_SAMPLE", mutid_col_name=NULL, 
	chr_col_name="CHROM", pos_col_name="POS", ref_col_name="REF", alt_col_name="ALT", 
	gene_col_name="GENE", AA_col_name="HGVS_P", CDNA_col_name="HGVS_C", 
	Z_col_name="ccf", diags_col_name="facetsLOHCall", boxes_col_name="clonalStatus", include_X=T, 
	anno_cols=NULL,sample_order=NULL, gene_order=NULL) {

	if(!sample_col_name %in% colnames(muts)) {stop("sample_col_name not a column name")}
	if (!is.null(mutid_col_name)){
		if(mutid_col_name %in% colnames(muts)) {
		cat ("Going to ignore chr_col_name, pos_col_name, ref_col_name and alt_col_name\n")
		muts$id <- muts[,mutid_col_name]
	}} else {
		if(!chr_col_name %in% colnames(muts)) {stop("chr_col_name not a column name")}
		if(!pos_col_name %in% colnames(muts)) {stop("pos_col_name not a column name")}
		if(!ref_col_name %in% colnames(muts)) {stop("ref_col_name not a column name")}
		if(!alt_col_name %in% colnames(muts)) {stop("alt_col_name not a column name")}
		if(any(is.na(muts[,chr_col_name]))|any(is.na(muts[,pos_col_name]))|any(is.na(muts[,ref_col_name]))|any(is.na(muts[,alt_col_name]))){
			stop("Something in chr_col_name, pos_col_name, ref_col_name, alt_col_name is NA - not allowed! Fix it!")
		}
		if(!include_X){ muts <- muts[which(!muts[,chr_col_name] %in% c("X", 23)),]}
		muts$id <- unlist(apply(muts[,c(chr_col_name, pos_col_name, ref_col_name, alt_col_name)], 1, toString))
	}
	if (!is.null(CDNA_col_name)){
		if (CDNA_col_name %in% colnames(muts)){
		if("." %in% muts[,AA_col_name]){
			cat ("Using cDNA position for non-coding positions\n")
			muts[which(muts[,AA_col_name]=="."),AA_col_name] <- muts[which(muts[,AA_col_name]=="."),CDNA_col_name]
		}}
	}
	annot <- unique(cbind(muts$id, paste(muts[, gene_col_name], " (", muts[, AA_col_name], ")", sep="")))
	if(!is.null(Z_col_name)){
		if( Z_col_name %in% colnames(muts)){
		if(length(grep("%", muts[,Z_col_name]))>0) {
			cat ("Converting Z_col_name from % to decimals\n")
			muts[,Z_col_name] <- as.numeric(gsub("%", "", muts[,Z_col_name]))/100
		}}
		muts[,Z_col_name] <- as.numeric(muts[,Z_col_name])
	}
	if(!is.null(diags_col_name)) {
		if( diags_col_name %in% colnames(muts)){
		if(!all(muts[,diags_col_name] %in% c(0,1))|!all(muts[,diags_col_name] %in% c(T,F))){
			cat ("Converting diags_col_name from % to 0s and 1s\n")
			muts[which(muts[,diags_col_name]==T),diags_col_name] <- 1
			muts[which(muts[,diags_col_name]=="true"),diags_col_name] <- 1
			muts[which(muts[,diags_col_name]==F),diags_col_name] <- 0
			muts[which(muts[,diags_col_name]=="false"),diags_col_name] <- 0
			if (any(muts[,diags_col_name]==".")) {
				cat ("There are .s in diags_col_name. Converting to 0.\n")
				muts[which(muts[,diags_col_name]=="."),diags_col_name] <- 0

			}
			if (any(is.na(muts[,diags_col_name]))) {
				cat ("There are NAs in diags_col_name. Converting to 0.\n")
				muts[which(is.na(muts[,diags_col_name])),diags_col_name] <- 0

			}

		}}
	}
	if(!is.null(boxes_col_name)){
		if( boxes_col_name %in% colnames(muts)){		
		if(!all(muts[,boxes_col_name] %in% c(0,1))|!all(muts[,boxes_col_name] %in% c(T,F))){
			cat ("Converting boxes_col_name from % to 0s and 1s\n")
			muts[which(muts[,boxes_col_name]=="clonal"),boxes_col_name] <- 1
			muts[which(muts[,boxes_col_name]!="clonal"),boxes_col_name] <- 0
		}}
	}






	if (is.null(sample_order)) { sample_order=rev(unique(muts[,sample_col_name]))}

	mat <- matrix(0, nrow=length(unique(muts$id)), ncol=length(sample_order))
	rownames(mat) <- unique(muts$id)
	colnames(mat) <- sample_order
	anno <- matrix(0, ncol=length(anno_cols), nrow=nrow(mat))
	rownames(anno) <- rownames(mat)
	colnames(anno) <- anno_cols

	res <- list(Z=mat, diags=mat, boxes=mat, anno=anno)
	names(res) <- c("Z", "diags", "boxes", "anno")
	
	for (i in 1:nrow(muts)) {
		if (!is.null(Z_col_name)) { res[["Z"]][muts$id[i], muts[i, sample_col_name]] <- muts[i, Z_col_name] }
		if (!is.null(diags_col_name)) { res[["diags"]][muts$id[i], muts[i, sample_col_name]] <- muts[i, diags_col_name] }
		if (!is.null(boxes_col_name)) { res[["boxes"]][muts$id[i], muts[i, sample_col_name]] <- muts[i, boxes_col_name] }
		if (!is.null(anno_cols)) { 
			for (j in anno_cols) {
				res[["anno"]][muts$id[i],j] <- muts[i,j]
			}
		}

		
	}
	res <- lapply(res, function(x) {
		rownames(x) <- annot[match(rownames(x), annot[,1]),2]
		x
	})
	if (is.null(Z_col_name)) { res[c("Z")]<- NULL }
	if (is.null(diags_col_name)) { res[c("diags")] <- NULL }
	if (is.null(boxes_col_name)) { res[c("boxes")] <- NULL }
	if (is.null(anno_cols)) { res[c("anno")] <- NULL }

	res
	

}



## Script to Pull out Variants Near Candidate SNPs ##
## December 1, 2015 ##
## Kristopher Standish ##

#############################################
## LOAD DATA ################################
#############################################

## TSCC Paths
DATE <- gsub("-","",Sys.Date())
PathToPrev <- "/projects/janssen/ASSOCIATION/20151201_Previous_CandGWAS/Previous_GWAS_rsIDs.Uniq.txt"
PathToP <- "/projects/janssen/ASSOCIATION/20150104_LT8_DEL_MNe_MN_DAS_BL_MN_PC1_PC2/20150104_LT8_DEL_MNe_MN_DAS_BL_MN_PC1_PC2.P"
PathToGenes <- "/home/kstandis/HandyStuff/GG-Gene_Names_DB.txt"
PathToProblem <- gsub("GG-Gene_Names_DB.txt","Problem_Genes.Filt.txt", PathToGenes )
PathToRaw <- "/projects/janssen/ASSOCIATION/20151201_Previous_CandGWAS/Nearby_SNPs.tpA.traw"
PathToRaw <- "/projects/janssen/ASSOCIATION/20151201_Previous_CandGWAS/20151202_CandidateGWAS/TAB-NearbySNPs.traw"
PathToPlot <- paste("/projects/janssen/ASSOCIATION/20151201_Previous_CandGWAS/",DATE,"_CandidateGWAS/",sep="")
if ( !file.exists(PathToPlot) ) { dir.create( PathToPlot ) }

## Mac Paths
DATE <- gsub("-","",Sys.Date())
PathToPrev <- "/Users/kstandis/Dropbox/Schork/JNJ11/Manuscripts/GWAS/Previous_GWAS_rsIDs.txt"
PathToP <- "/Users/kstandis/Downloads/20151026_Enrichment/20151016_Th1_Perm_Pheno.1_DAS_BL_MN_PC1_PC2.P"
PathToGenes <- "/Users/kstandis/Data/Genetics/Gene_Names/GG-Gene_Names_DB.txt"
PathToProblem <- "/Users/kstandis/Data/Genetics/Gene_Names/Problem_Genes.Filt.txt"
PathToPlot <- paste("/Users/kstandis/Data/Burn/Plots/",DATE,"_CandidateGWAS/",sep="")
if ( !file.exists(PathToPlot) ) { dir.create( PathToPlot ) }

## Load Data
CANDS.l <- read.table( PathToPrev, colClasses="character" )[,1]
P.l <- read.table( PathToP, sep="\t",header=T,colClasses=c("factor","character","numeric","numeric") )
RAW.l <- read.table( PathToRaw, header=T )
G.l <- read.table( PathToGenes, sep="\t",header=T, comment.char="", fill=T, quote="" )

#############################################
## ORGANIZE/FILTER DATA #####################
#############################################

## Get rid of Candidate SNPs that aren't in P file
CANDS.prev <- CANDS.l[ which(CANDS.l %in% P.l$SNP) ]

## Filter out P-Values that are "NA"
P <- P.l[ which(!is.na(P.l$P)), ]

## Pull out Genotype Info from loaded RAW file
RAW.i <- RAW.l[ ,1:6 ]
RAW <- RAW.l[ ,-(1:6) ]

## Rename/organize G Table a little
colnames(G.l) <- gsub("hg19.","", gsub("knownGene.","", gsub("kgXref.","", colnames(G.l) ) ) )
colnames(G.l)[1] <- "tx"
G.l$gtx <- paste( G.l$geneSymbol, G.l$tx, sep="_" )

## Filter G based on Chromosomes
 # Specify Chromosome Names
chroms <- c(1:22,"X","Y")
RM.chrom <- which( !(G.l$chrom %in% paste("chr",chroms,sep="")) )
if ( length(RM.chrom)>0 ) { G.1 <- G.l[ -RM.chrom, ] }else{ G.1 <- G.l }

## Problem Genes that appear on multiple Chromosomes
ALL_GENES <- as.character(unique( G.1$geneSymbol ))
if ( file.exists(PathToProblem) ) {
	GENES.prob <- as.character( read.table( PathToProblem,header=F )[,1] )
}else{
	start <- proc.time()
	GENES.prob.1 <- unlist(lapply( ALL_GENES, function(x) length(unique( G.1$chrom[which(G.1$geneSymbol==x)] )) ))
	round((proc.time()-start)[3],3)
	GENES.prob <- ALL_GENES[ which(GENES.prob.1>1) ]
}
 # Remove Bad Genes
RM.genes <- which( G.1$geneSymbol %in% GENES.prob )
if ( length(RM.genes)>0 ) { G.2 <- G.1[ -RM.genes, ] }else{ G.2 <- G.1 }
G <- G.2

#############################################
## CREATE FUNCTIONS #########################
#############################################

## FCT: Pull out SNPs near a Candidate SNP
Nearby <- function( P, Cand, Buffer ) {
	Which <- which(P$SNP==Cand)
	Cand_Chrom <- as.character( P$CHR[Which] )
	Cand_Pos <- as.numeric( P$BP[Which] )
	Pos_Rng <- Cand_Pos + c(-Buffer,Buffer)
	Which_Rows <- which( P$CHR==Cand_Chrom & P$BP>Pos_Rng[1]  & P$BP<Pos_Rng[2] )
	Nearby_SNPs <- P$SNP[ Which_Rows ]
	return(Nearby_SNPs)
}

## FCT: Pull out Genes near a Candidate SNP
Nearby_Gene <- function( P, Cand, Buffer ) {
	Which <- which(P$SNP==Cand)
	Cand_Chrom <- as.character( P$CHR[Which] )
	Cand_Pos <- as.numeric( P$BP[Which] )
		## Nearby Genes
	Which_Rows.G <- which( G$chrom==paste("chr",Cand_Chrom,sep="") & G$cdsStart<Cand_Pos+Buffer & G$cdsEnd>Cand_Pos-Buffer )
	G.which <- G[ Which_Rows.G, ]
	return(as.character(G.which$geneSymbol))
}

## FCT: Plot Gene Coordinates
Gene_Line <- function( G_Line, Rand_N ) {
	# print(G_Line)
	# print(G_Line[1,"cdsStart"])
	arrows(as.numeric(G_Line[1,"txStart"]),Rand_N,as.numeric(G_Line[1,"txEnd"]),Rand_N, angle=90,code=3,length=.05, col="steelblue4", lwd=2 )
	arrows(as.numeric(G_Line[1,"cdsStart"]),Rand_N,as.numeric(G_Line[1,"cdsEnd"]),Rand_N, angle=90,code=3,length=.05, col="steelblue1", lwd=2 )
	EX.starts <- as.numeric(strsplit(as.character(G_Line[1,"exonStarts"]),",")[[1]])
	EX.ends <- as.numeric(strsplit(as.character(G_Line[1,"exonEnds"]),",")[[1]])
	arrows( EX.starts, Rand_N, EX.ends, Rand_N, length=0, col="orange2", lwd=5 )
	LABEL <- paste(as.character(G_Line[1,"geneSymbol"]),as.character(G_Line[1,"tx"]),sep="_")
	text(as.numeric(G_Line[1,"txStart"]),Rand_N, label=LABEL, pos=2, col="steelblue3",cex=.4 )
}

## FCT: Manhattan Plot for SNPs near Candidate
Manhat <- function( P, Cand ) {
	## Candidate Info
	Which <- which(P$SNP==Cand)
	Cand_Chrom <- as.character( P$CHR[Which] )
	Cand_Pos <- as.numeric( P$BP[Which] )
	Buffer <- 1e5
	Pos_Rng <- Cand_Pos + c(-Buffer,Buffer)
	## Nearby SNP table info
	Nearby_SNPs <- Nearby.list[[paste("SNP",gsub(":",".",Cand),sep="_")]]
	Which_Rows <- which( P$SNP %in% Nearby_SNPs )
	P.which <- P[ Which_Rows, ]
	## Nearby Genes
	Which_Rows.G <- which( G$chrom==paste("chr",Cand_Chrom,sep="") & G$cdsStart>Pos_Rng[1]-1e4  & G$cdsStart<Pos_Rng[2]+1e4 )
	G.which <- G[ Which_Rows.G, ]
	## Plot it
	YLIM <- c( -2, max(4,-log10(P.which$P)) )
	jpeg( paste(PathToPlot,"1-Cand_Manhat.",Cand,".jpeg",sep=""), height=1400,width=2000,pointsize=34 )
	plot( -log10(P.which$P) ~ P.which$BP, pch="+", ylim=YLIM, xlab=paste("Chromosome",Cand_Chrom,"Position"), main=paste("Manhattan Plot Surrounding",Cand), col="slateblue4" )
	abline( h=seq(-10,10,1),lty=3,col="grey50" ) # abline( h=seq(-10,10,1),lty=c(rep(3,10),1,rep(3,10)),col="grey50" )
	abline( h=0,lty=1,lwd=2,col="black" )
	points( -log10(P$P[Which]) ~ P$BP[Which], cex=2,pch="+",col="firebrick1" )
	text( -log10(P$P[Which]) ~ P$BP[Which], label=Cand,pos=2,col="firebrick1" )
	GAP <- .09
	yval <- -GAP
	for ( i in 1:length(Which_Rows.G) ) {
		Gene_Line( G.which[i,], yval )
		yval <- ifelse( yval<=-2, -GAP, yval-GAP )
		# if ( yval <= -2 ) { yval <- -GAP
		# }else{
		# 	yval <- yval - GAP
		# }
	}
	# apply( G.which, 1, function(x) Gene_Line(x,runif(1,-2,0)) )
	dev.off()
}

Manhat.link <- function( P.which, Cand, Buffer ) {
	P.which <- P.which[ order(P.which$R2,decreasing=F), ]
	## Candidate Info
	Which <- which(P.which$SNP==Cand)
	Cand_Chrom <- as.character( P.which$CHR[Which] )
	Cand_Pos <- as.numeric( P.which$BP[Which] )
	Pos_Rng <- Cand_Pos + c(-Buffer,Buffer)
	## Nearby Genes
	Which_Rows.G.1 <- which( G$chrom==paste("chr",Cand_Chrom,sep="") & G$cdsStart>Pos_Rng[1]-1e4 & G$cdsStart<Pos_Rng[2]+1e4 )
	Which_Rows.G.2 <- which( G$chrom==paste("chr",Cand_Chrom,sep="") & G$cdsStart<Cand_Pos+1e4 & G$cdsEnd>Cand_Pos-1e4 )
	Which_Rows.G <- union( Which_Rows.G.1, Which_Rows.G.2 )
	G.which <- G[ Which_Rows.G, ]
	## Plot it
	 # Plotting Parameters
	YLIM <- c( -2, max(4,-log10(P.which$P)) )
	COLS.list <- c("black","slateblue4","steelblue3","springgreen3","gold2","chocolate1","firebrick1")
	COLS <- colorRampPalette(COLS.list)(101)
	 # Plot Data
	jpeg( paste(PathToPlot,"1-Cand_Manhat.",Cand,"b",Buffer,".jpeg",sep=""), height=1400,width=2000,pointsize=40 )
	plot( -log10(P.which$P) ~ P.which$BP, pch="+", ylim=YLIM, xlab=paste("Chromosome",Cand_Chrom,"Position"), main=paste("Manhattan Plot Surrounding",Cand), col=COLS[1+round(100*P.which$R2)] )
	abline( h=seq(-10,10,1),lty=3,col="grey50" ) # abline( h=seq(-10,10,1),lty=c(rep(3,10),1,rep(3,10)),col="grey50" )
	abline( h=0,lty=1,lwd=2,col="black" )
	points( -log10(P.which$P[Which]) ~ P.which$BP[Which], cex=2,pch=10,col="firebrick1", lwd=5 )
	text( -log10(P.which$P[Which]) ~ P.which$BP[Which], label=Cand,pos=3,col="firebrick1" )
	 # Genes
	GAP <- .09
	yval <- -GAP
	for ( i in 1:length(Which_Rows.G) ) {
		Gene_Line( G.which[i,], yval )
		yval <- ifelse( yval<=-2, -GAP, yval-GAP )
	}
	dev.off()
}

## FCT: Calculate LD between Candidate and Nearby SNPs
Linkage <- function( RAW, RAW.i, Nearby_SNPs, Cand ) {
	## Candidate Info
	Which <- which(RAW.i$SNP==Cand)
	RAW.cand <- c( RAW[ Which, ], recursive=T )
	## Nearby SNP RAW info
	# Nearby_SNPs <- Nearby.list[[paste("SNP",gsub(":",".",Cand),sep="_")]]
	RAW.which <- RAW[ which( RAW.i$SNP %in% Nearby_SNPs ), ]
	P.which <- P[ which( P$SNP %in% Nearby_SNPs ), ]
	## Calculate LD for each SNP
	R2 <- abs(apply( RAW.which, 1, function(x) cor(x,RAW.cand,method="spearman") ))
	P.which.2 <- data.frame( P.which, R2 )
	# D <- (x11)(x22) – (x12)(x21)
	return(P.which.2)
}

#############################################
## FINALIZE CANDIDATE SNP LIST ##############
#############################################

## Determine Custom List of SNPs to Plot
 # Get list of SNPs p<1e-6 (or other Threshold)
P.Thresh <- 5e-6
Cand_P <- P[ which(P$P <= P.Thresh), ]
To_Test <- as.character(Cand_P$SNP)
To_Keep <- c()
for ( r in 1:nrow(Cand_P) ) {
	snp <- as.character(Cand_P$SNP[r])
	if ( snp %in% To_Test ) {
		chr <- Cand_P$CHR[r]
		bp <- Cand_P$BP[r]
		which_nearby <- which(Cand_P$CHR==chr & Cand_P$BP>=bp-1e5 & Cand_P$BP<=bp+1e5 )
		snp_to_keep <- Cand_P$SNP[which_nearby][ which.min(Cand_P$P[which_nearby]) ]
		To_Keep <- c( To_Keep, snp_to_keep )
		To_Test <- setdiff( To_Test, as.character(Cand_P$SNP[which_nearby]) )
	}else{
		next
	}
}
Cand_P.keep <- Cand_P[ which(Cand_P$SNP %in% To_Keep), ]
CANDS.p <- as.character( Cand_P.keep$SNP )
 # Get list of SNPs in Candidate Genes
Cand_Gene.List <- c( "GRIN2B","CNTN5","CYP2E1" )
Cand_Gene.Which <- which( G$geneSymbol %in% Cand_Gene.List )
Cand_Gene.Rng.1 <- G[ Cand_Gene.Which, c("geneSymbol","chrom","txStart","txEnd") ]
TEMP.chrom <- aggregate( Cand_Gene.Rng.1[,"chrom"], by=list(Gene=Cand_Gene.Rng.1[,"geneSymbol"]), function(x)head(x,1) )
TEMP.start <- aggregate( Cand_Gene.Rng.1[,"txStart"], by=list(Gene=Cand_Gene.Rng.1[,"geneSymbol"]), min )
TEMP.end <- aggregate( Cand_Gene.Rng.1[,"txEnd"], by=list(Gene=Cand_Gene.Rng.1[,"geneSymbol"]), max )
Cand_Gene.Rng <- cbind( TEMP.chrom, TEMP.start[,2], TEMP.end[,2] ) ; colnames(Cand_Gene.Rng) <- colnames(Cand_Gene.Rng.1)
Cand_Gene.SNPs <- apply( Cand_Gene.Rng, 1, function(x) P[ which( P$CHR==gsub("chr","",x["chrom"]) & P$BP>=as.numeric(x["txStart"]) & P$BP<=as.numeric(x["txEnd"]) ), c("SNP","P") ] )
names(Cand_Gene.SNPs) <- as.character( Cand_Gene.Rng[,"geneSymbol"] )
CANDS.gene <- unlist(lapply( Cand_Gene.SNPs, function(x) x$SNP[ which.min(x$P) ] ))

## Compile Lists of Candidate SNPs
CANDS.custom <- union( CANDS.p, CANDS.gene )
CANDS <- union( CANDS.prev, CANDS.custom )
P.cands <- P[ which(P$SNP %in% CANDS), ]

#############################################
## Pull/Compile Nearby SNPs #################
#############################################

## Get SNP IDs for all SNPs near Candidate SNPs
Buffer <- 1e5
Nearby.list <- lapply( CANDS, function(x) Nearby(P,x,Buffer) )
names(Nearby.list) <- paste("SNP",gsub(":",".",CANDS),sep="_")

## Compile into 1 table to be written
SNP.all <- Reduce( union, Nearby.list )
write.table( SNP.all, paste(PathToPlot,"TAB-NearbyCompile.txt",sep=""), col.names=F,row.names=F,quote=F )

## Create "RAW" File in Output Directory, then load RAW file
tRAW_File <- paste(PathToPlot,"TAB-NearbySNPs",sep="")
COMMAND <- paste( "plink --bfile /projects/janssen/VCFs/PLINK/BED_FULL.ALL --extract",paste(PathToPlot,"TAB-NearbyCompile.txt",sep=""),"--recode A-transpose --out",tRAW_File )
system( COMMAND )
 # Load/Reorganize Genotype Info from loaded RAW file
RAW.l <- read.table( paste(tRAW_File,".traw",sep=""), header=T )
RAW.i <- RAW.l[ ,1:6 ]
RAW <- RAW.l[ ,-(1:6) ]

#############################################
## Get Linkage for all Candidate SNPs #######
#############################################

## Calculate Linkage w/ Candidate SNP
P.link <- lapply( CANDS, function(x) Linkage(RAW,RAW.i,Nearby.list[[paste("SNP",gsub(":",".",x),sep="_")]],x) )
names(P.link) <- paste("SNP",gsub(":",".",CANDS),sep="_")
 # Compile into 1 array
P.link.arr <- Reduce( rbind, P.link )
write.table( P.link.arr, paste(PathToPlot,"TAB-NearbyCompile.LD.all.txt",sep=""), col.names=F,row.names=F,quote=F )

## Save List of SNPs w/ LD >= X
for ( x in seq(.5,1,.1) ) {
	TEMP <- as.character( P.link.arr$SNP[ which(P.link.arr$R2>x) ] )
	write.table( TEMP, paste(PathToPlot,"TAB-NearbyCompile.LD",x,".txt",sep=""), col.names=F,row.names=F,quote=F )
}

#############################################
## Manhattan Plot Near Candidate SNPs #######
#############################################

## Manhattan Plots around Candidate SNPs
# Manhat.link( P.link[[paste("SNP",gsub(":",".",CANDS[5]),sep="_")]],CANDS[5] )
SCRAP <- lapply( CANDS, function(x) Manhat.link(P.link[[paste("SNP",gsub(":",".",x),sep="_")]],x,Buffer) )


#############################################
## Run Candidate SNP through Pipeline #######
#############################################

## Go for it...
x <- "rs1813443"
x <- "rs1532269"
x <- CANDS[20]
Xs <- c("rs1813443","rs1532269",CANDS[c(10,14,20)])
Xs <- CANDS.prev
Xs <- CANDS.gene
for ( x in Xs ) {
	print(x)
	for ( Buffer in c(1e5,1e6) ) {
		NEAR.x <- Nearby(P,x,Buffer)
		LINK.x <- Linkage(RAW,RAW.i,NEAR.x,x)
		SCRAP.x <- Manhat.link(LINK.x,x,Buffer)
		print("Next")
	}	
}
	



#############################################
## END OF DOC ###############################
#############################################

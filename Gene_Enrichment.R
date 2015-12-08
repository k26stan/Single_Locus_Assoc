## Get Enrichment Scores from GWAS P-Values for all Genes ##
## October 26, 2015 ##
## Kristopher Standish ##

################################################
## LOAD DATA ###################################
################################################

library(SuppDists)
# library(fBasics)
LINE <- commandArgs(trailingOnly = TRUE)
# LINE <- c("/projects/janssen/ASSOCIATION/20151016_Th1_Permute/20151016_Th1_Perm_Pheno.1_DAS_BL_MN_PC1_PC2/20151016_Th1_Perm_Pheno.1_DAS_BL_MN_PC1_PC2.P","/home/kstandis/HandyStuff/GG-Gene_Names_DB.txt","/projects/janssen/ASSOCIATION/20151102_Th1_Permute_Test/")
# LINE <- c("/projects/janssen/ASSOCIATION/20150104_LT8_DEL_MNe_MN_DAS_BL_MN_PC1_PC2/20150104_LT8_DEL_MNe_MN_DAS_BL_MN_PC1_PC2.P","/home/kstandis/HandyStuff/GG-Gene_Names_DB.txt","/projects/janssen/ASSOCIATION/20151102_EnrichCompile/")
PathToPFile <- LINE[1]
PathToGenes <- LINE[2]
PathToPlot <- LINE[3]
Script_or_Gene <- LINE[4]
Path_Split <- strsplit(PathToPFile,"/")[[1]]
File_Tag <- gsub(".P","",Path_Split[length(Path_Split)],fixed=T)
PathToProblem <- gsub("GG-Gene_Names_DB.txt","Problem_Genes.Filt.txt", PathToGenes )
# PathToPlot <- paste( c(Path_Split[-length(Path_Split)],"Plots",""), collapse="/")
if ( !file.exists(PathToPlot) ) { dir.create( PathToPlot ) }

## Set Date
if ( !exists("LINE") ) {
	## Mac Paths
	DATE <- gsub("-","",Sys.Date())
	PathToPFile <- "/Users/kstandis/Downloads/20151026_Enrichment/20151016_Th1_Perm_Pheno.1_DAS_BL_MN_PC1_PC2.P"
	PathToGenes <- "/Users/kstandis/Data/Genetics/Gene_Names/GG-Gene_Names_DB.txt"
	PathToProblem <- "/Users/kstandis/Data/Genetics/Gene_Names/Problem_Genes.Filt.txt"
	PathToPlot <- paste("/Users/kstandis/Data/Burn/Plots/",DATE,"_Enrichment/",sep="")
	Path_Split <- strsplit(PathToPFile,"/")[[1]]
	File_Tag <- gsub(".P","",Path_Split[length(Path_Split)],fixed=T)
	if ( !file.exists(PathToPlot) ) { dir.create( PathToPlot ) }
	Script_or_Gene <- "Gene"
}
print(paste("Path To GWAS Results:", PathToPFile ))
print(paste("Path To Gene Table:", PathToGenes ))
print(paste("Path To Output Directory:", PathToPlot ))
print(paste("File Tag:", File_Tag ))

## Load P-Value File & Gene File
print("Loading Files")
start <- proc.time()
P <- read.table( PathToPFile, sep="\t",header=T,colClasses=c("factor","character","numeric","numeric") )
G.l <- read.table( PathToGenes, sep="\t",header=T, comment.char="", fill=T, quote="" )
colnames(G.l) <- gsub("hg19.","", gsub("knownGene.","", gsub("kgXref.","", colnames(G.l) ) ) )
colnames(G.l)[1] <- "tx"
G.l$gtx <- paste( G.l$geneSymbol, G.l$tx, sep="_" )
print(paste( "Files Loaded:",round((proc.time()-start)[3],3),"seconds" ))

################################################
## ORGANIZE/FILTER DATA ########################
################################################

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

## Gene Areas from Actual GWAS
GOI <- list()
GOI[["Full"]] <- array(,c(2,2)) ; colnames(GOI[["Full"]]) <- c("BEST_P","AREA") ; rownames(GOI[["Full"]]) <- c("GRIN2B","CNTN5")
GOI[["Ex"]] <- GOI[["Full"]]
GOI[["Full"]][,"AREA"] <- c( 8.178442, 4.0526161 )
GOI[["Full"]][,"BEST_P"] <- c( 2.346e-08, 1.650e-06 )
GOI[["Ex"]][,"AREA"] <- c( 2.139639, 0.3252284 )
GOI[["Ex"]][,"BEST_P"] <- c( 7.554e-06, 2.136e-02 )

################################################
## GET GENOME-WIDE SUMMARY STATS ###############
################################################

## Sort P-Values (Low->High)
P.sort <- sort(P$P)

## Get P-Value at Nth best Variant
Order <- 5
Ranks <- sort(union( 1:100, c(1:9) %*% 10^t(0:Order) )) # sort(c(1:9) %*% 10^t(0:Order))
P.ranks <- P.sort[Ranks]
RANKS <- data.frame( Ranks, Value=P.ranks )

## Get Number of Variants beyond given P-Value
Orders <- 3:10
Thresh <- sort( c(1:9) %*% 10^-t(Orders) )
P.counts <- unlist(lapply( Thresh, function(x) length(which(P.sort<x)) ))
THRESH <- data.frame( Thresh, Count=P.counts )

## Write RANKS/THRESH Tables
write.table( RANKS, paste(PathToPlot,File_Tag,".TAB-RANKS.txt",sep=""),sep="\t",row.names=F,col.names=T,quote=F )
write.table( THRESH, paste(PathToPlot,File_Tag,".TAB-THRESH.txt",sep=""),sep="\t",row.names=F,col.names=T,quote=F )

## Plot These Results
COLS <- "tomato2"
jpeg( paste(PathToPlot,File_Tag,".P-ThreshRanks.jpeg",sep=""),height=900,width=1800,pointsize=30)
par(mfrow=c(1,2))
 # RANKS Table
plot( log10(RANKS$Ranks), -log10(RANKS$Value), main="P-Value vs Rank",xlab="Rank of Sorted P-Value (log)",ylab="-log10(P-Value)",pch="+",col=COLS[1],xaxt="n" )
# axis(1, at=0:10, label=paste("1e",0:10,sep="") )
axis(1, at=0:10, label=10^(0:10) )
abline(h=seq(0,20,1),v=seq(0,20,1),lty=3,col="grey50")
 # THRESH Table
plot( log10(THRESH$Thresh), log10(THRESH$Count), main="Number of P-Values beyond Threshold",xlab="Threshold (log)",ylab="# Variants beyond Threshold (log10)",pch="+",col=COLS[1],xaxt="n",ylim=c(0,max(log10(THRESH$Count),na.rm=T)) )
axis(1, at=-10:0, label=10^(-10:0) )
abline(h=seq(0,20,1),v=seq(-20,0,1),lty=3,col="grey50")
dev.off()

################################################
## FCT: COMPILE GENE-LEVEL STATS ###############
################################################

## Function to Calculate Area b/n Curve & Expected Line for QQ Plot of P-Values
Area_QQ <- function( Ps ) {
	NUM_VARS <- length(Ps)
	EXP <- -log10( 1:NUM_VARS / NUM_VARS )
	OBS <- -log10( sort( Ps ) )
	IND <- c(1:length(EXP),length(EXP),1)
	 # Calculate Area - Trapezoid Rule
	H_VALS <- EXP[2:NUM_VARS-1] - EXP[2:NUM_VARS]
	B_SUM <- ( OBS[2:NUM_VARS-1]-EXP[2:NUM_VARS-1] ) + ( OBS[2:NUM_VARS]-EXP[2:NUM_VARS] )
	AREA.traps <- .5*H_VALS*B_SUM
	AREA <- sum( AREA.traps )
	if ( length(OBS)!=length(EXP) ) { print( c(NUM_VARS,Ps[1],AREA) ) }
	E_VALS <- OBS - EXP
	MAX_E <- max( E_VALS )
	MAX_E.which <- which.max( E_VALS )
	 # Return Values of Interest
	OUTPUT <- list( AREA=AREA, MAX_E=MAX_E, MAX_E.which=MAX_E.which )
	return(OUTPUT)
}

## Function to Compile Stats (including N_VARS, BEST_P, AREA) for Genes
Gene_Stat <- function( P, G, Name, Script_or_Gene ) {

	## Pull out Coordinates ##
	if ( Script_or_Gene=="Gene" | Script_or_Gene=="gene" ) {
		
		## Which Rows have info about Gene?
		ROWS <- which( G$geneSymbol==Name )
		N.tx <- length(ROWS)

		## Check All chromosomes listed are the same?
		CHR.range.list <- as.character( G$chrom[ROWS] )
		if ( !all(CHR.range.list==CHR.range.list[1]) ) {
			## Compile Stats
			OUT.G_LEN <- rep("NA",2) ; names(OUT.G_LEN) <- c("Full","Ex")
			OUT.N_VAR <- rep("NA",2) ; names(OUT.N_VAR) <- c("Full","Ex")
			OUT.BEST_P <- rep("NA",2) ; names(OUT.BEST_P) <- c("Full","Ex")
			OUT.AREA <- rep("NA",2) ; names(OUT.AREA) <- c("Full","Ex")
			OUT.MAX_E <- rep("NA",2) ; names(OUT.MAX_E) <- c("Full","Ex")
			OUT.MAX_E.which <- rep("NA",2) ; names(OUT.MAX_E.which) <- c("Full","Ex")
			OUT <- list( G_LEN=OUT.G_LEN, N_VAR=OUT.N_VAR, BEST_P=OUT.BEST_P, AREA=OUT.AREA, MAX_E=OUT.MAX_E, MAX_E.which=OUT.MAX_E.which )
			# OUT <- list( G_LEN=OUT.G_LEN, N_VAR=OUT.N_VAR, BEST_P=OUT.BEST_P, AREA=OUT.AREA )
			return(OUT)
		}
		CHR.range <- CHR.range.list[1]
		 # Remove "chr" tag
		CHR <- gsub( "chr","", CHR.range )

		## TX Start & End Positions w/ Buffer
		TX_S.list <- as.numeric(G$txStart[ROWS]) # - 5000
		TX_E.list <- as.numeric(G$txEnd[ROWS]) # + 5000
		 # Find Biggest Possible Range of Values
		TX_S <- min(TX_S.list)
		TX_E <- max(TX_E.list)
		## Check coordinates for all transcripts are contiguous
		if ( N.tx>1 ) {
			Check <- unique( unlist(lapply(1:N.tx, function(x)TX_S.list[x]:TX_E.list[x])) )
			if ( length(Check)!=(1+TX_E-TX_S) ) {
				OUT.G_LEN <- rep("NA",2) ; names(OUT.G_LEN) <- c("Full","Ex")
				OUT.N_VAR <- rep("NA",2) ; names(OUT.N_VAR) <- c("Full","Ex")
				OUT.BEST_P <- rep("NA",2) ; names(OUT.BEST_P) <- c("Full","Ex")
				OUT.AREA <- rep("NA",2) ; names(OUT.AREA) <- c("Full","Ex")
				OUT.MAX_E <- rep("NA",2) ; names(OUT.MAX_E) <- c("Full","Ex")
				OUT.MAX_E.which <- rep("NA",2) ; names(OUT.MAX_E.which) <- c("Full","Ex")
				# OUT <- list( G_LEN=OUT.G_LEN, N_VAR=OUT.N_VAR, BEST_P=OUT.BEST_P, AREA=OUT.AREA )
				OUT <- list( G_LEN=OUT.G_LEN, N_VAR=OUT.N_VAR, BEST_P=OUT.BEST_P, AREA=OUT.AREA, MAX_E=OUT.MAX_E, MAX_E.which=OUT.MAX_E.which )
				return(OUT)
			}
		}

		## Exon Start & End Positions
		EX_B.list <- strsplit( as.character(G$exonStarts[ROWS]), "," )
		EX_E.list <- strsplit( as.character(G$exonEnds[ROWS]), "," )
		EX.list <- cbind( as.numeric(unlist(EX_B.list)),as.numeric(unlist(EX_E.list)) ) ; colnames(EX.list) <- c("EX_B","EX_E")
		if ( nrow(EX.list)>1 ) {
			RM.dups <- intersect( which(duplicated(EX.list[,1])), which(duplicated(EX.list[,2])) )
			if (length(RM.dups)>0 ) { EX.list <- EX.list[-RM.dups,] }
			EX.list <- EX.list[order(EX.list[,1],decreasing=F),]
			EX.list.fin <- aggregate( EX.list[,2], by=list(EX_B=EX.list[,1]), max )
		}else{ EX.list.fin <- EX.list }
		 # Final Start/End Positions
		EX_B <- EX.list.fin[,1]
		EX_E <- EX.list.fin[,2]

	}else{
		## Which Row has info about Transcript?
		ROW <- which( G$gtx==Name )

		## Get Chromosome
		CHR.range <- as.character( G$chrom[ROW] )
		 # Remove "chr" tag
		CHR <- gsub( "chr","", CHR.range )

		## TX Start & End Positions w/ Buffer
		TX_S <- as.numeric(G$txStart[ROW]) - 5000
		TX_E <- as.numeric(G$txEnd[ROW]) + 5000
		## CD Start & End Positions w/ Buffer
		CD_S <- as.numeric(G$cdsStart[ROW]) - 5000
		CD_E <- as.numeric(G$cdsEnd[ROW]) + 5000

		## Exon Start & End Positions
		EX_B <- as.numeric( strsplit( as.character(G$exonStarts[ROW]), "," )[[1]] )
		EX_E <- as.numeric( strsplit( as.character(G$exonEnds[ROW]), "," )[[1]] )
	}

	## Pull Chromosomal Positions from P file ##

	## Full Gene
	ROWS <- which( P$CHR==CHR & P$BP>=TX_S & P$BP<=TX_E )
	 # Chromosomal Positions
	POS <- P$BP[ROWS]

	## Exons
	ROWS.ex <- unique(ROWS[ unlist(lapply( 1:length(EX_B), function(x) which(POS>=EX_B[x] & POS<=EX_E[x]) )) ])
	 # Chromosomal Positions
	POS.ex <- P$BP[ROWS.ex]

	## Pull P-Values for Full Gene & Exons ##
	Ps <- P$P[ROWS]
	Ps.ex <- P$P[ROWS.ex]

	## Compile Stats about Gene ##
	 # Gene Length (end to end)
	G_LEN <- TX_E - TX_S
	G_LEN.ex <- length( Reduce( union, lapply( 1:length(EX_B), function(x) EX_B[x]:EX_E[x] ) ) )
	 # Number of Variants in Gene
	N_VAR <- length(ROWS)
	N_VAR.ex <- length(ROWS.ex)
	 # Best P-Value
	BEST_P <- ifelse( N_VAR>0, min(Ps,na.rm=T), "NA" ) # min( Ps, na.rm=T )
	BEST_P.ex <- ifelse( N_VAR.ex>0, min(Ps.ex,na.rm=T), "NA" ) # min( Ps.ex, na.rm=T )
	 # Area under QQ Plot (for Enrichment)
	# AREA <- ifelse( N_VAR>0, Area_QQ(Ps), "NA" ) # Area_QQ( Ps )
	# AREA.out <- ifelse( N_VAR>0, Area_QQ(Ps), "NA" ) # Area_QQ( Ps )
	if ( N_VAR>0 ) {
		AREA.out <- Area_QQ(Ps)
		AREA <- AREA.out$AREA
		MAX_E <- AREA.out$MAX_E
		MAX_E.which <- AREA.out$MAX_E.which
	}else{
		AREA <- MAX_E <- MAX_E.which <- "NA"
	}
	if ( N_VAR.ex>0 ) {
		AREA.out.ex <- Area_QQ(Ps.ex)
		AREA.ex <- AREA.out.ex$AREA
		MAX_E.ex <- AREA.out.ex$MAX_E
		MAX_E.which.ex <- AREA.out.ex$MAX_E.which
	}else{
		AREA.ex <- MAX_E.ex <- MAX_E.which.ex <- "NA"
	}
	# AREA.ex <- ifelse( N_VAR.ex>0, Area_QQ(Ps.ex), "NA" ) # Area_QQ( Ps.ex )
	# AREA.out.ex <- ifelse( N_VAR.ex>0, Area_QQ(Ps.ex), "NA" ) # Area_QQ( Ps.ex )
	# AREA.ex <- AREA.out.ex$AREA
	# MAX_E.ex <- AREA.out.ex$MAX_E
	# MAX_E.which.ex <- AREA.out.ex$MAX_E.which

	## Compile Stats
	OUT.G_LEN <- c( G_LEN, G_LEN.ex ) ; names(OUT.G_LEN) <- c("Full","Ex")
	OUT.N_VAR <- c( N_VAR, N_VAR.ex ) ; names(OUT.N_VAR) <- c("Full","Ex")
	OUT.BEST_P <- c( BEST_P, BEST_P.ex ) ; names(OUT.BEST_P) <- c("Full","Ex")
	OUT.AREA <- c( AREA, AREA.ex ) ; names(OUT.AREA) <- c("Full","Ex")
	OUT.MAX_E <- c( MAX_E, MAX_E.ex ) ; names(OUT.MAX_E) <- c("Full","Ex")
	OUT.MAX_E.which <- c( MAX_E.which, MAX_E.which.ex ) ; names(OUT.MAX_E.which) <- c("Full","Ex")
	OUT <- list( G_LEN=OUT.G_LEN, N_VAR=OUT.N_VAR, BEST_P=OUT.BEST_P, AREA=OUT.AREA, MAX_E=OUT.MAX_E, MAX_E.which=OUT.MAX_E.which )
	return(OUT)
}

################################################
## GET GENE SUMMARY STATS ######################
################################################
print("Calculating Gene/Transcript Summary Stats (by Chromosome)")

## Run on All Genes/Transcripts (by Chromosome)
 # Specify "Gene" [or] "Transcript"
if ( Script_or_Gene=="Gene" | Script_or_Gene=="gene" | Script_or_Gene=="GENE" ) {
	Name_Col <- "geneSymbol"
}else{ Name_Col <- "gtx" }
N_QUANT <- 20
 # Create Output Variable to compile Chromosome Results
OUT <- list()
TIMES <- numeric(length(chroms)*N_QUANT) ; names(TIMES) <- paste( rep(paste("chr",chroms,sep=""),rep(N_QUANT,length(chroms))), 1:N_QUANT,sep="_" )
N_GENES <- TIMES
 # Loop Thru Chromosomes
start <- proc.time()
for ( chr in 1:23 ) {
# for ( chr in 21:22 ) {
	chr_tag <- paste("chr",chroms[chr],sep="")	
	print(paste("Start",chr_tag,"Pvals -", round(proc.time()-start,3)[3] ))
	P.chr <- P[ which(P$CHR==chr), ]
	print(paste("Loaded",nrow(P.chr),"Pvals on",chr_tag," -", round(proc.time()-start,3)[3] ))
	if ( length(which(is.na(P.chr$P)))>0 ) { P.chr <- P.chr[ which(!is.na(P.chr$P)), ] }
	print(paste("Cut down to",nrow(P.chr),"Pvals on",chr_tag," -", round(proc.time()-start,3)[3] ))
	Quantiles <- quantile( range(P.chr$BP), c(0:N_QUANT)/N_QUANT )
	for ( q in 1:N_QUANT ) {
		chr_q_tag <- paste(chr_tag,q,sep="_")
		quant_rng <- Quantiles[c(q,q+1)]
		P.q <- P.chr[ which( P.chr$BP>=quant_rng[1] & P.chr$BP<quant_rng[2]+5e6 ), ]
		Names <- unique( G[ which( G$chrom==chr_tag & G$txStart>=quant_rng[1] & G$txStart<quant_rng[2]), Name_Col ] )
		Names <- setdiff( Names, GENES.prob )
		print(paste( length(Names),Script_or_Gene, "on", chr_q_tag )) # print(paste( length(Names),Script_or_Gene, "on", chr_tag,": Q",q ))
		OUT[[chr_q_tag]] <- lapply( Names, function(Name) Gene_Stat( P.q, G, Name, Script_or_Gene ) ) ; names(OUT[[chr_q_tag]]) <- Names
		TIMES[chr_q_tag] <- round((proc.time()-start)[3],3)
		N_GENES[chr_q_tag] <- length(Names)
		print(paste("Done with",chr_q_tag,"-", TIMES[chr_q_tag] ))
	}
	# print(paste("Done with",chr_tag,"-", TIMES ))
}
print(paste( "Total Runtime:",max(TIMES) ))

## Compile Results for all Transcripts
print("Compile Stats for Genome")
G_LENS <- matrix( as.numeric(unlist(lapply( OUT, function(x) lapply(x,function(y)y$G_LEN) ))), ncol=2,byrow=T )
N_VARS <- matrix( as.numeric(unlist(lapply( OUT, function(x) lapply(x,function(y)y$N_VAR) ))), ncol=2,byrow=T )
AREAS <- matrix( as.numeric(unlist(lapply( OUT, function(x) lapply(x,function(y)y$AREA) ))), ncol=2,byrow=T )
BEST_PS <- matrix( as.numeric(unlist(lapply( OUT, function(x) lapply(x,function(y)y$BEST_P) ))), ncol=2,byrow=T )
MAX_ES <- matrix( as.numeric(unlist(lapply( OUT, function(x) lapply(x,function(y)y$MAX_E) ))), ncol=2,byrow=T )
MAX_E.whichS <- matrix( as.numeric(unlist(lapply( OUT, function(x) lapply(x,function(y)y$MAX_E.which) ))), ncol=2,byrow=T )
colnames(AREAS) <- colnames(BEST_PS) <- colnames(N_VARS) <- colnames(G_LENS) <- colnames(MAX_ES) <- colnames(MAX_E.whichS) <- c("Full","Ex")
rownames(AREAS) <- rownames(BEST_PS) <- rownames(N_VARS) <- rownames(G_LENS) <- rownames(MAX_ES) <- rownames(MAX_E.whichS) <- unlist(lapply( OUT, names ))
# pairs(data.frame( N_VARS=log10(N_VARS[,1]), BEST_PS=-log10(BEST_PS[,1]), AREAS=AREAS[,1], MAX_ES=MAX_ES[,1] ), col=factor(MAX_E.whichS[,1]==N_VARS[,1]), cex=log10(N_VARS[,1]))
# pairs(data.frame( N_VARS=log10(N_VARS[,2]), BEST_PS=-log10(BEST_PS[,2]), AREAS=AREAS[,2], MAX_ES=MAX_ES[,2] ), col=factor(MAX_E.whichS[,2]<=N_VARS[,2]), cex=log10(N_VARS[,2]))
# plot( -log10(BEST_PS[,1]), MAX_ES[,1], col=as.numeric(factor(MAX_E.whichS[,1]==N_VARS[,1])), cex=log10(N_VARS[,1])) ; abline(0,1)

## Write Output Tables
print(paste("Writing Summary Tables to:",PathToPlot))
write.table( G_LENS, paste(PathToPlot,File_Tag,".TAB-G_LENS.txt",sep=""),sep="\t",row.names=T,col.names=T,quote=F )
write.table( N_VARS, paste(PathToPlot,File_Tag,".TAB-N_VARS.txt",sep=""),sep="\t",row.names=T,col.names=T,quote=F )
write.table( AREAS, paste(PathToPlot,File_Tag,".TAB-AREAS.txt",sep=""),sep="\t",row.names=T,col.names=T,quote=F )
write.table( BEST_PS, paste(PathToPlot,File_Tag,".TAB-BEST_PS.txt",sep=""),sep="\t",row.names=T,col.names=T,quote=F )
write.table( MAX_ES, paste(PathToPlot,File_Tag,".TAB-MAX_ES.txt",sep=""),sep="\t",row.names=T,col.names=T,quote=F )
write.table( MAX_E.whichS, paste(PathToPlot,File_Tag,".TAB-MAX_E.whichS.txt",sep=""),sep="\t",row.names=T,col.names=T,quote=F )

## Compile Stats into 2 Data Frames (Full & Exonic)
FULL <- cbind( G_LENS[,"Full"], N_VARS[,"Full"], BEST_PS[,"Full"], AREAS[,"Full"], MAX_ES[,"Full"], MAX_E.whichS[,"Full"], log10(G_LENS[,"Full"]), log10(N_VARS[,"Full"]), -log10(BEST_PS[,"Full"]) )
EXON <- cbind( G_LENS[,"Ex"], N_VARS[,"Ex"], BEST_PS[,"Ex"], AREAS[,"Ex"], MAX_ES[,"Ex"], MAX_E.whichS[,"Ex"], log10(G_LENS[,"Ex"]), log10(N_VARS[,"Ex"]), -log10(BEST_PS[,"Ex"]) )
# FULL <- cbind( G_LENS[,"Full"], N_VARS[,"Full"], BEST_PS[,"Full"], AREAS[,"Full"], log10(G_LENS[,"Full"]), log10(N_VARS[,"Full"]), -log10(BEST_PS[,"Full"]) )
# EXON <- cbind( G_LENS[,"Ex"], N_VARS[,"Ex"], BEST_PS[,"Ex"], AREAS[,"Ex"], log10(G_LENS[,"Ex"]), log10(N_VARS[,"Ex"]), -log10(BEST_PS[,"Ex"]) )
rownames(FULL) <- rownames(EXON) <- unlist(lapply( OUT, names ))
colnames(FULL) <- colnames(EXON) <- c("G_LEN","N_VAR","BEST_P","AREA","MAX_E","MAX_E.which","logG_LEN","logN_VAR","-logBEST_P")
write.table( FULL, paste(PathToPlot,File_Tag,".TAB-FULL.txt",sep=""),sep="\t",row.names=T,col.names=T,quote=F )
write.table( EXON, paste(PathToPlot,File_Tag,".TAB-EXON.txt",sep=""),sep="\t",row.names=T,col.names=T,quote=F )
 # Compile Stats on Genes of Interest
FULL.GOI <- cbind( GOI$Full, FULL[rownames(GOI$Full),c("G_LEN","N_VAR","logG_LEN","logN_VAR","MAX_E","MAX_E.which")], -log10(GOI$Full[,"BEST_P"]) ) ; colnames(FULL.GOI)[ncol(FULL.GOI)] <- "-logBEST_P"
EXON.GOI <- cbind( GOI$Ex, EXON[rownames(GOI$Ex),c("G_LEN","N_VAR","logG_LEN","logN_VAR","MAX_E","MAX_E.which")], -log10(GOI$Ex[,"BEST_P"]) ) ; colnames(EXON.GOI)[ncol(EXON.GOI)] <- "-logBEST_P"

################################################
## PLOT COMPILED RESULTS #######################
################################################
print("Creating Scatter Plots of Summary Stats")

## Plot Time to Calculate Stats
jpeg( paste(PathToPlot,File_Tag,".0-Runtime.Q",N_QUANT,".jpeg",sep=""),height=900,width=1800,pointsize=30)
par(mfrow=c(1,2))
plot( 1:length(TIMES), TIMES, type="l",col="slateblue3",lwd=2,main="Progress",xlab="Genomic Region",ylab="Runtime",xaxt="n")
abline(v=seq(0,1000,N_QUANT),lty=3,col="grey50",lwd=1 )
axis( 1, at=(1:length(TIMES))[seq(1,length(TIMES),10)],label=names(TIMES)[seq(1,length(TIMES),10)], las=2 )
Y_VALS <- c(TIMES[1],diff(TIMES))
plot( Y_VALS ~ N_GENES, pch="+",col="slateblue3",main="Runtime vs Chromosome",xlab="# Genes per Chrom",ylab="Runtime",ylim=c(0,max(Y_VALS)) )
text( N_GENES, Y_VALS, label=names(N_GENES),pos=2,col="tomato2" )
points( Y_VALS ~ N_GENES, pch="+",col="slateblue3" )
dev.off()

## Function to Plot Variables against one another
Plot_Variables <- function( X_Lab, Y_Lab, Color, tag ) {
	if ( tag=="Full" ) {
		DATA <- FULL
		DATA.GOI <- FULL.GOI
	}else{
		DATA <- EXON
		DATA.GOI <- EXON.GOI
	}
	X_Vals <- DATA[,X_Lab]
	Y_Vals <- DATA[,Y_Lab]
	SUBSET <- which( is.na(Y_Vals) | Y_Vals=="-Inf" | is.na(X_Vals) | X_Vals=="-Inf" )
	X_GOI <- DATA.GOI[,X_Lab]
	Y_GOI <- DATA.GOI[,Y_Lab]
	if ( length(SUBSET)>0 ) {
		X_LIM <- range(X_Vals[-SUBSET],X_GOI,na.rm=T)
		Y_LIM <- range(Y_Vals[-SUBSET],Y_GOI,na.rm=T)
		MOD <- lm( Y_Vals ~ X_Vals, subset=-SUBSET )
	}else{
		MOD <- lm( Y_Vals ~ X_Vals )
		X_LIM <- range(X_Vals,X_GOI,na.rm=T)
		Y_LIM <- range(Y_Vals,Y_GOI,na.rm=T)
	}
	plot( Y_Vals ~ X_Vals, xlab=X_Lab,ylab=Y_Lab,xlim=X_LIM,ylim=Y_LIM,main=paste(Y_Lab,"vs",X_Lab,"(",tag,")"), pch=20,col=Color )
	abline(MOD, col="black",lwd=2,lty=2 )
	P_Val <- summary(MOD)$coefficients["X_Vals","Pr(>|t|)"]
	if ( length(SUBSET)>0 ) {
		text( quantile(range(X_Vals[-SUBSET],na.rm=T),.01),quantile(range(Y_Vals[-SUBSET],na.rm=T),.95), label=paste("P =",formatC(P_Val,digits=2,format="e")),pos=4 )
	}else{ text( quantile(range(X_Vals,na.rm=T),.01),quantile(range(Y_Vals,na.rm=T),.95), label=paste("P =",formatC(P_Val,digits=2,format="e")),pos=4 ) }
	points( X_GOI, Y_GOI, col="tomato2",pch="+",lwd=2 )
	text( X_GOI, Y_GOI, label=rownames(GOI[[tag]]),col="tomato2",pos=4 )
}

## Plot Pairs for "Full" and "Exon"
COLS <- c("springgreen2","cadetblue2","steelblue2","slateblue3")
Tags <- c("Full","Ex")
for ( tag in Tags ) {
	## Raw Gene Lengths
	jpeg( paste(PathToPlot,File_Tag,".1-Scatter.",tag,".jpeg",sep=""),height=1200,width=1800,pointsize=30)
	par(mfrow=c(2,3))
	Plot_Variables( "G_LEN","N_VAR", COLS[1], tag )
	Plot_Variables( "G_LEN","-logBEST_P", COLS[2], tag )
	Plot_Variables( "G_LEN","AREA", COLS[3], tag )

	Plot_Variables( "-logBEST_P","AREA", COLS[4], tag )
	Plot_Variables( "N_VAR","-logBEST_P", COLS[2], tag )
	Plot_Variables( "N_VAR","AREA", COLS[3], tag )
	dev.off()

	## Log Transformed Gene Lengths
	jpeg( paste(PathToPlot,File_Tag,".1-Scatter.log.",tag,".jpeg",sep=""),height=1200,width=1800,pointsize=30)
	par(mfrow=c(2,3))
	Plot_Variables( "logG_LEN","logN_VAR", COLS[1], tag )
	Plot_Variables( "logG_LEN","-logBEST_P", COLS[2], tag )
	Plot_Variables( "logG_LEN","AREA", COLS[3], tag )

	Plot_Variables( "-logBEST_P","AREA", COLS[4], tag )
	Plot_Variables( "logN_VAR","-logBEST_P", COLS[2], tag )
	Plot_Variables( "logN_VAR","AREA", COLS[3], tag )
	dev.off()

}

################################################
## FIT DISTRIBUTION TO AREA DATA ###############
################################################
print("Plotting AREA Distribution")

## FCT: Null Distribution of AREA
Null_Area <- function( AREA.vec, GOI, Color, tag ) {
	## Normalize AREA Measurements
	# AREA.n <- ( AREA.vec - mean(AREA.vec,na.rm=T) ) / sd(AREA.vec,na.rm=T)
	AREA.n <- AREA.vec
	AREA.n <- AREA.n[which(!is.na(AREA.n))]

	## Johnson Fit of AREA Data
		# z=gamma+delta log(f(u)), with u=(x-xi)/lambda
		#      and where f( ) has four possible forms:
		#        SL:  f(u)=u the log normal
		#        SU:  f(u)=u+sqrt(1+u^2) an unbounded distribution
		#        SB:  f(u)=u/(1-u) a bounded distribution
		#        SN:  \exp(u) the normal
	Moments <- moments(AREA.n)
	O <- JohnsonFit( AREA.n, moment="quant" )
	XX <- seq(-20,20,.01)
	RNG <- extendrange(AREA.n)
	XLIM <- c( RNG[1], 1+max(RNG[2],GOI[[tag]][,"AREA"]) )
	BRKS <- seq( floor(RNG[1]), ceiling(RNG[2]), .1 )
	HIST <- hist(AREA.n,freq=F,breaks=BRKS,col=Color,main=paste("Distribution of Area Values -",tag),xlab="Area", xlim=XLIM )
	abline( v=seq(floor(RNG[1]),ceiling(XLIM[2]),1),h=seq(0,2,.2),lty=3,col="grey50",lwd=1 )
	HIST <- hist(AREA.n,freq=F,breaks=BRKS,col=Color,add=T)
	points( XX, dJohnson(XX,O), col="chocolate2",type="l",lwd=3) # plot(function(x)dJohnson(x,O), -20, 20, add=T, col="blue")
	text( quantile(AREA.n,.99,na.rm=T), max(HIST$density), label=paste(names(Moments),round(Moments,3),sep=" = ",collapse="\n" ),pos=1 )
	for ( g in 1:nrow(GOI[[tag]]) ) { gene_area <- GOI[[tag]][g,"AREA"]
		arrows( gene_area, .1, gene_area, 0, lwd=3, col="firebrick1" )
		P.dis <- pJohnson(gene_area,O,lower.tail=F)
		P.prm <- (1+length(which(AREA.n>gene_area))) / (length(AREA.n)+1)
		text( gene_area, .12, paste(rownames(GOI[[tag]])[g],"-",GOI[[tag]][g,"AREA"],"\np.dis =",formatC(P.dis,2,format="e"),"\np.prm =",formatC(P.prm,2,format="e")), srt=90, pos=4, col="firebrick1" )
	}
}

## Plot Null Distribution w/ GOIs
COLS <- c("chartreuse1","dodgerblue2","tomato2")
jpeg( paste(PathToPlot,File_Tag,".2-Area_Distrib.jpeg",sep=""),height=1200,width=1600,pointsize=24)
par(mfrow=c(2,1))
Null_Area( AREAS[,1], GOI, COLS[1], "Full" )
Null_Area( AREAS[,2], GOI, COLS[2], "Ex" )
dev.off()

## Plot Null Distribution w/ GOIs
COLS <- c("chartreuse1","dodgerblue2","tomato2")
jpeg( paste(PathToPlot,File_Tag,".2-MaxE_Distrib.jpeg",sep=""),height=1200,width=1600,pointsize=24)
par(mfrow=c(2,1))
Null_Area( MAX_ES[,1], GOI, COLS[1], "Full" )
Null_Area( MAX_ES[,2], GOI, COLS[2], "Ex" )
dev.off()

###########################################################################
###########################################################################
## END OF DOC #############################################################
###########################################################################
###########################################################################

# TEST <- OUT$chr1
# # ## Compile Results for all Transcripts
# AREAS <- matrix( as.numeric(unlist(lapply( TEST, function(x) x$AREA ))), byrow=T,ncol=2 ) ; colnames(AREAS) <- c("Full","Ex") ; rownames(AREAS) <- names(TEST)
# BEST_PS <- matrix( as.numeric(unlist(lapply( TEST, function(x) x$BEST_P ))), byrow=T,ncol=2 ) ; colnames(BEST_PS) <- c("Full","Ex") ; rownames(BEST_PS) <- names(TEST)
# N_VARS <- matrix( as.numeric(unlist(lapply( TEST, function(x) x$N_VAR ))), byrow=T,ncol=2 ) ; colnames(N_VARS) <- c("Full","Ex") ; rownames(N_VARS) <- names(TEST)

# ## 
# g.1 <- sum( (AREA.n-mean(AREA.n))^3 ) / ( length(AREA.n)*sd(AREA.n)^3 )
# g.2 <- sum( (AREA.n-mean(AREA.n))^4 ) / ( length(AREA.n)*sd(AREA.n)^4 ) - 3
# g.1 <- skewness(AREAS[,1],na.rm=T)
# g.2 <- kurtosis(AREAS[,1],na.rm=T)
# K <- ( g.1^2*(g.2+6)^2 ) / ( 4*(4*g.2-3*g.1^2+12)*(2*g.2-3*g.1^2) )

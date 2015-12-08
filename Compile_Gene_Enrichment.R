
################################################
## LOAD DATA ###################################
################################################
library(SuppDists)

## Set Paths
# PathToPermP <- "/projects/janssen/ASSOCIATION/20151016_Th"
# PathToPermP2 <- "/projects/janssen/ASSOCIATION/20151118_Th"
# PathToPermArea <- "/projects/janssen/ASSOCIATION/20151103_Th"
PathToPermP2 <- "/projects/janssen/ASSOCIATION/20151202_Th"
PathToPermArea <- "/projects/janssen/ASSOCIATION/20151202_Th"
PathToTrueP <- "/projects/janssen/ASSOCIATION/20150104_LT8_DEL_MNe_MN_DAS_BL_MN_PC1_PC2/20150104_LT8_DEL_MNe_MN_DAS_BL_MN_PC1_PC2.P"
# PathToTrueP2 <- "/projects/janssen/ASSOCIATION/20151102_EnrichCompile/"
# PathToTrueArea <- "/projects/janssen/ASSOCIATION/20151102_EnrichCompile/"
PathToTrueP2 <- "/projects/janssen/ASSOCIATION/20151201_GeneEnrichment_DAS/"
PathToTrueArea <- "/projects/janssen/ASSOCIATION/20151201_GeneEnrichment_DAS/"
PathToPlot <- "/projects/janssen/ASSOCIATION/20151201_GeneEnrichment_DAS/Compile_Plots/"
if ( !file.exists(PathToPlot) ) { dir.create( PathToPlot ) }
 # Get True File Tag
Split.1 <- strsplit( PathToTrueP, "/" )[[1]] 
True_Tag <- gsub( ".P","",Split.1[length(Split.1)],fixed=T )

## Load Data

## Load Compiled Summary Stats
Table_List <- c("FULL","EXON")
Metrics <- as.character( read.table( paste(PathToTrueArea,True_Tag,".TAB-FULL.txt",sep=""), header=F,nrow=1,sep="\t",colClasses="character" )[1,] )
Colnames.A <- c("Gene",Metrics)
 # True
A.f <- read.table( paste(PathToTrueArea,True_Tag,".TAB-FULL.txt",sep=""), header=F,skip=1,sep="\t" )
A.e <- read.table( paste(PathToTrueArea,True_Tag,".TAB-EXON.txt",sep=""), header=F,skip=1,sep="\t" )
colnames(A.f) <- colnames(A.e) <- Colnames.A
A <- data.frame( FULL=A.f, EXON=A.e )
 # Permuted
PERM.A <- list()
for ( tab in Table_List ) {
	PERM.A[[tab]] <- list()
	for ( t in 1:5 ) {
		Dir_Path <- paste(PathToPermArea,t,"_Permute/",sep="")
		N_Files <- length(grep( tab,list.files(Dir_Path), value=T))
		# print( N_Files )
		for ( i in 1:N_Files ) {
			tag <- paste("T",t,"I",i,sep="")
			PathToFile <- paste( Dir_Path,"20151016_Th",t,"_Perm_Pheno.",i,"_DAS_BL_MN_PC1_PC2.TAB-",tab,".txt",sep="")
			# print( PathToFile )
			PERM.A[[tab]][[tag]] <- read.table( PathToFile, header=F,skip=1,sep="\t" )
		}
		print(paste(tab,t))
	}
}

## Load Ranked/Thresholded P-Value Data
Table_List <- c("THRESH","RANKS")
 # True
P.t <- read.table( paste(PathToTrueP2,"20150104_LT8_DEL_MNe_MN_DAS_BL_MN_PC1_PC2.TAB-THRESH.txt",sep=""), header=T,sep="\t" )
P.r <- read.table( paste(PathToTrueP2,"20150104_LT8_DEL_MNe_MN_DAS_BL_MN_PC1_PC2.TAB-RANKS.txt",sep=""), header=T,sep="\t" )
 # Permuted
PERM.P <- list()
for ( tab in Table_List ) {
	PERM.P[[tab]] <- list()
	for ( t in 1:5 ) {
		Dir_Path <- paste(PathToPermP2,t,"_Permute/",sep="")
		N_Files <- length(grep( tab,list.files(Dir_Path), value=T))
		# print( N_Files )
		for ( i in 1:N_Files ) {
			tag <- paste("T",t,"I",i,sep="")
			PathToFile <- paste( Dir_Path,"20151016_Th",t,"_Perm_Pheno.",i,"_DAS_BL_MN_PC1_PC2.TAB-",tab,".txt",sep="")
			# print( PathToFile )
			PERM.P[[tab]][[tag]] <- read.table( PathToFile, header=T,sep="\t" )
		}
		print(paste(tab,t))
	}
}

# ## Load P-Values
#  # True
# P <- read.table( PathToTrueP, header=F,skip=1,sep="\t" )
#  # Permuted
# PERM.P <- list()
# for ( t in 1:5 ) {
# 	Dir_Path <- paste(PathToPermP,t,"_Permute/",sep="")
# 	N_Files <- length(grep( tab,list.files(Dir_Path), value=T))
# 	# print( N_Files )
# 	for ( i in 1:N_Files ) {
# 		tag <- paste("T",t,"I",i,sep="")
# 		File_Path.1 <- paste("20151118_Th",t,"_Perm_Pheno.",i,"_DAS_BL_MN_PC1_PC2")
# 		File_Path <- paste( File_Path.1,"/",File_Path.1,".P",sep="")
# 		PathToFile <- paste( Dir_Path,File_Path,sep="")
# 		# print( PathToFile )
# 		PERM.P[[tab]][[tag]] <- read.table( PathToFile, header=F,skip=1,sep="\t" )
# 	}
# 	print(paste(tab,t))
# }

## Compile Data
 # FULL/EXON Gene-Level Tables
PERM.f <- Reduce( rbind, PERM.A$FULL )
PERM.e <- Reduce( rbind, PERM.A$EXON )
colnames(PERM.f) <- colnames(PERM.e) <- Colnames.A
 # THRESH/RANK P-Value Tables
PERM.t <- Reduce( cbind, PERM.P$THRESH ) ; PERM.t <- PERM.t[,c(1,seq(2,ncol(PERM.t),2))]
PERM.t.2 <- data.frame( True=P.t[,-1], Perm=PERM.t[,-1] ) ; rownames(PERM.t.2) <- P.t[,1]
PERM.r <- Reduce( cbind, PERM.P$RANKS ) ; PERM.r <- PERM.r[,c(1,seq(2,ncol(PERM.r),2))]
PERM.r.2 <- data.frame( True=P.r[,-1], Perm=PERM.r[,-1] ) ; rownames(PERM.r.2) <- P.r[,1]

################################################
## NULL DISTRIBUTION of P-VALUES ###############
################################################

## Calculate Significance of Actual Data vs Permuted P-Values
SIG.t <- (apply(PERM.t.2, 1, function(x) length(x)-rank(x,ties.method="min")[1] ) + 1) / (ncol(PERM.t.2))
SIG.r <- (apply(PERM.r.2, 1, function(x) rank(x,ties.method="min")[1] ) + 1) / (ncol(PERM.r.2))
SIG.t <- (apply(PERM.t.2, 1, function(x) length(x)-rank(x,ties.method="average")[1] ) + 1) / (ncol(PERM.t.2))
SIG.r <- (apply(PERM.r.2, 1, function(x) rank(x,ties.method="average")[1] ) + 1) / (ncol(PERM.r.2))
data.frame( True.t=PERM.t.2[,"True"], SIG.t=round(SIG.t,3) )
data.frame( True.r=PERM.r.2[,"True"], SIG.r=round(SIG.r,3) )

## Plot Results
COLS <- "tomato2"
# jpeg( paste(PathToPlot,True_Tag,"_vsNull.P-ThreshRanks.jpeg",sep=""),height=1200,width=1600,pointsize=24)
# par(mfrow=c(2,1))
# par(mar=c(5,5,5,5))
 # RANKS Table
jpeg( paste(PathToPlot,True_Tag,"_vsNull.P-Ranks.jpeg",sep=""),height=1200,width=2000,pointsize=30)
XVALS <- as.numeric(rownames(PERM.r.2))
XLIM <- range( log10(XVALS) )
YLIM <- c( 0, max(-log10(PERM.r.2)) )
plot( 0,0,type="n",xlim=XLIM,ylim=YLIM, main="P-Value vs Rank",xlab="Rank of Sorted P-Value (log)",ylab="-log10(P-Value)",xaxt="n" )
axis(1, at=0:10, label=10^(0:10) )
abline(h=seq(0,20,1),v=seq(0,20,1),lty=3,col="grey50")
abline( h=YLIM[2]*seq(0,1,.1), lty=c(1,rep(3,4),1,rep(3,4),1),col="dodgerblue1",lwd=2)
apply( PERM.r.2[,-1], 2, function(x) points( log10(XVALS), -log10(x), pch=20,col="grey25",type="o" ) )
points( log10(XVALS), -log10(PERM.r.2[,1]), type="o",pch=20,col=COLS[1] )
axis( 4, at=seq(YLIM[1],YLIM[2],length.out=11),label=seq(0,1,.1),col="dodgerblue2" )
mtext("Permuted P-Value", side=4, line=3, cex.lab=1, col="dodgerblue3")
points( log10(XVALS),YLIM[2]*SIG.r,col="dodgerblue2",type="l",lwd=4 )
dev.off()
 # THRESH Table
jpeg( paste(PathToPlot,True_Tag,"_vsNull.P-Thresh.jpeg",sep=""),height=1200,width=2000,pointsize=30)
XVALS <- as.numeric(rownames(PERM.t.2))
XLIM <- range( log10(XVALS) )
YLIM <- range( log10(PERM.t.2[which(PERM.t.2!=0,arr.ind=T)]) )
plot( 0,0,type="n",xlim=XLIM,ylim=YLIM, main="Number of P-Values beyond Threshold",xlab="Threshold (log)",ylab="# Variants beyond Threshold (log10)",xaxt="n" )
axis(1, at=-10:0, label=10^(-10:0) )
abline(h=seq(0,20,1),v=seq(-20,0,1),lty=3,col="grey50")
abline( h=YLIM[2]*seq(0,1,.1), lty=c(1,rep(3,4),1,rep(3,4),1),col="dodgerblue1",lwd=2)
apply( PERM.t.2[,-1], 2, function(x) points( log10(XVALS), log10(x), pch=20,col="grey25",type="o" ) )
points( log10(XVALS), log10(PERM.t.2[,1]), type="o",pch=20,col=COLS[1] )
axis( 4, at=seq(YLIM[1],YLIM[2],length.out=11),label=seq(0,1,.1),col="dodgerblue2" )
mtext("Permuted P-Value", side=4, line=3, cex.lab=1, col="dodgerblue3")
points( log10(XVALS),YLIM[2]*SIG.t,col="dodgerblue2",type="l",lwd=4 )
dev.off()

################################################
## NULL DISTRIBUTION of AREAS ##################
################################################

################################################
## Pick out top N genes by Area from True Data
GOI.names.custom <- c("GRIN2B","CNTN5")#,"HLA-DRB1")
N.top <- 10
GOI.names <- list()
GOI.names$area.f <- as.character( A.f$Gene[order(A.f$AREA,decreasing=T)][1:N.top] )
GOI.names$area.e <- as.character( A.e$Gene[order(A.e$AREA,decreasing=T)][1:N.top] )
GOI.names$maxe.f <- as.character( A.f$Gene[order(A.f$MAX_E,decreasing=T)][1:N.top] )
GOI.names$maxe.e <- as.character( A.e$Gene[order(A.e$MAX_E,decreasing=T)][1:N.top] )
GOI.names$best.f <- as.character( A.f$Gene[order(A.f$BEST_P,decreasing=F)][1:N.top] )
GOI.names$best.e <- as.character( A.e$Gene[order(A.e$BEST_P,decreasing=F)][1:N.top] )
GOI.genes <- Reduce( union, GOI.names) # list(GOI.names$custom, GOI.names$area.f, GOI.names$area.e, GOI.names$maxe.f, GOI.names$maxe.e) )
GOI.rows <- match( GOI.genes, A.f$Gene )
GOI <- list()
GOI[["Full"]] <- A.f[ GOI.rows, c("BEST_P","AREA","MAX_E") ] ; rownames(GOI[["Full"]]) <- GOI.genes
GOI[["Ex"]] <- A.e[ GOI.rows, c("BEST_P","AREA","MAX_E") ] ; rownames(GOI[["Ex"]]) <- GOI.genes
 # Write Tables w/ Area, etc...
write.table( GOI$Full[ GOI.names$area.f, ], file=paste(PathToPlot,"TAB-GOI-AREA_F.txt",sep=""), sep="\t",row.names=T,col.names=T,quote=F )
write.table( GOI$Ex[ GOI.names$area.e, ], file=paste(PathToPlot,"TAB-GOI-AREA_E.txt",sep=""), sep="\t",row.names=T,col.names=T,quote=F )
write.table( GOI$Full[ GOI.names$maxe.f, ], file=paste(PathToPlot,"TAB-GOI-MAXE_F.txt",sep=""), sep="\t",row.names=T,col.names=T,quote=F )
write.table( GOI$Ex[ GOI.names$maxe.e, ], file=paste(PathToPlot,"TAB-GOI-MAXE_E.txt",sep=""), sep="\t",row.names=T,col.names=T,quote=F )
write.table( GOI$Full[ GOI.names$best.f, ], file=paste(PathToPlot,"TAB-GOI-BEST_F.txt",sep=""), sep="\t",row.names=T,col.names=T,quote=F )
write.table( GOI$Ex[ GOI.names$best.e, ], file=paste(PathToPlot,"TAB-GOI-BEST_E.txt",sep=""), sep="\t",row.names=T,col.names=T,quote=F )

# ## Gene Areas from Actual GWAS
# GOI <- list()
# GOI[["Full"]] <- array(,c(2,2))
# colnames(GOI[["Full"]]) <- c("BEST_P","AREA")
# rownames(GOI[["Full"]]) <- c("GRIN2B","CNTN5")
# GOI[["Ex"]] <- GOI[["Full"]]
# GOI[["Full"]][,"AREA"] <- c( 8.178442, 4.0526161 )
# GOI[["Full"]][,"BEST_P"] <- c( 2.346e-08, 1.650e-06 )
# GOI[["Ex"]][,"AREA"] <- c( 2.139639, 0.3252284 )
# GOI[["Ex"]][,"BEST_P"] <- c( 7.554e-06, 2.136e-02 )

## Compile Stats on Genes of Interest
FULL.GOI <- cbind( GOI$Full, A.f[match(rownames(GOI$Full),A.f$Gene),c("G_LEN","N_VAR","logG_LEN","logN_VAR")], -log10(GOI$Full[,"BEST_P"]) )
colnames(FULL.GOI)[ncol(FULL.GOI)] <- "-logBEST_P"
EXON.GOI <- cbind( GOI$Ex, A.e[match(rownames(GOI$Ex),A.e$Gene),c("G_LEN","N_VAR","logG_LEN","logN_VAR")], -log10(GOI$Ex[,"BEST_P"]) )
colnames(EXON.GOI)[ncol(EXON.GOI)] <- "-logBEST_P"

################################################
## FCT: Plot Metrics vs One Another
Plot_Variables <- function( X_Lab, Y_Lab, Color, tag ) {
	if ( tag=="Full" ) {
		DATA <- PERM.f
		DATA.GOI <- FULL.GOI
	}else{
		DATA <- PERM.e
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
	# ## Raw Gene Lengths
	# jpeg( paste(PathToPlot,True_Tag,"_vsNull.1-Scatter.",tag,".jpeg",sep=""),height=1200,width=1800,pointsize=30)
	# par(mfrow=c(2,3))
	# Plot_Variables( "G_LEN","N_VAR", COLS[1], tag )
	# Plot_Variables( "G_LEN","-logBEST_P", COLS[2], tag )
	# Plot_Variables( "G_LEN","AREA", COLS[3], tag )

	# Plot_Variables( "-logBEST_P","AREA", COLS[4], tag )
	# Plot_Variables( "N_VAR","-logBEST_P", COLS[2], tag )
	# Plot_Variables( "N_VAR","AREA", COLS[3], tag )
	# dev.off()

	## Log Transformed Gene Lengths
	jpeg( paste(PathToPlot,True_Tag,"_vsNull.1-Scatter.log.",tag,".jpeg",sep=""),height=1200,width=1800,pointsize=30)
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
## FCT: Null Distribution of AREA
Null_Area <- function( AREA.vec, GOI, Color, Metric, tag ) {
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
	XLIM <- c( RNG[1], 1+max(RNG[2],GOI[[tag]][,Metric],na.rm=T) )
	BRKS <- seq( floor(RNG[1]), ceiling(RNG[2]), .1 )
	HIST <- hist(AREA.n,freq=F,breaks=BRKS,col=Color,main=paste("Distribution of",Metric,"Values -",tag),xlab=Metric, xlim=XLIM )
	abline( v=seq(floor(RNG[1]),ceiling(XLIM[2]),1),h=seq(0,2,.2),lty=3,col="grey50",lwd=1 )
	HIST <- hist(AREA.n,freq=F,breaks=BRKS,col=Color,add=T)
	points( XX, dJohnson(XX,O), col="chocolate2",type="l",lwd=3) # plot(function(x)dJohnson(x,O), -20, 20, add=T, col="blue")
	text( quantile(AREA.n,.99,na.rm=T), max(HIST$density), label=paste(names(Moments),round(Moments,3),sep=" = ",collapse="\n" ),pos=1 )
	for ( g in 1:nrow(GOI[[tag]]) ) { gene_area <- GOI[[tag]][g,Metric]
	# for ( g in 1:nrow(GOI[[tag]]) ) { gene_area <- GOI[[tag]][g,"MAX_E"]
		if ( is.na(gene_area) ) { next }
		arrows( gene_area, .1, gene_area, 0, lwd=3, col="firebrick1" )
		P.prm <- (1+length(which(AREA.n>gene_area))) / (length(AREA.n)+1)
		# P.dis <- pJohnson(gene_area,O,lower.tail=F)
		# text( gene_area, .12, paste(rownames(GOI[[tag]])[g],"-",round(GOI[[tag]][g,"AREA"],4),"\np.dis =",formatC(P.dis,2,format="e"),"\np.prm =",formatC(P.prm,2,format="e")), srt=90, pos=4, col="firebrick1" )
		text( gene_area, .12, paste(rownames(GOI[[tag]])[g],"-",round(GOI[[tag]][g,Metric],4),"\np.prm =",formatC(P.prm,2,format="e")), srt=90, pos=4, col="firebrick1" )
	}
}

## AREA vs Null Distribution
GOI.temp <- lapply( GOI, function(x) x[c(GOI.names.custom,"CYP2E1"),] )
COLS <- c("chartreuse1","dodgerblue2","tomato2")
COLS <- c("steelblue3","slateblue3","tomato2")
jpeg( paste(PathToPlot,True_Tag,"_vsNull.2-Area_Distrib.jpeg",sep=""),height=1200,width=1600,pointsize=24)
par(mfrow=c(2,1))
Null_Area( PERM.f$AREA, GOI.temp, COLS[1], "AREA", "Full" )
Null_Area( PERM.e$AREA, GOI.temp, COLS[2], "AREA", "Ex" )
dev.off()

## MAX_E vs Null Distribution
GOI.temp <- lapply( GOI, function(x) x[c(GOI.names.custom,"CYP2E1"),] )
COLS <- c("steelblue3","slateblue3","tomato2")
jpeg( paste(PathToPlot,True_Tag,"_vsNull.2-MaxE_Distrib.jpeg",sep=""),height=1200,width=1600,pointsize=24)
par(mfrow=c(2,1))
Null_Area( PERM.f$MAX_E, GOI.temp, COLS[1], "MAX_E", "Full" )
Null_Area( PERM.e$MAX_E, GOI.temp, COLS[2], "MAX_E", "Ex" )
dev.off()


################################################
## Get Null Distribution for Top N Genes
 # (for EACH permutation)
N.top <- 1000
N.top <- 500
N.top <- 100
N.top <- 10

## Compile the Top Ranking Area Values for True & Permuted Data
Calc_TOPA <- function( N.top, Perm_or_True ) {
	which_col <- grep("AREA",colnames(PERM.f))
	if ( Perm_or_True=="Perm" | Perm_or_True=="Perm" ) {
		## Leave one out from Permuted Samples
		Which_Out <- sample(1:100,1)
		 # Calculate the AREA of the top N enriched Genes for ACTUAL data
		TOP.A.af <- sort(PERM.A$FULL[[Which_Out]][,which_col],decreasing=T)[1:N.top]
		TOP.A.ae <- sort(PERM.A$EXON[[Which_Out]][,which_col],decreasing=T)[1:N.top]
		 # Calculate the AREA of the top N enriched Genes for EACH PERMUTATION
		TOP.A.f <- matrix( unlist(lapply( PERM.A$FULL[-Which_Out], function(x) sort(x[,which_col],decreasing=T)[1:N.top] )), byrow=T,ncol=N.top )
		TOP.A.e <- matrix( unlist(lapply( PERM.A$EXON[-Which_Out], function(x) sort(x[,which_col],decreasing=T)[1:N.top] )), byrow=T,ncol=N.top )
	}else{
		## Use True Data for Example Plot
		 # Calculate the AREA of the top N enriched Genes for ACTUAL data
		TOP.A.af <- sort(A.f[,"AREA"],decreasing=T)[1:N.top]
		TOP.A.ae <- sort(A.e[,"AREA"],decreasing=T)[1:N.top]
		 # Calculate the AREA of the top N enriched Genes for EACH PERMUTATION
		TOP.A.f <- matrix( unlist(lapply( PERM.A$FULL, function(x) sort(x[,which_col],decreasing=T)[1:N.top] )), byrow=T,ncol=N.top )
		TOP.A.e <- matrix( unlist(lapply( PERM.A$EXON, function(x) sort(x[,which_col],decreasing=T)[1:N.top] )), byrow=T,ncol=N.top )
	}

	## Calculate Probability of Data at each Gene Rank
	TOP.A.pf <- unlist(lapply( 1:N.top, function(x) (1+length(which(TOP.A.af[x] < TOP.A.f[,x])))/(1+nrow(TOP.A.f)) ))
	TOP.A.pe <- unlist(lapply( 1:N.top, function(x) (1+length(which(TOP.A.ae[x] < TOP.A.e[,x])))/(1+nrow(TOP.A.e)) ))

	## Compile Outputs
	# TOP.A.arr <- data.frame( RANK=1:N.top, AF=TOP.A.af,PF=TOP.A.pf, AE=TOP.A.ae,PE=TOP.A.pe )
	# head( TOP.A.arr, min(N.top,30) )
	# head( TOP.A.arr[order(TOP.A.arr$PF),], min(N.top,30) )
	# head( TOP.A.arr[order(TOP.A.arr$PE),], min(N.top,30) )	
	TOP <- list(A.af=TOP.A.af,A.ae=TOP.A.ae,A.f=TOP.A.f,A.e=TOP.A.e,A.pf=TOP.A.pf,A.pe=TOP.A.pe )
	return( TOP )
}

## Plot the Top Ranking Area Values for True & Permuted Data
Plot_TOPA <- function( TOP.A.true, TOP.A.perm, TOP.A.pval ) {
	N.top <- length(TOP.A.true)
	XLIM <- c(1,N.top)
	YLIM <- c( 0, max(TOP.A.true,TOP.A.perm) )
	plot( 0,0,type="n",xlim=XLIM,ylim=YLIM,xlab="Rank",ylab="Area",main="Ranked Area: Actual vs Null")
	abline( h=seq(0,30,1),v=seq(0,N.top,N.top/10), lty=3,col="grey20",lwd=1)
	abline( h=YLIM[2]*seq(0,1,.1),v=seq(0,N.top,N.top/10), lty=c(1,rep(3,4),1,rep(3,4),1),col="dodgerblue1",lwd=2)
	lapply( 1:N.top, function(x) points(rep(x,nrow(TOP.A.perm)),TOP.A.perm[,x],pch=20,col="grey50") )
	points( 1:N.top,TOP.A.true,col="firebrick1",type="l",lwd=4 )
	axis( 4, at=seq(YLIM[1],YLIM[2],length.out=11),label=seq(0,1,.1),col="dodgerblue3" )
	points( 1:N.top,YLIM[2]*TOP.A.pval,col="dodgerblue2",type="l",lwd=4 )
}

## Run Through Various Top N Genes
N.tops <- c( 10, 50, 100, 500, 1000 )
for ( N.top in N.tops ) {
	TOP <- Calc_TOPA( N.top, "True" )
	jpeg( paste(PathToPlot,True_Tag,"_vsNull-Top",N.top,".jpeg",sep=""),height=1600,width=2400,pointsize=24)
	par(mfrow=c(2,1))
	Plot_TOPA( TOP$A.af, TOP$A.f, TOP$A.pf )
	Plot_TOPA( TOP$A.ae, TOP$A.e, TOP$A.pe )
	dev.off()	
}

################################################
## Correct for Gene Size/Number of Variants
MOD.G_LEN.f <- lm( AREA ~ logG_LEN, data=PERM.f )
MOD.N_VAR.f <- lm( AREA ~ logN_VAR, data=PERM.f )
MOD.G_LEN.e <- lm( AREA ~ logG_LEN, data=PERM.e )
MOD.N_VAR.e <- lm( AREA ~ logN_VAR, data=PERM.e )

GOI.res.nvar <- GOI
GOI.res.nvar$Full$AREA <- GOI.res.nvar$Full$AREA - predict( MOD.N_VAR.f, newdata=FULL.GOI )
GOI.res.nvar$Ex$AREA <- GOI.res.nvar$Ex$AREA - predict( MOD.N_VAR.e, newdata=EXON.GOI )
GOI.res.glen <- GOI
GOI.res.glen$Full$AREA <- GOI.res.glen$Full$AREA - predict( MOD.G_LEN.f, newdata=FULL.GOI )
GOI.res.glen$Ex$AREA <- GOI.res.glen$Ex$AREA - predict( MOD.G_LEN.e, newdata=EXON.GOI )

COLS <- c("chartreuse1","dodgerblue2","tomato2")
 # Null Distrib of Resids vs N_VARS
jpeg( paste(PathToPlot,True_Tag,"_vsNull.2-Area_Distrib.Resid.N_VAR.jpeg",sep=""),height=1200,width=1600,pointsize=24)
par(mfrow=c(2,1))
Null_Area( resid(MOD.N_VAR.f), GOI.res.nvar, COLS[1], "Full" )
Null_Area( resid(MOD.N_VAR.e), GOI.res.nvar, COLS[2], "Ex" )
dev.off()
 # Null Distrib of Resids vs G_LENS
jpeg( paste(PathToPlot,True_Tag,"_vsNull.2-Area_Distrib.Resid.G_LEN.jpeg",sep=""),height=1200,width=1600,pointsize=24)
par(mfrow=c(2,1))
Null_Area( resid(MOD.G_LEN.f), GOI.res.glen, COLS[1], "Full" )
Null_Area( resid(MOD.G_LEN.e), GOI.res.glen, COLS[2], "Ex" )
dev.off()


## Quantiles ##
QUANTS <- 10^(-8:0)
quantile( PERM.f$AREA, 1-QUANTS, na.rm=T )


QARRS <- list()

which_col <- grep("AREA",colnames(PERM.f))
QARRS$AREA <- matrix( unlist(lapply( PERM.d$FULL, function(x) quantile( x[,which_col], 1-QUANTS, na.rm=T ) )), byrow=T,ncol=length(QUANTS) )

which_col <- grep("BEST_P",colnames(PERM.f))
QARRS$BEST_P <- matrix( unlist(lapply( PERM.d$FULL, function(x) quantile( x[,which_col], QUANTS, na.rm=T ) )), byrow=T,ncol=length(QUANTS) )


lapply( QARRS, function(x) apply(x,2,mean) )
lapply( QARRS, function(x) apply(x,2,sd) )





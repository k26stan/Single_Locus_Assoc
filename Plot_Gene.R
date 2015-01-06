P <- read.table( "20150104_LT8_DEL_MNe_MN_DAS_BL_MN_PC1_PC2/20150104_LT8_DEL_MNe_MN_DAS_BL_MN_PC1_PC2.P", sep="\t", header=T )
TAB <- read.table( "/home/kstandis/RV/Tools/GG-Gene_Names_DB.txt", sep="\t",header=T, comment.char="",quote="")

GRIN <- which( TAB$hg19.kgXref.geneSymbol=="GRIN2B" )

P.GRIN.which <- which( P$CHR==12 & P$BP>13714409-5e3 & P$BP<14133022+5e3 )
P.GRIN <- P[P.GRIN.which,]



GG.g <- TAB[GRIN,]
MRG.g <- P.GRIN
t <- 1

## Compile Parameters
NUM_TRANS <- nrow(GG.g)
NUM_VARS <- nrow(MRG.g)
COLS <- c("cadetblue1","dodgerblue1","chartreuse1","mediumpurple2","sienna1","firebrick1") # Orig
COLS <- c("slateblue3","cadetblue1","steelblue2","firebrick2","goldenrod2","springgreen3")
X_LIM <- c(min(GG.g$hg19.knownGene.txStart),max(GG.g$hg19.knownGene.txEnd))
Y_LIM <- c( -2, 8 )

jpeg("/home/kstandis/20150106_TEMP.jpeg", height=1600,width=2400, pointsize=42)
plot(0,0,type="n", xlim=X_LIM,ylim=Y_LIM, xlab=paste("Chromosome 12 Position"),ylab="-log10(p)", main="GRIN2B" )
abline( h=-log10(5e-8), lty=2, lwd=3, col="firebrick2" )
abline( h=seq(0,10,1), lty=2, col="grey50" )
abline( v=seq(13e6,15e6,1e5), lty=2, col="grey50" )
## Get Outer Coordinates
TX_RNG <- c( GG.g$hg19.knownGene.txStart[t], GG.g$hg19.knownGene.txEnd[t] )
CD_RNG <- c( GG.g$hg19.knownGene.cdsStart[t], GG.g$hg19.knownGene.cdsEnd[t] )

## Get Exon Coordinates
EX_B <- as.numeric( strsplit( as.character(GG.g$hg19.knownGene.exonStarts), "," )[[t]] )
EX_E <- as.numeric( strsplit( as.character(GG.g$hg19.knownGene.exonEnds), "," )[[t]] )

## Plot this Shiz
Y <- -1
arrows( TX_RNG[1],Y,TX_RNG[2],Y, code=3,angle=90, lwd=5,col=COLS[1] )
abline( h=Y, col="grey50",lwd=1 )
polygon( CD_RNG[c(1,2,2,1)],c(-.1,-.1,.1,.1)+Y, col=COLS[2],lwd=1 )
for ( e in 1:length(EX_B) ) {
	X_COORDS <- c( EX_B[e],EX_E[e],EX_E[e],EX_B[e] )
	polygon( X_COORDS,c(-.2,-.2,.2,.2)+Y, col=COLS[3],lwd=1 )	
}
# arrows( TX_RNG[1],Y,TX_RNG[2],Y, code=2,angle=35, lwd=5,col=COLS[1] )
arrows( TX_RNG[1],Y,TX_RNG[2],Y, code=3,angle=90,length=0, lwd=5,col=COLS[1] )

points( P.GRIN$BP, -log10(P.GRIN$P), pch=20, col="deepskyblue2" )


dev.off()

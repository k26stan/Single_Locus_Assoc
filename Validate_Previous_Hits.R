## Plot GWAS Results for Previous GWAS Hits ##
## Attempt to Validate Previous Results ##
## September 22, 2015 ##
## Kristopher Standish ##

TAB.l <- read.table( "Data/Burn/Results/Single_Locus/Previous_GWAS_Candidates_Compiled.P",sep="\t",header=T, stringsAsFactors=F )
TAB <- TAB.l[ order(TAB.l$POS), ]
TAB <- TAB[ order(TAB$CHR), ]

COLS <- c("firebrick2","gold2","chartreuse2","dodgerblue2")
png( "/Users/kstandis/Data/Burn/Results/Single_Locus/20151001_PreviousHits/1_PrevHits.png",height=800,width=1600,pointsize=24)
par(mar=c(6,4,4,4))
plot( 0,0,type="n", xlim=c(1,nrow(TAB)), ylim=c(0,4), xaxt="n",xlab="",ylab="-log10(p)",main="ART3001 Response Association of Previous Candidate SNPs" )
axis(1, at=1:nrow(TAB), TAB$ID, las=2 )
abline( h=-log10(.05/(4*nrow(TAB))),lty=2,col="tomato2" )
abline( h=-log10(.05/nrow(TAB)),lty=2,col="chocolate1" )
abline( h=-log10(.05),lty=2,col="gold2" )
for ( col in 4:7 ) { points( 1:nrow(TAB), -log10(TAB[,col]), col=COLS[col-3],pch="+" ) }
legend( "topleft", pch="+",col=1:4,legend=colnames(TAB)[4:7], bg="white" )
dev.off()

png( "/Users/kstandis/Data/Burn/Results/Single_Locus/20151001_PreviousHits/1_PrevHits_QQ.png",height=1000,width=1000,pointsize=26)
plot( 0,0,type="n", xlim=c(0,4),ylim=c(0,4), xlab="-log10(Exp)",ylab="-log10(Obs)",main="ART3001 Response Association of Previous Candidate SNPs")
abline( 0,1 )
abline( h=seq(0,10,1),lty=3,col="grey50" )
abline( v=seq(0,10,1),lty=3,col="grey50" )
abline( h=-log10(.05/(4*nrow(TAB))),lty=2,col="tomato2" )
abline( h=-log10(.05/nrow(TAB)),lty=2,col="chocolate1" )
abline( h=-log10(.05),lty=2,col="gold2" )
for ( i in 1:100 ) { points( -log10(1:nrow(TAB)/nrow(TAB)), -log10(sort(runif(nrow(TAB),0,1))), pch=20,col="grey75" ) }
for ( col in 4:7 ) { points( -log10(1:nrow(TAB)/nrow(TAB)), -log10(sort(TAB[,col])), pch="+",col=COLS[col-3] ) }
legend( "topleft", pch="+",col=1:4,legend=colnames(TAB)[4:7], bg="white" )
dev.off()


## Rscript to Make Manhattan Plots for GWAS ##
## May 13, 2014 ##
## Kristopher Standish ##

## Name: Manhat_Plot.R

## Usage: Rscript Manhat_Plot.R ${P_FILE}

LINE <- commandArgs(trailingOnly = TRUE)
# LINE <- c("/projects/janssen/ASSOCIATION/20150104_LT8_DEL_MNe_MN_DAS_BL_MN_PC1_PC2/20150104_LT8_DEL_MNe_MN_DAS_BL_MN_PC1_PC2.P","LT8_DEL_MNe_MN.txt","DAS_BL_MN_PC1_PC2")
PathToHC <- LINE[1]
Split <- strsplit(PathToHC,"/")
PathToAssoc <- paste(Split[[1]][1:(length(Split[[1]])-1)],collapse="/")
Pheno <- LINE[2]
Pheno <- gsub(".txt","",Pheno)
Covars <- LINE[3]

#########################################################
## LOAD DATA ############################################
#########################################################

## Load Data
start <- proc.time()
HC <- read.table(PathToHC,sep="\t",header=T,colClass=c("factor","character","numeric","numeric")) ; proc.time()-start

## Chromosome Length Table
lens <- read.table(file="/projects/janssen/clinical/Chrom_Lengths.txt", sep="\t", header=T, stringsAsFactors=F)

####################################
## MANHATTAN PLOTS #################
####################################

## Set up variables for plotting/naming/etc...
SNP_Cnt <- nrow(HC)
Set_Range <- -log10( min(HC$P,na.rm=T) ) + 1
Set_Names <- "HC"
COLS <- rep(c("chartreuse","dodgerblue"),13)
PLOT_THRSH <- 2

## Loop it to plot multiple plots
jpeg( paste(PathToAssoc, "/Manhattan_",Pheno,"_",Covars,".jpeg", sep=""), height=1000, width=2000, pointsize=30)
Data <- HC[which(HC$P<10^-PLOT_THRSH),]
nrow(Data)
plot(0,0, type="n", xlim=c(0,3e9), ylim=c(PLOT_THRSH,ceiling(Set_Range)), xlab="Position", ylab="-log10(p-value)", main=paste("Single-Locus Association - ", Pheno, " (cov: ", Covars, ")", sep=""), xaxt="n" )
 # Lines
abline( h=seq(0,50,1), lty=3, col="grey50" )
abline(-log10(.05/mean(SNP_Cnt)),0, col="black", lty=2,lwd=2)
text(0,-log10(.05/mean(SNP_Cnt))+.13, labels="Bonferroni",col="black",pos=4)
abline(-log10(5e-8),0, col="firebrick3", lty=2,lwd=2)
text(0,-log10(5e-8)+.13, labels="p=5e-8", col="firebrick3",pos=4)
 # Axis
axis(1, at=apply(data.frame(lens$ADD[1:23],lens$ADD[2:24]),1,mean), labels=lens[1:23,1], tick=F)
 # Legend
legend(.9*3e9,ceiling(Set_Range), legend=c("Odd-Chr","Even-Chr"), pch=20, col=COLS[1:2] )
 # Data
for (chrom in 1:24) {
  points((lens[chrom,3]+Data$BP[which(Data$CHR==chrom)]), -log10(Data$P[which(Data$CHR==chrom)]), col=COLS[chrom], pch=20, cex=.5)
  arrows(lens[chrom,3], -4, lens[chrom,3], PLOT_THRSH, col="black", length=0, lwd=3 )
}
dev.off()

###########################
## Q-Q PLOTS ##############
###########################

# Set up variables for plotting/naming Q-Q Plots...
QQ_Cols <- c("black", "chartreuse")

jpeg( paste(PathToAssoc, "/QQ_Plot_",Pheno,"_",Covars,".jpeg", sep=""), height=1000, width=1000, pointsize=30 )
plot( 0,0, type="n", main=paste("QQ Plot -", Pheno, " (cov: ", Covars,")", sep=""), xlab="Expected (-logP)", ylab="Observed (-logP)", xlim=c(0,ceiling(Set_Range)), ylim=c(0,ceiling(Set_Range)), las=1, xaxs="i", yaxs="i", bty="l")
 # Lines
abline( h=seq(0,50,1), lty=3, col="grey50" )
abline( v=seq(0,50,1), lty=3, col="grey50" )
abline(-log10(.05/mean(SNP_Cnt)),0, col="black", lty=2,lwd=2)
text(0,-log10(.05/mean(SNP_Cnt))+.13, labels="Bonferroni",col="black",pos=4)
abline(-log10(5e-8),0, col="firebrick3", lty=2,lwd=2)
text(0,-log10(5e-8)+.13, labels="p=5e-8", col="firebrick3",pos=4)
abline( 0,1, lty=1, lwd=2, col=QQ_Cols[1] )
 # Specify Data
Data <- HC
lobs <- -log10(sort(Data$P))
lexp <- -log10(1:length(lobs)/(length(lobs)+1))
WHICH <- which(lobs>PLOT_THRSH)

 # Plot Data
# points(c(0,ceiling(max(lobs))), c(0,ceiling(max(lobs))), col=QQ_Cols[1], type="l", lwd=2)
points(lexp[WHICH], lobs[WHICH], col=QQ_Cols[2], pch=20) 
dev.off()

###########################
## COMPILE VARIANTS #######
###########################

# Which Variants Meet Genome-Wide Significance Thresholds
PRINT_THRSH <- 1e-4
HC_CND <- HC[ which(HC$P <= PRINT_THRSH),]

write.table(HC_CND, paste(PathToAssoc, "/CND_",Pheno,"_",Covars,".txt", sep=""), sep="\t", row.names=F, col.names=T, quote=F)


###########################
## END OF DOC #############
###########################

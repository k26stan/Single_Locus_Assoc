# Make Covariate File w/ EigenVectors #

# PCA_Plot.R <Path/To/EigenVectors> <Path/To/Covariate_Table> <Set>

LINE <- commandArgs(trailingOnly = TRUE)
# LINE <- c("/projects/janssen/scripts/AS-ASSOCIATION/20140514pca_HT_SEX_AGE/20140514pca-BGI_HT.eigenvec","/projects/janssen/scripts/AS-ASSOCIATION/COV_Tables/All_BigCall.txt")
# LINE <- c("/projects/janssen/scripts/AS-ASSOCIATION/20140514_HT_SEX_AGE_PC1_PC2_PC3/20140514-HC_FULL_HT.eigenvec","/projects/janssen/scripts/AS-ASSOCIATION/COV_Tables/All_BigCall.txt")
PathToHC <- LINE[1]
PathToCov <- LINE[2]
Split <- strsplit(PathToHC,"/")
PathToAssoc <- LINE[3] #paste(Split[[1]][1:(length(Split[[1]])-1)],collapse="/")

#################################################
## LOAD DATA ####################################

## Load Tables
HC <- read.table(PathToHC,header=T)
COV <- read.table(PathToCov,sep="\t",header=T)

## Merge Eigenvector Tables
HC <- HC[,c(1,3:ncol(HC))]
MRG_HC <- merge(x=COV,y=HC,by="FID")

#################################################
## PLOT #########################################

## Set up some Variables
RACE_COL <- c("tomato3","sienna1","steelblue2","goldenrod3","slateblue3","springgreen3")
ETHN_PCH <- c(17,19) #c("+","x")

## Plot PC1 vs PC2 vs PC3 for HC data sets (2x2, Legend in box 2)
LIM_H1 <- c(min(HC$PC1),max(HC$PC1))
LIM_H2 <- c(min(HC$PC2),max(HC$PC2))
LIM_H3 <- c(min(HC$PC3),max(HC$PC3))
jpeg(paste(PathToAssoc,"/PCA_1v2v3.jpeg",sep=""), width=2000, height=2000, pointsize=20)
par(mfrow=c(2,2))
plot(MRG_HC$PC1,MRG_HC$PC3, col=RACE_COL[factor(MRG_HC$RACE)],pch=ETHN_PCH[factor(MRG_HC$ETHN)],main="PC1 vs PC3 - HC Calls", xlab="PC1", ylab="PC3", xlim=LIM_H1, ylim=LIM_H3)
plot(0,0,type="n",xlim=c(0,1),ylim=c(0,1),xaxt="n",yaxt="n")
legend(0,.5, legend=c("Hisp/Lat","Not Hisp/Lat"), col="black", pch=ETHN_PCH, title="Ethnicity (Self)", cex=2.2)
legend(0,1, legend=levels(MRG_HC$RACE), col=RACE_COL, pch=20, title="Race (Self)", cex=2.2)
plot(MRG_HC$PC1,MRG_HC$PC2, col=RACE_COL[factor(MRG_HC$RACE)],pch=ETHN_PCH[factor(MRG_HC$ETHN)],main="PC1 vs PC2 - HC Calls", xlab="PC1", ylab="PC2", xlim=LIM_H1, ylim=LIM_H2)
plot(MRG_HC$PC3,MRG_HC$PC2, col=RACE_COL[factor(MRG_HC$RACE)],pch=ETHN_PCH[factor(MRG_HC$ETHN)],main="PC3 vs PC2 - HC Calls", xlab="PC3", ylab="PC2", xlim=LIM_H3, ylim=LIM_H2)
dev.off()

## Scatter Plot (3D) - PC1 vs PC2 vs PC3
library(scatterplot3d)

color=RACE_COL[factor(MRG_HC$RACE)]
jpeg(paste(PathToAssoc,"/PCA_HC_Scatter3D.jpeg",sep=""), width=2000, height=2000, pointsize=20)
scatterplot3d(MRG_HC$PC1,MRG_HC$PC2,MRG_HC$PC3, color, pch=ETHN_PCH[factor(MRG_HC$ETHN)])
dev.off()



###################################
## END OF DOC #####################
###################################



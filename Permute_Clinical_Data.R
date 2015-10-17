# Permute Clinical Data for Running GWAS Permutations #

# Perm_Clin_Data.R <Path/To/Home_Dir> <Path/To/Phenotype_Table> <Path/To/Covariate_Table> <Num_Perm>
# Rscript ${PERMUTE_CLIN_DAT} ${HOME_DIR} ${PHENO_FILE} ${COV_FILE} ${NUM_PERM}

LINE <- commandArgs(trailingOnly = TRUE)
# LINE <- c("/projects/janssen/ASSOCIATION","/projects/janssen/ASSOCIATION/PH-PHENOTYPES/LT8_rTJC_BL_MN.txt","/projects/janssen/ASSOCIATION/PH-PHENOTYPES/20151015_Full_Table.txt","RF_ACPA,logDIS_DUR,PC1,PC2","10")
PathToHome <- LINE[1]
PathToPheno <- LINE[2]
PathToCov <- LINE[3]
Covs_List <- LINE[4]
Num_Perm <- LINE[5]

## Load Tables
PH <- read.table(PathToPheno,header=T)
COV <- read.table(PathToCov,sep="\t",header=T)

## Separate Covariate List into Colnames
Covs_Colnames <- strsplit(Covs_List,"QQQ")[[1]]
if ( any(grepl("PC",Covs_Colnames)) ) { Covs_Colnames <- Covs_Colnames[-grep("PC",Covs_Colnames)] }
COV.temp <- COV[,c(1,which(colnames(COV)%in%Covs_Colnames))]

## Merge Tables
MRG <- merge(x=PH,y=COV.temp,by.x=colnames(PH)[1],by.y=colnames(COV.temp)[1] )
colnames(MRG)[1:2] <- c("FID","IID")

## Shuffle Sample Names
for ( n in 1:Num_Perm ) { 
	Shuffled <- sample( as.character(MRG[,"FID"]) )
	MRG[,"FID"] <- MRG[,"IID"] <- Shuffled

	## Write Tables
	PathToOut <- paste(PathToHome,"/Permuted_Phenos/Perm_Pheno.",n,".txt",sep="")
	write.table(MRG[,1:3],PathToOut,sep="\t",row.names=F,col.names=T,quote=F)
	PathToOut <- paste(PathToHome,"/Permuted_Phenos/Perm_Cov.",n,".txt",sep="")
	write.table(MRG[,-3],PathToOut,sep="\t",row.names=F,col.names=T,quote=F)

}









######### END OF DOC ###########
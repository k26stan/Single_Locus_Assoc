# Make Covariate File w/ EigenVectors #

# Make_Cov_Tab.R <Path/To/EigenVectors> <Path/To/Covariate_Table> <Covariate_List> <Path/To/New_Cov_Table>

LINE <- commandArgs(trailingOnly = TRUE)
# LINE <- c("/projects/janssen/ASSOCIATION/../ASSOCIATION/EIGEN/HC_FULL.eigenvec","/projects/janssen/ASSOCIATION/../ASSOCIATION/PH-PHENOTYPES/20151015_Full_Table.txt","RF_ACPA,logDIS_DUR","/projects/janssen/ASSOCIATION/20151015_LT8_DAS_BL_MN_RF_ACPA_logDIS_DUR/Cov_w_PCs.txt")
PathToVec <- LINE[1]
PathToCov <- LINE[2]
Covs_List <- LINE[3]
PathToNewCov <- LINE[4]

## Load Tables
VEC <- read.table(PathToVec,header=T)
COV <- read.table(PathToCov,sep="\t",header=T)

## Separate Covariate List into Colnames
Covs_Colnames <- strsplit(Covs_List,",")[[1]]
if ( any(grepl("PC",Covs_Colnames)) ) { Covs_Colnames <- Covs_Colnames[-grep("PC",Covs_Colnames)] }

COV.temp <- COV[,c(1,1,which(colnames(COV)%in%Covs_Colnames))]
for ( col in 3:ncol(COV.temp) ) {
	if ( length(unique(COV.temp[,col]))==2 ) {
		COV.temp[,col] <- as.numeric(as.factor(COV.temp[,col]))
	}
}
VEC.temp <- VEC[,c(1,3:ncol(VEC))]

## Merge Tables
# MRG <- merge(x=COV,y=VEC,by.x=colnames(COV)[1],by.y="FID" )
MRG <- merge(x=COV.temp,y=VEC.temp,by.x=colnames(COV.temp)[1],by.y=colnames(VEC.temp)[1] )
colnames(MRG)[1:2] <- c("FID","IID")

## Write Table
write.table(MRG,PathToNewCov,sep="\t",row.names=F,col.names=T,quote=F)

######### END OF DOC ###########

## Rscript to Compile Info for Janssen ##
## June 23, 2014 ##
## Kristopher Standish ##

## Name: Compile_Outs.R

## Usage: Rscript ${COMPILE_OUTS_R} ${CND_FILE} ${CND_FILE%%txt}hwe ${CND_FILE%%txt}adj ${CND_GENES}

LINE <- commandArgs(trailingOnly = TRUE)
# LINE <- paste("/projects/janssen/ASSOCIATION/20140623_LT8_DEL_MNe_MN_DAS_BL_MN_AGE_SEX_AGESEX_BMI_COUN_PC1_PC2/CND_LT8_DEL_MNe_MN_DAS_BL_MN_AGE_SEX_AGESEX_BMI_COUN_PC1_PC2.",c("txt","hwe","adj","Gene.txt","Annot.txt"),sep="")
# LINE <- paste("/Users/kstandis/SDSC_PROJ/ASSOCIATION/20140625_LT8_DEL_MNe_MN_DAS_BL_MN_AGESEX_logBMI_COUN_PC1_PC2/CND_LT8_DEL_MNe_MN_DAS_BL_MN_AGESEX_logBMI_COUN_PC1_PC2.",c("txt","hwe","adj","Gene.txt","Annot.txt"),sep="")
# LINE <- paste("/Users/kstandis/SDSC_PROJ/ASSOCIATION/20140625_LT8_DEL_28_BL_DAS_BL_AGESEX_logBMI_COUN_PC1_PC2/CND_LT8_DEL_28_BL_DAS_BL_AGESEX_logBMI_COUN_PC1_PC2.",c("txt","hwe","adj","Gene.txt","Annot.txt"),sep="")
# LINE <- paste("/projects/janssen/ASSOCIATION/20140625_LT8_DEL_28_BL_DAS_BL_AGESEX_logBMI_COUN_PC1_PC2/CND_LT8_DEL_28_BL_DAS_BL_AGESEX_logBMI_COUN_PC1_PC2.",c("txt","hwe","adj","frqx","pv","Gene.txt","Annot.txt"),sep="")
# LINE <- paste("/projects/janssen/ASSOCIATION/20140926_LT8_DEL_MNe_MN_DAS_BL_MN_AGESEX_logBMI_COUN_PC1_PC2/CND_LT8_DEL_MNe_MN_DAS_BL_MN_AGESEX_logBMI_COUN_PC1_PC2.",c("txt","hwe","adj","frqx","pv","Gene.txt","Annot.txt"),sep="")
# LINE <- paste("/projects/janssen/ASSOCIATION/20140926b_LT8_DEL_MNe_MN_DAS_BL_MN_AGESEX_logBMI_COUN_PC1_PC2/CND_LT8_DEL_MNe_MN_DAS_BL_MN_AGESEX_logBMI_COUN_PC1_PC2.",c("txt","hwe","adj","frqx","pv","Gene.txt","Annot.txt"),sep="")
PathToTxt <- LINE[1]
PathToHW <- LINE[2]
PathToADJ <- LINE[3]
PathToFRQ <- LINE[4]
PathToPV <- LINE[5]
PathToGene <- LINE[6]
PathToAnnot <- LINE[7]
Split <- strsplit(PathToTxt,"/")
PathToAssoc <- paste(Split[[1]][1:(length(Split[[1]])-1)],collapse="/")
PathToWrite <- gsub("txt","compiled",PathToTxt)

## Pull out Date, Pheno, Covariates
Dir_Name <- strsplit(PathToAssoc,"/")[[1]][length(strsplit(PathToAssoc,"/")[[1]])]
Date <- strsplit(Dir_Name,"_LT8_")[[1]][1]
Pheno <- paste( strsplit( strsplit(Dir_Name,"_LT8_")[[1]][2], "_" )[[1]][1:3], collapse="_" )
Covs <- paste( strsplit( strsplit(Dir_Name,"_LT8_")[[1]][2], "_" )[[1]][-(1:3)], collapse="_" )
Covs <- gsub("_","+",Covs) ; Covs <- gsub("DAS+BL+MN","DAS_BL_MN",Covs,fixed=T)  ; Covs <- gsub("DAS+BL","DAS_BL",Covs,fixed=T)

#########################################################
## LOAD DATA ############################################
#########################################################

## Load Data
TXT <- read.table(PathToTxt,sep="\t",header=T,colClass=c("factor","character","numeric","numeric"))
HW <- read.table(PathToHW,header=T)
ADJ <- read.table(PathToADJ,header=T)
FRQ <- read.table(PathToFRQ,header=T,sep="\t")
PV <- read.table(PathToPV,header=T)
GENE <- read.table(PathToGene,sep="\t",header=T, fill=T, quote="", comment.char="")
# ANNOT <- read.table(PathToAnnot,sep="\t",header=T,quote="",comment.char="",fill=T)

####################################
## PULL OUT COLUMNS OF INTEREST ####
####################################

## TXT Array
 # Use All
CHR_POS <- paste("chr",TXT$CHR,":",TXT$BP-1,"-",TXT$BP,sep="")
TXT_2 <- data.frame(TXT,CHR_POS)
names(TXT_2)[4] <- "P_Assoc"

## HW Array
HW_2 <- HW[,c("SNP","GENO","P")] # HW_2 <- HW[,c(2,4,5,6,9)]
names(HW_2)[which(names(HW_2)=="P")] <- "P_HW"

## ADJ Array
ADJ_2 <- ADJ[,c("SNP","QQ","FDR_BH","FDR_BY")]
names(ADJ_2)[which(names(ADJ_2)=="QQ")] <- "P_Exp"

## FRQ Array
FRQ_2 <- FRQ[,c("SNP","A1","A2","C.HOM.A1.","C.HET.","C.HOM.A2.")]

## PV Array
PV_2 <- PV[,c("SNP","BETA","STAT","P")]
names(PV_2)[which(names(PV_2)=="P")] <- "P_Assoc_2"

## Gene Array
CHR_POS_G <- paste(GENE$Chromosome,":",GENE$Begin,"-",GENE$End,sep="")
GENE_2 <- data.frame( GENE, CHR_POS=CHR_POS_G )
# GENE_2 <- GENE[which(GENE$VarType=="snp"),] # GENE_2 <- GENE[which(ANNOT$VarType=="snp"),]

####################################
## MERGE TABLES ####################
####################################

## Start Merging Tables
MG_1 <- merge( TXT_2,HW_2, by="SNP", all=T )
MG_2 <- merge( MG_1,ADJ_2, by="SNP", all=T )
MG_3 <- merge( MG_2,FRQ_2, by="SNP", all=T )
MG_4 <- merge( MG_3,PV_2, by="SNP", all=T )
MG_5 <- merge( MG_4,GENE_2, by="CHR_POS", all=T )
MG_6 <- MG_5

## Try to fix Merging Issues
NO_SNP <- which(is.na( MG_6$SNP ))
NO_VarType <- which(is.na( MG_6$VarType ))
 # Get CHR:BP for all vars missing Annotations but that have Plink Output (NO_VarType)
POS_NO_VarType <- paste("chr",MG_6$CHR[NO_VarType],":",MG_6$BP[NO_VarType],sep="")
 # Get Chromosome:Start & Chromosome:End for all vars missing Plink Output
Begin_NO_SNP <- paste(MG_6$Chromosome[NO_SNP],MG_6$Begin[NO_SNP],sep=":")
End_NO_SNP <- paste(MG_6$Chromosome[NO_SNP],MG_6$End[NO_SNP],sep=":")
 # Which are in POS* & Begin* [or] the POS* & End* lists??
In_POS_Begin <- intersect(POS_NO_VarType,Begin_NO_SNP)
In_POS_End <- intersect(POS_NO_VarType,End_NO_SNP)
MG_6[,1] <- as.character(MG_6[,1])
 # Loop Thru, b/c wtf not...
for ( b in 1:length(In_POS_Begin) ) {
	Beg_Pos <- In_POS_Begin[b]
	MG_Coord_NO_VarType <- NO_VarType[ which( POS_NO_VarType==Beg_Pos ) ]
	MG_Coord_NO_SNP <- NO_SNP[ which( Begin_NO_SNP==Beg_Pos ) ]
	MG_6[MG_Coord_NO_SNP,2:18] <- MG_6[MG_Coord_NO_VarType,2:18]
	MG_6[MG_Coord_NO_VarType,1] <- "XXXX"
} # length(which(is.na(MG_6))) ; length(which(is.na(MG_5)))
for ( e in 1:length(In_POS_End) ) {
	End_Pos <- In_POS_End[e]
	MG_Coord_NO_VarType <- NO_VarType[ which( POS_NO_VarType==End_Pos ) ]
	MG_Coord_NO_SNP <- NO_SNP[ which( End_NO_SNP==End_Pos ) ]
	# MG_6[MG_Coord_NO_SNP,19:ncol(MG_6)] <- MG_6[MG_Coord_NO_VarType,19:ncol(MG_6)]
	MG_6[MG_Coord_NO_SNP,2:18] <- MG_6[MG_Coord_NO_VarType,2:18]
	MG_6[MG_Coord_NO_VarType,1] <- "XXXX"
} # length(which(is.na(MG_6))) ; length(which(is.na(MG_5)))
# length(table( which(is.na(MG_6),arr.ind=T)[,1] )) ; length(table( which(is.na(MG_5),arr.ind=T)[,1] ))

 # Remove the Rows in "NO_SNP" that I matched above
XXXX_ROWS <- which( MG_6[,1]=="XXXX" )
MG <- MG_6[-XXXX_ROWS,] # length(table( which(is.na(MG),arr.ind=T)[,1] ))

####################################
## FORMAT/ORGANIZE TABLES ##########
####################################

## Make SVS Identifier
CHR_POS_SVS <- paste(MG$SNP,MG$VarType,sep="-")
CHR_POS_SVS <- gsub("snp","SNV",CHR_POS_SVS) ; CHR_POS_SVS <- gsub("del","DEL",CHR_POS_SVS) ; CHR_POS_SVS <- gsub("ins","INS",CHR_POS_SVS)

## Calculate Genotype Frequency for Cohort
 # Clarify Ref/Alt Alleles
WHICH_MATCH <- numeric(nrow(MG))
WHICH_MATCH[which(as.character(MG$A1)==as.character(MG$Reference))] <- 1
WHICH_MATCH[which(as.character(MG$A2)==as.character(MG$Reference))] <- 2
MISS <- which( WHICH_MATCH==0 )
for ( row in MISS ) {
	if ( nchar(gsub("-","",as.character(MG[row,"A1"])))>nchar(gsub("-","",as.character(MG[row,"A2"])))
		& nchar(gsub("-","",as.character(MG[row,"Reference"])))>nchar(gsub("-","",as.character(MG[row,"Allele"]))) ) {
		WHICH_MATCH[row] <- 1
	}
	if ( nchar(gsub("-","",as.character(MG[row,"A1"])))>nchar(gsub("-","",as.character(MG[row,"A2"])))
		& nchar(gsub("-","",as.character(MG[row,"Reference"])))<nchar(gsub("-","",as.character(MG[row,"Allele"]))) ) {
		WHICH_MATCH[row] <- 2
	}
	if ( nchar(gsub("-","",as.character(MG[row,"A1"])))<nchar(gsub("-","",as.character(MG[row,"A2"])))
		& nchar(gsub("-","",as.character(MG[row,"Reference"])))<nchar(gsub("-","",as.character(MG[row,"Allele"]))) ) {
		WHICH_MATCH[row] <- 1
	}
	if ( nchar(gsub("-","",as.character(MG[row,"A1"])))<nchar(gsub("-","",as.character(MG[row,"A2"])))
		& nchar(gsub("-","",as.character(MG[row,"Reference"])))>nchar(gsub("-","",as.character(MG[row,"Allele"]))) ) {
		WHICH_MATCH[row] <- 2
	}
}

 # Calculate Allele Frequencies
NUM_PATIENTS <- median(rowSums( MG[,c("C.HOM.A1.","C.HET.","C.HOM.A2.")] ), na.rm=T)
A1_FREQ <- ( 2*MG$C.HOM.A1. + MG$C.HET. ) / ( 2*NUM_PATIENTS )
A2_FREQ <- ( 2*MG$C.HOM.A2. + MG$C.HET. ) / ( 2*NUM_PATIENTS )
ALLELE_FREQS <- data.frame(A1_FREQ,A2_FREQ)
GENO_FREQS <- MG[,c("C.HOM.A1.","C.HOM.A2.")]/NUM_PATIENTS
FREQ_PRC <- array(, c(nrow(MG),5) ) ; colnames(FREQ_PRC) <- c("REF_ALL","ALT_ALL","REF_GENO","ALT_GENO","HET")
for ( row in 1:nrow(MG) ) {
	WHICH_REF <- WHICH_MATCH[row]
	if ( WHICH_REF > 0 ) {
		WHICH_ALT <- 3-WHICH_REF
		FREQ_PRC[row,"REF_ALL"] <- ALLELE_FREQS[row,WHICH_REF]
		FREQ_PRC[row,"ALT_ALL"] <- ALLELE_FREQS[row,WHICH_ALT]
		FREQ_PRC[row,"REF_GENO"] <- GENO_FREQS[row,WHICH_REF]
		FREQ_PRC[row,"ALT_GENO"] <- GENO_FREQS[row,WHICH_ALT]
		FREQ_PRC[row,"HET"] <- MG$C.HET.[row] / NUM_PATIENTS
	}
}
FREQS_2 <- data.frame(ALLELE_FREQS,FREQ_PRC)

## Pull out predicted functional impact for each variant
FUNC_PRED <- rep( "",nrow(MG) )
for ( row in 1:nrow(MG) ) {
	TEMP_STR <- strsplit( as.character(MG$Functional_Impact[row]), "Predicted " )
	FUNC_PRED[row] <- TEMP_STR[[1]][ length(TEMP_STR[[1]]) ]
}

## Put all columns together (including Pheno, Covs, Date)
MG_WRITE <- data.frame(MG, DATE=rep(Date,nrow(MG)), PHENO=rep(Pheno,nrow(MG)), COVS=rep(Covs,nrow(MG)), SVS_ID=CHR_POS_SVS, FREQS_2, FUNC_PRED )

## Re-order Chromosomes
MG_WRITE <- MG_WRITE[order(as.numeric(as.character(MG_WRITE$BP)),decreasing=F), ]
MG_WRITE <- MG_WRITE[order(as.numeric(as.character(MG_WRITE$CHR)),decreasing=F), ]
MG_WRITE <- MG_WRITE[c( which(MG_WRITE$CHR!=0),which(MG_WRITE$CHR==0),which(is.na(MG_WRITE$CHR)) ), ]

## Reorder ALL columns in more logical manner
COL_ORDER <- c("DATE","PHENO","COVS","CHR","BP","SVS_ID","SNP","VarType","CHR_POS","Chromosome","Begin","End",
	"A1","A2","Reference","Allele","GENO","REF_GENO","HET","ALT_GENO","A1_FREQ","A2_FREQ","REF_ALL","ALT_ALL","C.HOM.A1.","C.HET.","C.HOM.A2.","P_HW",
	"X1000genomes_AFR","X1000genomes_AMR","X1000genomes_ASN","X1000genomes_EUR","X1000GENOMES_AF","CG_WELLDERLY_AF",
	"BETA","STAT","P_Assoc","P_Exp","FDR_BY","FDR_BH",
	"Gene","Gene_Type","Location","Coding_Impact","Functional_Impact","FUNC_PRED",
	"Splice_Site_Pred","Protein_Impact_Prediction.Polyphen.","Protein_Impact_Prediction.SIFT.","Protein_Impact_Prediction.Condel.",
	"Gene_Ontology","Disease_Ontology","omimGene_ID.omimGene_association","Protein_Domain_Gene_Ontology","dbSNP_ID" )
setdiff(names(MG_WRITE),COL_ORDER)
 # Create Full Table
MG_WRITE_FULL <- MG_WRITE[,COL_ORDER] # t(MG_WRITE_FULL[875,])

## Remove redundant or superfluous columns from Table
COL_ORDER_SHORT <- c("DATE","PHENO","COVS","CHR","BP","SVS_ID","SNP",
	"Reference","Allele","A1","A2","GENO","REF_GENO","HET","ALT_GENO","REF_ALL","ALT_ALL","P_HW",
	"BETA","STAT","P_Assoc","P_Exp","FDR_BY","FDR_BH",
	"X1000genomes_AFR","X1000genomes_AMR","X1000genomes_ASN","X1000genomes_EUR","X1000GENOMES_AF","CG_WELLDERLY_AF",
	"Gene","Gene_Type","Location","Coding_Impact","FUNC_PRED","Functional_Impact",
	"Splice_Site_Pred","Protein_Impact_Prediction.Polyphen.","Protein_Impact_Prediction.SIFT.","Protein_Impact_Prediction.Condel.",
	"Gene_Ontology","Disease_Ontology" )
 # Create Short Table
MG_WRITE_SHORT <- MG_WRITE_FULL[,COL_ORDER_SHORT]

####################################
## WRITE TABLES ####################
####################################
write.table(MG_WRITE_FULL, gsub("compiled","compile.full",PathToWrite),sep="\t",row.names=F,col.names=T,quote=F)
write.table(MG_WRITE_SHORT, gsub("compiled","compile.short",PathToWrite),sep="\t",row.names=F,col.names=T,quote=F)

########################################################################
## COMPILE SUMMARY BY GENE/VARIANT #####################################
########################################################################
DAT <- MG_WRITE
## How many of each type of variant in each gene
Compile <- c()
# Parse Gene/Loc/Cod columns
for (var in 1:nrow(DAT)) {
	# Compile
	# print(paste("######",var))
	Chr_Pos <- as.character(DAT$CHR_POS[var])
	Gene_Split1 <- strsplit(as.character(DAT$Gene[var]),"///")[[1]]
	Gene_Split2 <- t(sapply(strsplit(as.character(Gene_Split1),"(",fixed=T),"[",1:2))
	Loc_Split <- strsplit(as.character(DAT$Location[var]),"///")[[1]]
	Cod_Split <- strsplit(as.character(DAT$Coding_Impact[var]),"///")[[1]]
	Pos_Vec <- rep(Chr_Pos,length(Loc_Split))
	# print(length(Pos_Vec))
	# print(nrow(Gene_Split2))
	# print(length(Cod_Split))
	Compile <- rbind( Compile,cbind(Pos_Vec,Gene_Split2,Loc_Split,Cod_Split,paste(Pos_Vec,Gene_Split2[,1],sep="_")) )
}
colnames(Compile) <- c("CHR_POS","Gene","Transcript","Location","Coding_Impact","Unq")
Compile <- data.frame(Compile)

## Get Location of Vars for each Variant/Gene Combo
Unq <- unique(Compile$Unq)
# Uniq_Genes <- unique(Compile[,"Gene"])
Loc_Comp <- array(0,dim=c(length(Unq),8))
colnames(Loc_Comp) <- c("CHR_POS","Gene","Upstream","UTR","Intron","Exon","Downstream","ncRNA")
for (u in 1:length(Unq)) {
	Which_Unq <- Unq[u]
	TEMP <- Compile[which(Compile$Unq==Which_Unq),]
	Loc_Comp[u,1] <- as.character(TEMP$CHR_POS[1])
	Loc_Comp[u,2] <- as.character(TEMP$Gene[1])
	# Loc_Comp
	if ( length(grep("Upstream",TEMP[,"Location"]))>0 ) { Loc_Comp[u,3] <- 1 }
	if ( length(grep("UTR",TEMP[,"Location"]))>0 ) { Loc_Comp[u,4] <- 1 }
	if ( length(grep("Intron",TEMP[,"Location"]))>0 ) { Loc_Comp[u,5] <- 1 }
	if ( length(grep("Exon",TEMP[,"Location"]))>0 ) { Loc_Comp[u,6] <- 1 }
	if ( length(grep("Downstream",TEMP[,"Location"]))>0 ) { Loc_Comp[u,7] <- 1 }
	if ( length(grep("noncoding_rna",TEMP[,"Location"]))>0 ) { Loc_Comp[u,8] <- 1 }
}
Loc_Comp <- data.frame(Loc_Comp)
	
# Compile Counts by Gene
Unq_G <- unique(Loc_Comp[,"Gene"])
GENE_SUMM <- array(,c(length(Unq_G),7))
colnames(GENE_SUMM) <- c("Total","Upstream","UTR","Intron","Exon","Downstream","ncRNA")
rownames(GENE_SUMM) <- Unq_G
for (g in 1:length(Unq_G)) {
	Gene <- Unq_G[g]
	TEMP <- Loc_Comp[which(Loc_Comp[,"Gene"]==Gene),]
	TEMP2 <- data.matrix(TEMP[,3:8])-1
	GENE_SUMM[g,1] <- sum(TEMP2)
# if (nrow(TEMP2)>1 ) {
		GENE_SUMM[g,2:7] <- colSums(TEMP2)
	# }
	# else { GENE_SUMM[g,2:6] <- TEMP2[3:7] }
}

####################################
## PLOT VARS BY FUNCTION ###########
####################################
COLS <- c("dodgerblue2","chocolate2","gold2","firebrick2","chartreuse2","purple2")
jpeg(gsub("compiled","plot.jpeg",PathToWrite),width=1000+30*nrow(GENE_SUMM),height=1400,pointsize=23)
barplot(t(GENE_SUMM[which(!is.na(rownames(GENE_SUMM))),2:7]),legend=T,col=COLS,las=2,main="Variants per Gene by Location",ylab="Counts")
dev.off()

####################################
## WRITE SUMMARY TABLES ############
####################################

# Organized by Transcript
write.table(Compile,gsub("compiled","transcript",PathToWrite),sep="\t",row.names=F,col.names=T,quote=F)

# Gene_SNP Identifier Table
write.table(Loc_Comp,gsub("compiled","location",PathToWrite),sep="\t",row.names=F,col.names=T,quote=F)

# Variant Counts by Gene
write.table(GENE_SUMM,gsub("compiled","plot.table",PathToWrite),sep="\t",row.names=T,col.names=T,quote=F)














# # MG_1 <- merge(TXT_2,GENE_2,by.x="CHR_POS",by.y="Chromosome.Begin.End",all=T)
# MG_1 <- merge(TXT_2,GENE_2,by="CHR_POS",all=T)
# MG_2 <- merge(HW_2,ADJ_2,by="SNP",all=T)
# MG_3 <- merge(MG_2,FRQ_2,by="SNP",all=T)
# MG <- merge(MG_1,MG_3,by="SNP",all=T) # MG <- merge(MG_1,MG_2,by="SNP",all=T)

## Perform GWAS on Variants in BGI and HC Data Sets ##
## Use Height or Soemthing as Phenotype (CRP, maybe?) ##
## This is for BigCall Paper to assess differences in BGI/HC Call Sets ##
## May 13, 2014 ##
## Unnecessary Comment
## Kristopher Standish ##

########################################################
## Input Parameters for Run ############################

# Specify some Input Parameters
DATE=$1 # Date/Identifier to be used in Filenames
PHENO=$2 # Path to Phenotype File
COVS=$3 # Path to Covarites File
PC_COUNT=$4
PHE_TYPE=$5
ASSOC_DIR=$6
ANNOTS=%7

########################################################
## Set Paths To Tools/Data #############################

# Specify list of Covariates to include (for command ad for filename)
if [ $PC_COUNT -eq 0 ]
then
COVS_COMMAND=`echo "${COVS}" | sed 's/QQQ/,/g'`
COVS_FILENAME=`echo "${COVS}" | sed 's/QQQ/_/g'`
else
PCS=`seq 1 ${PC_COUNT}`
PCS_COMMAND=`echo "PC"${PCS} | sed 's/ /QQQPC/g'`
COVS_COMMAND=`echo "${COVS}QQQ${PCS_COMMAND}" | sed 's/QQQ/,/g'`
COVS_FILENAME=`echo "${COVS}QQQ${PCS_COMMAND}" | sed 's/QQQ/_/g'`
fi

# Specify commands and extensions for Cont vs Bin Phenotype
if [ $PHE_TYPE = "C" ]
then
SUFFIX=linear
else
SUFFIX=logistic
fi

# Specify Some Paths (Constant)
HOME_DIR=/projects/janssen/ASSOCIATION
OUT_DIR=${ASSOC_DIR}/${DATE}_${PHENO}_${COVS_FILENAME}
BED_DIR=/projects/janssen/scripts/AS-ASSOCIATION/BED_FILES
mkdir ${OUT_DIR}

# Specify The Phenotype/Covariate Files
PHE_FILE=${HOME_DIR}/PH-PHENOTYPES/${PHENO}.txt
COV_FILE=${HOME_DIR}/PH-PHENOTYPES/COV.txt
COV_ANC_FILE=${HOME_DIR}/PH-PHENOTYPES/COV_ANC.txt
EIG_VEC=${HOME_DIR}/EIGEN/HC_FULL.eigenvec

# Specify Some Paths To Tools
R_MANH_PLOT=/projects/janssen/ASSOCIATION/SCRIPTS/Manhat_Plot.R
R_MAKE_COV_TAB=/projects/janssen/ASSOCIATION/SCRIPTS/Make_Cov_Tab.R
R_PCA_PLOT=/projects/janssen/ASSOCIATION/SCRIPTS/PCA_Plot.R
PLINK=/projects/janssen/Tools/plink_linux_x86_64/plink

SET=HC_FULL

# Make Path for Covariate and Output Files
NEW_COV_FILE=${OUT_DIR}/${SET}_COV_w_PC.txt
OUT_FILE_1=${OUT_DIR}/${DATE}-${SET}_${PHENO}_FULL
OUT_FILE_2=${OUT_DIR}/${DATE}-${SET}_${PHENO}_FULL.assoc.${SUFFIX}
P_FILE=${OUT_DIR}/${DATE}-${SET}_${PHENO}_FULL.P

echo $OUT_FILE_1
echo $OUT_FILE_2
echo $P_FILE

########################################################
## Use/Plot Principal Components #######################

# Make Covariate File for this Run (That includes PCs)
Rscript ${R_MAKE_COV_TAB} ${EIG_VEC} ${COV_FILE} ${SET} ${OUT_DIR}

# Make PCA Plots (colored by Self Report)
#Rscript ${R_PCA_PLOT} ${EIG_VEC} ${COV_ANC_FILE} ${OUT_DIR}

########################################################
## Perform Association Test ############################

# Run PLINK for Whole Variant Set
${PLINK} --bfile ${BED_DIR}/BED_${SET}_SNP \
--pheno ${PHE_FILE} \
--covar ${NEW_COV_FILE} --covar-name ${COVS_COMMAND} \
--${SUFFIX} hide-covar --adjust \
--allow-no-sex  \
--out ${OUT_FILE_1}

########################################################
## Compile Results #####################################

# Pull out just the p-values
cat ${OUT_FILE_2} | awk '{print $1"\t"$3"\t"$9}' > ${P_FILE}

########################################################
## Plot this Shiz ######################################

# Make manhattan plot
Rscript ${R_MANH_PLOT} ${P_FILE} ${PHENO} ${COVS_FILENAME}

########################################################
## Pull out dbSNP ######################################

########################################################
## Pull out Annots #####################################
if [ $ANNOTS != "F" ] ; then
python Get_Annot.py 
# var1_file - file with positions of variants you want annotations for
# annot_file - file with annotations at a bunch of variant positions
# output - outut file containing annotations for positions in var1



#############################
######## END OF DOC #########
#############################

## Single-Locus Tests ##

################################################################################
## RF_ACPA #####################################################################

## Names/Paths
echo \### Defining Set Variables and Paths at `date` \###
DATE=20150410
HOME_DIR=/projects/janssen/ASSOCIATION

## Files
VAR_FILE=${HOME_DIR}/../VCFs/PLINK/BED_MAF1.ALL.bed
ANNOTS=${HOME_DIR}/../ANNOTATE/JnJ_121613_all_annotations.txt.bgz
PHENO_FILE=${HOME_DIR}/../ASSOCIATION/PH-PHENOTYPES/LT8_RF.txt
PHENO_TYPE=B
COV_FILE=${HOME_DIR}/../ASSOCIATION/PH-PHENOTYPES/COV.txt
COVS=`echo SEX`
EIG_VEC=${HOME_DIR}/../ASSOCIATION/EIGEN/HC_FULL.eigenvec
PC_COUNT=2
START_STEP=1

COVS=`echo "$COVS" | sed 's/ /QQQ/g'`

########################################
#### Pheno = RF ########################
## Run The Script
/projects/janssen/ASSOCIATION/SCRIPTS/GWAS.sh \
${DATE} \
${HOME_DIR} \
${VAR_FILE} \
${ANNOTS} \
${PHENO_FILE} \
${PHENO_TYPE} \
${COV_FILE} \
${COVS} \
${EIG_VEC} \
${PC_COUNT} \
${START_STEP}

########################################
#### Pheno = RF_ACPA ###################
## Run The Script
PHENO_FILE=${HOME_DIR}/../ASSOCIATION/PH-PHENOTYPES/LT8_RF_ACPA.txt
/projects/janssen/ASSOCIATION/SCRIPTS/GWAS.sh \
${DATE} \
${HOME_DIR} \
${VAR_FILE} \
${ANNOTS} \
${PHENO_FILE} \
${PHENO_TYPE} \
${COV_FILE} \
${COVS} \
${EIG_VEC} \
${PC_COUNT} \
${START_STEP}

########################################
#### Pheno = ACPA ######################
## Run The Script
PHENO_FILE=${HOME_DIR}/../ASSOCIATION/PH-PHENOTYPES/LT8_ACPA.txt
/projects/janssen/ASSOCIATION/SCRIPTS/GWAS.sh \
${DATE} \
${HOME_DIR} \
${VAR_FILE} \
${ANNOTS} \
${PHENO_FILE} \
${PHENO_TYPE} \
${COV_FILE} \
${COVS} \
${EIG_VEC} \
${PC_COUNT} \
${START_STEP}




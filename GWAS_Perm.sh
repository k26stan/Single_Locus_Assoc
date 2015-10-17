## Perform GWAS on Variants for Population Data Sets ##
## April/June, 2014 ##
## Kristopher Standish ##
## UCSD/JCVI ##

###########################################################
## 1 ## Set up Paths ######################################
###########################################################
echo \### 1 - `date` \###
echo \### Defining Set Variables and Paths \###

###########################################################
## Manually Input Parameters ##

# Names/Paths
DATE=$1
HOME_DIR=$2

# Parameters and Files
VAR_FILE=$3
#IN_VCF=$4 # Path to VCF File with all the Variants
ANNOTS=$4 # Path to Annotation File
PHENO_FILE=$5 # Which Phenotype File are you using?
PHENO_TYPE=$6 # Is phenotype (B)inary or (C)ontinuous?
COV_FILE=$7 # Path to Covariate File or "F"
COVS=$8 # Which Covariates to Include?
EIG_VEC=$9 # Output from Plink's --pca command (MAF>1%) or "F"
PC_COUNT=${10} # How many PCs to Include as Covariates?
START_STEP=${11} # Which Step do you want to start on?

###########################################################
## Constant Paths ##

# Public Tools
GATK_JAR=/projects/janssen/Tools/gatk2.7-2/GenomeAnalysisTK.jar
REF_FA=/projects/janssen/ref/ref.fa
VCF_TOOLS=/projects/janssen/Tools/vcftools_0.1.11/bin/vcftools
PLINK=/projects/janssen/Tools/plink_linux_x86_64/plink 
GENE_TABLE=/home/kstandis/HandyStuff/GG-Gene_Names_DB.txt
# Custom Scripts
MANH_PLOT_R=/projects/janssen/ASSOCIATION/SCRIPTS/Manhat_Plot.R
MAKE_COV_TAB_R=/projects/janssen/ASSOCIATION/SCRIPTS/Make_Cov_Tab.R
PCA_PLOT_R=/projects/janssen/ASSOCIATION/SCRIPTS/PCA_Plot.R
GET_ANNOT_PY=/projects/janssen/ASSOCIATION/SCRIPTS/Get_Annot.py
GET_PLINK_OUT_PY=/projects/janssen/ASSOCIATION/SCRIPTS/Get_Plink_Out.py
COMPILE_OUTS_R=/projects/janssen/ASSOCIATION/SCRIPTS/Compile_Outs.R

###########################################################
## Pull some Info out of Parameters ##

# Get Names of Specific Files
DIRS=(${VAR_FILE//\// })
VAR_FILE_NAME=${DIRS[${#DIRS[@]} - 1]} # Get Name of Variant File
DIRS=(${PHENO_FILE//\// })
PHENO=${DIRS[${#DIRS[@]} - 1]} # Get Name of Phenotype File

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

# Incorporate Country/Site of Study as Binary Covariate (if Included)
if [[ $COVS == *COUN* ]]
then
COVS_COMMAND=`echo $COVS_COMMAND | sed 's/COUN/CN_ARG,CN_AUS,CN_COL,CN_HUN,CN_LTU,CN_MEX,CN_MYS,CN_NZL,CN_POL,CN_RUS,CN_UKR/g'`
fi

# Specify commands and extensions for Cont vs Bin Phenotype
if [ $PHENO_TYPE = "C" ]
then
SUFFIX=linear
else
SUFFIX=logistic
fi

# Set up Directory for today's adventure
echo \### Moving to Base Directory at `date`: \###
echo $ASSOC
ASSOC=${HOME_DIR}/${DATE}_${PHENO%%.txt}_${COVS_FILENAME}
mkdir ${ASSOC}
NEW_COV_FILE=${ASSOC}/Cov_w_PCs.txt
ASSOC_FILE=${ASSOC}/${DATE}_${PHENO%%.txt}_${COVS_FILENAME} # .assoc.${SUFFIX}
P_FILE=${ASSOC}/${DATE}_${PHENO%%.txt}_${COVS_FILENAME}.P
CND_FILE=${ASSOC}/CND_${PHENO%%.txt}_${COVS_FILENAME}.txt
CND_ANNOTS=${ASSOC}/CND_${PHENO%%.txt}_${COVS_FILENAME}.Annot.txt
CND_GENES=${ASSOC}/CND_${PHENO%%.txt}_${COVS_FILENAME}.Gene.txt
cd ${ASSOC}

# Specify a File to which to Write Updates
UPDATE_FILE=${ASSOC}/Update.txt

if [ "$START_STEP" -le 1 ]; then
echo `date` "1 - Got all the paths taken care of...I think" > ${UPDATE_FILE}
echo "This should be the GATK directory:" ${GATK_JAR} >> ${UPDATE_FILE}
printf "V\nV\nV\nV\nV\nV\nV\nV\n"

fi
###########################################################
## 2 ## Determine Variant Format and Adjust ###############
###########################################################
if [ "$START_STEP" -le 2 ]; then
 # If File is VCF, make a PED/MAP file, then BED file
 # If File is PED, make BED file
 # If File is BED, just skip ahead
echo \### 2 - `date` \###
echo \### Determining/Adjusting Variant File Formats \###
echo `date` "2 - Determining/Adjusting Variant File Formats" >> ${UPDATE_FILE}

# Determing File Type and Converty to .bed File
if [ ${VAR_FILE: -4} == ".vcf" ] ; then
	${VCF_TOOLS} --plink --vcf ${VAR_FILE} --out ${ASSOC}_${VAR_FILE_NAME%%.vcf}
	VAR_FILE=${ASSOC}_${VAR_FILE_NAME%%.vcf}.ped
fi
if [ ${VAR_FILE: -4} == ".ped" ] ; then
	${PLINK} --make-bed --file ${VAR_FILE} --out ${VAR_FILE%%.ped}
	VAR_FILE=${VAR_FILE%%.ped}.bed
fi

echo `date` "Variant File now in .bed Format" >> ${UPDATE_FILE}
printf "V\nV\nV\nV\nV\nV\nV\nV\n"

fi
###########################################################
## 3 ## Use/Calc/Plot Principal Components ################
###########################################################
if [ "$START_STEP" -le 3 ]; then
 # Determine if Principal Components have been calculated (given)
   # If not, Calculate them and save the Path
echo \### 3 - `date` \###
echo \### Calculate/Specify PCs \###
echo `date` "3 - Calculate/Specify PCs" >> ${UPDATE_FILE}

# If No Principal Components Exist, Make Them
if [ $EIG_VEC = "F" ] ; then
# Use BED to run PCA
${PLINK} --bfile ${VAR_FILE%%.bed} \
--pca header \
--allow-no-sex \
--out ${ASSOC}/${VAR_FILE_NAME%%.bed}
EIG_VEC=${ASSOC}/${VAR_FILE_NAME%%.bed}.eigenvec
fi

echo `date` "Principal Components Calculated/Specified" >> ${UPDATE_FILE}
printf "V\nV\nV\nV\nV\nV\nV\nV\n"

fi
###########################################################
## 4 ## Make New Covariate File Containing PCs ############
###########################################################
if [ "$START_STEP" -le 4 ]; then
	# Compile Covariate File with PCs
echo \### 4 - `date` \###
echo \### Compile Covariates with PCs \###
echo `date` "4 - Compile Covariates with PCs" >> ${UPDATE_FILE}

# Make new Covariate File
if [ $COV_FILE = "F" ] ; then
cp ${EIG_VEC} ${NEW_COV_FILE}
else
Rscript ${MAKE_COV_TAB_R} ${EIG_VEC} ${COV_FILE} ${COVS_COMMAND} ${NEW_COV_FILE}
fi

echo `date` "Principal Components Calculated/Specified" >> ${UPDATE_FILE}
printf "V\nV\nV\nV\nV\nV\nV\nV\n"

fi
###########################################################
## 5 ## Perform Single-Locus Association Test ############
###########################################################
if [ "$START_STEP" -le 5 ]; then
	# Perform Single-Locus Association Tests using Plink
echo \### 5 - `date` \###
echo \### Perform Single-Locus Association Test \###
echo `date` "5 - Perform Single-Locus Association Test" >> ${UPDATE_FILE}

# Perform Association Test
${PLINK} --bfile ${VAR_FILE%%.bed} \
--pheno ${PHENO_FILE} \
--covar ${NEW_COV_FILE} --covar-name ${COVS_COMMAND} \
--${SUFFIX} hide-covar \
--allow-no-sex \
--maf 0.01 \
--keep-allele-order \
--out ${ASSOC_FILE}

echo `date` "Single-Locus Tests Done" >> ${UPDATE_FILE}
printf "V\nV\nV\nV\nV\nV\nV\nV\n"

fi
###########################################################
## 6 ## Compile Results ###################################
###########################################################
if [ "$START_STEP" -le 6 ]; then
	# Pull out Just the P-Values and Identifiers for Easier Access
echo \### 6 - `date` \###
echo \### Compile/Consolidate Results \###
echo `date` "6 - Compile/Consolidate Results" >> ${UPDATE_FILE}

# Pull out just the p-values & identifying information
cat ${ASSOC_FILE}.assoc.${SUFFIX} | awk '{print $1"\t"$2"\t"$3"\t"$9}' > ${P_FILE}

# Delete Original Association File (to save space)
rm ${ASSOC_FILE}.assoc.${SUFFIX}

echo `date` "Results Compiled" >> ${UPDATE_FILE}
printf "V\nV\nV\nV\nV\nV\nV\nV\n"

fi
###########################################################
## X ## End of Doc ########################################
###########################################################

echo `date` "DONE DONE DONE DONE DONE DONE DONE DONE DONE" >> ${UPDATE_FILE}
echo `date`
echo "DONE DONE DONE DONE DONE DONE DONE DONE DONE"
printf "V\nV\nV\nV\nV\nV\nV\nV\n"



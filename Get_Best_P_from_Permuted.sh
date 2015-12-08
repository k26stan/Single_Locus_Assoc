## Find Minimum Value in Each Permuted Output File

BEST_P=/projects/janssen/ASSOCIATION/20151016_Permute_DEL_Best.P
Get_Best() {
	P_FILE=$1
	awk 'NR == 1 {line = $0; min = $4}
     NR > 1 && $4 < min {line = $0; min = $4}
     END{print line}' ${P_FILE} >> ${BEST_P}
}

for file in `ls 20151016_Th*/20*/*P`; do
	Get_Best ${file}
	echo $file Done
done


HOME_DIR=/projects/janssen/ASSOCIATION/20151016_Permute
BEST_P=${HOME_DIR}/Compile_Best.P
for n in `seq $NUM_PERM`; do
for n in `seq 1`; do
for t in `seq 5`; do
	P_DIR=${HOME_DIR}/20151016_Th${t}_Perm_Pheno.${n}_DAS_BL_MN_PC1_PC2
	P_FILE=${P_DIR}/20151016_Th${t}_Perm_Pheno.${n}_DAS_BL_MN_PC1_PC2.P
	Get_Best ${P_FILE}
done
done


#!/usr/bash
DIR_SCRIPTS=/Users/KANG/geneoscopy_dev/scripts
DIR_DATA=/Users/KANG/geneoscopy_dev/data/run_proj_b_c_d

NUM_SAMPLES=84
GROUP=N_vs_C
THLD_PVAL=0.005
THLD_FC=0
NORMALIZED_CHIPDATA=${DIR_DATA}/chipdata_rma.expression_console.txt
QC_TABLE=${DIR_DATA}/QC_table_combined_b_c_d.txt
SAMPLE_SHEET=${DIR_DATA}/sample_sheet_combined_b_c_d.txt
VALID_CHIPS=${DIR_DATA}/valid_chips.txt

if [ $GROUP = 'N_vs_C' ]; then
	SAMPLE_SHEET=${SAMPLE_SHEET%.txt}.normal_vs_cancer.txt
	VALID_CHIPS=${VALID_CHIPS%.txt}.normal_vs_cancer.txt
fi

echo "Running quality control ... "
python ${DIR_SCRIPTS}/quality_control.py -n $NUM_SAMPLES -g $GROUP -d $NORMALIZED_CHIPDATA -q $QC_TABLE -s $SAMPLE_SHEET -v foo -o ${DIR_DATA}/chipdata_geneset_x_valid_chips.txt

echo "Splitting training/testing sets ... "
mkdir -p ${DIR_DATA}/training
mkdir -p ${DIR_DATA}/testing
python ${DIR_SCRIPTS}/split_expr_train_vs_test.py -v $VALID_CHIPS -i ${DIR_DATA}/chipdata_geneset_x_valid_chips.txt -tr0 ${DIR_DATA}/training/chipdata.txt -tr1 ${DIR_DATA}/training/valid_chips.txt -te0 ${DIR_DATA}/testing/chipdata.txt -te1 ${DIR_DATA}/testing/valid_chips.txt

echo "Analyzing DE genes ... "
Rscript ${DIR_SCRIPTS}/de_analysis.r ${DIR_DATA}/training/chipdata.txt ${DIR_DATA}/training/valid_chips.txt $GROUP $THLD_PVAL $THLD_FC ${DIR_DATA}/training/top_de_genes.txt
NUM_TOP_GENES=$(( $( wc -l ${DIR_DATA}/training/top_de_genes.txt | awk '{print $1}') - 1 ))
echo $NUM_TOP_GENES "DE genes found"

ML_MODELS=(random_forest svm neural_net grad_boosting)
for ML_MODEL in "${ML_MODELS[@]}"; do
	echo "###" $ML_MODEL "###"
	echo "Training models ... "
	python ${DIR_SCRIPTS}/convert_expr_for_ml.py -t ${DIR_DATA}/training/top_de_genes.txt -i ${DIR_DATA}/training/chipdata.txt -o ${DIR_DATA}/training/training_set.txt
	python ${DIR_SCRIPTS}/convert_expr_for_ml.py -t ${DIR_DATA}/training/top_de_genes.txt -i ${DIR_DATA}/testing/chipdata.txt -o ${DIR_DATA}/testing/testing_set.txt
	rm -rf $ML_MODEL
	python ${DIR_SCRIPTS}/crc_training.py -i ${DIR_DATA}/training/training_set.txt -a $ML_MODEL -o ${DIR_DATA}/training/${ML_MODEL}

	echo "Testing prediction ... "
	python ${DIR_SCRIPTS}/crc_prediction.py -i ${DIR_DATA}/testing/testing_set.txt -a $ML_MODEL -m ${DIR_DATA}/training/${ML_MODEL}/${ML_MODEL}_model.pkl
done


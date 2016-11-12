#!/usr/bash
DIR_SCRIPTS=/Users/KANG/geneoscopy_dev/scripts

##### INPUT VARIABLES #####

DIR_DATA=/Users/KANG/geneoscopy_dev/data/run_proj_abcdefgh_1
NUM_SAMPLES=177
GROUP=N_vs_C
GENE_FILTER=nanostring
THLD_QC=0.75
THLD_PVAL=0.2
# THLD_PVAL=0.0025
THLD_FC=0
NORMALIZED_CHIPDATA_FULL=${DIR_DATA}/chipdata_rma.expression_console.txt
NORMALIZED_CHIPDATA=${DIR_DATA}/chipdata_rma.expression_console.${GENE_FILTER}.txt
# GENE_FILTER_LST=${DIR_DATA}/../external_data/nanostring/PanCancer_nanostring_genes_annotated.txt 
GENE_FILTER_LST=${DIR_DATA}/../external_data/Genecards_colon_cancer/GeneCards_Nanostring_genes_annotated.txt
QC_TABLE=${DIR_DATA}/QC_table_combined_abcdefgh.txt
SAMPLE_SHEET=${DIR_DATA}/sample_sheet_combined_abcdefgh.two_group_no_benign.txt
PATIENT_SHEET=${DIR_DATA}/patient_info_sheet.txt
VALID_CHIPS=${DIR_DATA}/valid_chips.two_group_no_benign.txt
CRC_PREDICTORS=${DIR_DATA}/../external_data/CIViC/civic_selected_genes_TCs.txt

##### END OF INPUT VARIABLES #####

echo "Num of total samples:" $NUM_SAMPLES
echo "Label groups:" $GROUP
echo "DE p-value treshold:" $THLD_PVAL
echo "DE fold change threshold:" $THLD_FC
echo "Normalized expr data:" $NORMALIZED_CHIPDATA
echo "QC table:" $QC_TABLE
echo "Sample sheet:" $SAMPLE_SHEET
echo ""

echo "Filtering gene set ... "
python ${DIR_SCRIPTS}/filter_gene_set.py -i $NORMALIZED_CHIPDATA_FULL -l $GENE_FILTER_LST -c 2 -o $NORMALIZED_CHIPDATA
# NORMALIZED_CHIPDATA=$NORMALIZED_CHIPDATA_FULL

echo "Running quality control ... "
python ${DIR_SCRIPTS}/quality_control.py -n $NUM_SAMPLES -g $GROUP -d $NORMALIZED_CHIPDATA -q $QC_TABLE -s $SAMPLE_SHEET -t $THLD_QC -v $VALID_CHIPS -o ${DIR_DATA}/chipdata_geneset_x_valid_chips.txt
python ${DIR_SCRIPTS}/quality_control.py -n $NUM_SAMPLES -g $GROUP -d $NORMALIZED_CHIPDATA_FULL -q $QC_TABLE -s $SAMPLE_SHEET -t $THLD_QC -v foo -o ${DIR_DATA}/chipdata_geneset_x_valid_chips_full.txt

echo "Splitting training/testing sets ... "
mkdir -p ${DIR_DATA}/training
mkdir -p ${DIR_DATA}/testing
python ${DIR_SCRIPTS}/split_expr_train_vs_test.py -v $VALID_CHIPS -i ${DIR_DATA}/chipdata_geneset_x_valid_chips.txt -tr0 ${DIR_DATA}/training/chipdata.txt -tr1 ${DIR_DATA}/training/valid_chips.txt -te0 ${DIR_DATA}/testing/chipdata.txt -te1 ${DIR_DATA}/testing/valid_chips.txt
python ${DIR_SCRIPTS}/split_expr_train_vs_test.py -v $VALID_CHIPS -i ${DIR_DATA}/chipdata_geneset_x_valid_chips_full.txt -tr0 ${DIR_DATA}/training/chipdata_full.txt -te0 ${DIR_DATA}/testing/chipdata_full.txt

THLD_PVALS=( 0.05 0.1 0.15 0.2 0.3 0.4 0.5 0.75 1 )
for THLD_PVAL in "${THLD_PVALS[@]}"; do
	echo "#################################"
	echo "P-value threshold -->" $THLD_PVAL
	echo "#################################"
	echo ""

	echo "Analyzing DE genes ... "
	Rscript ${DIR_SCRIPTS}/de_analysis.r ${DIR_DATA}/training/chipdata.txt ${DIR_DATA}/training/valid_chips.txt $GROUP $THLD_PVAL $THLD_FC ${DIR_DATA}/training/top_de_genes.txt

	### python ${DIR_SCRIPTS}/de_analysis_sd.py ${DIR_DATA}/training/chipdata.txt ${DIR_DATA}/training/valid_chips.txt ${DIR_DATA}/training/top_de_genes.txt
	NUM_TOP_GENES=$(( $( wc -l ${DIR_DATA}/training/top_de_genes.txt | awk '{print $1}') - 1 ))
	echo $NUM_TOP_GENES "DE genes found"
	echo ""

	# ML_MODELS=(random_forest svm neural_net grad_boosting gauss_process)
	ML_MODELS=(grad_boosting)
	for ML_MODEL in "${ML_MODELS[@]}"; do
		echo "###" $ML_MODEL "###"
		
		echo "Preparing training/testing datasets ... "
		python ${DIR_SCRIPTS}/convert_expr_for_ml.py -t ${DIR_DATA}/training/top_de_genes.txt -i ${DIR_DATA}/training/chipdata.txt -o ${DIR_DATA}/training/training_set.txt
		python ${DIR_SCRIPTS}/convert_expr_for_ml.py -t ${DIR_DATA}/training/top_de_genes.txt -i ${DIR_DATA}/testing/chipdata.txt -o ${DIR_DATA}/testing/testing_set.txt
		python ${DIR_SCRIPTS}/incorporate_patient_info.py -p ${PATIENT_SHEET} -d ${DIR_DATA}/training/training_set.txt
		python ${DIR_SCRIPTS}/incorporate_patient_info.py -p ${PATIENT_SHEET} -d ${DIR_DATA}/testing/testing_set.txt 

		python ${DIR_SCRIPTS}/convert_expr_for_ml.py -i ${DIR_DATA}/training/chipdata_full.txt -o ${DIR_DATA}/training/training_set_full.txt
		python ${DIR_SCRIPTS}/convert_expr_for_ml.py -i ${DIR_DATA}/testing/chipdata_full.txt -o ${DIR_DATA}/testing/testing_set_full.txt

		echo "Training models ... "
		rm -rf $ML_MODEL
		python ${DIR_SCRIPTS}/crc_training.py -i ${DIR_DATA}/training/training_set.txt -f ${DIR_DATA}/training/training_set_full.txt -a $ML_MODEL -p $CRC_PREDICTORS -o ${DIR_DATA}/training/${ML_MODEL} -s ${DIR_DATA}/training/predictor_normal_stats.txt

		echo "Testing prediction ... "
		python ${DIR_SCRIPTS}/crc_prediction.py -i ${DIR_DATA}/testing/testing_set.txt -f ${DIR_DATA}/testing/testing_set_full.txt -a $ML_MODEL -p $CRC_PREDICTORS -m ${DIR_DATA}/training/${ML_MODEL}/${ML_MODEL}_model.pkl -s ${DIR_DATA}/training/predictor_normal_stats.txt
		echo ""
	done
done

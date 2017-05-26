#!/usr/bash
DIR_SCRIPTS=/Users/KANG/geneoscopy_dev/scripts

##### INPUT VARIABLES #####

DIR_DATA=/Users/KANG/geneoscopy_dev/data/run_proj_batch1-17_1
GROUP=N_vs_P_vs_C
GENE_FILTER=genecards_nanostring_civic
NORMALIZED_CHIPDATA_FULL=${DIR_DATA}/chipdata_rma.expression_console.sst_rma.txt
# NORMALIZED_CHIPDATA=${DIR_DATA}/chipdata_rma.expression_console.${GENE_FILTER}.txt
# GENE_FILTER_LST=${DIR_DATA}/../external_data/nanostring/PanCancer_nanostring_genes_annotated.txt
# GENE_FILTER_LST=${DIR_DATA}/../external_data/Genecards_colon_cancer/GeneCards_Nanostring_genes_annotated.txt  
GENE_FILTER_LST=${DIR_DATA}/../external_data/Genecards_colon_cancer/GeneCards_Nanostring_CIViC_genes_annotated.txt 
# QC_TABLE=${DIR_DATA}/QC_table_combined_abcdefghi.txt
# SAMPLE_SHEET=${DIR_DATA}/sample_sheet_combined_abcdefghi.txt
PATIENT_SHEET=${DIR_DATA}/patient_info_sheet.txt
VALID_CHIPS=${DIR_DATA}/valid_chips.txt
# CRC_PREDICTORS=${DIR_DATA}/../external_data/CIViC/civic_selected_genes_TCs.txt
CV_FOLDS=10

##### END OF INPUT VARIABLES #####

# echo "Working directory:" $DIR_DATA
# echo "Label groups:" $GROUP
# # echo "Normalized expr data:" $NORMALIZED_CHIPDATA_FULL
# echo ""


# echo "Preparing normalized data ..."
# python ${DIR_SCRIPTS}/prepare_normalized_expr_data.py -i $NORMALIZED_CHIPDATA_FULL -l $GENE_FILTER_LST -c 2 -o1 ${DIR_DATA}/chipdata_geneset_x_valid_chips.txt -o2 ${DIR_DATA}/chipdata_geneset_x_valid_chips_full.txt

# #### Use full list of probesets
# # cp ${DIR_DATA}/chipdata_geneset_x_valid_chips_full.txt ${DIR_DATA}/chipdata_geneset_x_valid_chips.txt
# ####

# echo "Splitting training/testing sets ... "
# mkdir -p ${DIR_DATA}/training
# mkdir -p ${DIR_DATA}/testing
# python ${DIR_SCRIPTS}/split_expr_train_vs_test.py -v $VALID_CHIPS -i ${DIR_DATA}/chipdata_geneset_x_valid_chips.txt -tr0 ${DIR_DATA}/training/chipdata.txt -tr1 ${DIR_DATA}/training/valid_chips.txt -te0 ${DIR_DATA}/testing/chipdata.txt -te1 ${DIR_DATA}/testing/valid_chips.txt
# python ${DIR_SCRIPTS}/split_expr_train_vs_test.py -v $VALID_CHIPS -i ${DIR_DATA}/chipdata_geneset_x_valid_chips_full.txt -tr0 ${DIR_DATA}/training/chipdata_full.txt -te0 ${DIR_DATA}/testing/chipdata_full.txt


for i in $(seq 1 100); do
	echo "###########################"
	echo "Random selection run " $i
	echo "###########################"
	echo ""

	echo "Get random DE genes ... "
	python ${DIR_SCRIPTS}/random_gene_selection.py ${DIR_DATA}/training/chipdata.txt ${DIR_DATA}/training/top_de_genes_${i}.txt


	# ML_MODELS=(random_forest svm grad_boosting adaboost gauss_process)
	ML_MODELS=( nu_svm )
	for ML_MODEL in "${ML_MODELS[@]}"; do
		echo "###" $ML_MODEL "###"
		
		echo "Preparing training/testing datasets ... "
		python ${DIR_SCRIPTS}/convert_expr_for_ml.py -t ${DIR_DATA}/training/top_de_genes_${i}.txt -i ${DIR_DATA}/training/chipdata.txt -o ${DIR_DATA}/training/training_set.txt
		python ${DIR_SCRIPTS}/convert_expr_for_ml.py -t ${DIR_DATA}/training/top_de_genes_${i}.txt -i ${DIR_DATA}/testing/chipdata.txt -o ${DIR_DATA}/testing/testing_set.txt
		python ${DIR_SCRIPTS}/incorporate_patient_info.py -p ${PATIENT_SHEET} -d ${DIR_DATA}/training/training_set.txt
		python ${DIR_SCRIPTS}/incorporate_patient_info.py -p ${PATIENT_SHEET} -d ${DIR_DATA}/testing/testing_set.txt 

		python ${DIR_SCRIPTS}/convert_expr_for_ml.py -i ${DIR_DATA}/training/chipdata_full.txt -o ${DIR_DATA}/training/training_set_full.txt
		python ${DIR_SCRIPTS}/convert_expr_for_ml.py -i ${DIR_DATA}/testing/chipdata_full.txt -o ${DIR_DATA}/testing/testing_set_full.txt

		echo "Training models ... "
		rm -rf $ML_MODEL
		python ${DIR_SCRIPTS}/crc_training.py -i ${DIR_DATA}/training/training_set.txt -f ${DIR_DATA}/training/training_set_full.txt -a $ML_MODEL -o ${DIR_DATA}/training/${ML_MODEL} -s ${DIR_DATA}/training/predictor_normal_stats.txt
		# ## ONLY FOR ROC curve
		# python ${DIR_SCRIPTS}/crc_prediction.py -i ${DIR_DATA}/training/training_set.txt -f ${DIR_DATA}/training/training_set_full.txt -a $ML_MODEL -m ${DIR_DATA}/training/${ML_MODEL}/${ML_MODEL}_model.pkl -s ${DIR_DATA}/training/predictor_normal_stats.txt -v 1
		# ###

		echo "Testing prediction ... "
		python ${DIR_SCRIPTS}/crc_prediction.py -i ${DIR_DATA}/testing/testing_set.txt -f ${DIR_DATA}/testing/testing_set_full.txt -a $ML_MODEL -m ${DIR_DATA}/training/${ML_MODEL}/${ML_MODEL}_model.pkl -s ${DIR_DATA}/training/predictor_normal_stats.txt
		echo ""
	done


	# echo "Preparing training/testing datasets ... "
	# python ${DIR_SCRIPTS}/convert_expr_for_ml.py -t ${DIR_DATA}/training/top_de_genes_${i}.txt -i ${DIR_DATA}/training/chipdata.txt -o ${DIR_DATA}/training/training_set.txt
	# python ${DIR_SCRIPTS}/convert_expr_for_ml.py -t ${DIR_DATA}/training/top_de_genes_${i}.txt -i ${DIR_DATA}/testing/chipdata.txt -o ${DIR_DATA}/testing/testing_set.txt
	# python ${DIR_SCRIPTS}/incorporate_patient_info.py -p ${PATIENT_SHEET} -d ${DIR_DATA}/training/training_set.txt
	# python ${DIR_SCRIPTS}/incorporate_patient_info.py -p ${PATIENT_SHEET} -d ${DIR_DATA}/testing/testing_set.txt 

	# python ${DIR_SCRIPTS}/convert_expr_for_ml.py -i ${DIR_DATA}/training/chipdata_full.txt -o ${DIR_DATA}/training/training_set_full.txt
	# python ${DIR_SCRIPTS}/convert_expr_for_ml.py -i ${DIR_DATA}/testing/chipdata_full.txt -o ${DIR_DATA}/testing/testing_set_full.txt

	# echo "Training models ... "
	# rm -rf ${DIR_DATA}/training/ensemble
	# python ${DIR_SCRIPTS}/crc_training_ensemble.py -i ${DIR_DATA}/training/training_set.txt -f ${DIR_DATA}/training/training_set_full.txt -o ${DIR_DATA}/training/ensemble -s ${DIR_DATA}/training/predictor_normal_stats.txt

	# echo "Testing prediction ... "
	# python ${DIR_SCRIPTS}/crc_prediction_ensemble.py -i ${DIR_DATA}/testing/testing_set.txt -f ${DIR_DATA}/testing/testing_set_full.txt -m ${DIR_DATA}/training/ensemble -s ${DIR_DATA}/training/predictor_normal_stats.txt
	# echo ""

done

#!/usr/bash
DIR_SCRIPTS=/Users/KANG/geneoscopy_dev/scripts

##### INPUT VARIABLES #####

DIR_DATA=/Users/KANG/geneoscopy_dev/data/run_proj_batch1-17_3
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

echo "Working directory:" $DIR_DATA
echo "Label groups:" $GROUP
# echo "Normalized expr data:" $NORMALIZED_CHIPDATA_FULL
echo ""


echo "Preparing normalized data ..."
python ${DIR_SCRIPTS}/prepare_normalized_expr_data.py -i $NORMALIZED_CHIPDATA_FULL -l $GENE_FILTER_LST -c 2 -o1 ${DIR_DATA}/chipdata_geneset_x_valid_chips.txt -o2 ${DIR_DATA}/chipdata_geneset_x_valid_chips_full.txt
#### Use full list of probesets
# cp ${DIR_DATA}/chipdata_geneset_x_valid_chips_full.txt ${DIR_DATA}/chipdata_geneset_x_valid_chips.txt
####

echo "Splitting training/testing sets ... "
mkdir -p ${DIR_DATA}/training
mkdir -p ${DIR_DATA}/testing
python ${DIR_SCRIPTS}/split_expr_train_vs_test.py -v $VALID_CHIPS -i ${DIR_DATA}/chipdata_geneset_x_valid_chips.txt -tr0 ${DIR_DATA}/training/chipdata.txt -tr1 ${DIR_DATA}/training/valid_chips.txt -te0 ${DIR_DATA}/testing/chipdata.txt -te1 ${DIR_DATA}/testing/valid_chips.txt
python ${DIR_SCRIPTS}/split_expr_train_vs_test.py -v $VALID_CHIPS -i ${DIR_DATA}/chipdata_geneset_x_valid_chips_full.txt -tr0 ${DIR_DATA}/training/chipdata_full.txt -te0 ${DIR_DATA}/testing/chipdata_full.txt


# echo "Analyzing DE genes ... "
# 	Rscript ${DIR_SCRIPTS}/de_analysis.r ${DIR_DATA}/training/chipdata.txt ${DIR_DATA}/training/valid_chips.txt $GROUP 1 0 ${DIR_DATA}/training/top_de_genes.txt

# python ${DIR_SCRIPTS}/combine_de_genes.py ${DIR_DATA}/training/top_de_genes.constrained_gene_set.txt ${DIR_DATA}/training/top_de_genes.all_gene_symbol.txt 0.005 0.0005 ${DIR_DATA}/training/top_de_genes.txt

# echo "Analyzing CV DE genes ... "
# mkdir -p ${DIR_DATA}/training/cv_folds
# python ${DIR_SCRIPTS}/split_cv_folds.py -i ${DIR_DATA}/training/chipdata.txt -v ${DIR_DATA}/training/valid_chips.txt -f $CV_FOLDS -g $GROUP -o ${DIR_DATA}/training/cv_folds
# for ((i=1; i<=$CV_FOLDS; i++)); do
# 	echo "    cv fold "$i
# 	Rscript ${DIR_SCRIPTS}/de_analysis.r ${DIR_DATA}/training/cv_folds/chipdata_random_${i}.txt ${DIR_DATA}/training/cv_folds/valid_chips_random_${i}.txt $GROUP 1 0 ${DIR_DATA}/training/cv_folds/top_de_genes_${i}.txt
# done


# echo "Analyzing PvsN, CvsN DE genes ... "
# Rscript ${DIR_SCRIPTS}/de_analysis.r ${DIR_DATA}/training/chipdata.txt ${DIR_DATA}/training/valid_chips.p_vs_n.txt N_vs_C 1 0 ${DIR_DATA}/training/top_de_genes.p_vs_n.txt
# Rscript ${DIR_SCRIPTS}/de_analysis.r ${DIR_DATA}/training/chipdata.txt ${DIR_DATA}/training/valid_chips.c_vs_n.txt N_vs_C 1 0 ${DIR_DATA}/training/top_de_genes.c_vs_n.txt


NUM_TOP_GENES_LIST=( 50 75 100 125 150 175 200 250 300 400 500 1000 )
for NUM_TOP_GENES in "${NUM_TOP_GENES_LIST[@]}"; do
	echo "#################################"
	echo "top genes -->" $NUM_TOP_GENES
	echo "#################################"
	echo ""
	# python ${DIR_SCRIPTS}/summarize_cv_de_genes.py -i ${DIR_DATA}/training/cv_folds -f $CV_FOLDS -t $NUM_TOP_GENES -o ${DIR_DATA}/training/top_de_genes.txt

	# head -$((NUM_TOP_GENES/2+1)) ${DIR_DATA}/training/top_de_genes.p_vs_n.txt > ${DIR_DATA}/training/top_de_genes.txt
	# head -$((NUM_TOP_GENES/2+1)) ${DIR_DATA}/training/top_de_genes.c_vs_n.txt | tail -$((NUM_TOP_GENES/2)) >> ${DIR_DATA}/training/top_de_genes.txt


	echo "Analyzing DE genes ... "
	Rscript ${DIR_SCRIPTS}/de_analysis.r ${DIR_DATA}/training/chipdata.txt ${DIR_DATA}/training/valid_chips.txt $GROUP 1 0 ${DIR_DATA}/training/top_de_genes.txt
	cp ${DIR_DATA}/training/top_de_genes.txt ${DIR_DATA}/training/tmp.txt
	head -$((NUM_TOP_GENES+1)) ${DIR_DATA}/training/tmp.txt > ${DIR_DATA}/training/top_de_genes.txt
	rm ${DIR_DATA}/training/tmp.txt


	ML_MODELS=(random_forest svm grad_boosting adaboost)
	# ML_MODELS=( svm )
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
		python ${DIR_SCRIPTS}/crc_training.py -i ${DIR_DATA}/training/training_set.txt -f ${DIR_DATA}/training/training_set_full.txt -a $ML_MODEL -o ${DIR_DATA}/training/${ML_MODEL} -s ${DIR_DATA}/training/predictor_normal_stats.txt

		echo "Testing prediction ... "
		python ${DIR_SCRIPTS}/crc_prediction.py -i ${DIR_DATA}/testing/testing_set.txt -f ${DIR_DATA}/testing/testing_set_full.txt -a $ML_MODEL -m ${DIR_DATA}/training/${ML_MODEL}/${ML_MODEL}_model.pkl -s ${DIR_DATA}/training/predictor_normal_stats.txt
		echo ""
	done


	# echo "Preparing training/testing datasets ... "
	# python ${DIR_SCRIPTS}/convert_expr_for_ml.py -t ${DIR_DATA}/training/top_de_genes.txt -i ${DIR_DATA}/training/chipdata.txt -o ${DIR_DATA}/training/training_set.txt
	# python ${DIR_SCRIPTS}/convert_expr_for_ml.py -t ${DIR_DATA}/training/top_de_genes.txt -i ${DIR_DATA}/testing/chipdata.txt -o ${DIR_DATA}/testing/testing_set.txt
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

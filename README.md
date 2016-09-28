# geneoscopy_dev

This project applies feature selection to highly confident differentailly expressed (DE) transcripts, and employs supervised machine learning algorithms (e.g. Random Forest, Support Vector Machine, Neural Network, Gradient Boosting) to train model and predict patient labels (i.e. cancer, polys, normal).

##### Git URL https://github.com/dragontt/geneoscopy_dev.git

### Requirement

##### Intall Python: [dowload]

https://www.python.org/ftp/python/2.7.10/python-2.7.10-macosx10.6.pkg

##### Install pip: [command line]
```python
python get-pip.py
```
##### Intall Python packages: [command line]
```shell
pip install numpy
pip install scipy
pip install matplotlib
pip install scikit-learn
```
If permission denied, using the following command:
```shell
sudo pip install <pkg>
```

### Download Resouces:

Download EE datasets to 
```shell 
data/<project_name>
```

### Train Machine Learning Models

This stage trains a supervised model using normalized microarray data. The following example shows how to train a Random Forest model from 191 samples at p-value of 0.005 and outputs the pkl formatted model. This pkl file can be used directly for EE prediction without retraining. 

##### Example:
```python
python scripts/train_ee.py -a random_forest -i resources/ee_data/Exp\ 1\ modeling\ 191N\ p005.txt -o output/ee_data/Exp_1_modeling_191N_p005_random_forest/
```
##### For additional help: 
```python
python scripts/train_ee.py -h
```

### Apply Machine Learning Models For Prediction

This stage applies a supervised model that uses microarray data to make EE category preidction. The following example shows how to use a trained Random Forest model to predict EE categories from 50 samples of p-value of 0.005. 

##### Example:
```python
python scripts/predict_ee.py -a random_forest -i resources/ee_data/Exp\ 1\ validation\ 50N\ for\ modeling\ 191N\ p005.txt -m output/ee_data/Exp_1_modeling_191N_p005_random_forest/rf_model.pkl
```
##### For additional help: 
```python
python scripts/predict_ee.py -h
```

### Clustering Differential Expression (log2 fold change)

This module computes the differential expression (DE) against the median of gene expressions of the healhty samples, and computes the heirachical clustering of the DE, along with dendrogram display.

##### Example:
```python
python scripts/display_expr.py -o clustering -i resources/ee_data/Exp\ 1\ modeling\ 191N\ p005.txt
```

### Perform Principal Component Analysis (PCA)

This module computes PCA with whitening and display first 3 components.

##### Example:
```python
python scripts/display_expr.py -o pca -i resources/ee_data/Exp\ 1\ modeling\ 191N\ p005.txt
```

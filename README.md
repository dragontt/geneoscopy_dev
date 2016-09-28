# geneoscopy_dev

This project applies feature selection to highly confident differentailly expressed (DE) transcripts, and employs supervised machine learning algorithms (e.g. Random Forest, Support Vector Machine, Neural Network, Gradient Boosting) to train model and predict patient labels (i.e. cancer, polys, normal).

### Requirement

##### Install pip:
```python
python get-pip.py
```
##### Intall Python packages:
```shell
sudo pip install numpy
sudo pip install scipy
sudo pip install matplotlib
sudo pip install scikit-learn

### Data Resouces

Put normalized expression matrix, quality control file, sample sheet file into 
```shell 
data/<project_name>
```

### Run pipeline

This pipeline includes assessment of sample quality, split of training/testing sets, analysis of significantly DE genes in training set, training of ML models, and testing of prediction quality. Input arguments are scripts directory, data resource directory, normalized expression matrix file, quality control table, sample sheet. Make changes of those arguments in scripts/run_pipeline.sh. 

##### Example:
```shell
bash scripts/run_pipeline.sh
```


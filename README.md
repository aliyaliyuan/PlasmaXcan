## Environment Setup

```
git clone https://github.com/aliyaliyuan/PlasmaXcan.git
cd PlasmaXcan
conda env create -f env.yml
conda activate scPrediXcan
```

## Data Acquisition 

You must start with a Seurat object from single-cell data. I used liver data acquired from single-cell eQTL experiments. The dataset I used is available on GEO at https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE289173 .

## Pre-processing Data
After acquiring the Seurat object from your single-cell data, epigenetic features from Enformer (https://github.com/hakyimlab/scPrediXcan/blob/master/Scripts/Enformer_epigenomic_features/Geuvadis_individuals_epigenome.txt), it is
time to pre-process the data. 

Start with preprocess.R to generate the correct tsv file that has the observed gene expression data from the single-cell data you are using. This tsv file is The "IbG_" tsv file produced by this script contains the following columns: gene_name, individual_1, individual_2, individual_3, etc. It is the individual-by-gene matrix that shows the expression value of each gene per individual. 

Then, go to the ctPred directory and run preprocess_IbG.py to convert that tsv into the appropriate csv for ctPred. preprocess.py will add the epigenetic feature information from Enformer and will calculate 
mean expression and add column ("percentile"), which is rank-based percentiles. The output csv will have the following columns: gene_name, chromo, feature_1, ..., feature_5313, mean_expression, percentile. 

## ctPred

ctPred.py will generate the weights needed to generate the gene expression predictions. This will also produce a scatterplot that compares the observed genetic expression versus the predicted genetic expression of the model. If the data is of good quality, expect to see a positive linear association between observed and predicted expression. 

Below is an example of the scatterplot of the expression comparison of the periportal hepatoctye model I generated: 
<img width="2074" height="1638" alt="scPred_HG00096_periportal" src="https://github.com/user-attachments/assets/4b3f2f10-41a1-444d-b21a-53514811c1be" />

## References

https://github.com/hakyimlab/scPrediXcan
  ctPred.py is based on Tutorial.ipynb from this repository. Also, the env.yml was sourced from this repository 
  
https://github.com/hakyimlab/scPrediXcan/blob/master/Scripts/Enformer_epigenomic_features/Geuvadis_individuals_epigenome.txt
  Links to the Zenodo databases with the epigenetic features derived from deep learning model Enformer. They come in the form of .h5 files. 

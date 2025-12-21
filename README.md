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

Start with preprocess.R to generate the correct tsv file that has the observed gene expression data from the single-cell data you are using. 

Then, go to the ctPred directory and run preprocess.py to convert that tsv into the appropriate csv for ctPred. preprocess.py will add the epigenetic feature information from Enformer and will convert
mean expression into rank-based percentiles. 

## ctPred

ctPred.py will generate the weights needed to generate the gene expression predictions. 

## References

https://github.com/hakyimlab/scPrediXcan
https://github.com/hakyimlab/scPrediXcan/blob/master/Scripts/Enformer_epigenomic_features/Geuvadis_individuals_epigenome.txt

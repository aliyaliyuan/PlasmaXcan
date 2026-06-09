## Environment Setup

```
git clone https://github.com/aliyaliyuan/PlasmaXcan.git
cd PlasmaXcan
conda env create -f env.yml
conda activate scPrediXcan
```
You will need to activate the conda env to run any script in this pipeline. Do not forget (but you probably will, so if you get an error, just check if it's because you forgot to activate your conda env). 

## Data Acquisition 

I started with a Seurat object https://zenodo.org/records/14586466 . I used liver data acquired from single-cell eQTL experiments. The Seurat object I used was dervied from single-cell liver data available on GEO at https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE289173 . For this pipeline, you need a Seurat object (.Rds) format, but .h5ad are also common single-cell data formats (especially if you are using data from the Human Cell Atlas). Conversion script pending. 

## Pre-processing Data

**Obtain the following data files:**
- metadata.csv from Geuvadis Enformer-dervied epigenetic features: https://zenodo.org/records/15477910
- Seurat object (.rds)

**Run the following scripts in this order**
- pseudobulk_sample.Rmd (requires editing script to change path files, each part is clearly indicated)
- make_ctPred_training_file.py (requires editing script to change path files, just 3 lines at beginning)

pseudobulk_sample.Rmd converts the single-cell data into pseudobulk data, which means it generates the gene expression for each gene in each cell population per individual sample in the data you chose. In my analysis, I only used controls so that I could generate models based on healthy cells. make_ctPred_training_file.py generates a ctPred training file (includes gene_name, chromosome, 5313 epigenomic features derived from deep learning model Enformer, mean_expression, and percentile). Percentile ranks the measured expression of each gene relative to the others. The features from Enformer are predicted locations of enhancers and promoters, which helps train the ctPred model by giving it more information to "learn" from. 

## ctPred

ctPred.py will generate the weights needed to generate the gene expression predictions. This will also produce a scatterplot that compares the observed genetic expression versus the predicted genetic expression of the model. If the data is of good quality, expect to see a positive linear association between observed and predicted expression. I got an average Spearman correlation of 0.80 for all 17 liver cell neural network models I generated. 

Below is an example of the scatterplot of the expression comparison of the periportal hepatoctye model I generated: 
<img width="2074" height="1638" alt="scPred_HG00096_periportal" src="https://github.com/user-attachments/assets/4b3f2f10-41a1-444d-b21a-53514811c1be" />

## l-ctPred



## References

https://github.com/hakyimlab/scPrediXcan
  ctPred.py is based on Tutorial.ipynb from this repository. Also, the env.yml was sourced from this repository 
  
https://github.com/hakyimlab/scPrediXcan/blob/master/Scripts/Enformer_epigenomic_features/Geuvadis_individuals_epigenome.txt
  Links to the Zenodo databases with the epigenetic features derived from deep learning model Enformer. They come in the form of .h5 files. 

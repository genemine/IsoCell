# IsoCell
## 1. IsoCell
An approach to enhance single cell clustering by integrating isoform-level expression through orthogonal projection
### 1.1 Description
IsoCell is an approach to enhance single cell clustering by integrating isoform-level expression through orthogonal projection. First, the orthogonal projection matrix of genes was constructed. And then some isoforms were selected and projected orthogonally to genes to integrate with genes. Based on fifteen real world scRNA-seq datasets. It was found that the integration of alternative splicing information led to better clustering performances in most datasets, compared with using only gene-level expression data.

## 2. Download & Install
IsoCell is implemented as an R package, which is freely available for non-commercial use. You can download from github directly:
```
> install.packages("devtools")
> 
> library(devtools)
> devtools::install_github("genemine/IsoCell")
```

## 3. Usage
Notes: IsoCell was was tested on linux and Windows.

Using IsoCell is very simple. Just follow the steps below:

Step 1: open your R or Rstudio

Step 2: in the R command window, run the following command to load the R package
```
> library(IsoCell)
```

Step 3: in R command window, run the following command to see the help document for running IsoCell. Then, you should be able to see a help page.
```
> ?IsoCell
```

At the end of the help page, there is an example code. Copy these codes to command to run as follows:  
Step 4: load demo data containing gene-level expression, isoform-level expression and the corresponding cell label. The cells in the demo data come from 3 different types.
```
> data(demo)
> # gene_expr (323*16): a expression data matrix, genes in rows and samples in columns
> # isoform_expr (2316*16): a expression data matrix, isoforms in rows and samples in columns
> # label (1*16): corresponding cell label
```

Step 5: Running IsoCell function  
Parameters:
```
> # gene_expr:  A gene-level expression data matrix, genes in rows and samples in columns.
> # isoform_expr:  An isoform-level expression data matrix, isoforms in rows and samples in columns.
> # t (default: 0.1): The retention ratio of isoforms, isoforms whose residual vector length is ranked in the top t percent are retained.
> # theta (default: 0.9): The threshold of the cumulative variance contribution rate of the principal components: keep the principal components whose cumulative variance contribution rate is greater than theta.
> # s (default: FALSE): An optional parameter that can be used to control the number of principal components retained, instead of theta. Note that when the parameter s is selected, the parameter theta will no longer work. (If s is selected, then s must be an integer and less than min(dim(gene_expr))).
```
Return:
```
> # P: The orthogonal projection matrix of genes.
> # len_residual_vector: The length of the residual vector of normalized isoform.
> # selisoform: The isoforms IsoCell select for the final combining.
> # combined_data: The final combined data for clustering.
```
Run code:
```
> result = IsoCell(gene_expr,isoform_expr)
```

## 4. Performance evaluation
Step 1: Clustering test (An example with SIMLR is as follows):
```
> library(SIMLR)
> # clustering data with only gene-level expression
> example1 = SIMLR(gene_expr,4)
> cluster1 = example1$y$cluster
> 
> # clustering data with combined data
> example2 = SIMLR(result$combined_data,4)
> cluster2 = example2$y$cluster
```
Step 2: Evaluate clustering performance [evaluate.R](https://github.com/genemine/IsoCell/blob/main/code/evaluate.R)

Parameters:
```
> # truelabel: A numeric vector of true labels of each sample.
> # prelabel: A numeric vector of predicted labels of each sample.
```
Return:
```
> # NMI: Value of normalized mutual information.
> # ARI: Value of adjusted rand index.
```
Run code:
```
> evaluate(truelabel = label, prelabel = cluster1)
> NMI      ARI 
> 0.55209  0.31373
>
> evaluate(truelabel = label, prelabel = cluster2)
> NMI      ARI 
> 0.58609  0.37500
```

## 5. Contact
If any questions, please do not hesitate to contact us at: 

Hongdong Li, hongdong@csu.edu.cn

Jianxin Wang, jxwang@csu.edu.cn


## 6. How to cite?
If you use this tool, please cite the following work.

Yingyi Liu, Hongdong Li, Yunpei Xu, Jianxin Wang, IsoCell: An approach to enhance single cell clustering by integrating isoform-level expression through orthogonal projection, 2020, submitted

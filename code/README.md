# Analysis
After getting the results of IsoCell (result = IsoCell(gene_expr,isoform_expr)), relevant evaluation and analysis can be performed on the results.

## 1. Fold Change
Based on the known cell labels of each dataset, we calculated the fold change (FC) value of expression between every pair of cell types. The calculation of fold change requires gene/isoform expression data (path/to/expression), each row is gene/isoform, each column is a cell sample, the true label of the dataset (path/to/label) and the number of cell clusters in the dataset (clusters).
[fc.py](https://github.com/genemine/IsoCell/blob/main/code/fc.py)
```
 python fc.py path/to/expression path/to/label clusters

 eg:python fc.py E:/isocell/gene.csv E:/isocell/labels.csv 4
```
(seliso)

## 2. p-value
To statistically tested whether gene/isoform is differentially expressed between cell types, the Mann Whitney U test is used for statistical testing and the p-value is obtained. The calculation of fold change also requires gene/isoform expression data (path/to/expression), each row is gene/isoform, each column is a cell sample, the true label of the dataset (path/to/label) and the number of cell clusters in the dataset (clusters).
[pvalue.py](https://github.com/genemine/IsoCell/blob/main/code/pvalue.py)
```
 python pvalue.py path/to/expression path/to/label clusters

 eg:python pvalue.py E:/isocell/gene.csv E:/isocell/labels.csv 4
```
(fdr)


## 3. Performance comparison between different isoform retention rate 
Step 1: Run IsoCell with different isoform retention rate (t), for example:
```
> library(IsoCell)
> data(demo)
> result_0.5 = IsoCell(gene_expr,isoform_expr,t=0.5)
> result_0.7 = IsoCell(gene_expr,isoform_expr,t=0.7)
```
Step 2: Clustering test (An example with SIMLR is as follows):
```
> library(SIMLR)
> # clustering data with only gene-level expression
> example1 = SIMLR(result_0.5$combined_data,4)
> cluster1 = example1$y$cluster
> 
> # clustering data with combined data
> example2 = SIMLR(result_0.7$combined_data,4)
> cluster2 = example2$y$cluster
```
Step 4: Evaluate clustering performance [evaluate.R](https://github.com/genemine/IsoCell/blob/main/code/evaluate.R)
```
> evaluate(truelabel = label[,1], prelabel = cluster1)
> NMI      ARI 
> 0.55209  0.31373
>
> evaluate(truelabel = label[,1], prelabel = cluster2)
> NMI      ARI 
> 0.58609  0.37500
```

## 4. Performance comparison between using the residual matrix and the original expression of the selected 10% isoforms
Step 1: Run IsoCell
```
> library(IsoCell)
> data(demo)
> result = IsoCell(gene_expr,isoform_expr)
```
Step 2: Combine the gene expression matrix to both the selected 10% isoforms and the corresponding residual matrix
```
> # Combine the gene expression matrix to the original expression of the selected 10% isoforms
> seliso_oriexpr=isoform_expr[result$selisoform,]
> gene_seliso_oriexpr=bind_rows((gene_expr),(seliso_oriexpr))\

> # Combine the gene expression matrix to the residual matrix of the selected 10% isoforms
> seliso_resexpr=result$combined_data[result$selisoform,]
> gene_seliso_oriexpr=bind_rows((gene_expr),(seliso_resexpr))
```
Step 3: Clustering test (An example with SIMLR is as follows):
```
> library(SIMLR)
> # clustering data with only gene-level expression
> example1 = SIMLR(gene_seliso_oriexpr,4)
> cluster1 = example1$y$cluster
> 
> # clustering data with combined data
> example2 = SIMLR(gene_seliso_oriexpr,4)
> cluster2 = example2$y$cluster
```
Step 4: Evaluate clustering performance [evaluate.R](https://github.com/genemine/IsoCell/blob/main/code/evaluate.R)
```
> evaluate(truelabel = label[,1], prelabel = cluster1)
> NMI      ARI 
> 0.55209  0.31373
>
> evaluate(truelabel = label[,1], prelabel = cluster2)
> NMI      ARI 
> 0.58609  0.37500
```

## 5. Performance comparison between using the combined data and gene data
Step 1: Run IsoCell
```
> library(IsoCell)
> data(demo)
> result = IsoCell(gene_expr,isoform_expr)
```
Step 2: Clustering test (An example with SIMLR is as follows):
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
Step 3: Evaluate clustering performance [evaluate.R](https://github.com/genemine/IsoCell/blob/main/code/evaluate.R)
```
> evaluate(truelabel = label[,1], prelabel = cluster1)
> NMI      ARI 
> 0.55209  0.31373
>
> evaluate(truelabel = label[,1], prelabel = cluster2)
> NMI      ARI 
> 0.58609  0.37500
```

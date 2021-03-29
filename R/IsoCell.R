#' Isocell
#'
#' intergrating isoform-level expression through orthogonal projection
#' @param gene_expr  A gene-level expression data matrix, genes in rows and samples in columns.
#' @param isoform_expr  An isoform-level expression data matrix, isoforms in rows and samples in columns.
#' @param t The retention ratio of isoforms, isoforms whose residual vector length is ranked in the top t percent are retained.
#' @return P: The orthogonal projection matrix of genes.
#' @return len_residual_vector: The length of the residual vector for each normalized isoform vector.
#' @return selisoform: The isoforms IsoCell select for the final combining.
#' @return combined_data: The final combined data for clustering.
#' @keywords orthogonal projection
#' @export
#' @author Hongdong Li, lhdcsu@gmail.com, Central South University.
#' @examples
#' data(demo)
#' result = IsoCell(gene_expr,isoform_expr)



IsoCell<-function(gene_expr,isoform_expr,t=0.1){
  
  if(ncol(gene_expr) == ncol(isoform_expr)){
    
    nsample = ncol(gene_expr)
    ngene = nrow(gene_expr)
    nisoform = nrow(isoform_expr)
    
    #step 1: construct orthogonal projection matrix
    cat("Construct orthogonal projection matrix.\n")
    
    #run PCA for gene_expr
    #cat("Run PCA for gene-level expression.\n")
    
    pca_prcomp <- prcomp(t(gene_expr))
    pca_info = summary(pca_prcomp)$importance
    
    estimate = which(pca_info['Cumulative Proportion',]>0.9)[1]
    #PCA1=pca1$loadings[,1]
    #screeplot(pca1,type="line") 
    
    pca_data <- predict(pca_prcomp)
    gene_pc = pca_data[,1:estimate]
    
    #orthogonal projection matrix
    library(MASS)
    projm = (gene_pc)%*%(ginv(t(gene_pc)%*%gene_pc))%*%(t(gene_pc))
    
    
    #step 2: isoform selection
    cat("Perform isoform selection.\n")
    
    #isoform normlazition, the length of each normalized isoform is 1
    isoform_expr = as.matrix(t(isoform_expr))
    isoform_len = sqrt(apply(isoform_expr*isoform_expr,2,sum))#colSums
    norm_isofrom = as.matrix(t(t(isoform_expr)/isoform_len))
    
    #select isoform according to its residual vector
    projnormiso = projm %*% norm_isofrom
    r_normiso = norm_isofrom - projnormiso
    
    len_r_normiso = sort (sqrt(apply(r_normiso*r_normiso,2,sum)),decreasing = T)
    #select isoforms ranked in the top t percent 
    selisoform = names(len_r_normiso[1:(nisoform*t)])
    isoform_expr = isoform_expr[,selisoform]
    
    #step3: data combination
    
    projiso = projm %*% isoform_expr
    r_isoform = isoform_expr - projiso
    
    cat("Combine data.\n")
    gene_expr = as.matrix(gene_expr)
    combined_data = cbind(t(gene_expr),r_isoform)
    combined_data = t(combined_data)
    
    nfeatrue=nrow(data)
    
  }else{cat("Input error.\n")}
  
  return(list(P = projm,
              len_residual_vector = len_r_normiso,
              selisoform = selisoform,
              combined_data = combined_data))
}
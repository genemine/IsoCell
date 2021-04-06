#' Isocell
#'
#' intergrating isoform-level expression through orthogonal projection
#' @param gene_expr  A gene-level expression data matrix, genes in rows and samples in columns.
#' @param isoform_expr  An isoform-level expression data matrix, isoforms in rows and samples in columns.
#' @param t The retention ratio of isoforms, isoforms whose residual vector length is ranked in the top t percent are retained. The default is 0.1.
#' @param theta The threshold of the cumulative variance contribution rate of the principal components: keep the principal components whose cumulative variance contribution rate is greater than \code{theta}. The default is 0.9.
#' @param s An optional parameter that can be used to control the number of principal components retained, instead of \code{theta}. Note that when the parameter \code{s} is selected, the parameter \code{theta} will no longer work. (If \code{s} is selected, then \code{s} must be an integer and less than min(dim(gene_expr)).) The default is FASLE.
#' @return \code{P}: The orthogonal projection matrix of genes.
#' @return \code{len_residual_vector}: The length of the residual vector for each normalized isoform vector.
#' @return \code{selisoform}: The isoforms IsoCell select for the final combining.
#' @return \code{combined_data}: The final combined data for clustering.
#' @keywords orthogonal projection
#' @export
#' @author Hongdong Li, lhdcsu@gmail.com, Central South University.
#' @examples
#' data(demo)
#' result = IsoCell(gene_expr,isoform_expr)



IsoCell<-function(gene_expr,isoform_expr,t=0.1,theta=0.9,s=FALSE){

    if(ncol(gene_expr) == ncol(isoform_expr)){

        nsample = ncol(gene_expr)
        ngene = nrow(gene_expr)
        nisoform = nrow(isoform_expr)

        #step 1: construct orthogonal projection matrix
        cat("Construct orthogonal projection matrix.\n")

        #run PCA for gene_expr
        #cat("Run PCA for gene-level expression.\n")
        gene_pc = princomponent(gene_expr,theta,s)

        #orthogonal projection matrix
        projm=(gene_pc)%*%(t(gene_pc))
        #residual matrix
        isoform_expr = as.matrix(t(isoform_expr))
        projiso = projm %*% isoform_expr
        r_isoform = isoform_expr - projiso


        #step 2: isoform selection
        cat("Perform isoform selection.\n")

        #isoform normlazition, the length of each normalized isoform is 1
        isoform_len = sqrt(apply(isoform_expr*isoform_expr,2,sum))#colSums
        norm_isofrom = as.matrix(t(t(isoform_expr)/isoform_len))

        #select isoform according to its residual vector
        r_normiso=t(t(r_isoform)/isoform_len)

        len_r_normiso = sort (sqrt(apply(r_normiso*r_normiso,2,sum)),decreasing = T)
        #select isoforms ranked in the top t percent
        selisoform = names(len_r_normiso[1:(nisoform*t)])
        r_isoform=r_isoform[,selisoform]


        #step3: data combination
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


princomponent<-function(gene_expr,theta,s){
    #get the principal components of gene-level expression
    #gene_expr: A gene-level expression data matrix, genes in rows and samples in columns.
    #theta: The parameter \code{theta} is used to control the number of principal components retained. If \code{theta} is greater than 1, then keep \code{theta} principal components; if \code{theta} is less than 1, keep the principal components whose cumulative variance contribution rate reaches \code{theta}. The default is 0.9，which means that the principal components whose cumulative variance contribution rate reaches 90% are retained.

    if(s == FALSE){

        pca_prcomp <- prcomp_u(t(gene_expr))
        pca_info = summary(pca_prcomp)$importance

        estimate = which(pca_info['Cumulative Proportion',]>theta)[1]
        gene_pc = pca_prcomp$u[,1:estimate]

    }else{

        library(irlba)
        library(MASS)
        pca_prcomp <- prcomp_irlba(t(gene_expr),s)
        pca_info = summary(pca_prcomp)$importance
        u=pca_prcomp$x %*% ginv(diag(pca_prcomp$sdev * sqrt(ncol(gene_expr)-1)))

        gene_pc = u
    }

    return(gene_pc)
}


# Modify the prcomp function to get the u value in svd
prcomp_u<-function (x, retx = TRUE, center = TRUE, scale. = FALSE, tol = NULL,
    rank. = NULL, ...)
{
    chkDots(...)
    x <- as.matrix(x)
    x <- scale(x, center = center, scale = scale.)
    cen <- attr(x, "scaled:center")
    sc <- attr(x, "scaled:scale")
    if (any(sc == 0))
    stop("cannot rescale a constant/zero column to unit variance")
    n <- nrow(x)
    p <- ncol(x)
    k <- if (!is.null(rank.)) {
        stopifnot(length(rank.) == 1, is.finite(rank.), as.integer(rank.) >
            0)
        min(as.integer(rank.), n, p)
    }
    else min(n, p)
    s <- svd(x, nv = k)
    j <- seq_len(k)
    s$d <- s$d/sqrt(max(1, n - 1))
    if (!is.null(tol)) {
        rank <- sum(s$d > (s$d[1L] * tol))
        if (rank < k) {
            j <- seq_len(k <- rank)
            s$v <- s$v[, j, drop = FALSE]
        }
    }
    dimnames(s$v) <- list(colnames(x), paste0("PC", j))
    r <- list(sdev = s$d, rotation = s$v, center = if (is.null(cen)) FALSE else cen,
        scale = if (is.null(sc)) FALSE else sc)
    if (retx)
    r$u <- s$u
    class(r) <- "prcomp"
    r
}

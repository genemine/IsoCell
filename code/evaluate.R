#' evaluate
#'
#' using NMI and ARi to evaluate clustering performance
#' @param gene_expr  A gene-level expression data matrix, genes in rows and samples in columns.
#' @param isoform_expr  An isoform-level expression data matrix, isoforms in rows and samples in columns.
#' @param t The retention ratio of isoforms, isoforms whose residual vector length is ranked in the top t percent are retained.
#' @return NMI: Normalized mutual information value.
#' @return ARI: Adjusted rand index value.
#' @keywords clustering performance evaluation
#' @export
#' @examples
#' truelabel = sample(1:3, size=10, replace=TRUE)
#' predlabel = sample(1:3, size=10, replace=TRUE)
#' evaluate(truelabel,predlabel)

evaluate<-function(truelabel,prelabel){
  if(length(truelabel)!=length(prelabel))
    stop("truelabel and prelabel must have the same length")
  total = length(truelabel)
  x_ids = unique(truelabel)
  y_ids = unique(prelabel)
  #Mutual information
  MI = 0.0
  for(idx in x_ids){
    for(idy in y_ids){
      idxOccur = which(truelabel==idx)
      idyOccur = which(prelabel==idy)
      idxyOccur = intersect(idxOccur,idyOccur)
      if(length(idxyOccur)>0){
        MI = MI + (length(idxyOccur)/total)*log2((length(idxyOccur)*total)/(length(idxOccur)*length(idyOccur)));
      }
    }
  }
  
  #Normalized Mutual information
  Hx = 0; #Entropies
  for(idx in x_ids){
    idxOccurCount = length(which(truelabel==idx));
    Hx = Hx - (idxOccurCount/total) * log2(idxOccurCount/total);
  }
  Hy = 0;#Entropies
  for(idy in y_ids){
    idyOccurCount = length(which(prelabel==idy));
    Hy = Hy - (idyOccurCount/total) * log2(idyOccurCount/total);
  }
  nmi = 2 * MI / (Hx+Hy)
  
  #(adjusted) Rand Index
  tab = table(truelabel,prelabel)
  conv_df = as.data.frame.matrix(tab)
  n <- sum(tab)
  ni <- apply(tab, 1, sum)
  nj <- apply(tab, 2, sum)
  n2 <- choose(n, 2)
  nis2 <- sum(choose(ni[ni > 1], 2))
  njs2 <- sum(choose(nj[nj > 1], 2))
  ri = 1 + (sum(tab^2) - (sum(ni^2) + sum(nj^2))/2)/n2
  ari=c(sum(choose(tab[tab > 1], 2)) - (nis2 * njs2)/n2)/((nis2 + njs2)/2 - (nis2 * njs2)/n2)
  
  out = c(nmi,ari)
  names(out)=c("NMI","ARI")
  return(out)
}
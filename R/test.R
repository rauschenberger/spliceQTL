#' @export
#' @title
#' Conduct single test
#' 
#' @description
#' This function tests for alternative splicing.
#' 
#' @param Y
#' exon expression\strong{:}
#' matrix with \eqn{n} rows (samples) and \eqn{p} columns (exons)
#' 
#' @param X
#' SNP genotype\strong{:}
#' matrix with \eqn{n} rows (samples) and \eqn{q} columns (SNPs)
#' 
#' @param W 
#'  numeric matrix\strong{:}
#'  a square matrix with as many rows as there are covariates in the independent data set
#'    It represents the correlation structure expected in the independent data 
#'    
#' @param w.type
#'  string\strong{:}
#'    if \code{W=NULL} and w.type is a given string, take the type: "cov" is the only one allowed. 
#'    Then the inner product of the indep data \code{X} is taken
#'    
#' @param map
#' list with names "genes", "exons", and "snps"
#' (output from \code{\link{map.genes}}, \code{\link{map.exons}},
#' and \code{\link{map.snps}})
#' 
#' @param i
#' gene index\strong{:}
#' integer between \eqn{1} and \code{nrow(map$genes)}
#' 
#' @param limit
#' cutoff for rounding \code{p}-values
#' 
#' @param steps
#' size of permutation chunks\strong{:}
#' integer vector
#' 
#' @details
#' The maximum number of permutations equals \code{sum(steps)}. Permutations is
#' interrupted if at least \code{limit} test statistics for the permuted data
#' are larger than the test statistic for the observed data.
#' 
#' @examples
#' NA
#' 
test.single <- function(Y, X, map, i, limit=NULL, steps=NULL, W=NULL, w.type = NULL){
  
  if(is.null(limit)){limit <- 5}
  if(is.null(steps)){steps <- c(10,20,20,50)}
  
  # check input
  if(!is.numeric(limit)){
    stop("Argument \"limit\" is not numeric.",call.=FALSE)
  }
  if(limit<1){
    stop("Argument \"limit\" is below one.",call.=FALSE)
  }
  if(!is.numeric(steps)|!is.vector(steps)){
    stop("Argument \"steps\" is no numeric vector.",call.=FALSE)
  }
  if(sum(steps)<2){
    stop("Too few permutations \"sum(steps)\".",call.=FALSE)
  }
  
  # extract data
  ys <- map$exons[[i]]
  y <- Y[, ys, drop=FALSE]
  xs <- seq(from=map$snps$from[i], to=map$snps$to[i], by=1)
  x <- X[, xs, drop=FALSE]
  rm(Y,X); silent <- gc()
  
  # Compute W if required
  if(is.null(W)) { if(is.null(w.type)) { W <- diag(ncol(x)) } else {
    if(w.type=="cov") { W <- crossprod(x) }
  }  }
  
  # test effects

    tstat <- spliceQTL:::G2.multin(
      dep.data=y, indep.data=x, nperm=steps[1]-1, W=W)$Sg
    for(nperm in steps[-1]){
      tstat <- c(tstat, 
                 spliceQTL:::G2.multin(
                  dep.data=y, indep.data=x, nperm=nperm, W=W )$Sg[-1])
      if(sum(tstat >= tstat[1]) >= limit){break}
    }
    pvalue <- mean(tstat >= tstat[1], na.rm=TRUE)
  
  
  return(pvalue)
}


#' @export
#' @title
#' Conduct multiple tests
#' 
#' @description
#' This function tests for alternative splicing.
#'  
#' @param spec
#' number of cores\strong{:}
#' positive integer
#' 
#' @param min
#' minim chunk size\strong{:}
#' positive integer
#' 
#' @param steps
#' number of iteration chunks\strong{:}
#' positive integer
#' 
#' @inheritParams test.single
#' 
#' @details
#' Automatic adjustment of the number of permutations
#' such that Bonferroni-significant p-values are possible.
#' 
#' @examples
#' NA
#' 
test.multiple <- function(Y, X, map, w.type = NULL, spec=1, min=100, steps=20){
  
  p <- nrow(map$genes)
  
  # permutations
  if(FALSE){ # old
    min <- 5
    max <- p/0.05+1
    limit <- ceiling(0.05*max/p)
    base <- 1.5 # adjust sequence
    from <- log(min,base=base)
    to <- log(max,base=base)
    steps <- c(min,diff(unique(round(base^(seq(from=from,to=to,length.out=20))))))
  }
  
  if(FALSE){ # new
    max <- p/0.05+1
    limit <- ceiling(0.05*max/p)
    steps <- diff(limit^seq(from=1,to=log(max)/log(limit),length.out=pmin(p,steps))) # was (p,20)
    steps <- c(limit,round(steps)) # Or replace "limit" by "minimum # of permutations"!
    steps[steps==0] <- 1
    steps[length(steps)] <- max-sum(steps[-length(steps)])
  }
  
  if(TRUE){
    max <- p/0.05+1
    limit <- ceiling(0.05*max/p)
    steps <- diff(limit^seq(from=log(min),to=log(max)/log(limit),length.out=steps)) # was pmin(p,steps)
    steps[steps<min] <- min
    #for(i in 1:10){
    #    cond <- steps>10^i & steps<10^(i+1)
    #    steps[cond] <- ceiling(steps[cond]/10^i)*10^i 
    #}
    steps = signif(steps,digits=1)
    steps <- steps[cumsum(steps)<=max]
    steps[length(steps)+1] <- max-sum(steps)
  }
  
  if(any(steps<0)){stop("negative step",call.=FALSE)}
  if(max != sum(steps)){stop("invalid step",call.=FALSE)}
  
  if(spec==1){
    set.seed(1)
    pvalue <- lapply(X=seq_len(p),
                     FUN=function(i) spliceQTL::test.single(Y=Y, X=X, map=map, i=i, limit=limit, steps=steps, w.type=w.type))
  } else {
    type <- ifelse(test=.Platform$OS.type=="windows", yes="PSOCK", no="FORK")
    cluster <- parallel::makeCluster(spec=spec, type=type)
    parallel::clusterSetRNGStream(cl=cluster, iseed=1)
    #parallel::clusterExport(cl=cluster,varlist=c("Y","X","map","limit","steps","rho"),envir=environment())
    #parallel::clusterEvalQ(cl=cluster,library(spliceQTL,lib.loc="/virdir/Scratch/arauschenberger/library"))
    pvalue <- parallel::parLapply(cl=cluster, X=seq_len(p), fun=function(i) test.single(Y=Y, X=X, map=map, i=i, limit=limit, steps=steps, w.type=w.type))
    #pvalue <- parallel::parLapply(cl=cluster,X=seq_len(p),fun=function(i) test.trial(y=Y[,map$exons[[i]],drop=FALSE],x=X[,seq(from=map$snps$from[i],to=map$snps$to[i],by=1),drop=FALSE],limit=limit,steps=steps,rho=rho))
    parallel::stopCluster(cluster)
    rm(cluster)
  }
  
  # tyding up
  pvalue <- do.call(what=rbind,args=pvalue)
#  colnames(pvalue) <- paste0("rho=",rho)
  rownames(pvalue) <- map$genes$gene_id
  
  return(pvalue)
}



#--- spliceQTL test functions --------------------------------------------------

#
# Function: get.g2stat.multin
# Computes the G2 test statistic given two data matrices, under a multinomial distribution
# Here we write the test statistic as a function of W, with no particular structure given
# We also compute it in a more efficient way than we used to - using a trace property
# This is used internally by the G2 function
# Inputs: 
#  U = U1 U1^t, where U1 = (I-H)Y, a n*K matrix where n=number obs and K=number multinomial responses possible
#  tau.mat = X' W X, a n*n matrix 
# Output: test statistic (single value)
# 
get.g2stat.multin <- function(U, tau.mat)
{
  g2tstat <- sum( U * tau.mat )
  g2tstat
}


# Function: G2.multin.rho
# This is to compute the G2 test statistic under the assumption that the response follows a multinomial distribution
### Input 
### dep data and indep data with samples on the rows and genes on the columns
### grouping: Either a logical value = F or a matrix with a single column and same number of rows as samples. 
###         Column name should be defined.
###         Contains clinical information of the samples. 
###         Should have two groups only. 
### nperm : number of permutations 
### rho: the null correlation between SNPs
### mu: the null correlation between observations corresponding to different exons and different individuals
#
### Output
### A list containing G2 p.values and G2 test statistics
#
### Example : G2T = G2(dep.data = cgh, indep.data = expr, grouping=F, stand=TRUE, nperm=1000)
### G2 p.values : G2T$G2p
### G2 TS : G2T$$Sg
#
G2.multin.rho <- function(dep.data,indep.data,stand=TRUE,nperm=100,grouping=F,rho=0,mu=0){
  
  nperm = nperm
  ## check for the number of samples in dep and indep data
  
  
  if (nrow(dep.data)!=nrow(indep.data)){
    cat("number of samples not same in dep and indep data","\n")
  }
  
  if(any(abs(rho)>1)){
    cat("correlations rho larger than abs(1) are not allowed")
  }
  
  nresponses <- ncol(dep.data)
  ncovariates <- ncol(indep.data)
  ### centering and standardizing the data are not done in this case
  
  #  dep.data = scale(dep.data,center=T,scale=stand)
  #  indep.data = scale(indep.data,center=T,scale=stand)
  
  #### No  grouping of the samples.
  
  ## Calculate U=(I-H)Y and UU', where Y has observations on rows; also tau.mat=X*W.rho*X', 
  ##   where X has observations on rows and variables on columns
  ##  and W.rho = I + rho*(J-I), a square matrix with as many rows as columns in X
  ## NOTE: this formulation uses X with n obs on the rows and m covariates no the columns, so it is the transpose of the first calculations
  nsamples <- nrow(dep.data)
  n.persample <- rowSums(dep.data)
  n.all <- sum(dep.data)
  H <- (1/n.all)*matrix( rep(n.persample,each=nsamples),nrow=nsamples,byrow=T)
  U <- (diag(rep(1,nsamples)) - H) %*% dep.data
  ## Now we may have a vector of values for rho - so we define tau.mat as an array, with the 3rd index corresponding to the value of rho
  tau.mat <- array(0,dim=c(nsamples,nsamples,length(rho)))
  for(xk in 1:length(rho))  
  {  
    if (rho[xk]==0) { tau.mat[,,xk] <- tcrossprod(indep.data) } 
    else { w.rho <- diag(rep(1,ncovariates)) + rho[xk]*(tcrossprod(rep(1,ncovariates)) -diag(rep(1,ncovariates))  )
    tau.mat[,,xk] <- indep.data %*% w.rho %*% t(indep.data)}
    
  }
  ######################################
  ### NOTES ARMIN START ################
  # all(X %*% t(X) == tau.mat[,,1]) # rho = 0 -> TRUE
  # all(X %*% (t(X) %*% X) %*% t(X) == tau.mat[,,1]) # rho = 1
  # plot(as.numeric(X %*% (t(X) %*% X) %*% t(X)),as.numeric(tau.mat[,,1]))
  ### NOTES ARMIN END ##################
  ######################################
  samp_names = 1:nsamples ## this was rownames(indep.data), but I now do this so that rownames do not have to be added to the array tau.mat
  Sg = get.g2stat.multin(U,mu=mu,rho=rho,tau.mat)
  ### now we will have a vector as result, with one value per combination of values of rho and mu
  #
  ### G2 
  ### Permutations
  # When using permutations: only the rows of tau.mat are permuted
  # To check how the permutations can be efficiently applied, see tests_permutation_g2_multin.R
  
  
  perm_samp = matrix(0, nrow=nrow(indep.data), ncol=nperm)   ## generate the permutation matrix
  for(i in 1:ncol(perm_samp)){
    perm_samp[,i] = samp_names[sample(1:length(samp_names),length(samp_names))]
  }
  
  ## permutation starts - recompute tau.mat  (or recompute U each time)
  for (perm in 1:nperm){
    tau.mat.perm = tau.mat[perm_samp[,perm],,,drop=FALSE]          # permute rows
    tau.mat.perm = tau.mat.perm[,perm_samp[,perm],,drop=FALSE]     # permute columns
    
    Sg = c(Sg,spliceQTL:::get.g2stat.multin(U, mu=mu,rho=rho,tau.mat.perm) )
  }
  
  ########################################################################
  
  #### G2 test statistic
  # *** recompute for a vector of values for each case - just reformat the result with as many rows as permutations + 1,
  # and as many columns as combinations of values of rho and mu
  Sg = matrix(Sg,nrow=nperm+1,ncol=length(mu)*length(rho))
  colnames(Sg) <- paste(rep("rho",ncol(Sg)),rep(1:length(rho),each=length(mu)),rep("mu",ncol(Sg)),rep(1:length(mu),length(rho)) )
  
  ### Calculte G2 pval
  G2p =  apply(Sg,2,spliceQTL:::get.pval.percol) 
  
  return (list(perm = perm_samp,G2p = G2p,Sg = Sg))
}

#
# Function: G2.multin
# This is to compute the G2 test statistic under the assumption that the response follows a multinomial distribution
### Input 
### dep.data and indep.data with samples on the rows and genes on the columns
### nperm : number of permutations 
### W : a matrix that represents the (expected) correlation sttructure between effects of different SNPs on the (same or different) exons
###    If given, it must be a square, symmetric matrix with as many rows as the number of genes/covariates/columns in indep.data
###
### Output
### A list containing G2 p.values and G2 test statistics

### Example : G2T = G2(dep.data = cgh, indep.data = expr, grouping=F, stand=TRUE, nperm=1000)
### G2 p.values : G2T$G2p
### G2 TS : G2T$$Sg


G2.multin <- function(dep.data, indep.data, stand = TRUE, nperm=100, W=NULL){
  
  ## check for the number of samples in dep and indep data
  
  if (nrow(dep.data)!=nrow(indep.data)){
    stop("number of samples not the same in dep and indep data","\n")
  }
  
  ## check if W has compatible dimensions
  if( !is.null(W) ){ if(nrow(W) != ncol(W)) { stop("W must be a square, symmetric matrix","\n") }
    if( nrow(W) != ncol(indep.data) ) { stop("W must have as many rows as the number of columns in indep.data","\n") }
  } else { W <- diag(rep(1,ncol(indep.data))) }
  
  nresponses <- ncol(dep.data)
  ncovariates <- ncol(indep.data)
  ### centering and standardizing the data are not done in this case
  
  #  dep.data = scale(dep.data,center=T,scale=stand)
  #  indep.data = scale(indep.data,center=T,scale=stand)
  
  #### No  grouping of the samples.
  
  ## Calculate U=(I-H)Y and UU', where Y has observations on rows;  
  ##  also tau.mat=X*W*X', where X has observations on rows and variables on columns
  ## NOTE: this formulation uses X with n obs on the rows and m covariates on the columns, so it is the transpose of the first calculations
  nsamples <- nrow(dep.data)
  n.persample <- rowSums(dep.data)
  n.all <- sum(dep.data)
  H <- (1/n.all)*matrix( rep(n.persample,each=nsamples),nrow=nsamples,byrow=T)
  U1 <- (diag(rep(1,nsamples)) - H) %*% dep.data
  U <- tcrossprod(U1)
  tau.mat <- indep.data %*% W %*% t(indep.data)
  samp_names = 1:nsamples ## this was rownames(indep.data), but I now do this so that rownames do not have to be added to the array tau.mat
  Sg = get.g2stat.multin(U,tau.mat)
  #
  ### G2 
  ### Permutations
  # When using permutations: only the rows of tau.mat are permuted
  # To check how the permutations can be efficiently applied, see tests_permutation_g2_multin.R
  
  
  perm_samp = matrix(0, nrow=nsamples, ncol=nperm)   ## generate the permutation matrix
  for(i in 1:ncol(perm_samp)){
    perm_samp[,i] = samp_names[sample(1:length(samp_names),length(samp_names))]
  }
  
  ## permutation starts - recompute tau.mat  (or recompute U each time)
  for (perm in 1:nperm){
    tau.mat.perm = tau.mat[perm_samp[,perm],,drop=FALSE]          # permute rows
    tau.mat.perm = tau.mat.perm[,perm_samp[,perm],drop=FALSE]     # permute columns
    
    Sg = c(Sg, get.g2stat.multin(U,tau.mat.perm) )
  }
  
  
  ########################################################################
  
  #### G2 test statistic
  # *** recompute for a vector of values for each case - just reformat the result with as many rows as permutations + 1,
  # and as many columns as combinations of values of rho and mu
  Sg = t(as.matrix(Sg)) 
  
  ### Calculte G2 pval
  G2p = mean(Sg[1]<= c(Inf , Sg[2:(nperm+1)]))
  names(G2p) = "G2 pval"
  
  return (list(perm = perm_samp,G2p = G2p,Sg = Sg))
}


# Function: get.g2stat.multin
# Computes the G2 test statistic given two data matrices, under a multinomial distribution
# This is used internally by the G2 function
# Inputs: 
#  U = (I-H)Y, a n*K matrix where n=number obs and K=number multinomial responses possible
#  tau.mat = X' W.rho X, a n*n matrix : both square, symmetric matrices with an equal number of rows
# Output: test statistic (single value)
# 
get.g2stat.multin <- function(U, mu, rho, tau.mat){
  g2tstat <- NULL
  for(xk in 1:length(rho))
  {
    for(xj in 1:length(mu))
    {
      if(mu[xj]==0) { g2tstat <- c(g2tstat, sum( diag( tcrossprod(U) %*% tau.mat[,,xk] ) ) )
      } else {
        g2tstat <- c(g2tstat, (1-mu[xj])*sum(diag( tcrossprod(U) %*% tau.mat[,,xk] ) ) + mu[xj]*sum( t(U) %*% tau.mat[,,xk] %*% U )  )
      }
      
    }
  }
  g2tstat
}


# Function: get.pval.percol
# This function takes a vector containing the observed test stat as the first entry, followed by values generated by permutation,
# and computed the estimated p-value
# Input
# x: a vector with length nperm+1
# Output
# the pvalue computed
get.pval.percol <- function(x){
  pval = mean(x[1]<= c(Inf , x[2:length(x)]))
  pval
}

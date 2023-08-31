# The 'full.ase' function carries out an adjacency spectral embedding on a given adjacency matrix 'A'.
# It returns the singular values and the embedded coordinates of the graph.
# Parameters:
#   A:        The adjacency matrix.
#   d:        Dimension of the embedded space.
#   diagaug:  If TRUE, diagonal augmentation is performed on 'A'.
#   approx:   If TRUE, an approximate singular value decomposition method is used.
full.ase <- function(A, d, diagaug=TRUE, approx=TRUE) {
  n <- nrow(A)
  
  # If diagaug is TRUE, perform diagonal augmentation
  if (diagaug) {
    diag(A) <- rowSums(A) / (nrow(A)-1)
  }
  
  # Choose SVD method based on the 'approx' flag
  if (approx) {
    A.svd <- irlba(A, d, maxit=10000, tol=1e-10)
  } else {
    A.svd <- svd(A, d)
  }
  
  # Compute the embedded coordinates
  Xhat <- as.matrix(A.svd$u) %*% diag(sqrt(A.svd$d), nrow=d, ncol=d)
  return(list(eval=A.svd$d, Xhat=Matrix(Xhat)))
}


# This function calculates latent positions for multiple graphs via COSIE model.
# It takes in a list of adjacency matrices 'A' and returns their latent positions.
# Parameters:
#   A:             List of adjacency matrices.
#   dSVD:          Dimension for joint singular value decomposition.
#   dASE:          Dimension for adjacency spectral embedding.
#   center:        If TRUE, mean centering is applied to adjacency matrices.
#   verbose:       If TRUE, additional messages are printed.
#   python:        Flag (unused in the current function but could be for potential future extensions).
#   latent.form:   Form of the latent representation. Can take values 1, 2, or 3.
jrdpg.latent <- function(A, dSVD=NA, dASE, center=FALSE, verbose=FALSE, python=FALSE, latent.form=2) {
  
  m <- length(A)
  A <- lapply(A, function(x) (x) )
  
  # If 'center' is TRUE, perform mean centering on the adjacency matrices
  if (center) {
    if (verbose) cat("do centering...\n")
    Abar <- Reduce('+', A) / length(A)
    A <- lapply(A, function(x) as.matrix(x-Abar))
  }   
  
  # Compute the joint embedding using the 'mase' function
  jrdpg <- mase(A, d = dSVD, d_vec=dASE)
  d2 <- dim(jrdpg$R[[1]])[1]
  
  # Obtain the latent representation based on the 'latent.form' value
  if(latent.form==3){
    B.svd <- lapply(jrdpg$R, function (x) svd(x) )
    Xhat <- lapply(B.svd, function(x) (jrdpg$V) %*% as.matrix(x$u) %*% diag(sqrt(x$d), nrow=d2, ncol=d2) %*% diag(sign(x$d)) %*% t(as.matrix(x$u) )  ) 
  }else if (latent.form==2){
    Xhat <- lapply(jrdpg$R, function(x) as.matrix((jrdpg$V) %*% x) )
  }else{
    P <- lapply(jrdpg$R, function (x) as.matrix((jrdpg$V) %*% x %*% t(jrdpg$V)) )
    Xhat <- lapply(P, function(x) as.matrix( full.ase(x, d2, diagaug=FALSE, approx=FALSE )$Xhat))  
  }
  
  return(Xhat)
}

diagAug <- function(A)
{
  diag(A) <- rowSums(A) / (nrow(A)-1)
  A
}

# Alist: a list of n x n adjacency matrices or igraph objects
buildOmni <- function(Alist, diagaug=FALSE, center=FALSE, verbose=FALSE)
{
  require(igraph)
  require(Matrix)
  
  if (class(Alist[[1]]) == "igraph") {
    Alist <- lapply(Alist, get.adjacency)
  }
  
  if (center) {
    if (verbose) cat("do centering...\n")
    Abar <- Reduce('+', Alist) / length(Alist)
    Alist <- lapply(Alist, function(x) as.matrix(x-Abar))
  }   
  
  if (diagaug) {
    Alist <- lapply(Alist, diagAug)
  }   
  
  m <- length(Alist)
  nvec <- sapply(Alist, nrow)
  nsum <- c(0, cumsum(nvec))
  
  omni <- as.matrix(bdiag(Alist))
  for (i in 1:(m-1)) {
    irng <- (nsum[i]+1):nsum[i+1] 
    for (j in (i+1):m) {
      jrng <- (nsum[j]+1):nsum[j+1]
      omni[irng,jrng] <- (as.matrix(Alist[[i]]) + as.matrix(Alist[[j]])) / 2
      omni[jrng,irng] <- omni[irng,jrng]
    }
  }
  return(omni)
}

# Sample undirected non-weighted adjacency graph G~RDPG(X)
rdpg.sample <- function(X) {
  P <- X %*% t(X)
  n <-  nrow(P)
  U <- matrix(0, nrow = n, ncol = n)
  U[col(U) > row(U)] <- runif(n*(n-1)/2)
  U <- (U + t(U))
  diag(U) <- runif(n)
  A <- (U < P) + 0 ;
  diag(A) <- 0
  return(graph.adjacency(A,"undirected"))
}

## Given X a nxd matrix
## and Y a mxd matrix
## This function compute the nxm matrix whose ij element is 
## the Euclidean distance between the i-th row of X 
## and the j-th row of Y.
pdistXY <- function(X,Y) {
  D <- rowSums((as.matrix(X) - as.matrix(Y))^2)
  return(sqrt(D))
}

# Generate a time series graph with a specific perturbation at time points 6 and 7 as Scenario 1
# @param n Number of nodes.
# @param nperturb Number of nodes to perturb.
# @param cperturb Constant for perturbation.
# @param rmin Minimum range value for uniform distribution.
# @param rmax Maximum range value for uniform distribution.
# @param tmax Total number of time points.
# @return A list of graphs.
genTSG <- function(n, nperturb, cperturb=NULL, rmin, rmax, tmax)
{    
  X1 <- runif(n, rmin, rmax)
  pert1 <- ifelse(is.null(cperturb), min(X1), cperturb) 
  pert2 <- ifelse(is.null(cperturb), 1-max(X1), cperturb) 
  Xlist <- rep(list(X1), tmax)
  Xlist[[6]] <- Xlist[[1]] + c(rep(c(pert2, -pert1), each=nperturb/2), rep(0, n-nperturb))
  Xlist[[7]] <- Xlist[[1]] + c(rep(c(-pert1, pert2), each=nperturb/2), rep(0, n-nperturb))
  sapply(Xlist, range)
  
  glist <- lapply(Xlist, function(x) rdpg.sample(matrix(x,ncol=1)))
  return(glist)
}

# Generate a time series graph with a specific perturbation at time points 6 and 7 as Scenario 2 and 3.
# @param n Number of nodes.
# @param nperturb Number of nodes to perturb.
# @param cperturb Constant for perturbation.
# @param rmin Minimum range value for uniform distribution.
# @param rmax Maximum range value for uniform distribution.
# @param tmax Total number of time points.
# @param d_true Dimension for true latent positions
# @param alpha Parameter in Dirichlet distribution(1_{d_true}*alpha)
# @return A list of graphs.
genTSGdd <- function(n, nperturb, cperturb=NULL, rmin, rmax, tmax, d_true, alpha)
{    
  k <- d_true
  Z <- matrix(0, n, d_true)
  initialZ <- matrix(0, n, d_true)
  tempinitialZ <- matrix(0, n, d_true)
  if(alpha<.1^4){
    for (i in 1:k) {
      Z[((n/k)*(i-1)+1):((n/k)*(i)), i] = 1
    }
  }else{
    tempinitialZ <- na.omit(rdirichlet(n, alpha * rep(1, d_true) ))
    while(dim(tempinitialZ)[1]!=n){
      i <- n - dim(tempinitialZ)[1]
      tempinitialZ <- rbind(tempinitialZ, na.omit(rdirichlet(i, alpha * rep(1, d_true) )))
    }
    initialZ <- tempinitialZ
    Z <- initialZ
  }
  
  
  q <- .5 * runif(1, rmin, rmax)
  B <- (  .5 * runif(1, rmin, rmax) - q + .5 )* diag(rep(1,d_true) ) +  q* rep(1, d_true) %*% t(rep(1, d_true))
  
  temp <-  svd(Z%*%B %*% t(Z))
  X1 <- temp$u[,1:d_true] %*%  diag( sqrt(temp$d[1:d_true]) ) 
  order_knn <- order(apply( X1, MARGIN = 1, function(x) sum((x-X1[1,])^2) ))
  temp <- X1
  temp[(1:(n/k)),] <- X1[order_knn[(1:(n/k))],]
  temp[-(1:(n/k)),] <- X1[-order_knn[(1:(n/k))],]
  X1 <- temp
  Xlist <- rep(list(X1), tmax)
  Z_gamma <- matrix(0, n, k)
  Z_pert <- t(matrix(rep(.6*rdirichlet(1, rep(1,  d_true) )+.2, times = n/k), d_true,n/k) )
  Z_gamma <- rbind(Z_pert, matrix(0, n-n/k, d_true)) 
  
  Xlist[[6]] <- Xlist[[1]] + cperturb * Z_gamma
  Xlist[[7]] <- Xlist[[1]] - cperturb * Z_gamma
  sapply(Xlist, range)
  glist <- lapply(Xlist, function(x) rdpg.sample(matrix(x,ncol= d_true)))
  return(glist)
}

# Generate a time series graph with a specific perturbation at time points 6 and 7 similar to Scenario 2 and 3.
# Fix perturbation to avoid confounding effects of randomness in perturbation, but results are similar for random perturbations.
# @param n Number of nodes.
# @param nperturb Number of nodes to perturb.
# @param cperturb Constant for perturbation.
# @param rmin Minimum range value for uniform distribution.
# @param rmax Maximum range value for uniform distribution.
# @param tmax Total number of time points.
# @param d_true Dimension for true latent positions
# @param alpha Parameter in Dirichlet distribution(1_{d_true}*alpha)
# @return A list of graphs.
genTSGonepar <- function(n, nperturb, cperturb=NULL, rmin, rmax, tmax, d_true, alpha){
  
  require(gtools)
  k <- d_true
  Z <- matrix(0, n, d_true)
  initialZ <- matrix(0, n, d_true)
  tempinitialZ <- matrix(0, n, d_true)
  if(alpha<.1^4){
    for (i in 1:k) {
      Z[((n/k)*(i-1)+1):((n/k)*(i)), i] = 1
    }
  }else{
    tempinitialZ <- na.omit(rdirichlet(n, alpha * rep(1, d_true) ))
    while(dim(tempinitialZ)[1]!=n){
      i <- n - dim(tempinitialZ)[1]
      tempinitialZ <- rbind(tempinitialZ, na.omit(rdirichlet(i, alpha * rep(1, d_true) )))
    }
    initialZ <- tempinitialZ
    Z <- initialZ
  }
  
  q <- .3
  B <- (  .3  - q + .5 )* diag(rep(1,d_true) ) +  q* rep(1, d_true) %*% t(rep(1, d_true))
  
  temp <-  svd(Z%*%B %*% t(Z))
  X1 <- temp$u[,1:d_true] %*%  diag( sqrt(temp$d[1:d_true]) )
  order_knn <- order(apply( X1, MARGIN = 1, function(x) sum((x-X1[1,])^2) ))
  temp <- X1
  temp[(1:(n/k)),] <- X1[order_knn[(1:(n/k))],]
  temp[-(1:(n/k)),] <- X1[-order_knn[(1:(n/k))],]
  X1 <- temp
  Xlist <- rep(list(X1), tmax)
  Z_gamma <- matrix(0, n, k)
  Z_pert <- t(matrix(rep(.3, times = n/k), d_true,n/k) )
  Z_gamma <- rbind(Z_pert, matrix(0, n-n/k, d_true))
  
  
  Xlist[[6]] <- Xlist[[1]] + cperturb * Z_gamma
  Xlist[[7]] <- Xlist[[1]] - cperturb * Z_gamma
  sapply(Xlist, range)
  glist <- lapply(Xlist, function(x) rdpg.sample(matrix(x,ncol= d_true)))
  return(glist)
}

# Variant of genTSG, but perturbs at time points 16 and 17.
conchartgenTSG <- function(n, nperturb, cperturb=NULL, rmin, rmax, tmax)
{    
  X1 <- runif(n, rmin, rmax)
  pert1 <- ifelse(is.null(cperturb), min(X1), cperturb) 
  pert2 <- ifelse(is.null(cperturb), 1-max(X1), cperturb) 
  Xlist <- rep(list(X1), tmax)
  Xlist[[10+6]] <- Xlist[[1]] + c(rep(c(pert2, -pert1), each=nperturb/2), rep(0, n-nperturb))
  Xlist[[10+7]] <- Xlist[[1]] + c(rep(c(-pert1, pert2), each=nperturb/2), rep(0, n-nperturb))
  sapply(Xlist, range)
  
  glist <- lapply(Xlist, function(x) rdpg.sample(matrix(x,ncol=1)))
  return(glist)
}

# Calculate the true Frobenius norm and pairwise distances between consecutive time points in the time series graph.
# This function perturbs at time points 16 and 17.
truenorm.aux <- function(n, nperturb, cperturb=NULL, rmin, rmax, tmax){
  
  X1 <- runif(n, rmin, rmax)
  pert1 <- ifelse(is.null(cperturb), min(X1), cperturb) 
  pert2 <- ifelse(is.null(cperturb), 1-max(X1), cperturb) 
  Xlist <- rep(list(X1), tmax)
  Xlist[[10+6]] <- Xlist[[1]] + c(rep(c(pert2, -pert1), each=nperturb/2), rep(0, n-nperturb))
  Xlist[[10+7]] <- Xlist[[1]] + c(rep(c(-pert1, pert2), each=nperturb/2), rep(0, n-nperturb))
  
  norm <- sapply(1:(tmax-1), function(x) norm(Xlist[[x]]-Xlist[[x+1]], "2"))
  pdist <- sapply(1:(tmax-1), function(x) pdistXY(Xlist[[x]], Xlist[[x+1]]))
  return(list(tnorm=norm, pdist=pdist))
  
}

# Variant of truenorm.aux, but perturbs at time points 6 and 7.
truenorm <- function(n, nperturb, cperturb=NULL, rmin, rmax, tmax){
  
  X1 <- runif(n, rmin, rmax)
  pert1 <- ifelse(is.null(cperturb), min(X1), cperturb) 
  pert2 <- ifelse(is.null(cperturb), 1-max(X1), cperturb) 
  Xlist <- rep(list(X1), tmax)
  Xlist[[6]] <- Xlist[[1]] + c(rep(c(pert2, -pert1), each=nperturb/2), rep(0, n-nperturb))
  Xlist[[7]] <- Xlist[[1]] + c(rep(c(-pert1, pert2), each=nperturb/2), rep(0, n-nperturb))
  #X1 is a vector so 2-norm is fine here
  norm <- sapply(1:(tmax-1), function(x) norm(Xlist[[x]]-Xlist[[x+1]], "2"))
  pdist <- sapply(1:(tmax-1), function(x) pdistXY(Xlist[[x]], Xlist[[x+1]]))
  return(list(tnorm=norm, pdist=pdist))
  
}

# This function computes norms and pairwise distances for embeddings obtained from the omnibus matrix.
# Arguments:
#   glist: A list of igraph objects.
#   nomni: Number of omnibus embeddings to use.
#   dmax: The maximum embedding dimension.
#   center: Whether to center the adjacency matrices.
#   approx: Whether to use approximation methods.
doOmni <- function(glist, nomni=2, dmax=NULL, center=FALSE, approx=TRUE)
{
  n <- vcount(glist[[1]])

  tmax <- length(glist)
  if (nomni == 2) {
    omni <- lapply(1:(tmax-1), function(x) buildOmni(glist[x:(x+1)], center=center, diagaug=TRUE))
    ase <- lapply(omni, function(x) full.ase(x[], dmax, diagaug=FALSE, approx=approx))
    norm <- sapply(ase, function(x) norm(x$Xhat[1:n,] - x$Xhat[-(1:n),], "F"))#"2"))
    pdist <- sapply(ase, function(x) pdistXY(x$Xhat[1:n,],x$Xhat[-(1:n),])) #abs(x$Xhat[1:n,] - x$Xhat[-(1:n),]))
  } else {
    omni <- buildOmni(glist, center=center, diagaug=TRUE)
    ase <- full.ase(omni, dmax, diagaug=FALSE, approx=approx)
    chunk <- running(1:nrow(omni), width=n, by=n, fun=function(x) x)
    norm <- sapply(1:(tmax-1), function(x) norm(ase$Xhat[chunk[,x],] - ase$Xhat[chunk[,x+1],], "F"))
    pdist <- sapply(1:(tmax-1), function(x) pdistXY(ase$Xhat[chunk[,x],], ase$Xhat[chunk[,x+1],]))
  }
  #    sapply(omni, function(x) is.connected(graph_from_adjacency_matrix(x, "undirected")))
  return(list(tnorm=norm, pdist=pdist))
}

# This function computes norms and pairwise distances for the Multiple Adjacency Spectral Embedding (MASE) of a list of graphs.
# Arguments:
#   glist: A list of igraph objects.
#   latpos.list: A list of latent position matrices.
#   nmase: The s of MASE to perform.
#   dSVD, dASE: Dimension parameters. dSVD is for joint SVD, dASE is for individual ASE.
#   center: Whether to center the adjacency matrices.
#   approx: Whether to use approximation method.
#   python: A flag (unused).
#   latent.form: Determines the method to construct latent positions.
doMase <- function(glist, nmase=2, dSVD=NA, dASE=NULL, center=FALSE, approx=TRUE, python=TRUE,latent.form=3)
{
  n <- vcount(glist[[1]])
  
  tmax <- length(glist)
  #Use Frobenius norm instead
  if (nmase == 2) {
    adj <- lapply(1:(tmax-1), function(x) lapply(glist[x:(x+1)], get.adjacency))
    Xhat <- lapply(adj, function(x) jrdpg.latent(as.matrix(x), dSVD, dASE, center=center,python = python, latent.form = latent.form))
    norm <- sapply(1:(tmax-1), function(x) norm((Xhat[[x]])[[1]] - (Xhat[[x]])[[2]], "F"))
    pdist <- sapply(1:(tmax-1), function(x) pdistXY((Xhat[[x]])[[1]] , (Xhat[[x]])[[2]] ))
  } else {
    adj <- lapply(lapply(glist, get.adjacency), as.matrix )
    Xhat <- jrdpg.latent(adj,dSVD, dASE, center=center,python = python, latent.form = latent.form)
    
    norm <- sapply(1:(tmax-1), function(x) norm(Xhat[[x]]-Xhat[[x+1]], "F"))
    pdist <- sapply(1:(tmax-1), function(x) pdistXY(Xhat[[x]], Xhat[[x+1]] ) )
     
  }
  
  return(list(tnorm=norm, pdist=pdist))
}

# This function computes norms and pairwise distances for the scan statistics of a list of graphs.
doScan <- function(glist)
{
  n <- vcount(glist[[1]])
  tmax <- length(glist)
  sstat <- lapply(1:(tmax-1), function(x) lapply(glist[(x):(x+1)], local.scan))
  
  norm <- sapply(1:(tmax-1), function(x) norm(max(sstat[[x]][[1]]) - max(sstat[[x]][[2]]), "2"))
  pdist <- sapply(1:(tmax-1), function(x) pdistXY(as.matrix(sstat[[x]])[[1]] , as.matrix(sstat[[x]])[[2]]))
  return(list(tnorm=norm, pdist=pdist))
}

# This function computes Frob norms and pairwise Euclidean distances for a list of graphs.
doAdj <- function(glist)
{
  n <- vcount(glist[[1]])
  tmax <- length(glist)
  adjstat <- lapply(1:(tmax-1), function(x) lapply(glist[(x):(x+1)], function(y) get.adjacency(y) ))
  
  norm <- sapply(1:(tmax-1), function(x) norm(adjstat[[x]][[1]] - adjstat[[x]][[2]], "F"))
  pdist <- sapply(1:(tmax-1), function(x) pdistXY(as.matrix(adjstat[[x]])[[1]] , as.matrix(adjstat[[x]])[[2]]))
  return(list(tnorm=norm, pdist=pdist))
}

#' Calculate local scan statistics for each vertex in TSG with weighted argument included
#'
#' @param g igraph object or matrix.
#' @param gp igraph object or matrix. It's typically at t=t-1, for "them" statistics.
#' @param k radius of the neighborhood.
#' @param mode one of "in", "out", "total".
#' @param FUN scan function. It's one of "ecount" or "vcount".
#' @param weighted whether or not the edges are weighted in locality stat calculation.
#' @return a vector of scan statistics, one for each vertex in the graph.
#' @export 
local.scan <- function(g,gp=NULL,k=1,mode="out",weighted=FALSE) {
  if (k < 0) stop("Error: k should be a non-negative integer!\n")
  require(igraph)
  if (is.matrix(g) | is.matrix(gp)) {
    gmode <- ifelse((mode == "out" | mode == "in"), "directed", "undirected")
    g <- simplify(graph.adjacency(g, mode = gmode))
    if (!is.null(gp)) 
      gp <- simplify(graph.adjacency(gp, mode = gmode))
  }
  
  n <- vcount(g)
  local.stat <- numeric(n)
  
  find.common.edge0<-function(nbrhd,g,gp,x,mode){
    nbrhd=nbrhd[nbrhd!=x]
    l=length(nbrhd)
    vps=numeric(4*l)
    vps[c(seq(1,(l*2),2),seq(2*l+2,4*l,2))]=x
    vps[c(seq(2,(l*2),2),seq(2*l+1,4*l,2))]=nbrhd
    if(mode=="all") {
      e.list=get.edge.ids(gp,vps,directed=T)       
    } else if(mode=="out") {
      e.list=get.edge.ids(gp,vps[1:(l*2)],directed=T)        
    } else if(mode=="in") {
      e.list=get.edge.ids(gp,vps[(2*l+1):(l*4)],directed=T)         
    }
    e.list=e.list[e.list!=0]
    return(unique(e.list))
  }
  
  # unweighted graph
  if(!isTRUE(weighted)) {
    if (is.null(gp)) {
      if (k == 0) local.stat=degree(g,mode=mode)
      else local.stat=sapply(graph.neighborhood(g,k,1:n,mode),ecount)
    }
    else {
      if (k==0) 
        local.stat= sapply(1:n,function(x) {
          com.nbrh=intersect(unlist(neighborhood(g,k+1,x,mode)), unlist(neighborhood(gp,k+1,x,mode)))
          if(length(com.nbrh)==1) return(0)
          else {
            e.list=find.common.edge0(com.nbrh,g,gp,x,mode)
            length(e.list)
          }
        })
      else local.stat = sapply(neighborhood(g,k,1:n,mode),function(x) {
        #GDC modification
        #ecount(induced.subgraph(gp,x))
        ecount(induced.subgraph(gp,V(gp)[as_ids(x)]))
      })
    }
  }
  
  # weighted graph
  else{
    if (is.null(gp)) {
      if (k == 0) local.stat=sapply(1:n, function(x) {
        incident.edges = incident (g, x, mode=mode )
        sum(E(g)[incident.edges]$weight)
      })
      else local.stat= sapply(graph.neighborhood(g, k, V(g),mode=mode), function(x) {
        sum(E(x)$weight)
      })
    }
    else {
      if (k==0) 
        local.stat=sapply(1:n,function(x) {
          com.nbrh=intersect(unlist(neighborhood(g,k+1,x,mode)), unlist(neighborhood(gp,k+1,x,mode)))
          if(length(com.nbrh)==1) return(0)
          else {
            e.list=find.common.edge0(com.nbrh,g,gp,x,mode=mode)
            sum(E(gp)[e.list]$weight)
          }
        })
      else local.stat = sapply(neighborhood(g,k,1:n,mode),function(x) {
        #GDC modification
        #sum(E(induced.subgraph(gp,x))$weight)
        sum(E(induced.subgraph(gp,V(gp)[as_ids(x)]))$weight)
      })
    }
  }
  return(local.stat)
}



plot.qcc.vertex <- function(x, m2,add.stats = TRUE, chart.all = TRUE, 
                            label.limits = "UCL",#c("LCL ", "UCL"),
                            title, xlab, ylab, ylim, axes.las = 0,
                            digits =  getOption("digits"),
                            restore.par = TRUE, ...) 
{
  l.length <- length(x)
  df <- c() #c(2,6)
  for (i in c(2,6)) {
    object <- x[[i]]  # Argh.  Really want to use 'object' anyway
    if ((missing(object)) | (!inherits(object, "qcc")))
      stop("an object of class `qcc' is required")
    
    # collect info from object
    type <- object$type
    std.dev <- object$std.dev
    data.name <- object$data.name
    center <- object$center
    stats <- object$statistics
    limits <- object$limits 
    lcl <- limits[,1]
    ucl <- limits[,2]
    newstats <- object$newstats
    newdata.name <- object$newdata.name
    violations <- object$violations
    if(chart.all) 
    { statistics <- c(stats, newstats)
    indices <- 1:length(statistics) }
    
    library(ggplot2)
    tmax <- l.length+1
    tvec <- 1:tmax
    minx <- names(running(tvec, width=2))
    vioinx <- rep(0,length(indices))
    runinx <- rep(0,length(indices))
    vioinx[violations$beyond.limits] <- 1
    runinx[violations$violating.runs] <- 1
    idf <- data.frame(vertex=indices, y= statistics, timepoints=factor(minx[i], levels = minx),lcl=lcl,ucl=ucl,center=center,lim=vioinx,run=runinx)
    df <- rbind(df,idf)
    
  }
  
  p <- ggplot(df,aes(x=vertex, y=y))+#geom_hline(aes(yintercept=lcl), linetype="dashed")+
    geom_hline(aes(yintercept=ucl), linetype="dashed")+
    geom_vline(aes(xintercept=20), linetype="dashed",color="red")+
    geom_hline(aes(yintercept=center), linetype="solid")+facet_wrap(~timepoints,nrow = 1,strip.position = "left", scales = "fixed") +
    geom_point(data = df%>% filter(lim!=1), alpha=1, color="grey70",shape=20,size=2)+ylab(TeX("$y_{i}^{(t)}$"))+
    #geom_point(data = df %>% filter(run==1), alpha=1, color="yellow")+
    geom_point(data = df %>% filter(lim==1), alpha=1, color="red",shape=17,size=2)+
    #geom_line() + 
    labs(x="vertex")+theme_bw()+
    theme(axis.title.x = element_text(size = 14), strip.text.x = element_text(size = 14),legend.position = "none",plot.title = element_text(hjust = 0.5,size=14, face='bold'),plot.subtitle = element_text(hjust = 0.5,size=14, face='bold'),text = element_text(size=14, face='bold' )) +ggtitle(title)
  
  return(p)
  invisible()
}



plot.qcc.true <- function(x, m2,add.stats = TRUE, chart.all = TRUE, 
                          label.limits = "UCL",#c("LCL ", "UCL"), ylim, axes.las = 0,
                          digits =  getOption("digits"),
                          restore.par = TRUE, ...) 
{
  object <- x  # Argh.  Really want to use 'object' anyway
  if ((missing(object)) | (!inherits(object, "qcc")))
    stop("an object of class `qcc' is required")
  # collect info from object
  type <- object$type
  std.dev <- object$std.dev
  data.name <- object$data.name
  center <- object$center
  stats <- object$statistics
  limits <- object$limits 
  lcl <- limits[,1]
  ucl <- limits[,2]
  newstats <- object$newstats
  newdata.name <- object$newdata.name
  violations <- object$violations
  if(chart.all) 
  { statistics <- c(stats, newstats)
  indices <- 1:length(statistics) }
  
  library(ggplot2)
  tmax <- length(indices)+1
  tvec <- 1:tmax
  minx <- names(running(tvec, width=2))
  vioinx <- rep(0,length(indices))
  runinx <- rep(0,length(indices))
  vioinx[violations$beyond.limits] <- 1
  runinx[violations$violating.runs] <- 1
  df <- data.frame(time=indices, y= statistics,lcl=lcl,ucl=ucl,center=center,lim=vioinx,run=runinx)
  
  p <- ggplot(df,aes(x=time, y=y))+
    geom_step(aes(x=time, y=ucl), linetype="dashed")+
    geom_step(aes(x=time, y=center), linetype="solid") +
    geom_point(data = df %>% filter(lim!=1), alpha=1, color="grey70", size=2,shape=20)+
    geom_line(aes(x=time, y=y),color="grey70")+
    # geom_point(data = df %>% filter(run==1), alpha=1, color="yellow", size=2)+
    geom_point(data = df %>% filter(lim==1), alpha=1, color="red", size=2,shape=17) + 
    scale_x_discrete(name ="time points", 
                     limits=c(m2))+theme_bw()+
    annotate("text", label = "UCL", 
             x =  tmax-.7, y = rev(ucl)[1] )+ylab(TeX("$||X^{(t-1)}-X^{(t)}||$"))+
    annotate("text", label = "CL", 
             x =  tmax-.7, y = rev(center)[1] )+
    theme(legend.position = "none",plot.title = element_text(hjust = 0.5,size=25, face='bold'),plot.subtitle = element_text(hjust = 0.5),text = element_text(size=rel(4.5), face='bold' ) ) 
  
  return(p)
}



plot.qcc <- function(x, m2,add.stats = TRUE, chart.all = TRUE, 
                     label.limits = "UCL",#c("LCL ", "UCL"),
                     title, xlab, ylab, ylim, axes.las = 0,
                     digits =  getOption("digits"),
                     restore.par = TRUE, ...) 
{
  object <- x  # Argh.  Really want to use 'object' anyway
  if ((missing(object)) | (!inherits(object, "qcc")))
    stop("an object of class `qcc' is required")
  # collect info from object
  type <- object$type
  std.dev <- object$std.dev
  data.name <- object$data.name
  center <- object$center
  stats <- object$statistics
  limits <- object$limits 
  lcl <- limits[,1]
  ucl <- limits[,2]
  newstats <- object$newstats
  newdata.name <- object$newdata.name
  violations <- object$violations
  if(chart.all) 
  { statistics <- c(stats, newstats)
  indices <- 1:length(statistics) }
  
  library(ggplot2)
  tmax <- length(indices)+1
  tvec <- 1:tmax
  minx <- names(running(tvec, width=2))
  vioinx <- rep(0,length(indices))
  runinx <- rep(0,length(indices))
  vioinx[violations$beyond.limits] <- 1
  runinx[violations$violating.runs] <- 1
  df <- data.frame(time=indices, y= statistics,lcl=lcl,ucl=ucl,center=center,lim=vioinx,run=runinx)
  
  p <- ggplot(df,aes(x=time, y=y))+
    geom_step(aes(x=time, y=ucl), linetype="dashed")+
    geom_step(aes(x=time, y=center), linetype="solid") +
    geom_point(data = df%>% filter(lim!=1), alpha=1, color="grey70", size=1,shape=20)+
    geom_line(aes(x=time, y=y),color="grey70")+
    # geom_point(data = df %>% filter(run==1), alpha=1, color="yellow", size=2)+
    geom_point(data = df %>% filter(lim==1), alpha=1, shape=17,color="red", size=2) + 
    scale_x_discrete(name ="time points", 
                     limits=c(m2))+theme_bw()+
    annotate("text", label = "UCL", 
             x =  tmax-.7, y = rev(ucl)[1] )+ylab(TeX("$y^{(t)}$"))+
    annotate("text", label = "CL", 
             x =  tmax-.7, y = rev(center)[1] )+
    theme(legend.position = "none",plot.title = element_text(hjust = 0.5,size=25, face='bold'),plot.subtitle = element_text(hjust = 0.5),text = element_text(size=rel(4.5), face='bold' ) ) +ggtitle(title)
  
  return(p)
  invisible()
}


#duplicated from /Users/guodongchen/Documents/R/getElbows.R on sep 24, 2022
#'
#' Given a decreasingly sorted vector, return the given number of elbows
#'
#' @param dat a input vector (e.g. a vector of standard deviations), or a input feature matrix
#' @param n the number of returned elbows
#' @param threshold either FALSE or a number. If threshold is a number, then all the elements in d that are not larger than the threshold will be ignored.
#' @param plot logical. When T, it depicts a scree plot with highlighted elbows
#' @param main a string of the plot title
#'
#' @return a vector of length \eqn{n}
#'
#' @references Zhu, Mu and Ghodsi, Ali (2006), Automatic dimensionality selection from
#' the scree plot via the use of profile likelihood, Computational
#' Statistics & Data Analysis, Volume 51 Issue 2, pp 918-930, November, 2006.
#'
#' @author Youngser Park <youngser@jhu.edu>
#' @export
#' @importFrom graphics plot points
#' @importFrom stats sd dnorm
#'
getElbows <- function(dat, n = 3, threshold = FALSE, plot = TRUE, main="", ...) {
  ## Given a decreasingly sorted vector, return the given number of elbows
  ##
  ## Args:
  ##   dat: a input vector (e.g. a vector of standard deviations), or a input feature matrix.
  ##   n: the number of returned elbows.
  ##   threshold: either FALSE or a number. If threshold is a number, then all
  ##   the elements in d that are not larger than the threshold will be ignored.
  ##   plot: logical. When T, it depicts a scree plot with highlighted elbows.
  ##
  ## Return:
  ##   q: a vector of length n.
  ##
  ## Reference:
  ##   Zhu, Mu and Ghodsi, Ali (2006), "Automatic dimensionality selection from
  ##   the scree plot via the use of profile likelihood", Computational
  ##   Statistics & Data Analysis, Volume 51 Issue 2, pp 918-930, November, 2006.
  
  #  if (is.unsorted(-d))
  
  if (is.matrix(dat)) {
    d <- sort(apply(dat,2,sd), decreasing=TRUE)
  } else {
    d <- sort(dat,decreasing=TRUE)
  }
  
  if (!is.logical(threshold))
    d <- d[d > threshold]
  
  p <- length(d)
  if (p == 0)
    stop(paste("d must have elements that are larger than the threshold ",
               threshold), "!", sep="")
  
  lq <- rep(0.0, p)                     # log likelihood, function of q
  for (q in 1:p) {
    mu1 <- mean(d[1:q])
    mu2 <- mean(d[-(1:q)])              # = NaN when q = p
    sigma2 <- (sum((d[1:q] - mu1)^2) + sum((d[-(1:q)] - mu2)^2)) /
      (p - 1 - (q < p))
    lq[q] <- sum( dnorm(  d[1:q ], mu1, sqrt(sigma2), log=TRUE) ) +
      sum( dnorm(d[-(1:q)], mu2, sqrt(sigma2), log=TRUE) )
  }
  
  q <- which.max(lq)
  # print(n)
  # print(q)
  # print(p)
  if (n > 1 && q < (p-1)) {
    q <- c(q, q + getElbows(d[(q+1):p], n-1, plot=FALSE))
  }
  
  if (plot==TRUE) {
    if (is.matrix(dat)) {
      sdv <- d # apply(dat,2,sd)
      plot(sdv,type="b",xlab="dim",ylab="stdev",main=main,...)
      points(q,sdv[q],col=2,pch=19)
    } else {
      plot(dat, type="b",main=main,...)
      points(q,dat[q],col=2,pch=19)
    }
  }
  
  return(q)
}

## Asssembed code supporting the text "Bayesian classification, anomaly detection, and survival analysis using network inputs with application to the microbiome"
# see https://github.com/KolaczykResearch/GP-Networks/tree/master
## https://projecteuclid.org/journals/annals-of-applied-statistics/volume-17/issue-1/Bayesian-classification-anomaly-detection-and-survival-analysis-using-network-inputs/10.1214/22-AOAS1623.short

## Citation

## Josephs, Nathaniel, Lizhen Lin, Steven Rosenberg, and Eric D. Kolaczyk. "Bayesian classification, anomaly detection, and survival analysis using network inputs with application to the microbiome." (2023). _The Annals of Applied Statistics_.

dist.frobenius <- function(G, X = NULL) {
  
  p.G <- ifelse(missing(G), 0, length(G))
  p.X <- ifelse(is.null(ncol(X)), 0, ncol(X))
  
  if (!missing(G)) {
    N <- length(G[[1]])
    D <- array(0, dim = c(N, N, p.G + p.X))
    n <- nrow(G[[1]][[1]])
    
    for (pp in 1:p.G) {
      V <- sapply(G[[pp]], as.vector)
      D.frobenius <- as.matrix(dist(t(V))^2) / (n*(n-1))
      D[, , pp] <- D.frobenius
    }
  }
  
  if (!is.null(X)) {
    N <- nrow(X); p <- ncol(X)
    
    if (missing(G)) {
      D <- array(0, dim = c(N, N, p.X))
    }
    
    for (pp in 1:p.X) {
      D[, , p.G + pp] <- as.matrix(dist(X[, pp])^2)
    }
  }
  
  return(D)
}


# Functions ----
kern.fun <- function(D, theta) {
  # computes squared-exponential kernel given distance matrix (array) D
  # NB: does not square distance
  
  N <- nrow(D)
  M <- ncol(D)
  p <- dim(D)[3]
  K0 <- matrix(0, N, M)
  sigma <- theta[1]
  
  for (pp in 1:p) {
    ell <- exp(theta[pp + 1])
    
    K0 <- K0 + exp(-ell * D[, , pp])
  }
  
  K <- sigma * K0  # sigma is overall signal variance
  
  return(K)
}

log1pe <- function (x) { # vectorized version: `x` can be a vector 
  l <- ifelse(x > 0, x, 0) # shift 
  x <- ifelse(x > 0, -x, x) # range reduction: `x = -abs(x)` 
  ifelse(x < log(.Machine$double.eps), l, l + log(1 + exp(x))) }

ilogit <- function (x) {
  if (x >= 0) {
    1 / (1 + exp(-x))  
  } else {
    z <- exp(x)
    z / (1 + z) 
  }
}

rmvnorm <- function(n = 1, mu, C){
  p <- length(mu)
  
  if (p == 1) {
    X <- rnorm(1, mu, C)
  } else{
    Z <- matrix(rnorm(p*n), p, n)
    
    X <- crossprod(C, Z)
    X <- sweep(X, 1, mu, FUN = `+`)
  }
  
  return(X)
}

chol.fun <- function(x, jit = 1e-6) {
  
  n <- nrow(x)
  L <- chol(x + diag(n) * jit^2)
  
  return(L)
}

slice <- function(theta, f, alpha, Y, D, a, b, sigma = 10) {
  
  log.lik <- function(u, f, C, g, theta, a, b, Y) {
    log(u) +                                                 # u
      -sum(log1pe(-Y*f)) +                                   # L(f)
      -.5 * (crossprod(backsolve(C, g, transpose = TRUE))) + # N(g; 0, Sigma + S)
      sum(-a * theta - b / exp(theta))                       # p(theta)
  }
  
  n <- length(f)
  p <- length(theta)
  Sigma <- kern.fun(D, c(alpha, theta))
  S <- diag(n) * alpha # auxillary noise
  S.inv <- diag(n) / alpha
  
  # 1. draw surrogate data
  S.chol <- S / sqrt(alpha)
  g <- rmvnorm(1, f, S.chol)
  
  # 2. compute implied latent variants
  C <- chol.fun(Sigma + S)
  R <- S - S %*% chol2inv(C) %*% S
  m <- R %*% S.inv %*% g
  
  L <- chol.fun(R)
  eta <- backsolve(L, f - m, transpose = TRUE)
  
  # 3. randomly center a bracket
  v <- runif(p, 0, sigma)
  theta.min <- theta - v
  theta.max <- theta.min + sigma
  
  # 4. draw u
  u <- runif(1)
  
  # 5. determine threshold
  log.y <- log.lik(u, f, C, g, theta, a, b, Y)
  
  # 6. draw proposal
  theta.s <- runif(p, theta.min, theta.max)
  
  # 7. compute function
  Sigma.s <- kern.fun(D, c(alpha, theta.s))
  C.s <- chol.fun(Sigma.s + S)
  R.s <- S - S %*% chol2inv(C.s) %*% S
  m.s <- R.s %*% S.inv %*% g
  
  L.s <- chol.fun(R.s)
  f.s <- crossprod(L.s, eta) + m.s
  
  log.yp <- log.lik(1, f.s, C.s, g, theta.s, a, b, Y)
  
  while (log.yp <= log.y) {
    for (pp in 1:p) {
      if (theta.s[pp] < theta[pp]) theta.min[pp] <- theta.s[pp] else theta.max[pp] <- theta.s[pp]
    }
    
    theta.s <- runif(p, theta.min, theta.max)
    
    Sigma.s <- kern.fun(D, c(alpha, theta.s))
    C.s <- chol.fun(Sigma.s + S)
    R.s <- S - S %*% chol2inv(C.s) %*% S
    m.s <- R.s %*% S.inv %*% g
    
    L.s <- chol.fun(R.s)
    f.s <- crossprod(L.s, eta) + m.s
    
    log.yp <- log.lik(1, f.s, C.s, g, theta.s, a, b, Y)
  }
  
  return(list("f" = f.s, "theta" = theta.s, "lp" = log.yp))
  
}

elliptical <- function(f0, C, Y) {
  # Samples f1 using elliptical slice sampler
  #
  # Args:
  #   f0       : current latent function
  #   C        : C'C = Sigma, f ~ GP(0, Sigma)
  #
  # Returns:
  #   New latent function f1
  
  log.L <- function(f) -sum(log1pe(-Y*f))
  
  n <- nrow(C)
  nu <- rmvnorm(1, rep(0, n), C)
  
  u <- runif(1)
  log.y <- log.L(f0) + log(u)
  
  theta <- runif(1, 0, 2*pi)
  theta.min <- theta - 2*pi; theta.max <- theta
  f1 <- f0 * cos(theta) + nu * sin(theta)
  
  log.L.f1 <- log.L(f1)
  
  while (log.L.f1 <= log.y) {
    if (theta < 0) theta.min <- theta else theta.max <- theta
    theta <- runif(1, theta.min, theta.max)
    f1 <- f0 * cos(theta) + nu * sin(theta)
    log.L.f1 <- log.L(f1)
  }
  
  return(list("f" = f1, "log.L" = log.L.f1))
  
}

gp.class <- function(dist, Y, split_index, a, b, ns = 1000, monitor = TRUE) {
  # Samples from p(f | X, Y)
  #
  # Args:
  #   dist        : distance array D s.t kernel K(D) with f ~ GP(0, K)
  #   Y           : Binary response vector
  #   split_index : Train/validate assignment
  #   a, b        : Hyperpriors for theta ~ InvGamma(a, b)
  #   ns          : number of samples
  #
  # Returns:
  #   Sample from latent posterior
  
  N.train <- sum(split_index == "train")
  train <- which(split_index == "train")
  D.train <- dist[train, train, , drop = FALSE]
  Y.train <- Y[train]
  
  # Initialize
  ## p(f | X, y)
  p <- length(a)
  theta <- matrix(nrow = p, ncol = ns)      # theta = (sigma, l_1, ..., l_p)
  theta[1, 1] <- 1                                      # sigma ~ IG(0, 0)
  theta[-1, 1] <- b[-1] / (a[-1] + 1)                   # ell ~ IG(a, b)
  f <- matrix(nrow = N.train, ncol = ns); f[, 1] <- 0   # prior mean
  lp.f <- numeric(ns); lp.f[1] <- 0
  
  for (t in 2:ns) {
    if (monitor && t %% (ns / 10) == 0) print(paste(t, "out of", ns))
    
    # jointly sample f and length scale(s)
    samp <- slice(theta = theta[-1, t - 1]
                  , f = f[, t - 1]
                  , alpha = theta[1, t - 1]
                  , Y = Y.train
                  , D = D.train
                  , a = a[-1], b = b[-1])
    
    f[, t] <- samp[["f"]]
    theta[-1, t] <- samp[["theta"]]
    lp.f[t] <- samp[["lp"]]
    
    # cheap updates of f
    K0 <- kern.fun(D.train, c(1, theta[-1, t])) # K0 = K / sigma
    K0.chol <- chol.fun(K0)
    K.chol <- K0.chol * sqrt(theta[1, t - 1])
    for (ii in 1:10) {
      slice <- elliptical(f[, t], K.chol, Y.train)
      f[, t] <- slice[["f"]]
    }
    lp.f[t] <- lp.f[t] + slice[["log.L"]]
    
    # sigma
    theta[1, t] <- rinvgamma(1, shape = a[1] + N.train / 2
                             , rate = crossprod(backsolve(K0.chol
                                                          , f[, t]
                                                          , transpose = TRUE)
                             ) / 2 + b[1])
    
    lp.f[t] <- lp.f[t] + dinvgamma(theta[1, t], shape = a[1] + N.train / 2
                                   , rate = crossprod(backsolve(K0.chol
                                                                , f[, t]
                                                                , transpose = TRUE)
                                   ) / 2 + b[1]
                                   , log = TRUE)
    K.chol <- K0.chol * sqrt(theta[1, t])
    
  }
  return(rbind(f, theta, t(as.matrix(lp.f))))
}

gp.class_pred <- function(gpc, D, split_index, p, burn = .2, thin = 10, avg = TRUE) {
  ns <- dim(gpc)[1]
  samp <- seq(burn*ns, ns, thin)
  
  N.train <- sum(split_index == "train")
  
  if (avg) {
    # posterior mean
    
    f.avg <- apply(gpc[samp, 1:N.train], 2, mean)
    
    # average hyperparameters
    theta.avg <- apply(gpc[samp, (N.train+1):(N.train+p+1)], 2, mean)
    K.avg <- kern.fun(D, replace(theta.avg, 1, 1)) # note that signal variance cancels
    fstar.avg <- crossprod(K.avg[split_index == "train"
                                 , split_index == "validate"]
                           , chol2inv(chol.fun(K.avg[split_index == "train"
                                                     , split_index == "train"]))) %*% f.avg
    y.avg <- sapply(fstar.avg, ilogit)
    
    return(y.avg)
    
  } else {
    # p(f* | X, y, x*)
    samp.fstar <- sapply(samp, function(ind) {
      theta.samp <- gpc[ind, (N.train+1):(N.train+p+1)]
      K <- kern.fun(D, replace(theta.samp, 1, 1))
      K.inv <- chol2inv(chol.fun(K[split_index == "train", split_index == "train"]))
      fstar <- crossprod(K[split_index == "train", split_index == "validate"]
                         , K.inv) %*% gpc[ind, 1:N.train]
    })
    # p(y* = 1 | f*, y, x*)
    if (is.matrix(samp.fstar)) {
      ystar <- apply(samp.fstar, 2, function(i) sapply(i, ilogit))
      y.map <- rowMeans(ystar)
      y.cred <- apply(ystar, 1, quantile, probs = c(.025, .975))
    } else {
      ystar <- sapply(samp.fstar, ilogit)
      y.map <- mean(ystar)
      y.cred <- quantile(ystar, probs = c(.025, .975))
    }
    
    return(c(y.map, y.cred))
  }
}

mcmc_array <- function (ns, nchains, params) {
  nparams <- length(params)
  array(dim = c(ns, nchains, nparams),
        dimnames = list(iterations = NULL,
                        chains = paste0("chain", 1:nchains),
                        parameters = params))
}

mcc <- function(y, y.pred) {
  
  y.pred <- factor(y.pred, levels = c(-1, 1)) # in case only one class predicted
  
  cm <- table(y, y.pred)
  
  TP <- cm["1", "1"]
  TN <- cm["-1", "-1"]
  FP <- cm["-1", "1"]
  FN <- cm["1", "-1"]
  
  num <- TP*TN - FP*FN
  den <- (TP+FP)*(TP+FN)*(TN+FP)*(TN+FN)
  den <- ifelse(den == 0, 1, den)
  
  return(num/sqrt(den)) 
}

roc <- function(labels, scores, plot = FALSE) {
  
  if (min(labels) == -1) labels <- (labels + 1) / 2
  
  Labels <- labels[order(scores, decreasing = TRUE)]
  FPR <- cumsum(!Labels)/sum(!Labels)
  TPR <- cumsum(Labels)/sum(Labels)
  
  if (plot) {
    plot(TPR ~ FPR, type = "o")
    abline(c(0, 1))
  }
  
  N <- length(labels)
  AUC <- sum((FPR[2:N] - FPR[1:N-1]) * TPR[2:N])
  
  return(AUC)
}

f1 <- function(y, y.pred) {
  
  y.pred <- factor(y.pred, levels = c(-1, 1)) # in case only one class predicted
  
  cm <- table(y, y.pred)
  
  TP <- cm["1", "1"]
  TN <- cm["-1", "-1"]
  FP <- cm["-1", "1"]
  FN <- cm["1", "-1"]
  
  precision <- TP / (TP + FP)
  recall <- TP / (TP + FN)
  
  F1 <- 2 * (precision * recall) / (precision + recall)
  
  return(F1) 
}

gp.occ <- function(gpc, D, burn = .2, thin = 10) {
  ns <- dim(gpc)[1]
  samp <- seq(burn*ns, ns, thin)
  
  # p(f* | X, y, x*)
  samp.fstar <- sapply(samp, function(ind) {
    theta.samp <- gpc[ind, (m.train+1):(m.train+p+1)]
    K <- kern.fun(D, replace(theta.samp, 1, 1))
    K.inv <- chol2inv(chol.fun(K[split_index == "train", split_index == "train"]))
    fstar <- crossprod(K[split_index == "train", split_index == "validate"]
                       , K.inv) %*% gpc[ind, 1:m.train]
  })
  
  # p(y* = 1 | f*, y, x*)
  ystar <- apply(samp.fstar, 2, function(i) sapply(i, ilogit))
  y.map <- rowMeans(ystar)
  
  # average hyperparameters
  f.avg <- apply(gpc[samp, 1:m.train], 2, mean)
  
  theta.avg <- apply(gpc[samp, (m.train+1):(m.train+p+1)], 2, mean)
  K.avg <- kern.fun(D, replace(theta.avg, 1, 1)) # note that signal variance cancels
  fstar.avg <- crossprod(K.avg[split_index == "train"
                               , split_index == "validate"]
                         , chol2inv(chol.fun(K.avg[split_index == "train"
                                                   , split_index == "train"]))) %*% f.avg
  y.avg <- sapply(fstar.avg, ilogit)
  
  # variance
  sigma_star <- apply(ystar, 1, sd)
  
  return(list(mu = fstar.avg, pi = y.avg, sigma = sigma_star, H = fstar.avg / sqrt(sigma_star)))
  
}


# Duplicated /Users/guodongchen/Documents/Graduate_project/carey/AD/code/masetsg/mase.R sep 24, 2022
#' 
#' Function to perform multiple adjacency spectral embedding
#' 
#' @param Adj_list a list of adjacency matrices with the same size n x n
#' @param d number of joint embedding dimensions. If NA, dimension is chosen automatically
#' @param d_vec vector of individual embedding dimensions for each graph. If NA, it is the same as d.
#' @param scaled.ASE logical. When true, the scaled adjacency spectral embedding is used.
#' @param diag.augment logical. When true, the diagonal augmentation method for ASE is performed.
#' @param elbow_graph number of elbow selected in Zhu & Ghodsi method for the scree plot of each individual graph eigenvalues.
#' @param elbow_mase number of elbow selected in Zhu & Ghodsi method for the scree plot of the singular values of the concatenated spectral embeddings of MASE.
#' @param show.scree.results when TRUE, the histogram of the estimated d for each graph, and the scree plot of the singular values of  the graph is shown if d is not specified.
#' @param par whether to run in parallel (TRUE/FALSE)
#' @param numpar number of clusters for parallel implmentation
#' 
#' @return A list containing a matrix V of size n x d, with the 
#' estimated invariant subspace, and a list R with the individual score parameters for each graph (size d x d each).
#' 
#' @references 
#' 
#' @author Jes\'us Arroyo <jesus.arroyo@jhu.edu>
mase <- function(Adj_list, d = NA, d_vec = NA,
                 scaled.ASE = TRUE, diag.augment = TRUE, 
                 elbow_graph = 1, elbow_mase = 2,
                 show.scree.results = FALSE,
                 par = FALSE, numpar = 12) {
  if(is.na(d_vec)) {
    d_vec = rep(d, length(Adj_list))
  }
  # modify by GDC
  if(length(d_vec)<length(Adj_list)) {
    d_vec = rbind(d_vec,rep(d_vec[1], length(Adj_list)-length(d_vec)))
  }
  # running in parallel
  if(par) {
    require(parallel)
    cl <- makeCluster(numpar)
    clusterEvalQ(cl, source("R/loadAll.R"))
    clusterExport(cl = cl, varlist = list("ase", "eig_embedding", "getElbows", "Adj_list",
                                          "elbow_graph", "d_vec", "diag.augment"), envir = environment())
    if(scaled.ASE) {
      latpos.list <- parLapply(cl = cl, 1:length(Adj_list), function(i) 
        ase(Adj_list[[i]], d = d_vec[i], diag.augment = diag.augment, elbow = elbow_graph))
    }else{
      latpos.list <- parLapply(cl = cl, 1:length(Adj_list), function(i) 
        eig_embedding(Adj_list[[i]], d = d_vec[i], diag.augment = diag.augment, elbow = elbow_graph))
    }
    stopCluster(cl)
  } else {
    if(scaled.ASE) {
      latpos.list <- lapply(1:length(Adj_list), function(i) 
        ase(Adj_list[[i]], d = d_vec[i], diag.augment = diag.augment, elbow = elbow_graph))
    }else{
      latpos.list <- lapply(1:length(Adj_list), function(i) 
        eig_embedding(Adj_list[[i]], d = d_vec[i], diag.augment = diag.augment, elbow = elbow_graph))
    }
  }
  V_all  <- Reduce(cbind, latpos.list)
  require(rARPACK)
  jointsvd <- svd(V_all)
  if(is.na(d)) {
    if(show.scree.results) {
      hist(sapply(latpos.list, ncol), main = "Estimated d for each graph")
    }
    d = getElbows(jointsvd$d, plot = show.scree.results)[elbow_mase]
  }
  V = jointsvd$u[, 1:d, drop = FALSE]
  R <- project_networks(Adj_list, V)
  return(list(V = V, sigma = jointsvd$d, R = R))
}



#' 
#' Function to perform graph adjacency spectral embedding (ASE)
#' 
#' @param A adjacency matrix
#' @param d number of joint embedding dimensions. If NA, dimension is chosen automatically
#' @param d.max maximum number of embedding dimensions to try when d is not provided. Default is sqrt(ncol(A)).
#' @param diag.augment whether to do diagonal augmentation (TRUE/FALSE)
#' @param elbow number of elbow selected in Zhu & Ghodsi method. Default is 1.
#' @return A matrix with n rows and d columns containing the estimated latent positions
#' 
#' @author Jes\'us Arroyo <jesus.arroyo@jhu.edu>
#' 
ase <- function(A, d = NA, d.max = sqrt(ncol(A)), diag.augment = TRUE, elbow = 1) {
  require(rARPACK)
  # Diagonal augmentation
  if(diag.augment & sum(abs(diag(A))) == 0) {
    deg = colSums(A)
    n = ncol(A)
    diag(A) = deg / (n-1)
  }
  if(is.na(d)) {
    eig <- eigs(as(A, "dgeMatrix"), d.max)
    vals <- sort(x =  abs(eig$values), decreasing = TRUE)
    d = getElbows(vals, plot = F)[elbow]
    selected.eigs <- which(abs(eig$values) >= vals[d])
    V <- eig$vectors[,selected.eigs, drop = F]
    D <- diag(sqrt(abs(eig$values[selected.eigs])), nrow = d)
    X <- V %*% D
    return(X)
  } else {
    eig <- eigs(as(A, "dgeMatrix"), k = d)
    X <- eig$vectors %*% diag(sqrt(abs(eig$values)), nrow = d)
    return(X) 
  }
}

#' 
#' Function to compute the graph unscaled adjacency spectral embedding (top eigenvectors)
#' 
#' @param A adjacency matrix
#' @param d number of joint embedding dimensions. If NULL, dimension is chosen automatically
#' @param d.max maximum number of embedding dimensions to choose (if d is given, this number is ignored)
#' @param diag.augment whether to do diagonal augmentation (TRUE/FALSE)
#' @param elbow number of elbow selected in Zhu & Ghodsi method. Default is 1.
#' 
#' @return A list containing a matrix with n rows and d columns representing the estimated latent positions, and the estimated
#' indefinite orthogonal projection matrix
#' 
#' @author Jes\'us Arroyo <jesus.arroyo@jhu.edu>
eig_embedding <- function(A, d = NA, d.max = ncol(A), diag.augment = FALSE, elbow = 1) {
  require(rARPACK)
  n <- ncol(A)
  if(diag.augment & sum(abs(diag(A))) == 0) {
    deg = colSums(A)
    n = ncol(A)
    diag(A) = deg / (n-1)
  }
  if(is.na(d)) {
    eig <- eigs(as(A, "dgeMatrix"), d.max)
    vals <- sort(x =  abs(eig$values), decreasing = TRUE)#[1:sqrt(n)]
    d = getElbows(vals, plot = F)[elbow]
    selected.eigs <- which(abs(eig$values) >= vals[d])
    eig <- eig$vectors[,selected.eigs, drop = FALSE]
  }else {
    eig <- eigs(as(A, "dgeMatrix"), k = d)$vectors 
  }
  return(eig)
}


#' 
#' Function to estimated the score matrices of a list of graphs given the common invariant subspace V
#' 
#' @param Adj_list list of adjacency matrices, of size n x n
#' @param V common invariant subspace. A matrix of size n x d.
#' @return A list containing the score matrices
#' 
#' @author Jes\'us Arroyo <jesus.arroyo@jhu.edu>
project_networks <- function(Adj_list, V) {
  require(Matrix)
  lapply(Adj_list, function(A) crossprod(crossprod(A, V), V))
}
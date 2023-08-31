shewhart.rules <- function(object, limits = object$limits, run.length = qcc.options("run.length"))
{
  # Return a list of cases beyond limits and violating runs
  bl <- beyond.limits(object, limits = limits)
  vr <- violating.runs(object, run.length = run.length)
  list(beyond.limits = bl, violating.runs = vr)
}

beyond.limits <- function(object, limits = object$limits)
{
  # Return cases beyond limits
  statistics <- c(object$statistics, object$newstats) 
  lcl <- limits[,1]
  ucl <- limits[,2]
  index.above.ucl <- seq(along = statistics)[statistics > ucl]
  #index.below.lcl <- seq(along = statistics)[statistics < lcl]
  return(c(index.above.ucl))#, index.below.lcl))
}


rdpg.sample.weight <- function(X) {
  P <- X %*% t(X)
  n <-  nrow(P)
  U <- matrix(0, nrow = n, ncol = n)
  U[col(U) > row(U)] <- runif(n*(n-1)/2)
  U <- (U + t(U))
  diag(U) <- runif(n)
  A <- (U < P) + 0 ;
  diag(A) <- 0
  return(graph.adjacency(A,"undirected", weighted = TRUE))
}


# This function creates a customized quality control chart (QCC) visualization.
# It relies on ggplot2 for rendering the chart, and takes into consideration
# control limits, statistics, and various design elements.
#
# Parameters:
#   x: The primary QCC object to be plotted.
#   m2: Not used in the function, deprecated.
#   add.stats: Not used in the function, deprecated.
#   chart.all: Not used in the function, deprecated.
#   s: Not used in the function, deprecated.
#   label.limits: Labeling for the control limits.
#   title: The title for the plot.
#   xlab: X-axis label (unused in this function).
#   ylab: Y-axis label (unused in this function).
#   ylim: Limits for the y-axis (unused in this function).
#   axes.las: Orientation of axis labels (unused in this function).
#   digits: Precision for numbers.
#   minx: Custom labels or tick marks for the x-axis.
#   restore.par: Restore graphical parameters (unused in this function).
#   ... : Additional arguments passed on.

plot.qcc <- function(x, m2, add.stats = TRUE, chart.all = TRUE, s=s,
                     label.limits = "UCL", 
                     title, xlab, ylab, ylim, axes.las = 0, minx=FALSE,
                     digits =  getOption("digits"),
                     restore.par = TRUE, ...) 
{
  # Extracting necessary information from the qcc object
  object <- x
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
  
  # Logic for which indices to include based on 'chart.all' argument
  if(chart.all) {
    statistics <- c(stats, newstats)
    indices <- 1:length(statistics) 
  }
  
  tmax <- length(indices)+1
  
  # Logic to generate x-axis labels
  if (length(minx) == 1) {
    # Some custom logic to format month abbreviations
    tvec <- (4+s):12
    n1 = running(tvec, width=2, fun = function(x) paste0(as.character(month.abb[x[1]]), ":", as.character(month.abb[x[2]])))
    tvec <- 1:4
    n3 = running(tvec, width=2, fun = function(x) paste0(as.character(month.abb[x[1]]), ":", as.character(month.abb[x[2]])))
    minx <- as.vector(c(n1,"Dec:Jan",n3))
  }
  
  # Identifying data points violating the control limits
  vioinx <- rep(0, length(indices))
  runinx <- rep(0, length(indices))
  vioinx[violations$beyond.limits] <- 1
  runinx[violations$violating.runs] <- 1
  
  # Creating a data frame to be used in ggplot
  df <- data.frame(time=indices, y=statistics, lcl=lcl, ucl=ucl, center=center, lim=vioinx, run=runinx)
  
  # Building the ggplot visualization
  p <- ggplot(df, aes(x=time, y=y)) +
    # Lines and points to represent the data
    geom_step(aes(x=time, y=ucl), linetype="dashed") +
    geom_step(aes(x=time, y=center), linetype="solid") +
    geom_point(data = df, alpha=1, color="grey70", shape=20) +
    geom_line(aes(x=time, y=y), color="grey70") +
    geom_point(data = df %>% filter(lim==1), alpha=1, color="red", shape=17) + 
    scale_x_discrete(name ="time points", limits=c(minx)) +
    theme_bw() +
    annotate("text", label = "UCL", x = tmax-1.6, y = rev(ucl)[1]+3) +
    annotate("text", label = "CL", x = tmax-1.6, y = rev(center)[1]+1 ) +
    theme_classic(base_size = 18) +
    theme(axis.title = element_text(size = 10), axis.text.x = element_text(size = 10, angle=90), legend.position = "none", plot.title = element_text(hjust = 0.5, size=12, face='bold'), plot.subtitle = element_text(hjust = 0.5)) +
    ggtitle(title)
  
  return(p)
  invisible()
}


# This function creates a customized quality control chart (QCC) visualization for a 'vertex' perspective.
# It relies on ggplot2 for rendering the chart and takes into consideration
# control limits, statistics, and various design elements.
# 
# Parameters:
#   x: List of QCC objects to be plotted.
#   m2: Not used in the function, deprecated.
#   add.stats: Not used in the function, deprecated.
#   chart.all: Not used in the function, deprecated.
#   s: Not used in the function, deprecated.
#   label.limits: Labeling for the control limits.
#   title: The title for the plot.
#   xlab, ylab, ylim: Axis labels and limits (unused in this function).
#   axes.las: Orientation of axis labels (unused in this function).
#   minx: Custom labels or tick marks for the x-axis.
#   digits: Precision for numbers (unused in this function).
#   restore.par: Restore graphical parameters (unused in this function).
#   artpoints: Parameter for potential customization (unused in this function).
#   ... : Additional arguments passed on.

plot.qcc.vertex <- function(x, m2, add.stats = TRUE, chart.all = TRUE, s=s,
                            label.limits = "UCL", 
                            title, xlab, ylab, ylim, axes.las = 0, minx=FALSE,
                            digits =  getOption("digits"),
                            restore.par = TRUE, artpoints=NA, ...) 
{
  l.length <- length(x)
  df <- c() # Initialize the dataframe
  
  # Loop through each QCC object in the list x
  for (i in seq(l.length)) {
    object <- x[[i]]
    
    # Check if the object is of class 'qcc'
    if ((missing(object)) | (!inherits(object, "qcc")))
      stop("an object of class `qcc' is required")
    
    # Extract necessary information from the qcc object
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
    
    # Logic for which indices to include based on 'chart.all' argument
    if(chart.all) {
      statistics <- c(stats, newstats)
      indices <- 1:length(statistics) 
    }
    
    # Identifying data points violating the control limits
    vioinx <- rep(0, length(indices))
    runinx <- rep(0, length(indices))
    vioinx[violations$beyond.limits] <- 1
    runinx[violations$violating.runs] <- 1
    
    # Creating a data frame and binding it with the main dataframe
    idf <- data.frame(vertex=indices, y=statistics, timepoints=factor(minx[i], levels = minx), lcl=lcl, ucl=ucl, center=center, lim=vioinx, run=runinx)
    df <- rbind(df, idf)
  }
  
  # Building the ggplot visualization
  p <- ggplot(df, aes(x=vertex, y=y)) +
    geom_point(data=df, alpha=.1, color="grey70", shape=20, size=.5) +
    geom_point(data=df %>% filter(lim==1), alpha=.3, color="red", shape=17, size=1) +
    geom_hline(aes(yintercept=center), linetype="solid") +
    geom_hline(aes(yintercept=ucl), linetype="dashed") +
    labs(x="vertex") +
    facet_wrap(~timepoints, nrow=3, switch="x", scales="fixed") +
    theme_bw() +
    ylab(TeX("$y_{i}^{(t)}$")) +
    theme(axis.text.x=element_text(angle=45), legend.position="none", plot.title=element_text(hjust=0.5, size=10, face='bold'), plot.subtitle=element_text(hjust=0.5)) +
    ggtitle(title)
  
  return(p)
  invisible()
}


# This function calculates latent positions for multiple graphs via COSIE model, given individual ASE.
# Arguments:
#   A: A list of adjacency matrices.
#   latpos.list: A list of latent position matrices obtained by individual ASE.
#   d1, d2, d3: Dimensions parameters. d1 is for joint SVD, d3 is for individual ASE. d2 = d1 is fixed.
#   center: Whether to center the adjacency matrices.
#   verbose: Whether to print progress information.
#   python: A flag (unused).
#   SVD: Determines the method to construct latent positions in the COSIE models. 1: ASE(VRV^T) 2: VR 3: V|R|^{1/2}V^{T}

jrdpg.latent <- function(A, latpos.list, d1=2,d2=1, d3=1,  center=FALSE, verbose=FALSE,python=TRUE, SVD=3){
  d2 = d1
  
  m <- length(A)
  A <- lapply(A, function(x) (x) )
  
  if (center) {
    if (verbose) cat("do centering...\n")
    Abar <- Reduce('+', A) / length(A)
    A <- lapply(A, function(x) as.matrix(x-Abar))
  }   
  jrdpg <- mase(A, latpos.list,d1, d_vec=d3, scaled.ASE = TRUE,par = FALSE,show.scree.results = FALSE)
  
  if(SVD==3){
    d2 <- dim(jrdpg$R[[1]])[1]
    if (is.na(d2)){
      B.svd <- lapply(jrdpg$R, function (x) svd(x) )
      Xhat <- list()
      for (i in 1:m) {
        Xhat[[i]] <-  (jrdpg$V) %*% as.matrix(B.svd[[i]]$u) %*% diag(sqrt(B.svd[[i]]$d)) %*% diag(sign((B.svd[[i]]$d)) ) %*% t(as.matrix((B.svd[[i]]$u)) )  
      }
    }else{
      B.svd <- lapply(jrdpg$R, function (x) svd(x) )
      Xhat <- lapply(B.svd, function(x) (jrdpg$V) %*% as.matrix(x$u) %*% diag(sqrt(x$d), nrow=d2, ncol=d2) %*% diag(sign(x$d)) %*% t(as.matrix(x$u) )  ) 
    }
    
  }else if (SVD==2){
    Xhat <- lapply(jrdpg$R, function(x) as.matrix((jrdpg$V) %*% x) )
  }else{
    #(1)
    P <- lapply(jrdpg$R, function (x) as.matrix((jrdpg$V) %*% x %*% t(jrdpg$V)) )
    Xhat <- lapply(P, function(x) as.matrix( truncated.ase(x, d2,diagaug=FALSE, approx=FALSE )$Xhat))  
  }
    

  return(list(Xhat=Xhat))
}

# This function computes the test stats given the Adjacency Spectral Embeddings (MASE) of a list of graphs.
# Arguments:
#   glist: A list of igraph objects.
#   latpos.list: A list of latent position matrices.
#   nmase: time span you used for embeddings. 2 or others.
#   dmax, dsvd, d3: Dimension parameters. dmax is for joint SVD, d3 is for individual ASE.
#   center: Whether to center the adjacency matrices.
#   approx: Whether to use approximation method.
#   python: A flag (unused).
#   SVD: Determines the method to construct latent positions.
#   attrweigth: Specifies the attribute to be used for the adjacency matrix if the graph has weighted edges.
#newmase
doMase <- function(glist, latpos.list, nmase=2, dmax, dsvd,d3=1, center=FALSE, approx=TRUE, python=TRUE,SVD=3,attrweigth=NULL)
{
  n <- vcount(glist[[1]])
  tmax <- length(glist)
  
  if (nmase == 2) {
    Xhat <- lapply(1:(tmax-1), function(x) jrdpg.latent(lapply(glist[x:(x+1)], function(y) get.adjacency(y,attr=attrweigth) ), latpos.list[x:(x+1)] , dmax, dsvd,d3, center=center,python = python, SVD = SVD))
    norm <- sapply(1:(tmax-1), function(x) norm((Xhat[[x]]$Xhat)[[1]] - (Xhat[[x]]$Xhat)[[2]], "2"))
    pdist <- sapply(1:(tmax-1), function(x) pdistXY((Xhat[[x]]$Xhat)[[1]] , (Xhat[[x]]$Xhat)[[2]]))
    dsvd <- sapply(Xhat, function(x) x$d)
  } else {
    adj <- lapply(glist, function(y) get.adjacency(y, attr = attrweigth) )
    Xhat <- jrdpg.latent(adj, latpos.list,dmax,dsvd,d3, center=center,python = python, SVD = SVD)
    norm <- sapply(1:(tmax-1), function(x) norm(Xhat$Xhat[[x]]-Xhat$Xhat[[x+1]], "2"))
    pdist <- sapply(1:(tmax-1), function(x) pdistXY(Xhat$Xhat[[x]], Xhat$Xhat[[x+1]]))
    dsvd <- Xhat$d
  }
  
  return(list(tnorm=norm, pdist=pdist,d=dsvd))
}
# This `mase` function applies the Multiple Adjacency Spectral Embedding (MASE) method 
# to a list of adjacency matrices to embed them into a common latent position space via SVD, 
# give the individual Adjacency Spectral Embeddings (ASE) of a list of graphs.
# 
# Parameters:
#   Adj_list: A list of adjacency matrices.
#   latpos.list: A list of estimated individual latent positions for each graph.
#   d: Number of dimensions for the embedding. Default is NA.
#   d_vec: Not used.
#   scaled.ASE: Not used
#   diag.augment: Not used
#   elbow_graph: Not used
#   elbow_mase: Number of elbows to consider for joint SVD in MASE.
#   show.scree.results: Whether to display the scree plot results.
#   par: Not used
#   numpar: Not used

mase <- function(Adj_list, latpos.list, d = NA, d_vec = NA,
                 scaled.ASE = TRUE, diag.augment = TRUE, 
                 elbow_graph = 1, elbow_mase = 2,
                 show.scree.results = FALSE,
                 par = FALSE, numpar = 12) {
  
  # Replicating the dimensions vector for all adjacency matrices in the list
  d_vec = rep(d_vec, length(Adj_list))
  
  # Combining all estimated latent positions
  V_all  <- Reduce(cbind, latpos.list)
  
  # Performing Singular Value Decomposition (SVD) on the combined latent positions
  require(rARPACK)
  jointsvd <- svd(V_all)
  
  # If the number of dimensions 'd' is not provided, determine it using the elbow method
  if(is.na(d)) {
    # Optionally, show the estimated dimensions histogram for each graph
    if(show.scree.results) {
      hist(sapply(latpos.list, ncol), main = "Estimated d for each graph")
    }
    
    # Apply the elbow method on the singular values from the SVD
    result = getElbows(jointsvd$d, plot = show.scree.results)
    
    # Assign the result based on the outcome of the elbow method
    if(length(result) == 1){
      d = result
    } else {
      d = result[elbow_mase]
    }
  }
  
  # Extract the leading eigenvectors corresponding to the determined number of dimensions
  V = jointsvd$u[, 1:d, drop = FALSE]
  
  # Project the graphs onto the common latent position space
  R <- project_networks(Adj_list, V)
  
  # Return the embedded positions, singular values, and projections
  return(list(V = V, sigma = jointsvd$d, R = R))
}


# This function computes the truncated Adjacency Spectral Embedding (ASE) for an adjacency matrix, choose elbow for first 200 dimensions.
# Arguments:
#   A: The adjacency matrix.
#   d: The dimension for the embedding.
#   diagaug: Whether to augment the diagonal of the adjacency matrix.
#   approx: Whether to use approximation methods.
#   elbow: The number of elbows to consider for dimension selection.
truncated.ase <- function(A, d=NULL, diagaug=TRUE, approx=TRUE, elbow=2) {
  #    set.seed(123)
  n <- nrow(A)
  dmax <- min(nrow(A), ncol(A)-1, 200)
  show.scree.results <- FALSE

  if (diagaug) {
    diag(A) <- rowSums(A) / (nrow(A)-1)
  }
  
  if (approx) {

    if (is.null(d)){
      A.svd <- irlba::irlba(A, dmax, maxit=10000, tol=1e-10)

      result = getElbows(A.svd$d, plot = show.scree.results)#[elbow_mase]
      if(length(result)==1){
        d = result
      }else{
        d = result[elbow]
      }
    }
    A.svd <- irlba::irlba(A,d, maxit=10000, tol=1e-10)
  } else {
    if (is.null(d)){
      A.svd <- svd(A,dmax)
      d <- dim_select(A.svd$d)
    }
    A.svd <- svd(A,d)
  }
  
  Xhat <- as.matrix(A.svd$u) %*% diag(sqrt(A.svd$d), nrow=d, ncol=d)
  return(list(eval=A.svd$d, Xhat=Matrix(Xhat), d=d))
}

# This function constructs the omnibus matrix for 2 graphs.
# Alist: a list of n x n adjacency matrices or igraph objects
myfast2buildOmni <- function(Alist, diagaug=FALSE, center=FALSE, verbose=FALSE, attrweigth=NULL )
{
  require(igraph)
  require(Matrix)
  
  if (class(Alist[[1]]) == "igraph") {
    Alist <- lapply(Alist, function(x) get.adjacency(x, attr = attrweigth ) )
  }
  
  if (center) {
    if (verbose) cat("do centering...\n")
    Abar <- Reduce('+', Alist) / length(Alist)
    Alist <- lapply(Alist, function(x) as.matrix(x-Abar))
  }   
  
  if (diagaug) {
    Alist <- lapply(Alist, diagAug)
  }   
  
  omni <- cbind(Alist[[1]], (Alist[[1]]+Alist[[2]])/2 )
  omni <- rbind(omni, cbind((Alist[[1]]+Alist[[2]])/2, Alist[[2]]))
  return(omni)
}

# This function constructs the omnibus matrix for more than 2 graphs.
# Arguments are similar to 'myfast2buildOmni'.
myfast12buildOmni <- function(Alist, diagaug=FALSE, center=FALSE, verbose=FALSE, attrweigth=NULL)
{
  require(igraph)
  require(Matrix)
  
  if (class(Alist[[1]]) == "igraph") {
    Alist <- lapply(Alist, function(x) get.adjacency(x, attr = attrweigth ) )
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
  
  
  n = ncol(Alist[[1]])
  Omni1 <- lapply(Alist, function(A) (Reduce(function(U,j) cbind(U, A), c(list(A), lapply(1:(m-1), function(i) i) ))))
  Omni2 <- lapply(Alist, function(A) (Reduce(function(U,j) rbind(U, A), c(list(A), lapply(1:(m-1), function(i) i) ))))
  Omni = (Reduce(rbind, Omni1) + Reduce(cbind, Omni2))/2
  return(Omni)
}

# This function computes norms and pairwise distances for embeddings obtained from the omnibus matrix.
# Arguments:
#   glist: A list of igraph objects.
#   omniase: An object containing the results of the omnibus embedding.
#   nomni: time span you used for embeddings. 2 or others.
#   dmax: The maximum embedding dimension.
#   center: Whether to center the adjacency matrices.
#   approx: Whether to use approximation methods.
doOmni <- function(glist, omniase=a,nomni=2, dmax=NULL, center=FALSE, approx=TRUE)
{
  n <- vcount(glist[[1]])
  tmax <- length(glist)
  
  if (nomni == 2) {
    norm <- sapply(omniase, function(x) norm((x[1:n,] - x[-(1:n),]),"2"))
    pdist <- sapply(omniase, function(x) pdistXY(x[1:n,],x[-(1:n),])) 
    
  } else {
    chunk <- running(1:nrow(omniase$Xhat), width=n, by=n, fun=function(x) x)
    norm <- sapply(1:(tmax-1), function(x) norm(omniase$Xhat[chunk[,x],] - omniase$Xhat[chunk[,x+1],], "2"))
    pdist <- sapply(1:(tmax-1), function(x) pdistXY(omniase$Xhat[chunk[,x],], omniase$Xhat[chunk[,x+1],]))
  }
  return(list(tnorm=norm, pdist=pdist))
}

# realdata_latpos_list_wrapper function: Computes latent positions for a list of graphs via doing individual ASE.
# glist: List of igraph objects
# fixedd: Specified embedding dimensions. Default is NULL.
# elbow_graph: Parameter for the elbow method to choose the number of dimensions, i.e. number of elbow, default 1.
# graph_attr: Graph attribute (e.g., "weight")
realdata_latpos_list_wrapper <- function(glist, fixedd=NULL, elbow_graph=1, graph_attr="weight") {
  attrweigth = graph_attr
  # approx <- TRUE
  center <- FALSE
  # #1
  diag.augment = TRUE
  if (is.null(fixedd)) {
    d_list <- list()
    for (i in 1:length(glist)) {
      latpos <- truncated.ase( get.adjacency(glist[[i]], attr=attrweigth )  , d = NULL, approx = FALSE, elbow = elbow_graph)
      d_list <- c(d_list, latpos$d)
      
    }
    d <- round(median(unlist(d_list)))
  }else{
    d <- fixedd
  }
  
  latpos.list <- list()
  for (i in 1:length(glist)) {
    latpos <- ase( get.adjacency(glist[[i]], attr=attrweigth )  , d = d, d.max=vcount(glist[[1]]), diag.augment = diag.augment, elbow = elbow_graph)
    latpos.list <- c(latpos.list, list(latpos))
    
  }
  return(latpos.list)
}

# realdata_doMase_wrapper function: Performs MASE (Multiple Adjacency Spectral Embedding) on a list of graphs to obtain test stats with  VR to get latent positions.
# glist: List of igraph objects
# latpos.list: List of latent positions. Default is NULL.
# fixedd: Specified embedding dimensions. Default is NULL.
# elbow_graph: Parameter for the elbow method to choose the number of dimensions, i.e. number of elbow, default 1.
# embed_span: Embedding span. Default is 2.
# graph_attr: Graph attribute (e.g., "weight")
realdata_doMase_wrapper <- function(glist, latpos.list=NULL, fixedd=NULL, elbow_graph=1, embed_span=2, graph_attr="weight") {
  attrweigth = graph_attr
  py = FALSE
  method3 = 2  # latent position estimates are VR
  # approx <- TRUE
  center <- FALSE
  # #1
  diag.augment = TRUE
  if (is.null(latpos.list)) {
    latpos.list <- realdata_latpos_list_wrapper(glist, fixedd=fixedd, elbow_graph=elbow_graph, graph_attr=graph_attr)
  }
  # no used of approx
  out <- doMase(glist,latpos.list, embed_span, dmax=NA, dsvd=NA,d3=NA, center, approx=FALSE, py, method3,attrweigth)
  return(out)
}

# realdata_doMase_wrapper function: Performs MASE (Multiple Adjacency Spectral Embedding) on a list of graphs to obtain test stats with  V|R|^1/2V^T to get latent positions.
# glist: List of igraph objects
# latpos.list: List of latent positions. Default is NULL.
# fixedd: Specified embedding dimensions. Default is NULL.
# elbow_graph: Parameter for the elbow method to choose the number of dimensions, i.e. number of elbow, default 1.
# embed_span: Embedding span. Default is 2.
# graph_attr: Graph attribute (e.g., "weight")
realdata_doMase_wrapper3 <- function(glist, latpos.list=NULL, fixedd=NULL, elbow_graph=1, embed_span=2, graph_attr="weight") {
  attrweigth = graph_attr
  py = FALSE
  # latent position estimates are V|R|^1/2V^T
  method3 = 3
  center <- FALSE
  # #1
  diag.augment = TRUE
  if (is.null(latpos.list)) {
    latpos.list <- realdata_latpos_list_wrapper(glist, fixedd=fixedd,elbow_graph=elbow_graph, graph_attr=graph_attr)
  }
  # no used of approx
  out <- doMase(glist,latpos.list, embed_span, dmax=NA, dsvd=NA,d3=NA, center, approx=FALSE, py, method3,attrweigth)
  return(out)
}


# This function wraps around the Omni embedding process for a list of graphs to obtain test stats with time span 2.
# Parameters:
#   glist: A list of graph objects
#   latpos.list: (Not used in the current function, placeholder)
#   fixedd: Specified embedding dimensions. Default is NULL.
#   elbow_graph: Parameter for the elbow method to choose the number of dimensions, i.e. number of elbow, default 1.
#   embed_span: (Not used in the current function, placeholder)
#   graph_attr: Graph attribute to use. Default is "weight".
realdata_doOmni2_wrapper <- function(glist, latpos.list=NULL, fixedd=NULL, elbow_graph=1, embed_span=2, graph_attr="weight") {
  attrweigth = graph_attr
  
  # Configuration parameters for the function
  center <- FALSE
  diag.augment = TRUE
  
  # Determine the maximum dimensions based on the graph size or provided value
  if (is.null(fixedd)) {
    dmax = vcount(glist[[1]]) * 2
  } else {
    dmax = fixedd
  }
  
  allomni = list()
  # Process pairs of consecutive graphs from the list to do OMNI embeddings
  for(i in 1:(length(glist)-1)) {
    O = myfast2buildOmni(glist[i:(i+1)], center=center, diagaug=TRUE, attrweigth=attrweigth)
    a = ase(O, d = NA, diag.augment = FALSE, d.max=dmax, elbow = elbow_graph)
    allomni = c(allomni, list(a))
  }
  
  # Finalize the Omni process to calculate the adjacent difference
  out = doOmni(glist, allomni, nomni=2, dmax, center=center, approx=FALSE)
  
  return(out)
}


# This function wraps around the Omni embedding process for a list of graphs to obtain test stats with time span all.
# Parameters:
#   glist: A list of graph objects
#   latpos.list: (Not used in the current function, placeholder)
#   fixedd: Specified embedding dimensions. Default is NULL.
#   elbow_graph: Parameter for the elbow method to choose the number of dimensions, i.e. number of elbow, default 1.
#   embed_span: (Not used in the current function, placeholder)
#   graph_attr: Graph attribute to use. Default is "weight".
realdata_doOmni12_wrapper <- function(glist, latpos.list=NULL, fixedd=NULL, elbow_graph=1, embed_span=2, graph_attr="weight") {
  attrweigth = graph_attr
  
  # Configuration parameters for the function
  center <- FALSE
  diag.augment = TRUE
  dmax = fixedd
  
  # Building and applying the Omni process
  O = myfast12buildOmni(glist, center=center, diagaug=TRUE, attrweigth=attrweigth)
  a = truncated.ase(O, d=dmax, diagaug=FALSE, approx=TRUE, elbow=elbow_graph)
  out = doOmni(glist, a, nomni=12, dmax, center=center, approx=FALSE)
  
  return(out)
}

# This function wraps around the process for a list of graphs to obtain test stats using scan methods.
# Parameters:
#   glist: A list of graph objects
#   latpos.list: Not used, placeholder.
#   fixedd: Not used, placeholder.
#   elbow_graph: Not used, placeholder.
#   embed_span: Not used, placeholder.
#   graph_attr: Not used, placeholder.
realdata_doScan_wrapper <- function(glist, latpos.list=NULL, fixedd=NULL, elbow_graph=1, embed_span=2, graph_attr="weight") {
  out <- doScan(glist)
  
  return(out)
}

# get_graph_adj_pvals function: Computes adjusted p-values for given graphs on graphAD.
# out: Object containing test statistics for graphAD and vertexAD.
# t_list: List of times.
# latpos.list: List of latent positions for bootstrap. Default is NULL for no performing bootstraps.
# bootstrap_d: Dimensions for performing wrapping embedding functions on bootstrapped graphs. Default is NULL.
# elbow_graph: Parameter for the elbow method to choose the number of dimensions, i.e. number of elbow, default 1.
# realdata_wrapper: A function used as a wrapper to calculate test stat using specific methods, default is 'realdata_doMase_wrapper'.
# embed_span: Embedding span. Default is 2.
# xlab: x-axis label for plotting. Default is "time".
# title: Title for the plot. Default is NULL.
# minx: x-axis value for plotting. Default is NULL.
# t_window_size: Size of the time window to obtain null distribution. Default is 11.
# method: Method for p-value adjustment. Default is "BH", Benjamini-Hochberg Procedure.
# return_plot: Boolean to decide if a plot should be returned. Default is TRUE. Otherwise return a qcc object.
# bootstrap_verbose: Boolean to indicate if bootstrapping is used for obtaining null distribution, . Default is FALSE, only use previous graphs in moving time window to obtain null distribution
# graph_attr: Graph attribute. Default is NULL.
get_graph_adj_pvals <- function(out, t_list, latpos.list, bootstrap_d=NULL, elbow_graph=1,
                                realdata_wrapper=realdata_doMase_wrapper, embed_span=2, xlab="time", 
                                title=NULL, minx=NULL, t_window_size=11, method="BH", return_plot=TRUE, 
                                bootstrap_verbose=FALSE, graph_attr=NULL) {
  
  # Length of the time-normalized data 
  tmax <- length(out$tnorm)
  
  # Alpha level for significance testing
  alpha_level <- 0.05
  
  # Method for adjusting p-values (default is Benjamini-Hochberg)
  adjust.method <- method
  
  # Initialize vectors for raw and adjusted p-values
  pvalues <- rep(0, tmax - t_window_size)
  p.adj.pvalues <- rep(0, tmax - t_window_size)
  
  # Number of bootstrapped samples to generate null distribution
  nmB <- 20 
  
  bootstrap_attr = graph_attr
  
  # Loop through each time point
  for (i in 1:(tmax - t_window_size)) {
    
    # If bootstrap_verbose is TRUE, compute p-values using bootstrapping
    if (bootstrap_verbose) {
      null.dist <- c() # Placeholder for the null distribution
      for (b in 1:nmB) {
        # Create random graphs based on the latent positions
        glist <- lapply(latpos.list[(i):(i+t_window_size)], function(x) rdpg.sample(x))
        
        # Get test stats from the random graphs
        temp <- realdata_wrapper(glist, fixedd=bootstrap_d, elbow_graph=elbow_graph,
                                 embed_span=embed_span, graph_attr=NULL)
        null.dist <- c(null.dist, temp$tnorm)
      }
      # Compute empirical p-value based on the bootstrapped null distribution
      pvalues[i] <- sum(out$tnorm[i+t_window_size] < null.dist) / (t_window_size * nmB)
    }
    # Otherwise, compute p-values based on observed data
    else {
      pvalues[i] <- sum(out$tnorm[i+t_window_size] < out$tnorm[(i):(i+t_window_size-1)]) / t_window_size
    }
  }
  
  # Adjust p-values for multiple testing
  p.adj.pvalues <- p.adjust(pvalues, adjust.method)
  
  # If 'return_plot' is TRUE, generate a plot of the adjusted p-values
  if (return_plot) {
    indices <- 1:length(pvalues)
    
    # Format time labels
    if (is.null(minx)) {
      minx <- lapply(t_list, function(s) substr(s, 3, nchar(s)))
      minx <- paste(minx, lead(minx), sep = ":")
      minx <- minx[(t_window_size+1):tmax]
    }
    lcl <- alpha_level
    dfm2 <- data.frame(time=indices, y=p.adj.pvalues, lcl=lcl)
    
    # Plot adjusted p-values with significance threshold
    pm <- ggplot(dfm2, aes(x=time, y=y)) +
      ylab(TeX("Adjusted $p$-value")) +
      geom_line(aes(x=time, y=y), alpha=1, color="grey70", shape=20) +
      geom_point(data = dfm2 %>% filter(y < alpha_level), alpha=1, color="red", shape=17) +
      scale_x_discrete(name = xlab, limits=c(minx)) +
      theme_bw() +
      geom_hline(data=dfm2, aes(yintercept=alpha_level), color="red") +
      theme_classic(base_size = 18) +
      theme(axis.title = element_text(size = 12), axis.text.x = element_text(size = 10, angle=90),
            legend.position = "none", plot.title = element_text(hjust = 0.5, size=12, face='bold'),
            plot.subtitle = element_text(hjust = 0.5)) +
      ggtitle(title)
    
    return(pm)
  }
  else {
    return(p.adj.pvalues)
  }
}

# get_graph_qcc function: Computes QCC (Quality Control Chart) for a given set of graphs.
# gip: A list of igraph objects.
# t_list: List of times.
# latpos.list: List of latent positions. Default is NULL.
# fixedd: Specified embedding dimensions. Default is NULL.
# elbow_graph: Parameter for the elbow method to choose the number of dimensions, i.e. number of elbow, default 1.
# realdata_wrapper: Function to obtain test stats.
# embed_span: Embedding span. Default is 2.
# xlab: x-axis label for plotting. Default is "time points"
get_graph_qcc <- function(gip, t_list, latpos.list=NULL, fixedd=NULL, elbow_graph=1, 
                          realdata_wrapper=realdata_doMase_wrapper, embed_span=2, xlab="time points", 
                          title=NULL, minx=NULL, t_window_size=11, return_plot=TRUE) {
  
  # Setting the start index based on window size
  s <- t_window_size + 1
  
  # Total number of time points
  tmax <- length(gip)
  
  # Initialize vectors to store mean and standard deviation
  mean2 <- rep(1,(tmax-1)-(s-1))
  std2 <- rep(1,(tmax-1)-(s-1))
  
  # Loop through each window
  for (w in 1:(tmax-1-(s-1))) {
    
    # Check if 'latpos.list' is provided and create a sub-list accordingly
    if (is.null(latpos.list)) {
      temp.latpos.list <- NULL
    } else {
      temp.latpos.list <- latpos.list[((w):(w+s-1))]
    }
    
    # Call the data wrapper function to calculate test stats
    out2 <- realdata_wrapper(gip[(w:(w+s-1))], latpos.list=temp.latpos.list, fixedd=fixedd, 
                             elbow_graph=elbow_graph, embed_span=embed_span, graph_attr="weight")
    
    # Calculate the difference in normalized values
    diff_norm <- matrix(out2$tnorm, s-1, 1)
    
    # Store mean and standard deviation for the window
    mean2[w] <- mean(diff_norm)
    std2[w] <- sd.xbar.one(diff_norm, rep(1,s-1), "MR")
  }
  
  # Get normalized values for the entire graph
  out2 <- realdata_wrapper(gip, latpos.list=latpos.list, fixedd=fixedd, elbow_graph=elbow_graph, 
                           embed_span=embed_span, graph_attr="weight")
  df <- out2$tnorm
  diff_norm <- matrix(df, tmax-1, 1)
  
  # Generate the QCC
  c1 <- qcc(diff_norm[s:(tmax-1)], type="xbar.one", center=mean2, std.dev=std2, nsigmas=3, plot=FALSE)
  
  # If 'return_plot' is TRUE, create and return the plot
  if (return_plot) {
    
    # Format time labels if 'minx' is not provided
    if (is.null(minx)) {
      minx <- lapply(t_list, function(s) substr(s, 3, nchar(s)))
      minx <- paste(minx, lead(minx), sep=":")
      minx <- minx[(t_window_size+1):(tmax-1)]
    }
    
    # Plot the QCC
    pm <- plot.qcc(c1, add.stats=FALSE, chart.all=TRUE, s=s, minx=minx,
                   label.limits=c("LCL", "UCL"), title=paste("Control Chart", title), xlab=xlab, m2=m2, ylab="y",
                   axes.las=0, digits=getOption("digits"), restore.par=FALSE)
    
    return(pm)
  } else {
    return(c1)
  }
}
# This function computes the vertex Quality Control Chart (QCC) for a given time-series graph.
# Parameters:
#   gip - A list of igraph objects representing the time series of graphs.
#   t_list - List of time points. 
#   middle.max.inx - An optional index for the "middle max".
#   latpos.list - List of latent positions (default is NULL).
#   fixedd - Specified dimensions (default is NULL).
#   elbow_graph - Numerical value to control elbow graph plotting (default is 1).
#   realdata_wrapper - A function used as a wrapper to calculate test stat using specific methods, default is 'realdata_doMase_wrapper'.
#   embed_span - Span for the embedding (default is 2).
#   xlab - Label for the x-axis (default is "time points").
#   title - Optional title for the output plot.
#   minx - x-axis text in the plot (default is NULL).
#   t_window_size - Size of the time window (default is 11).
#   return_plot - Boolean indicating whether to return a plot (default is TRUE).
get_vertex_qcc <- function(gip, t_list, middle.max.inx=NULL, latpos.list=NULL, fixedd=NULL, elbow_graph=1, 
                           realdata_wrapper=realdata_doMase_wrapper, embed_span=2, xlab="time points", 
                           title=NULL, minx=NULL, t_window_size=11, return_plot=TRUE) {
  
  # Get the number of vertices from the first graph
  n = vcount(gip[[1]])
  
  # Setting the start index based on window size
  s <- t_window_size + 1
  
  # Total number of time points
  tmax <- length(gip)
  
  # Initialize vectors to store mean and standard deviation for each vertex
  mean2v <- rep(1,(tmax-1)-(s-1))
  std2v <- rep(1,(tmax-1)-(s-1))
  
  # Loop through each window
  for (w in 1:(tmax-1-(s-1))) {
    
    # Check if 'latpos.list' is provided and create a sub-list accordingly
    if (is.null(latpos.list)) {
      temp.latpos.list <- NULL
    } else {
      temp.latpos.list <- latpos.list[((w):(w+s-1))]
    }
    
    # Call the data wrapper function to calculate test stats
    out2 <- realdata_wrapper(gip[(w:(w+s-1))], latpos.list=temp.latpos.list, fixedd=fixedd, 
                             elbow_graph=elbow_graph, embed_span=embed_span, graph_attr="weight")
    
    # Extract pairwise distances for the current window
    df <- out2$pdist
    
    # Convert distances to a matrix format
    diff_distv <- matrix(df, s-1, 1*n, byrow=TRUE)
    
    # Store mean and standard deviation for the window for each vertex
    mean2v[w] <- mean(diff_distv)
    std2v[w] <- sd.xbar(diff_distv, rep(1*n, s-1), "UWAVE-SD")
  }
  
  # Get pairwise distances for the entire graph
  out2 <- realdata_wrapper(gip, latpos.list=latpos.list, fixedd=fixedd, elbow_graph=elbow_graph, 
                           embed_span=embed_span, graph_attr="weight")
  df <- out2$pdist
  diff_distv <- matrix(df, tmax-1, 1*n, byrow=TRUE)
  
  # Generate a list of QCCs for each vertex
  c1m2v <- list()
  for (i in s:(tmax-1)) {
    c1m2v[[i-s+1]] <- qcc(diff_distv[i,], type="xbar.one", center=mean2v[i-s+1], std.dev=std2v[i-s+1], nsigmas=3, plot=FALSE)
  }
  
  # If 'return_plot' is TRUE, create and return the vertex plot
  if (return_plot) {
    
    # Format time labels if 'minx' is not provided
    if (is.null(minx)) {
      minx <- lapply(t_list, function(s) substr(s, 3, nchar(s)))
      minx <- paste(minx, lead(minx), sep=":")
      minx <- minx[(t_window_size+1):(tmax-1)]
    }
    
    # Plot the QCC for vertices
    pm <- plot.qcc.vertex(c1m2v, add.stats=FALSE, chart.all=TRUE, s=s, minx=minx,
                          label.limits=c("LCL", "UCL"), title=paste("Control Chart", title), xlab="time points", m2=seq(n), ylab="y", 
                          axes.las=0, digits=getOption("digits"), artpoints=middle.max.inx, restore.par=FALSE)
    
    return(pm)
  } else {
    return(c1m2v)
  }
}

# preprocess
#' Run pass-to-rank on a weighted graph.
#'
#' It extracts (non-zero) edge weight vector \eqn{W} from a graph and replaces it with \eqn{2*R / (|E|+1)} where \eqn{R} is the rank of \eqn{W} and \eqn{|E|} is the number of edges. This does 'no-op' for an unweighted graph.
#'
#' @param g a graph in \code{igraph} format or an n x 2 edge list or an n x n adjacency matrix
#'
#' @author Youngser Park <youngser@jhu.edu>
#' @export
#' @import igraph
ptr <- function(g)
{
  if (class(g) != "igraph") {
    if (!is.matrix(g)) stop("the input has to be either an igraph object or a matrix!")
    else {
      if (ncol(g)==2) g <- igraph::graph_from_edgelist(g)
      else if (nrow(g)==ncol(g)) g <- igraph::graph_from_adjacency_matrix(g, weighted = TRUE)
      else stop("the input matrix is not a graph format!")
    }
  }
  
  if (igraph::is.weighted(g)) {
    W <- E(g)$weight
  } else { # no-op!
    W <- rep(1,igraph::ecount(g))
  }
  
  E(g)$weight <- rank(W)*2 / (igraph::ecount(g)+1)
  return(g)
}
#' find largest connected component in a graph
#' It extracts (non-zero) largest connected subgraph .
#'
#' @param graph a graph in \code{igraph} format
#'
#' @author Guodong Chen <gchen35@jhu.edu>
#' @export
#' @import igraph
giant.component <- function(graph, ...) {
  #require(igraph)
  cl <- igraph::clusters(graph, ...)
  igraph::induced.subgraph(graph, which(cl$membership == which.max(cl$csize)))
}
#' remove edges which has zero weights for all graphs (if any) and find jointly largest connected component in graphs. Finally it removes all self-loops,
#' It extracts (non-zero) igraph list \eqn{gip} and removes all edges with zero edge weights and return a list of jointly largest connected component in graphs without self-loops .
#'
#' @param gip a list of graphs in \code{igraph} format
#'
#' @author Guodong Chen <gchen35@jhu.edu>
#' @export
#' @import igraph
#' @return A list of graphs in igraph format
jlcc <- function(gip){
  l.length <- length(gip)
  for (i in 1:l.length) {
    gip[[i]] <- igraph::delete.edges(gip[[i]], which(E(gip[[i]])$weight==0))
  }
  #find joint largest connected component
  df1 <- igraph::as_data_frame(gip[[1]])[,1:2]
  df2 <- igraph::as_data_frame(gip[[2]])[,1:2]
  g <- igraph::graph.intersection(gip[[1]],gip[[2]],keep.all.vertices = TRUE)
  if(l.length>2){
    for (i in 3:l.length) {
      g <- igraph::graph.intersection(g,gip[[i]],keep.all.vertices = TRUE)
    }
  }
  lcc <- giant.component(g)
  glist <- gip
  gip <- list()
  for (i in 1:l.length) {
    gip[[i]] <- igraph::induced_subgraph(glist[[i]], V(lcc)$name);
    gip[[i]] <- igraph::permute.vertices(gip[[i]], match(V(gip[[i]])$name, V(gip[[1]])$name));
  }
  for (i in 1:l.length) {
    gip[[i]] <- igraph::simplify(gip[[i]])
  }
  return(gip)
}


# Function to add missing vertices to a graph
add_missing_vertices <- function(g, all_vertices) {
  missing_vertices <- setdiff(all_vertices, V(g)$name)
  g <- add_vertices(g, length(missing_vertices), name = missing_vertices)
  g
}


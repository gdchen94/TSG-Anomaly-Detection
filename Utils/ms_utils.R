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
# modified the one in utils.R so it works on sparse matrices.
ase <- function(A, d = NA, d.max = sqrt(ncol(A)), diag.augment = TRUE, elbow = 1) {
  require(rARPACK)
  # Diagonal augmentation
  if(diag.augment & sum(abs(diag(A))) == 0) {
    deg = colSums(A)
    n = ncol(A)
    diag(A) = deg / (n-1)
  }
  if(is.na(d)) {
    eig <- eigs(A, d.max)
    vals <- sort(x =  abs(eig$values), decreasing = TRUE)
    d = getElbows(vals, plot = F)[elbow]
    selected.eigs <- which(abs(eig$values) >= vals[d])
    V <- eig$vectors[,selected.eigs, drop = F]
    D <- diag(sqrt(abs(eig$values[selected.eigs])), nrow = d)
    X <- V %*% D
    return(X)
  } else {
    eig <- eigs(A, k = d)
    X <- eig$vectors %*% diag(sqrt(abs(eig$values)), nrow = d)
    return(X) 
  }
}
# Function to create a control chart for vertexAD across time points
plot.qcc.vertex <- function(x, m2, add.stats = TRUE, chart.all = TRUE, s=s,
                            label.limits = "UCL", #c("LCL ", "UCL"),
                            title, xlab, ylab, ylim, axes.las = 0,
                            digits =  getOption("digits"),
                            restore.par = TRUE, artpoints=NA, ...) 
{
  # Initialize the length of input list and an empty dataframe
  l.length <- length(x)
  df <- c() #c(2,6)
  
  # Loop through each item in the input list
  for (i in seq(l.length)) {
    object <- x[[i]]
    
    # Validate input object class
    if ((missing(object)) | (!inherits(object, "qcc")))
      stop("an object of class `qcc' is required")
    
    # Extract relevant data from the object
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
    
    # Combine stats if 'chart.all' option is TRUE
    if(chart.all) {
      statistics <- c(stats, newstats)
      indices <- 1:length(statistics)
    }
    
    # Construct timepoint labels for the plot
    tvec <- (4+s):12
    n1=running(tvec, width=2, fun = function(x) paste0(as.character(month.abb[x[1]]), ":", as.character(month.abb[x[2]])))
    tvec <- 1:4
    n3=running(tvec, width=2, fun = function(x) paste0(as.character(month.abb[x[1]]), ":", as.character(month.abb[x[2]])))
    minx <- as.vector(c(n1,"Dec:Jan",n3))
    
    # Identify violations and runs
    vioinx <- rep(0,length(indices))
    runinx <- rep(0,length(indices))
    vioinx[violations$beyond.limits] <- 1
    runinx[violations$violating.runs] <- 1
    
    # Compile all data into a dataframe
    idf <- data.frame(vertex=indices, y=statistics, timepoints=factor(minx[i], levels=minx), lcl=lcl, ucl=ucl, center=center, lim=vioinx, run=runinx)
    df <- rbind(df, idf)
  }
  
  # Plot using ggplot2
  p <- ggplot(df, aes(x=vertex, y=y)) +
    geom_point(data=df, alpha=.1, color="grey70", shape=20, size=.5) +
    geom_point(data=df %>% filter(lim==1), alpha=1, color="red", shape=17, size=1) +
    geom_point(data=df %>% filter(vertex %in% artpoints), color="blue", shape=0, size=0.7) +
    geom_hline(aes(yintercept=center), linetype="solid") +
    facet_wrap(~timepoints, nrow=3, switch="x", scales="fixed") +
    geom_hline(aes(yintercept=ucl), linetype="dashed") +
    labs(x="vertex") +
    theme_bw() +
    ylab(TeX("$\\log(y_{i}^{(t)})$")) +
    theme(axis.title.x=element_text(size=16),
          axis.title.y=element_text(size=16),
          axis.text.x=element_text(size=12, angle=90),
          axis.text.y=element_text(size=12),
          strip.text=element_text(size=14),
          legend.position="none",
          plot.title=element_text(hjust=0.5, size=16, face='bold'),
          plot.subtitle=element_text(hjust=0.5)) +
    ggtitle(title)
  
  # Return the ggplot object
  return(p)
  invisible() # Ensures that function returns only the ggplot object and not any additional information
}


# Function to create a zoomed-in control chart for vertexAD across specific time points
plot.qcc.vertex.zoom <- function(x, m2, add.stats = TRUE, chart.all = TRUE, s=s,
                                 label.limits = "UCL", #c("LCL ", "UCL"),
                                 title, xlab, ylab, ylim, axes.las = 0,
                                 digits =  getOption("digits"),
                                 restore.par = TRUE, artpoints=NA, ...) 
{
  # Get the length of input list
  l.length <- length(x)
  df <- c() 
  
  # Loop through specific indices in the input list (i.e., 2, 3, 4, 5)
  for (i in c(2,3,4,5)) {
    object <- x[[i]]
    
    # Validate if the object is of class `qcc`
    if ((missing(object)) | (!inherits(object, "qcc")))
      stop("an object of class `qcc' is required")
    
    # Extract relevant information from the object
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
    
    # If 'chart.all' is TRUE, combine stats
    if(chart.all) {
      statistics <- c(stats, newstats)
      indices <- 1:length(statistics)
    }
    
    # Construct timepoint labels for the plot
    tvec <- (4+s):12
    n1 = running(tvec, width=2, fun = function(x) paste0(as.character(month.abb[x[1]]), ":", as.character(month.abb[x[2]])))
    tvec <- 1:4
    n3 = running(tvec, width=2, fun = function(x) paste0(as.character(month.abb[x[1]]), ":", as.character(month.abb[x[2]])))
    minx <- as.vector(c(n1, "Dec:Jan", n3))
    
    # Identify violations and runs
    vioinx <- rep(0,length(indices))
    runinx <- rep(0,length(indices))
    vioinx[violations$beyond.limits] <- 1
    runinx[violations$violating.runs] <- 1
    
    # Construct a dataframe with all the extracted and processed data
    idf <- data.frame(vertex=indices, y=statistics, timepoints=factor(minx[i], levels=minx), lcl=lcl, ucl=ucl, center=center, lim=vioinx, run=runinx)
    df <- rbind(df, idf)
  }
  
  # Use ggplot2 to create the visualization
  p <- ggplot(df, aes(x=vertex, y=y)) +
    geom_point(data=df, alpha=.1, color="grey70", shape=20, size=.5) +
    geom_point(data=df %>% filter(lim==1), color="red", shape=17, size=1) +
    geom_point(data=df %>% filter(vertex %in% artpoints), color="blue", shape=0, size=0.7) +
    geom_hline(aes(yintercept=center), linetype="solid") +
    facet_wrap(~timepoints, nrow=3, switch="x", scales="fixed") +
    geom_hline(aes(yintercept=ucl), linetype="dashed") +
    labs(x="vertex") +
    theme_bw() +
    ylab(TeX("$\\log(y_{i}^{(t)})$")) +
    theme(axis.title.x=element_text(size=16),
          axis.title.y=element_text(size=16),
          axis.text.x=element_text(size=12, angle=90),
          axis.text.y=element_text(size=12),
          strip.text=element_text(size=14),
          legend.position="none",
          plot.title=element_text(hjust=0.5, size=16, face='bold'),
          plot.subtitle=element_text(hjust=0.5)) +
    ggtitle(title)
  
  # Return the ggplot object
  return(p)
  invisible() # Ensures that function returns only the ggplot object and not any additional information
}





options(warn = -1)  # suppress warnings
# options(warn = 0)  # show warnings

# Install and load the packages if not already installed
packages <- unique(c("ggrepel", "readr", "anytime", "tidyverse", "lubridate", "pracma",
                     "igraph", "gridExtra", "ggplotify", "rstudioapi", "changepoints",
                     "gtools", "rARPACK", "metafolio", "latex2exp", "qcc", "rstudioapi",
                     "doParallel", "irlba", "pdp", "reshape2","mclust","factoextra",
                     "pROC", "kernlab", "graphkernels","foreach"))

for(pkg in packages) {
  if(!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg)
  }
  library(pkg, character.only = TRUE)
}


# Get the path of the currently active document in RStudio
current_script_path <- getSourceEditorContext()$path

# Set the working directory to the directory containing the current script
setwd(dirname(current_script_path))

source("../Utils/utils.R")
source("../Utils/qcc.R")
source("../Utils/enron_utils.R")
load("../Data/df.tsg-extinction-98.RData")

# R version 4.3.1
# This code will take about 16 hours to run on a regular PC.
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
    if (minx[i] != "") {
      idf <- data.frame(vertex=indices, y=statistics, timepoints=minx[i], lcl=lcl, ucl=ucl, center=center, lim=vioinx, run=runinx)
      df <- rbind(df, idf) 
    }
  }
  desired_order <- c("15:16", "19:20", "79:80", "139:140")
  
  # Convert the timepoints to a factor with the specified levels
  df$timepoints <- factor(df$timepoints, levels = desired_order)
  # Building the ggplot visualization
  p <- ggplot(df, aes(x=vertex, y=y)) +
    geom_point(data=df, alpha=.1, color="grey70", shape=20, size=.5) +
    geom_point(data=df %>% filter(lim==1), alpha=1, color="red", shape=17, size=1) +
    geom_hline(aes(yintercept=center), linetype="solid") +
    geom_point(data=df %>% filter(vertex %in% artpoints), color="blue", shape=0, size=0.7) +
    geom_hline(aes(yintercept=ucl), linetype="dashed") +
    labs(x="vertex") +
    facet_wrap(~timepoints, nrow=2, switch="x", scales="fixed") +
    theme_bw() +
    ylab(ylab) +
    theme(axis.text.x=element_text(angle=45), legend.position="none", plot.title=element_text(hjust=0.5, size=10, face='bold'), plot.subtitle=element_text(hjust=0.5)) +
    ggtitle(title)
  
  return(p)
  invisible()
}

plot.qcc.all <- function(x, add.stats = TRUE, chart.all = TRUE, s=s,
                         label.limits = "UCL", 
                         axes.las = 0, minx=FALSE,
                         digits =  getOption("digits"),
                         restore.par = TRUE, last_df=NULL, last_inx=0, ...) 
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
  
  
  # Identifying data points violating the control limits
  vioinx <- rep(0, length(indices))
  runinx <- rep(0, length(indices))
  vioinx[violations$beyond.limits] <- 1
  runinx[violations$violating.runs] <- 1
  
  # Creating a data frame to be used in ggplot
  
  df <- data.frame(time = indices, y = statistics, lcl = lcl, ucl = ucl, center = center, lim = vioinx, run = runinx)
  if (!is.null(last_df)) {
    df$inx <- last_inx + 1
    df <- rbind(last_df, df)
  } else {
    df$inx <- 0
  }
  return (df)
  
}

#print(load(url("https://www.cis.jhu.edu/~parky/MBTSG/df.tsg-extinction-98.RData")))

n_layers = length(df.tsg$tg)

# Specify the directory where you want to save the file
output_dir <- paste0(dirname(current_script_path), "_realdata_plots_", format(Sys.time(), "%Y-%m-%d"))

# Check if the directory exists, if not, create it
if (!dir.exists(output_dir)) {
  dir.create(output_dir)
}

adj_pval_pm2_list <- list()
adj_pval_pm12_list <- list()
adj_pval_po2_list <- list()
adj_pval_po12_list <- list()
adj_pval_ps2_list <- list()
cpt_hats <- list()
for (i in seq(n_layers)) {
  # check simple graph
  unlist(lapply(df.tsg$tg[[i]], function(x) is.simple(x)))
  
  # check undirected graph
  unlist(lapply(df.tsg$tg[[i]], function(x) is.directed(x)))
  
  A_list <- lapply(df.tsg$tg[[i]], function(g) as_adjacency_matrix(as.undirected(g), attr = "weight", sparse = FALSE))
  gip <- lapply(df.tsg$tg[[i]], function(g) as.undirected(g))
  
  aip <- lapply(gip, function(g) as_adjacency_matrix(g, attr = "weight", sparse = FALSE))
  tmax <- length(A_list)
  
  # WBS
  data_mat_list <- lapply(aip, function(x) x[lower.tri(x, diag = FALSE)])
  data_mat <- do.call(cbind, data_mat_list)
  M = 10 # number of random intervals for WBS
  latpos.list.wrapped <- realdata_latpos_list_wrapper(gip, dASE=NA, d.max="sqrt", approx=FALSE, elbow_graph=2, attr_weight="weight")
  
  d = median(do.call(cbind, lapply(latpos.list.wrapped, function(x) ncol(x)))) # parameter for scaled PCA algorithm
  delta = 2
  intervals = WBS.intervals(M = M, lower = 1, upper = ncol(data_mat))
  if (i!=97){
    WBS_result = WBS.nonpar.RDPG(data_mat, lowerdiag = TRUE, d,
                                 Alpha = intervals$Alpha, Beta = intervals$Beta, delta)
    cpt_hat = tuneBSnonparRDPG(WBS_result, data_mat, lowerdiag = TRUE, d)
    cpt_hats <- append(cpt_hats, list(cpt_hat))
    
  }
  
  # Generate all pairs using the running function
  # all_pairs <- lapply(seq(tmax), function(x) as.character(x))
  
  # Define a function to check if a pair should be excluded
  should_exclude <- function(pair) {
    # Extract the numbers from the pair string
    nums <- as.numeric(pair)
    # Check if the pair matches the pattern x2:x3, x3:x4, ..., x9:x10
    return(nums[1] %% 10 != 1)
  }
  
  # Replace the pairs that match the exclusion pattern with an empty string
  t_list <- seq(tmax)
  minx_all <- rep("", tmax-1)
  minx_all[16] <- "15:16"
  minx_all[20] <- "19:20"
  minx_all[80] <- "79:80"
  minx_all[140] <- "139:140"
  
  # d.max <- "sqrt"
  attr_weight <- "weight"
  # Wrapper function to retrieve the latent position list based on the given graphs and attributes
  latpos.list.wrapped <- realdata_latpos_list_wrapper(gip, dASE=NA, d.max="sqrt", approx=FALSE, elbow_graph=2, attr_weight="weight")
  
  
  # Iterate over different time window sizes for analysis
  for (t_window_size in c(12)) {#seq(10, 12)) {
    
    minx <- minx_all[(t_window_size+1):(tmax-1)]
    
    # Apply MASE, OMNI, and SCAN wrappers to process and prepare the graph to get test stats
    out2wrapper <- realdata_doMase_wrapper(gip, embed_span=2, attr_weight=attr_weight)
    out12wrapper <- realdata_doMase_wrapper(gip, embed_span=12,  attr_weight=attr_weight)
    # out2omniwrapper <- realdata_doOmni2_wrapper(gip,  attr_weight=attr_weight)
    # out12omniwrapper <- realdata_doOmni12_wrapper(gip, attr_weight=attr_weight)
    # outscanwrapper <- realdata_doScan_wrapper(gip, attr_weight=attr_weight)
    
    # Compute the adjusted p-values for various methods
    # adj_pval_pm2 <- get_graph_adj_pvals(out2wrapper, t_list, latpos.list=latpos.list.wrapped, realdata_wrapper=realdata_doMase_wrapper, embed_span=2, xlab="time (yy/mm)", title="MASE(2)",minx=minx,  t_window_size=t_window_size, return_plot=TRUE, number_bootstrap=50)
    # adj_pval_pm12 <- get_graph_adj_pvals(out12wrapper, t_list, latpos.list=latpos.list.wrapped, realdata_wrapper=realdata_doMase_wrapper, embed_span=12, xlab="time (yy/mm)", title="MASE(12)",minx=minx,  t_window_size=t_window_size, return_plot=TRUE, number_bootstrap=50)
    # adj_pval_po2 <- get_graph_adj_pvals(out2omniwrapper, t_list, latpos.list=latpos.list.wrapped, realdata_wrapper=realdata_doOmni2_wrapper, embed_span=2, xlab="time (yy/mm)", title="OMNI(2)",minx=minx,  t_window_size=t_window_size, return_plot=TRUE, number_bootstrap=50)
    # adj_pval_po12 <- get_graph_adj_pvals(out12omniwrapper, t_list, latpos.list=latpos.list.wrapped, realdata_wrapper=realdata_doOmni12_wrapper, embed_span=12, xlab="time (yy/mm)", title="OMNI(12)",minx=minx,  t_window_size=t_window_size, return_plot=TRUE, number_bootstrap=50)
    # adj_pval_ps2 <- get_graph_adj_pvals(outscanwrapper, t_list, latpos.list=latpos.list.wrapped, realdata_wrapper=realdata_doScan_wrapper, embed_span=2, xlab="time (yy/mm)", title="SCAN",minx=minx,  t_window_size=t_window_size, return_plot=TRUE, number_bootstrap=50)
    
    # adj_pval_pm2 <- get_graph_adj_pvals(out2wrapper, t_list, realdata_wrapper=realdata_doMase_wrapper, embed_span=2, xlab="time (yy/mm)", title="MASE(2)", minx=minx, t_window_size=t_window_size, return_plot=FALSE, number_bootstrap=0, attr_weight=attr_weight)
    # adj_pval_pm12 <- get_graph_adj_pvals(out12wrapper, t_list, realdata_wrapper=realdata_doMase_wrapper, embed_span=12, xlab="time (yy/mm)", title="MASE(12)", minx=minx, t_window_size=t_window_size, return_plot=FALSE, number_bootstrap=0, attr_weight=attr_weight)
    # adj_pval_po2 <- get_graph_adj_pvals(out2omniwrapper, t_list, realdata_wrapper=realdata_doOmni2_wrapper, embed_span=2, xlab="time (yy/mm)", title="OMNI(2)", minx=minx, t_window_size=t_window_size, return_plot=FALSE, number_bootstrap=0, attr_weight=attr_weight)
    # adj_pval_po12 <- get_graph_adj_pvals(out12omniwrapper, t_list, realdata_wrapper=realdata_doOmni12_wrapper, embed_span=12, xlab="time (yy/mm)", title="OMNI(12)", minx=minx, t_window_size=t_window_size, return_plot=FALSE, number_bootstrap=0, attr_weight=attr_weight)
    # adj_pval_ps2 <- get_graph_adj_pvals(outscanwrapper, t_list, realdata_wrapper=realdata_doScan_wrapper, embed_span=2, xlab="time (yy/mm)", title="SCAN", minx=minx, t_window_size=t_window_size, return_plot=FALSE, number_bootstrap=0, attr_weight=attr_weight)
    # 
    # # Append the results to the lists
    # adj_pval_pm2_list[[i]] <- adj_pval_pm2
    # adj_pval_pm12_list[[i]] <- adj_pval_pm12
    # adj_pval_po2_list[[i]] <- adj_pval_po2
    # adj_pval_po12_list[[i]] <- adj_pval_po12
    # adj_pval_ps2_list[[i]] <- adj_pval_ps2
    
    MC_tsg_pm2 <- get_graph_qcc(gip, t_list, latpos.list=latpos.list.wrapped, realdata_wrapper=realdata_doMase_wrapper, embed_span=2, xlab="time points", title="MASE(2)", minx=minx, t_window_size=t_window_size, return_plot=FALSE)
    MC_tsg_pm12 <- get_graph_qcc(gip, t_list, latpos.list=latpos.list.wrapped, realdata_wrapper=realdata_doMase_wrapper, embed_span=12, xlab="time points", title="MASE(12)", minx=minx, t_window_size=t_window_size, return_plot=FALSE)
    # MC_tsg_po2 <- get_graph_qcc(gip, t_list, realdata_wrapper=realdata_doOmni2_wrapper, embed_span=2, xlab="time points", title="OMNI(2)", minx=minx, t_window_size=t_window_size, return_plot=FALSE)
    # MC_tsg_po12 <- get_graph_qcc(gip, t_list, realdata_wrapper=realdata_doOmni12_wrapper, embed_span=12, xlab="time points", title="OMNI(12)", minx=minx, t_window_size=t_window_size, return_plot=FALSE)
    # MC_tsg_scan <- get_graph_qcc(gip, t_list,realdata_wrapper=realdata_doScan_wrapper,embed_span=2,xlab="time points", title="SCAN", minx=minx, t_window_size=t_window_size, return_plot=FALSE)
    
    # Append the results to the lists
    if (i == 1) {
      MC_tsg_pm2_qcc <- plot.qcc.all(MC_tsg_pm2, add.stats = FALSE, chart.all = TRUE, s = s, minx = minx, label.limits = c("LCL", "UCL"), axes.las = 0, digits = getOption("digits"), restore.par = FALSE)
      MC_tsg_pm12_qcc <- plot.qcc.all(MC_tsg_pm12, add.stats = FALSE, chart.all = TRUE, s = s, minx = minx, label.limits = c("LCL", "UCL"), axes.las = 0, digits = getOption("digits"), restore.par = FALSE)
      # MC_tsg_po2_qcc <- plot.qcc.all(MC_tsg_po2, add.stats = FALSE, chart.all = TRUE, s = s, minx = minx, label.limits = c("LCL", "UCL"), axes.las = 0, digits = getOption("digits"), restore.par = FALSE)
      # MC_tsg_po12_qcc <- plot.qcc.all(MC_tsg_po12, add.stats = FALSE, chart.all = TRUE, s = s, minx = minx, label.limits = c("LCL", "UCL"), axes.las = 0, digits = getOption("digits"), restore.par = FALSE)
      # MC_tsg_pscan_qcc <- plot.qcc.all(MC_tsg_pscan, add.stats = FALSE, chart.all = TRUE, s = s, minx = minx, label.limits = c("LCL", "UCL"), axes.las = 0, digits = getOption("digits"), restore.par = FALSE)
    } else {
      MC_tsg_pm2_qcc <- plot.qcc.all(MC_tsg_pm2, add.stats = FALSE, chart.all = TRUE, s = s, minx = minx, label.limits = c("LCL", "UCL"), axes.las = 0, digits = getOption("digits"), restore.par = FALSE, last_df = MC_tsg_pm2_qcc, last_inx = i)
      MC_tsg_pm12_qcc <- plot.qcc.all(MC_tsg_pm12, add.stats = FALSE, chart.all = TRUE, s = s, minx = minx, label.limits = c("LCL", "UCL"), axes.las = 0, digits = getOption("digits"), restore.par = FALSE, last_df = MC_tsg_pm12_qcc, last_inx = i)
      # MC_tsg_po2_qcc <- plot.qcc.all(MC_tsg_po2, add.stats = FALSE, chart.all = TRUE, s = s, minx = minx, label.limits = c("LCL", "UCL"), axes.las = 0, digits = getOption("digits"), restore.par = FALSE, last_df = MC_tsg_po2_qcc, last_inx = i)
      # MC_tsg_po12_qcc <- plot.qcc.all(MC_tsg_po12, add.stats = FALSE, chart.all = TRUE, s = s, minx = minx, label.limits = c("LCL", "UCL"), axes.las = 0, digits = getOption("digits"), restore.par = FALSE, last_df = MC_tsg_po12_qcc, last_inx = i)
      # MC_tsg_pscan_qcc <- plot.qcc.all(MC_tsg_pscan, add.stats = FALSE, chart.all = TRUE, s = s, minx = minx, label.limits = c("LCL", "UCL"), axes.las = 0, digits = getOption("digits"), restore.par = FALSE, last_df = MC_tsg_pscan_qcc, last_inx = i)
    }
    
    # adj_pval_pm2 <- get_graph_adj_pvals(out2wrapper, t_list, realdata_wrapper=realdata_doMase_wrapper, embed_span=2, xlab="time (yy/mm)", title="MASE(2)", minx=minx, t_window_size=t_window_size, return_plot=TRUE, number_bootstrap=0, attr_weight=attr_weight)
    # adj_pval_pm12 <- get_graph_adj_pvals(out12wrapper, t_list, realdata_wrapper=realdata_doMase_wrapper, embed_span=12, xlab="time (yy/mm)", title="MASE(12)", minx=minx, t_window_size=t_window_size, return_plot=TRUE, number_bootstrap=0, attr_weight=attr_weight)
    # adj_pval_po2 <- get_graph_adj_pvals(out2omniwrapper, t_list, realdata_wrapper=realdata_doOmni2_wrapper, embed_span=2, xlab="time (yy/mm)", title="OMNI(2)", minx=minx, t_window_size=t_window_size, return_plot=TRUE, number_bootstrap=0, attr_weight=attr_weight)
    # adj_pval_po12 <- get_graph_adj_pvals(out12omniwrapper, t_list, realdata_wrapper=realdata_doOmni12_wrapper, embed_span=12, xlab="time (yy/mm)", title="OMNI(12)", minx=minx, t_window_size=t_window_size, return_plot=TRUE, number_bootstrap=0, attr_weight=attr_weight)
    # adj_pval_ps2 <- get_graph_adj_pvals(outscanwrapper, t_list, realdata_wrapper=realdata_doScan_wrapper, embed_span=2, xlab="time (yy/mm)", title="SCAN", minx=minx, t_window_size=t_window_size, return_plot=TRUE, number_bootstrap=0, attr_weight=attr_weight)
    
    # MC_tsg_pm2 <- get_graph_qcc(gip, t_list, latpos.list=latpos.list.wrapped, realdata_wrapper=realdata_doMase_wrapper, embed_span=2, xlab="time points", title="MASE(2)", minx=minx, t_window_size=t_window_size, return_plot=TRUE)
    # MC_tsg_pm12 <- get_graph_qcc(gip, t_list, latpos.list=latpos.list.wrapped, realdata_wrapper=realdata_doMase_wrapper, embed_span=12, xlab="time points", title="MASE(12)", minx=minx, t_window_size=t_window_size, return_plot=TRUE)
    # MC_tsg_po2 <- get_graph_qcc(gip, t_list, realdata_wrapper=realdata_doOmni2_wrapper, embed_span=2, xlab="time points", title="OMNI(2)", minx=minx, t_window_size=t_window_size, return_plot=TRUE)
    # MC_tsg_po12 <- get_graph_qcc(gip, t_list, realdata_wrapper=realdata_doOmni12_wrapper, embed_span=12, xlab="time points", title="OMNI(12)", minx=minx, t_window_size=t_window_size, return_plot=TRUE)
    # MC_tsg_ps <- get_graph_qcc(gip, t_list,realdata_wrapper=realdata_doScan_wrapper,embed_span=2,xlab="time points", title="SCAN", minx=minx, t_window_size=t_window_size, return_plot=TRUE)
    
    
    
    # # Figures 11, 12, 13
    # plots <- list(adj_pval_pm2, adj_pval_pm12, adj_pval_po2, adj_pval_po12, adj_pval_ps2)
    # 
    # for (plot in plots) {
    #   print(plot)
    # }
    # png(paste0(output_dir, paste0("/MC_tsg_g_pval_t_window_size",as.character(t_window_size),"_", "layer_", as.character(i), "_",format(Sys.time(), "%Y_%m_%d")),".png"), width = 10, height = 9, units="in", res=400)
    # grid.arrange(adj_pval_pm2,adj_pval_pm12,adj_pval_po2,adj_pval_po12,adj_pval_ps2, nrow=3, heights=c(4,4,4))
    # dev.off()
    # 
    # # Figures 8, 9, 10
    # plots <- list(MC_tsg_pm2,MC_tsg_pm12,MC_tsg_po2,MC_tsg_po12,MC_tsg_ps)
    # 
    # for (plot in plots) {
    #   print(plot)
    # }
    # png(paste0(output_dir, paste0("/MC_tsg_g_cc_t_window_size",as.character(t_window_size),"_", "layer_", as.character(i), "_",format(Sys.time(), "%Y_%m_%d")),".png"), width = 10, height = 9, units="in", res=400)
    # grid.arrange(MC_tsg_pm2,MC_tsg_pm12,MC_tsg_po2,MC_tsg_po12,MC_tsg_ps, nrow=3, heights=c(4,4,4))
    # dev.off()
  }
}

MC_tsg_pm2_qcc_f <- MC_tsg_pm2_qcc
MC_tsg_pm12_qcc_f <- MC_tsg_pm12_qcc
# MC_tsg_po2_qcc_f <- MC_tsg_po2_qcc
# MC_tsg_po12_qcc_f <- MC_tsg_po12_qcc

MC_tsg_pm2_qcc_f$inx <- factor(MC_tsg_pm2_qcc$inx, levels = unique(MC_tsg_pm2_qcc$inx))
MC_tsg_pm12_qcc_f$inx <- factor(MC_tsg_pm12_qcc$inx, levels = unique(MC_tsg_pm12_qcc$inx))
# MC_tsg_po2_qcc_f$inx <- factor(MC_tsg_po2_qcc$inx, levels = unique(MC_tsg_po2_qcc$inx))
# MC_tsg_po12_qcc_f$inx <- factor(MC_tsg_po12_qcc$inx, levels = unique(MC_tsg_po12_qcc$inx))


# valid_times <- MC_tsg_pm2_qcc_f %>%
#   group_by(time) %>%
#   summarise(all_lim_1 = all(lim == 1)) %>%
#   filter(all_lim_1) %>%
#   pull(time)

valid_times <- MC_tsg_pm2_qcc_f %>%
  group_by(time) %>%
  summarise(proportion_lim_1 = mean(lim == 1)) %>%
  filter(proportion_lim_1 >= 0.95) %>%
  pull(time)
# 
# ggplot(MC_tsg_pm2_qcc_f, aes(x=time, y=y, group=inx)) +
#   # Lines and points to represent the data
#   # geom_step(aes(x=time, y=ucl), linetype="dashed", alpha=.17) +
#   geom_line(aes(x=time, y=ucl), linetype="solid", color="grey70", alpha=1) +
#   # geom_step(aes(x=time, y=center), linetype="solid", alpha=.17) +
#   geom_point(data = MC_tsg_pm2_qcc_f, alpha=.17, color="blue", shape=20) +
#   # geom_line(aes(x=time, y=y), color="grey70", alpha=.17) +
#   geom_point(data = MC_tsg_pm2_qcc_f %>% filter(lim == 1 & time %in% valid_times), alpha = 1, color = "red", shape = 17)+ 
#   geom_vline(xintercept = c(16, 80, 140, 20)-t_window_size, linetype = "dotted") +
#   scale_x_discrete(name ="time points", limits=minx)+#limits=c(minx)) +
#   theme_bw() +
#   theme_classic(base_size = 18) +
#   theme(axis.title = element_text(size = 10), axis.text.x = element_text(size = 10, angle=90), legend.position = "none", plot.title = element_text(hjust = 0.5, size=12, face='bold'), plot.subtitle = element_text(hjust = 0.5)) +
#   ggtitle("Control Chart MASE(2)")


pm2 <-  ggplot(MC_tsg_pm2_qcc_f, aes(x=time, y=y, group=inx)) +
  # Lines and points to represent the data
  # geom_step(aes(x=time, y=ucl), linetype="dashed", alpha=.17) +
  geom_line(aes(x=time, y=ucl), linetype="solid", color="grey70", alpha=1) +
  # geom_step(aes(x=time, y=center), linetype="solid", alpha=.17) +
  geom_point(data = MC_tsg_pm12_qcc_f %>% filter(lim != 1), alpha=.17, color="blue", shape=20) +
  # geom_line(aes(x=time, y=y), color="grey70", alpha=.17) +
  geom_point(data = MC_tsg_pm2_qcc_f %>% filter(lim == 1 & !time %in% valid_times), alpha = .17, color = "blue", shape = 17)+ 
  geom_point(data = MC_tsg_pm2_qcc_f %>% filter(lim == 1 & time %in% valid_times), alpha = 1, color = "red", shape = 17)+ 
  geom_vline(xintercept = c(16, 80, 140, 20)-t_window_size, linetype = "dotted") +
  scale_x_discrete(name ="time points", limits=minx)+#limits=c(minx)) +
  theme_bw() +
  theme_classic(base_size = 18) +
  theme(axis.title = element_text(size = 10), axis.text.x = element_text(size = 10, angle=90), legend.position = "none", plot.title = element_text(hjust = 0.5, size=12, face='bold'), plot.subtitle = element_text(hjust = 0.5)) +
  ggtitle("Control Chart MASE(2)")

png(paste0(output_dir, paste0("/MC_mase2_tsg_g_cc_t_window_size",as.character(t_window_size),"_", "all_layer_thres_", "_",format(Sys.time(), "%Y_%m_%d")),".png"), width = 10, height = 9, units="in", res=400)
pm2
dev.off()

# valid_times <- MC_tsg_pm12_qcc_f %>%
#   group_by(time) %>%
#   summarise(all_lim_1 = all(lim == 1)) %>%
#   filter(all_lim_1) %>%
#   pull(time)

valid_times <- MC_tsg_pm12_qcc_f %>%
  group_by(time) %>%
  summarise(proportion_lim_1 = mean(lim == 1)) %>%
  filter(proportion_lim_1 >= 0.95) %>%
  pull(time)


pm12 <- ggplot(MC_tsg_pm12_qcc_f, aes(x=time, y=y, group=inx))+
  # Lines and points to represent the data
  geom_line(aes(x=time, y=ucl), linetype="solid", color="grey70", alpha=1)+
  # geom_step(aes(x=time, y=center), linetype="solid", alpha=.17) +
  geom_point(data = MC_tsg_pm12_qcc_f %>% filter(lim != 1), alpha=.17, color="blue", shape=20) +
  # geom_line(aes(x=time, y=y), color="grey70", alpha=.17) +
  geom_point(data = MC_tsg_pm2_qcc_f %>% filter(lim == 1 & !time %in% valid_times), alpha = .17, color = "blue", shape = 17)+ 
  geom_point(data = MC_tsg_pm12_qcc_f %>% filter(lim == 1 & time %in% valid_times), alpha=1, color="red", shape=17) +
  geom_vline(xintercept = c(16, 80, 140, 20)-t_window_size, linetype = "dotted") +
  scale_x_discrete(name ="time points", limits=minx)+#limits=c(minx)) +
  theme_bw() +
  theme_classic(base_size = 18) +
  theme(axis.title = element_text(size = 10), axis.text.x = element_text(size = 10, angle=90), legend.position = "none", plot.title = element_text(hjust = 0.5, size=12, face='bold'), plot.subtitle = element_text(hjust = 0.5)) +
  ggtitle("Control Chart MASE(12)")

# png(paste0(output_dir, paste0("/MC_mase12_tsg_g_cc_t_window_size",as.character(t_window_size),"_", "all_layer_thres_", "_",format(Sys.time(), "%Y_%m_%d")),".png"), width = 10, height = 9, units="in", res=400)
# pm12
# dev.off()
# 
# 
# valid_times <- MC_tsg_po2_qcc_f %>%
#   group_by(time) %>%
#   summarise(all_lim_1 = all(lim == 1)) %>%
#   filter(all_lim_1) %>%
#   pull(time)
# 
# 
# po2 <- ggplot(MC_tsg_po2_qcc_f, aes(x=time, y=y, group=inx))+
#   # Lines and points to represent the data
#   geom_line(aes(x=time, y=ucl), linetype="solid", color="grey70", alpha=1)+
#   # geom_step(aes(x=time, y=center), linetype="solid", alpha=.17) +
#   geom_point(data = MC_tsg_po2_qcc_f, alpha=.17, color="blue", shape=20) +
#   # geom_line(aes(x=time, y=y), color="grey70", alpha=.17) +
#   geom_point(data = MC_tsg_po2_qcc_f %>% filter(time %in% valid_times), alpha=1, color="red", shape=17) + 
#   geom_vline(xintercept = c(16, 80, 140, 20)-t_window_size, linetype = "dotted") +
#   scale_x_discrete(name ="time points", limits=minx)+#limits=c(minx)) +
#   theme_bw() +
#   theme_classic(base_size = 18) +
#   theme(axis.title = element_text(size = 10), axis.text.x = element_text(size = 10, angle=90), legend.position = "none", plot.title = element_text(hjust = 0.5, size=12, face='bold'), plot.subtitle = element_text(hjust = 0.5)) +
#   ggtitle("Control Chart OMNI(2)")
# png(paste0(output_dir, paste0("/MC_omni2_tsg_g_cc_t_window_size",as.character(t_window_size),"_", "all_layer_thres_", "_",format(Sys.time(), "%Y_%m_%d")),".png"), width = 10, height = 9, units="in", res=400)
# po2
# dev.off()
# 
# valid_times <- MC_tsg_po12_qcc_f %>%
#   group_by(time) %>%
#   summarise(all_lim_1 = all(lim == 1)) %>%
#   filter(all_lim_1) %>%
#   pull(time)
# 
# 
# po12 <- ggplot(MC_tsg_po12_qcc_f, aes(x=time, y=y, group=inx))+
#   # Lines and points to represent the data
#   geom_line(aes(x=time, y=ucl), linetype="solid", color="grey70", alpha=1)+
#   # geom_step(aes(x=time, y=center), linetype="solid", alpha=.17) +
#   geom_point(data = MC_tsg_po12_qcc_f, alpha=.17, color="blue", shape=20) +
#   # geom_line(aes(x=time, y=y), color="grey70", alpha=.17) +
#   geom_point(data = MC_tsg_po12_qcc_f %>% filter(time %in% valid_times), alpha=1, color="red", shape=17) + 
#   geom_vline(xintercept = c(16, 80, 140, 20)-t_window_size, linetype = "dotted") +
#   scale_x_discrete(name ="time points", limits=minx)+#limits=c(minx)) +
#   theme_bw() +
#   theme_classic(base_size = 18) +
#   theme(axis.title = element_text(size = 10), axis.text.x = element_text(size = 10, angle=90), legend.position = "none", plot.title = element_text(hjust = 0.5, size=12, face='bold'), plot.subtitle = element_text(hjust = 0.5)) +
#   ggtitle("Control Chart OMNI(12)")
# png(paste0(output_dir, paste0("/MC_omni12_tsg_g_cc_t_window_size",as.character(t_window_size),"_", "all_layer_thres_", "_",format(Sys.time(), "%Y_%m_%d")),".png"), width = 10, height = 9, units="in", res=400)
# po12
# dev.off()
# 
# png(paste0(output_dir, paste0("/MC_tsg_pts_g_cc_t_window_size",as.character(t_window_size),"_", "all_layer_thres_", "_",format(Sys.time(), "%Y_%m_%d")),".png"), width = 10, height = 9, units="in", res=400)
# grid.arrange(pm2, pm12, po2, po12, nrow=2, heights=c(2,2))
# dev.off()
# 
# png(paste0(output_dir, paste0("/MC_tsg_pts_pm_g_cc_t_window_size",as.character(t_window_size),"_", "all_layer_thres_", "_",format(Sys.time(), "%Y_%m_%d")),".png"), width = 10, height = 9, units="in", res=400)
# grid.arrange(pm2, pm12, nrow=1, heights=c(2))
# dev.off()

# valid_times <- MC_tsg_pscan_qcc %>%
#   group_by(time) %>%
#   summarise(all_lim_1 = all(lim == 1)) %>%
#   filter(all_lim_1) %>%
#   pull(time)
# 
# 
# ps <- ggplot(MC_tsg_pscan_qcc, aes(x=time, y=y)) +
#   # Lines and points to represent the data
#   geom_step(aes(x=time, y=ucl), linetype="solid", alpha=.17) +
#   # geom_step(aes(x=time, y=center), linetype="solid", alpha=.17) +
#   # geom_point(data = MC_tsg_pm2_qcc, alpha=.17, color="blue", shape=20) +
#   # geom_line(aes(x=time, y=y), color="grey70", alpha=.17) +
#   geom_point(data = MC_tsg_pm2_qcc %>% filter(time %in% valid_times), alpha=1, color="red", shape=17) + 
#   geom_vline(xintercept = c(16, 80, 140, 20)-t_window_size, linetype = "dotted") +
#   scale_x_discrete(name ="time points", limits=minx)+#limits=c(minx)) +
#   theme_bw() +
#   theme_classic(base_size = 18) +
#   theme(axis.title = element_text(size = 10), axis.text.x = element_text(size = 10, angle=90), legend.position = "none", plot.title = element_text(hjust = 0.5, size=12, face='bold'), plot.subtitle = element_text(hjust = 0.5)) +
#   ggtitle("Control Chart OMNI(12)")
# png(paste0(output_dir, paste0("/MC_omni12_tsg_g_cc_t_window_size",as.character(t_window_size),"_", "all_layer_thres_", "_",format(Sys.time(), "%Y_%m_%d")),".png"), width = 10, height = 9, units="in", res=400)
# ps
# dev.off()

c1m12v_list <- list()
anomalous_v_list_1 <- list()
anomalous_v_list_2 <- list()
anomalous_v_list_3 <- list()
anomalous_v_list_4 <- list()
for (i in seq(n_layers)) {# n_layers
  # check simple graph
  unlist(lapply(df.tsg$tg[[i]], function(x) is.simple(x)))
  
  # check undirected graph
  unlist(lapply(df.tsg$tg[[i]], function(x) is.directed(x)))
  
  A_list <- lapply(df.tsg$tg[[i]], function(g) as_adjacency_matrix(as.undirected(g), attr = "weight", sparse = FALSE))
  gip <- lapply(df.tsg$tg[[i]], function(g) as.undirected(g))
  
  aip <- lapply(gip, function(g) as_adjacency_matrix(g, attr = "weight", sparse = FALSE))
  tmax <- length(A_list)
  
  # WBS
  # data_mat_list <- lapply(aip, function(x) x[lower.tri(x, diag = FALSE)])
  # data_mat <- do.call(cbind, data_mat_list)
  # M = 10 # number of random intervals for WBS
  # d = 4 # parameter for scaled PCA algorithm
  # delta = 2
  # intervals = WBS.intervals(M = M, lower = 1, upper = ncol(data_mat))
  # WBS_result = WBS.nonpar.RDPG(data_mat, lowerdiag = TRUE, d,
  #                              Alpha = intervals$Alpha, Beta = intervals$Beta, delta)
  # cpt_hat = tuneBSnonparRDPG(WBS_result, data_mat, lowerdiag = TRUE, d)
  
  # Generate all pairs using the running function
  # all_pairs <- lapply(seq(tmax), function(x) as.character(x))
  
  # Define a function to check if a pair should be excluded
  should_exclude <- function(pair) {
    # Extract the numbers from the pair string
    nums <- as.numeric(pair)
    # Check if the pair matches the pattern x2:x3, x3:x4, ..., x9:x10
    return(nums[1] %% 10 != 1)
  }
  
  # Replace the pairs that match the exclusion pattern with an empty string
  t_list <- seq(tmax)
  minx_all <- rep("", tmax-1)
  minx_all[16] <- "15:16"
  minx_all[20] <- "19:20"
  minx_all[80] <- "79:80"
  minx_all[140] <- "139:140"
  
  # d.max <- "sqrt"
  attr_weight <- "weight"
  # Wrapper function to retrieve the latent position list based on the given graphs and attributes
  latpos.list.wrapped <- realdata_latpos_list_wrapper(gip, dASE=NA, d.max="sqrt", approx=FALSE, elbow_graph=2, attr_weight="weight")
  
  
  # Iterate over different time window sizes for analysis
  for (t_window_size in c(12)) {#seq(10, 12)) {
    
    minx <- minx_all[(t_window_size+1):(tmax-1)]
    #c1m2v <- get_vertex_qcc(gip, t_list, middle.max.inx=NULL, latpos.list=latpos.list.wrapped, d.max=d.max, elbow_graph=elbow_graph, realdata_wrapper=realdata_doMase_wrapper3,embed_span=2,xlab="time points", title="MASE(2)", minx=minx,t_window_size=t_window_size, return_plot=FALSE)
    c1m12v <- get_vertex_qcc(gip, t_list, middle.max.inx=NULL, latpos.list=latpos.list.wrapped, d.max=d.max, elbow_graph=elbow_graph, realdata_wrapper=realdata_doMase_wrapper3,embed_span=12,xlab="time points", title="MASE(2)", minx=minx,t_window_size=t_window_size, return_plot=FALSE)
    
    
  }
  c1m12v_list[[length(c1m12v_list) + 1]] <- c1m12v
  anomalous_v_list_1[[length(anomalous_v_list_1) + 1]] <- c1m12v[[16-t_window_size ]]$violations$beyond.limits
  anomalous_v_list_2[[length(anomalous_v_list_2) + 1]] <- c1m12v[[20-t_window_size ]]$violations$beyond.limits
  anomalous_v_list_3[[length(anomalous_v_list_3) + 1]] <- c1m12v[[80-t_window_size ]]$violations$beyond.limits
  anomalous_v_list_4[[length(anomalous_v_list_4) + 1]] <- c1m12v[[140-t_window_size ]]$violations$beyond.limits
  
}


intersection_result <- Reduce(intersect, anomalous_v_list_2)

# Combine all lists into a single vector
all_points_1 <- unlist(anomalous_v_list_1)

# Count occurrences of each point
point_counts_1 <- table(all_points_1)

# Number of lists
num_lists_1 <- length(anomalous_v_list_1)

# Filter points that appear in at least 50% of the lists
points_in_95_percent_1 <- names(point_counts_1[point_counts_1 >= (num_lists_1 * 0.95)])

# Convert the names back to their original type if necessary
intersection_result_1 <- as.numeric(points_in_95_percent_1)

# Display the result
print(intersection_result_1)

# Combine all lists into a single vector
all_points_3 <- unlist(anomalous_v_list_3)

# Count occurrences of each point
point_counts_3 <- table(all_points_3)

# Number of lists
num_lists_3 <- length(anomalous_v_list_3)

# Filter points that appear in at least 50% of the lists
points_in_95_percent_3 <- names(point_counts_3[point_counts_3 >= (num_lists_3 * 0.95)])

# Convert the names back to their original type if necessary
intersection_result_3 <- as.numeric(points_in_95_percent_3)

# Display the result
print(intersection_result_3)

# Combine all lists into a single vector
all_points_4 <- unlist(anomalous_v_list_4)

# Count occurrences of each point
point_counts_4 <- table(all_points_4)

# Number of lists
num_lists_4 <- length(anomalous_v_list_4)

# Filter points that appear in at least 50% of the lists
points_in_95_percent_4 <- names(point_counts_4[point_counts_4 >= (num_lists_4 * 0.95)])

# Convert the names back to their original type if necessary
intersection_result_4 <- as.numeric(points_in_95_percent_4)

# Display the result
print(intersection_result_4)



# Combine all lists into a single vector
all_points <- unlist(anomalous_v_list_2)

# Count occurrences of each point
point_counts <- table(all_points)

# Number of lists
num_lists <- length(anomalous_v_list_2)

# Filter points that appear in at least 50% of the lists
points_in_95_percent <- names(point_counts[point_counts >= (num_lists * 0.95)])

# Convert the names back to their original type if necessary
intersection_result_2 <- as.numeric(points_in_95_percent)

# Display the result
print(intersection_result_2)

# Combine all lists into a single vector
all_points <- unlist(anomalous_v_list_2)

# Count occurrences of each point
point_counts <- table(all_points)

# Number of lists
num_lists <- length(anomalous_v_list_2)

# Filter points that appear in at least 50% of the lists
points_in_95_percent <- names(point_counts[point_counts >= (num_lists * 0.95)])

# Convert the names back to their original type if necessary
intersection_result_2 <- as.numeric(points_in_95_percent)

# Display the result
print(intersection_result_2)

intersection_result <- Reduce(intersect, list(intersection_result_1, intersection_result_3, intersection_result_4))

# anomalous_v_list <- list()
# c1m12v <- get_vertex_qcc(gip, t_list, middle.max.inx=NULL, latpos.list=latpos.list.wrapped, d.max=d.max, elbow_graph=elbow_graph, realdata_wrapper=realdata_doMase_wrapper3,embed_span=2,xlab="time points", title="MASE(2)", minx=minx,t_window_size=t_window_size, return_plot=FALSE)


d.max <- "sqrt"
elbow_graph <- 2
# Extracting MASE(2) vertex quality control chart data and creating a zoomed plot of it
# c1m12v <- get_vertex_qcc(gip, t_list, middle.max.inx=NULL, latpos.list=latpos.list.wrapped, d.max=d.max, elbow_graph=elbow_graph, realdata_wrapper=realdata_doMase_wrapper3,embed_span=2,xlab="time points", title="MASE(2)", minx=minx,t_window_size=t_window_size, return_plot=TRUE)
MC_tsg_cc_vertex_m12 <- plot.qcc.vertex(c1m12v, add.stats = FALSE, chart.all = TRUE, s=s,minx = minx,
                                         label.limits = c("LCL ", "UCL"), title=paste("Control Chart","MASE(12)"), xlab="time points",m2=seq(vcount(gip[[1]])), ylab=TeX("$\\log(y_{i}^{(t)})$"),
                                         axes.las = 0, digits = getOption("digits"),artpoints = intersection_result_2,
                                         restore.par = FALSE)
MC_tsg_cc_vertex_m12 <- MC_tsg_cc_vertex_m12+ scale_y_log10()

intersection_result <- Reduce(intersect, cpt_hats)

install.packages("clue")  # for solving the assignment problem (Hungarian algorithm)
library(clue)

# Define the true and estimated change points
true_change_points <- c(16, 20, 80, 140)
# Function to calculate the MSE for a single set of estimated change points
calculate_mse <- function(true_points, estimated_points) {
  num_true <- length(true_points)
  num_estimated <- length(estimated_points)
  
  if (num_estimated < num_true) {
    # Pad the estimated points to make the matrix square
    padding_length <- num_true - num_estimated
    padded_estimated_points <- c(estimated_points, rep(NA, padding_length))
    distance_matrix <- outer(true_points, padded_estimated_points, FUN = function(t, e) ifelse(is.na(e), 1e6, (t - e)^2))
  } else {
    distance_matrix <- outer(true_points, estimated_points, FUN = function(t, e) (t - e)^2)
  }
  
  # Solve the assignment problem to find the optimal pairing
  assignment <- solve_LSAP(distance_matrix)
  
  # Extract valid assignments
  valid_assignments <- assignment[assignment <= num_estimated]
  
  # Create pairs
  paired_true <- true_points[seq_along(valid_assignments)]
  paired_estimated <- estimated_points[valid_assignments]
  
  # Calculate squared errors for paired points
  squared_errors <- (paired_true - paired_estimated)^2
  
  # Handle any unmatched true change points
  if (num_true > num_estimated) {
    # If there are more true change points than estimated points
    # add penalty for unmatched points
    unmatched_count <- num_true - num_estimated
    penalty <- max(squared_errors, na.rm = TRUE)  # Choose a penalty, here max squared error
    squared_errors <- c(squared_errors, rep(penalty, unmatched_count))
  }
  
  # Calculate the Mean Squared Error (MSE)
  mse <- mean(squared_errors, na.rm = TRUE)
  
  return(mse)
}
# Calculate the MSE for each set of estimated change points and average them
mse_values <- sapply(cpt_hats, function(estimated_change_points) {
  calculate_mse(true_change_points, estimated_change_points)
})

# Calculate the average MSE
average_mse <- mean(mse_values)

# Output the results
print("MSE Values for Each Detection:")
print(mse_values)
print("Average Mean Squared Error (MSE):")
print(average_mse)


# Function to calculate the Hausdorff distance for a single set of estimated change points
calculate_hausdorff <- function(true_points, estimated_points) {
  # Calculate the Hausdorff distance
  distance <- hausdorff_dist(true_points, estimated_points)
  return(distance)
}

# Calculate the Hausdorff distance for each set of estimated change points and average them
hausdorff_values <- sapply(cpt_hats, function(estimated_change_points) {
  calculate_hausdorff(true_change_points, estimated_change_points)
})

# Calculate the average Hausdorff distance
average_hausdorff <- mean(hausdorff_values)

# Output the results
print("Hausdorff Distance Values for Each Detection:")
print(hausdorff_values)
print("Average Hausdorff Distance:")
print(average_hausdorff)

# Create a data frame for plotting
hausdorff_data <- data.frame(HausdorffDistance = hausdorff_values)

# Plot the kernel density with a vertical dashed line for the mean
ggplot(hausdorff_data, aes(x = HausdorffDistance)) +
  geom_density(fill = "blue", alpha = 0.3) +
  geom_vline(aes(xintercept = average_hausdorff), color = "red", linetype = "dashed", size = 1) +
  labs(title = "Kernel Density Plot of Hausdorff Distances",
       x = "Hausdorff Distance",
       y = "Density") +
  theme_minimal()

cpt_pm2 <- MC_tsg_pm2_qcc_f %>%
  filter(lim == 1) %>%
  split(.$inx) %>%
  lapply(function(sub_df) sub_df$time+t_window_size)

cpt_pm12 <- MC_tsg_pm12_qcc_f %>%
  filter(lim == 1) %>%
  split(.$inx) %>%
  lapply(function(sub_df) sub_df$time+t_window_size)

cpt_pm2 <- append(cpt_pm2[1:96],cpt_pm2[98])
cpt_pm12 <- append(cpt_pm12[1:96],cpt_pm12[98])#cpt_pm12[1:96]


calculate_hausdorff <- function(estimated_points, true_points) {
  distance <- hausdorff_dist(estimated_points, true_points)
  return(distance)
}

# Calculate the Hausdorff distance for each set of estimated change points
hausdorff_cpt_hats <- sapply(cpt_hats, function(estimated_change_points) {
  calculate_hausdorff(estimated_change_points, true_change_points)
})

hausdorff_cpt_pm2 <- sapply(cpt_pm2, function(estimated_change_points) {
  calculate_hausdorff(estimated_change_points, true_change_points)
})

hausdorff_cpt_pm12 <- sapply(cpt_pm12, function(estimated_change_points) {
  calculate_hausdorff(estimated_change_points, true_change_points)
})

# Assuming true_change_points, cpt_hats, cpt_pm2, and cpt_pm12 are defined

# True number of change points
true_num_cpts <- length(true_change_points)

# Calculate the absolute errors
num_cpt_hats_errors <- sapply(cpt_hats, function(estimated_change_points) {
  abs(length(estimated_change_points) - true_num_cpts)
})

num_cpt_pm2_errors <- sapply(cpt_pm2, function(estimated_change_points) {
  abs(length(estimated_change_points) - true_num_cpts)
})

num_cpt_pm12_errors <- sapply(cpt_pm12, function(estimated_change_points) {
  abs(length(estimated_change_points) - true_num_cpts)
})

# Combine the absolute errors into a single data frame
combined_num_cpts_errors <- data.frame(
  Method = rep(c("RDPG-CPD", "MASE(2)", "MASE(12)"), each = length(num_cpt_hats_errors)),
  AbsoluteError = c(num_cpt_hats_errors, num_cpt_pm2_errors, num_cpt_pm12_errors)
)

# Calculate the mean absolute errors
mean_rdpg_error <- mean(num_cpt_hats_errors)
mean_mase2_error <- mean(num_cpt_pm2_errors)
mean_mase12_error <- mean(num_cpt_pm12_errors)

# Plot the kernel density with vertical dashed lines for the means
density_plot_num_cpts_errors <- ggplot(combined_num_cpts_errors, aes(x = AbsoluteError, fill = Method, color = Method)) +
  geom_density(alpha = 0.3) +
  geom_vline(aes(xintercept = mean_rdpg_error, color = "RDPG-CPD"), linetype = "dashed", size = 1, show.legend = TRUE) +
  geom_vline(aes(xintercept = mean_mase2_error, color = "MASE(2)"), linetype = "dashed", size = 1, show.legend = TRUE) +
  geom_vline(aes(xintercept = mean_mase12_error, color = "MASE(12)"), linetype = "dashed", size = 1, show.legend = TRUE) +
  labs(
    x = "Absolute Error of Number of Change Points",
    y = "Density",
    fill = "Method",
    color = "Method") +
  scale_fill_manual(values = c("RDPG-CPD" = "red", "MASE(2)" = "blue", "MASE(12)" = "green")) +
  scale_color_manual(values = c("RDPG-CPD" = "red", "MASE(2)" = "blue", "MASE(12)" = "green")) +
  theme_bw() +
  theme(axis.title = element_text(size = 10), axis.text.x = element_text(size = 10, angle=90), legend.position = c(0.8, 0.8), plot.title = element_text(hjust = 0.5, size=12, face='bold'), plot.subtitle = element_text(hjust = 0.5)) +
  ggtitle("Kernel Density Plot of Absolute Error of Number of Change Points") +
  theme_classic(base_size = 18)

# Print the plot
print(density_plot_num_cpts_errors)

png(paste0(output_dir, paste0("/MC_num_cpts_density",as.character(t_window_size),"_", "all_layer_thres_", "_",format(Sys.time(), "%Y_%m_%d")),".png"), width = 10, height = 9, units="in", res=400)
density_plot_num_cpts_errors
dev.off()

# Combine the Hausdorff distances into a single data frame
combined_hausdorff <- data.frame(
  Method = rep(c("RDPG-CPD", "MASE(2)", "MASE(12)"), each = length(hausdorff_cpt_hats)),
  HausdorffDistance = c(hausdorff_cpt_hats, hausdorff_cpt_pm2, hausdorff_cpt_pm12)
)

# Calculate the mean Hausdorff distances
mean_rdpg <- mean(hausdorff_cpt_hats)
mean_mase2 <- mean(hausdorff_cpt_pm2)
mean_mase12 <- mean(hausdorff_cpt_pm12)

# Plot the kernel density with vertical dashed lines for the means
density_plot_0 <- ggplot(combined_hausdorff, aes(x = HausdorffDistance, fill = Method, color = Method)) +
  geom_density(alpha = 0.3) +
  geom_vline(aes(xintercept = mean_rdpg, color = "RDPG-CPD"), linetype = "dashed", size = 1, show.legend = TRUE) +
  geom_vline(aes(xintercept = mean_mase2, color = "MASE(2)"), linetype = "dashed", size = 1, show.legend = TRUE) +
  geom_vline(aes(xintercept = mean_mase12, color = "MASE(12)"), linetype = "dashed", size = 1, show.legend = TRUE) +
  labs(
    x = "Hausdorff Distance",
    y = "Density",
    fill = "Method",
    color = "Method") +
  scale_fill_manual(values = c("RDPG-CPD" = "red", "MASE(2)" = "blue", "MASE(12)" = "green")) +
  scale_color_manual(values = c("RDPG-CPD" = "red", "MASE(2)" = "blue", "MASE(12)" = "green")) +
  theme_bw() +
  theme(axis.title = element_text(size = 10), axis.text.x = element_text(size = 10, angle=90), legend.position = c(0.8, 0.8), plot.title = element_text(hjust = 0.5, size=12, face='bold'), plot.subtitle = element_text(hjust = 0.5)) +
  ggtitle("Kernel Density Plot of Hausdorff Distances")
theme_classic(base_size = 18)


png(paste0(output_dir, paste0("/MC_tsg_pts_pm12_g_v_cc_density_t_window_size",as.character(t_window_size),"_", "all_layer_thres_", "_",format(Sys.time(), "%Y_%m_%d")),".png"), width = 10, height = 9, units="in", res=400)
grid.arrange(pm2, pm12, MC_tsg_cc_vertex_m12, density_plot_0, nrow=2, heights=c(2,2))
dev.off()



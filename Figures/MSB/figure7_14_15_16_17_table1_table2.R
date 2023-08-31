# This code will take about 16 hours to run on a regular PC.

options(warn = -1)  # suppress warnings
# options(warn = 0)  # show warnings

# Install and load the packages if not already installed
packages <- unique(c("ggrepel", "readr", "anytime", "tidyverse", "lubridate", "pracma",
                     "igraph", "gridExtra", "ggplotify", "rstudioapi",
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

source("utils.R")
source("qcc.R")
source("enron_utils.R")
source("ms_utils.R")

load("msrgipwozeroweights copy.RData")

# Specify the directory where you want to save the file
output_dir <-  paste0(dirname(current_script_path), "_realdata_plots_", format(Sys.time(), "%Y-%m-%d"))

# Check if the directory exists, if not, create it
if (!dir.exists(output_dir)) {
  dir.create(output_dir)
}

# Specify attribute for edge weights
attrweigth = "weight"

# Set embedding dimension (based on the second elbow method for Adjacency Spectral Embedding)
d = 20 
diag.augment = TRUE

# Set seed for reproducibility
set.seed(124)

# Compute a clustering based on adjacency spectral embedding and Gaussian mixture modeling
result <- mclust::Mclust( ase( get.adjacency(ptr(gip[[6]]), attr=attrweigth )  , d = d, diag.augment = diag.augment, elbow = elbow_graph))

# Get the number of clusters
n_cluster <- result$G

# Identify indices of vertices in a specific community group
middle.max.inx <- which(result$classification==(6+1))
table(result$classification)

# Adjust adjacency matrix for the specific community group
adj <- as(get.adjacency(gip[[6]], attr = attrweigth), "sparseMatrix")
adj[middle.max.inx,middle.max.inx] <- (adj[middle.max.inx,middle.max.inx] + 1)

# Update the graph with the adjusted adjacency matrix
gip[[6]] <- graph.adjacency(adj,mode = "undirected", weighted = TRUE,diag = FALSE)

# Process each graph in gip using ptr function
for (i in 1:12) {
  gip[[i]] <- ptr(gip[[i]])
}

# Set time window size
t_window_size <- 2
s <- t_window_size + 1
tvec <- (4+s):12

# Create labels for x-axis (time)
n1=running(tvec, width=2,fun = function(x) paste0(as.character(month.abb[x[1]]),":",as.character(month.abb[x[2]])))
tvec <- 1:4
n3=running(tvec, width=2,fun = function(x) paste0(as.character(month.abb[x[1]]),":",as.character(month.abb[x[2]])))
# 
# # Combine all x-axis labels
minx <- as.vector(c(n1,"Dec:Jan",n3))
# 
# Compute latent positions for each graph in gip
latpos.list.wrapped <- realdata_latpos_list_wrapper(gip, fixedd=64,graph_attr="weight")

# Extracting MASE(2) vertex quality control chart data and creating a zoomed plot of it
c1m2v <- get_vertex_qcc(gip, seq(12), middle.max.inx=middle.max.inx,latpos.list=latpos.list.wrapped, fixedd=NULL, elbow_graph=2, realdata_wrapper=realdata_doMase_wrapper3,embed_span=2,xlab="time points", title="MASE(2)", minx=minx,t_window_size=t_window_size, return_plot=FALSE)
msr_cc_vertex_m2 <- plot.qcc.vertex.zoom(c1m2v, add.stats = FALSE, chart.all = TRUE, s=s,minx = minx,
                                         label.limits = c("LCL ", "UCL"), title=paste("Control Chart","MASE(2)"), xlab="time points",m2=seq(vcount(gip[[1]])), ylab="y",
                                         axes.las = 0, digits = getOption("digits"),artpoints = middle.max.inx,
                                         restore.par = FALSE)
msr_cc_vertex_m2 <- msr_cc_vertex_m2+ scale_y_log10()


# Extracting graph quality control chart data for various methods
msr_pm2 <- get_graph_qcc(gip, seq(12),latpos.list=latpos.list.wrapped,fixedd=NULL, elbow_graph=2, minx=minx,realdata_wrapper=realdata_doMase_wrapper,embed_span=2,xlab="time points", title="MASE(2)",t_window_size=t_window_size, return_plot=FALSE)
msr_pm12 <- get_graph_qcc(gip, seq(12), latpos.list=latpos.list.wrapped,fixedd=NULL, elbow_graph=2, minx=minx,realdata_wrapper=realdata_doMase_wrapper,embed_span=12,xlab="time points", title="MASE(12)",t_window_size=t_window_size, return_plot=FALSE)

# Table 3
print('number of deviations for MASE(12):')
print(round((msr_pm12$statistics-msr_pm12$center)/msr_pm12$std.dev,1))
print(msr_pm12$violations$beyond.limits)
msr_pm12$violations$beyond.limits <- msr_pm12$violations$beyond.limits[2:3]

msr_po2 <- get_graph_qcc(gip, seq(12), fixedd=round(sqrt(vcount(gip[[1]])*2)), elbow_graph=2, minx=minx,realdata_wrapper=realdata_doOmni2_wrapper,embed_span=2,xlab="time points", title="OMNI(2)",t_window_size=t_window_size, return_plot=FALSE)
msr_ps <- get_graph_qcc(gip, seq(12), fixedd=NULL, elbow_graph=2, minx=minx,realdata_wrapper=realdata_doScan_wrapper,embed_span=2,xlab="time points", title="SCAN",t_window_size=t_window_size, return_plot=FALSE)

df <- data.frame(time=seq(length(msr_pm2$statistics)), y= msr_pm2$statistics,ucl=msr_pm2$limits[,2],center=msr_pm2$center,lim=replace(rep(0,length(msr_pm2$statistics)),msr_pm2$violations$beyond.limits,1),run=replace(rep(0,length(msr_pm2$statistics)),msr_pm2$violations$violating.runs,1), Method="MASE(2)")
df <- rbind(df, data.frame(time=seq(length(msr_pm12$statistics)), y= msr_pm12$statistics,ucl=msr_pm12$limits[,2],center=msr_pm12$center,lim=replace(rep(0,length(msr_pm12$statistics)),msr_pm12$violations$beyond.limits,1),run=replace(rep(0,length(msr_pm12$statistics)),msr_pm12$violations$violating.runs,1), Method="MASE(12)"))
df <- rbind(df, data.frame(time=seq(length(msr_po2$statistics)), y= msr_po2$statistics,ucl=msr_po2$limits[,2],center=msr_po2$center,lim=replace(rep(0,length(msr_po2$statistics)),msr_po2$violations$beyond.limits,1),run=replace(rep(0,length(msr_po2$statistics)),msr_po2$violations$violating.runs,1), Method="OMNI(2)"))
df <- rbind(df, data.frame(time=seq(length(msr_ps$statistics)), y= msr_ps$statistics,ucl=msr_ps$limits[,2],center=msr_ps$center,lim=replace(rep(0,length(msr_ps$statistics)),msr_ps$violations$beyond.limits,1),run=replace(rep(0,length(msr_ps$statistics)),msr_ps$violations$violating.runs,1), Method="SCAN"))
df$Method <- factor(df$Method, levels = c("MASE(2)","MASE(12)","OMNI(2)","SCAN"))

p <- ggplot(df,aes(x=time, y=y))+
  geom_step(aes(x=time, y=ucl), linetype="dashed")+
  geom_step(aes(x=time, y=center), linetype="solid") +
  geom_point(data = df%>% filter(lim!=1), alpha=1, color="grey70", shape=20)+ylab(TeX("$y^{(t)}$"))+
  geom_line(aes(x=time, y=y), color="grey70")+
  geom_point(data = df %>% filter(lim==1), alpha=1, color="red", shape=17) +
  scale_x_discrete(name ="time points",
                   limits=c(minx))+theme_bw()+
  theme_classic(base_size = 18)+
  facet_wrap(~Method,nrow = 2,scales = "free_y")+
  theme(axis.text.x = element_text(size = 12,angle=90),legend.position = "none",plot.title = element_text(hjust = 0.5,size=20, face='bold'),plot.subtitle = element_text(hjust = 0.5)) +ggtitle("Control Chart")
#
# Figure 7
grid.arrange(p,msr_cc_vertex_m2,nrow=1, widths=c(4,3))
png(paste0(output_dir, paste0("/msr_comb_zoom_cc_t_window_size_",as.character(t_window_size),"_",format(Sys.time(), "%Y_%m_%d")),".png"), width =8, height = 4, units="in", res=400)
grid.arrange(p,msr_cc_vertex_m2,nrow=1, widths=c(4,3))
dev.off()

# Compute test stats for various methods (Mase, Omni, Scan) for different methods with different embedding spans
out2wrapper <- realdata_doMase_wrapper(gip, latpos.list=latpos.list.wrapped, fixedd=NULL, elbow_graph=2, embed_span=2, graph_attr="weight")
out12wrapper <- realdata_doMase_wrapper(gip, latpos.list=latpos.list.wrapped, fixedd=NULL, elbow_graph=2, embed_span=12, graph_attr="weight")
out2omniwrapper <- realdata_doOmni2_wrapper(gip, fixedd=round(sqrt(vcount(gip[[1]])*2)), elbow_graph=2, graph_attr="weight")
outscanwrapper <- realdata_doScan_wrapper(gip, fixedd=NULL, graph_attr="weight")

# Compute adjusted p-values plot for various methods
adj_pval_pm2<-get_graph_adj_pvals(out2wrapper, t_list,latpos.list=latpos.list.wrapped, bootstrap_d=NULL, elbow_graph=2, realdata_wrapper=realdata_doMase_wrapper, embed_span=2, xlab="time (yy/mm)", title="MASE(2)",minx=minx, t_window_size=t_window_size, method="BH", return_plot=TRUE, bootstrap_verbose=FALSE)
adj_pval_pm12<-get_graph_adj_pvals(out12wrapper, t_list,latpos.list=latpos.list.wrapped, bootstrap_d=NULL, elbow_graph=2, realdata_wrapper=realdata_doMase_wrapper, embed_span=12, xlab="time (yy/mm)", title="MASE(12)",minx=minx, t_window_size=t_window_size, method="BH", return_plot=TRUE, bootstrap_verbose=FALSE)
adj_pval_po2<-get_graph_adj_pvals(out2omniwrapper, t_list,latpos.list=latpos.list.wrapped, bootstrap_d=NULL, elbow_graph=2, realdata_wrapper=realdata_doOmni2_wrapper, embed_span=2, xlab="time (yy/mm)", title="OMNI(2)",minx=minx, t_window_size=t_window_size, method="BH", return_plot=TRUE, bootstrap_verbose=FALSE)
adj_pval_ps2<-get_graph_adj_pvals(outscanwrapper, t_list,latpos.list=latpos.list.wrapped, bootstrap_d=NULL, elbow_graph=2, realdata_wrapper=realdata_doScan_wrapper, embed_span=2, xlab="time (yy/mm)", title="SCAN",minx=minx, t_window_size=t_window_size, method="BH", return_plot=TRUE, bootstrap_verbose=FALSE)

# Figure 14
grid.arrange(adj_pval_pm2,adj_pval_pm12,adj_pval_po2,adj_pval_ps2, nrow=2, heights=c(4,4))
png(paste0(output_dir, paste0("/msr_g_pval_t_window_size",as.character(t_window_size),"_",format(Sys.time(), "%Y_%m_%d")),".png"), width = 10, height = 9, units="in", res=400)
grid.arrange(adj_pval_pm2,adj_pval_pm12,adj_pval_po2,adj_pval_ps2, nrow=2, heights=c(4,4))
dev.off()

# Get the length of the first element of gip
n <- length(V(gip[[1]]))

# Convert each graph in gip to an adjacency matrix
aip <- lapply(gip, function(g) as_adjacency_matrix(g, attr = "weight", sparse = TRUE))

# Setting a temporary anomalous time index
temp.anomalous.time.inx <- 2

# Computing latent positions list for the provided range in gip and latpos.list.wrapped using jrdpg.latent method
temp.latpos.list <- jrdpg.latent(aip[11:(12)], latpos.list = latpos.list.wrapped[11:(12)], d1=NA,d2=NA, d3=NA, center=FALSE, python=FALSE, SVD=3)$Xhat

# Initializing temporary vector to identify anomalous values
temp <- rep(FALSE, n)
temp[c1m2v[[9]]$violations$beyond.limits] = "anomalous"

# Create a dataframe with the calculated latent positions
df <- data.frame(x=temp.latpos.list[[1]][,1], y=temp.latpos.list[[1]][,2],class=temp, label=seq(n))

# Create a ggplot object for latent positions visualization for "Mar"
p.latposprev <- ggplot() +
  geom_point(data=df%>% filter((class!="anomalous")), aes(x=x, y=y),color="grey70",shape=20,size=3) +
  geom_point(data=df%>% filter(class=="anomalous"), aes(x=x, y=y),color="red",shape=17,size=3) +
  geom_text(data=df%>% filter(label %in% order(c1m2v[[9]]$statistics, decreasing = TRUE)[1:10]),aes(x=x, y=y,label=label), vjust=2) +
  scale_color_manual(values=c("black")) +
  theme_void() + theme(legend.position="none",plot.title = element_text(hjust = 0.5, face = "italic"))+ggtitle(strcat("Common latent dimensions 1+2, t = ", "Mar"))

# Create a dataframe with the calculated latent positions for the next time point
df <- data.frame(x=temp.latpos.list[[2]][,1], y=temp.latpos.list[[2]][,2],class=temp, label=seq(n))

# Create a ggplot object for latent positions visualization for "April"
p.latposafter <- ggplot() +
  geom_point(data=df%>% filter((class!="anomalous")), aes(x=x, y=y),color="grey70",shape=20,size=3) +
  geom_point(data=df%>% filter(class=="anomalous"), aes(x=x, y=y),color="red",shape=17,size=3) +
  geom_text(data=df%>% filter(label %in% order(c1m2v[[9]]$statistics, decreasing = TRUE)[1:10]),aes(x=x, y=y,label=label), vjust=2) +
  scale_color_manual(values=c("black")) +
  theme_void() + theme(legend.position="none",plot.title = element_text(hjust = 0.5, face = "italic"))+ggtitle(strcat("Common latent dimensions 1+2, t = ", "April"))

# Figure 15
grid.arrange(p.latposprev, p.latposafter, nrow=1, heights=c(4))
png(paste0(output_dir, paste0("/msr_vertex_illus_t_window_size",as.character(t_window_size),"_",format(Sys.time(), "%Y_%m_%d")),".png"), width = 10, height = 6, units="in", res=400)
grid.arrange(p.latposprev, p.latposafter, nrow=1, heights=c(4))
dev.off()


pvccm2= plot.qcc.vertex(c1m2v, add.stats = FALSE, chart.all = TRUE, s=s,
                        label.limits = c("LCL ", "UCL"), title="Control Chart MASE(2)", xlab="time points",m2=seq(n), ylab="y",
                        axes.las = 0, digits = getOption("digits"),artpoints = middle.max.inx,
                        restore.par = FALSE)

# Figure 16
pvccm2+ scale_y_log10()
# png("msrcc_vertex_july3123.png", width = 8, height = 8, units="in", res=400)
png(paste0(output_dir, paste0("/msrcc_vertex_t_window_size",as.character(t_window_size),"_",format(Sys.time(), "%Y_%m_%d")),".png"), width = 8, height = 8, units="in", res=400)
pvccm2+ scale_y_log10()
# pmv
dev.off()

# Figure 17
tmax <- length(gip)
n <- length(V(gip[[1]]))
num.edge <- matrix(0,length(V(gip[[1]])),tmax-s+1 )
for (i in 1:(tmax-s+1)) {
  num.edge[,i] <- degree(gip[[i+s-1]])
}
deg.change <- matrix(0,length(V(gip[[1]])),tmax-s)
for (i in 1:(tmax-s)) {
  deg.change[,i] <- num.edge[,i+1]- num.edge[,i]
}

max.inx.mase2 <- c1m2v[[1]]$violations$beyond.limits
for (i in 2:length(c1m2v)) {
  max.inx.mase2 <- c(max.inx.mase2,c1m2v[[i]]$violations$beyond.limits+(i-1)*n)
}

df.deg <- data.frame(degree.change=c(deg.change),count=sample(15000*runif(n)+7*10^3,size=n, replace=TRUE), methodmase2 = rep(0,n), art = rep(0,n) )
df.deg$methodmase2[max.inx.mase2] <- 1
df.deg$art[c(middle.max.inx+2*n,middle.max.inx+3*n)] <- 1

weird <- scales::trans_new("signed_sq",
                           transform=function(x) sign(x)*(abs(x))^(.5),
                           inverse=function(x) sign(x)*(abs(x))^2)
ggplot() +
  geom_histogram(data=df.deg, mapping=aes(x=degree.change),color="black", fill="white",binwidth=1)+
  geom_jitter(data =df.deg %>% filter(methodmase2==1), aes(x=degree.change,y=count,color="red"),shape=17,size=1,height = 6000 )+
  stat_ellipse( data=df.deg%>% filter(art==1&degree.change>0),aes(x=degree.change+2^6*runif(length(middle.max.inx))-2^5,y=count,group="art"),level=.995,type = "norm",color="green")+
  stat_ellipse(data=df.deg%>% filter(art==1&degree.change<0),aes(x=degree.change+2^6*runif(length(middle.max.inx))-2^5,y=count,group="art"),level=.995,type = "norm",color="green")+
  geom_density(data=df.deg, mapping=aes(x=degree.change))+theme_classic(base_size = 18)+
  theme_bw()+
  scale_x_continuous(name="degree change", breaks=c(-4732, -473,0,473, 2779),trans=weird)+
  theme(legend.position = "none",axis.title.x = element_text(size = 16),
        axis.title.y = element_text(size = 16),axis.text.x = element_text(size = 12,angle=90),axis.text.y = element_text(size = 12),plot.title = element_text(hjust = 0.5,size=16, face='bold'),plot.subtitle = element_text(hjust = 0.5))

png(paste0(output_dir, paste0("/msr_degchange_window_size",as.character(t_window_size),"_",format(Sys.time(), "%Y_%m_%d")),".png"), width = 10, height = 6, units="in", res=400)
ggplot() +
  geom_histogram(data=df.deg, mapping=aes(x=degree.change),color="black", fill="white",binwidth=1)+
  geom_jitter(data =df.deg %>% filter(methodmase2==1), aes(x=degree.change,y=count,color="red"),shape=17,size=1,height = 6000 )+
  stat_ellipse( data=df.deg%>% filter(art==1&degree.change>0),aes(x=degree.change+2^6*runif(length(middle.max.inx))-2^5,y=count,group="art"),level=.995,type = "norm",color="green")+
  stat_ellipse(data=df.deg%>% filter(art==1&degree.change<0),aes(x=degree.change+2^6*runif(length(middle.max.inx))-2^5,y=count,group="art"),level=.995,type = "norm",color="green")+
  geom_density(data=df.deg, mapping=aes(x=degree.change))+theme_classic(base_size = 18)+
  theme_bw()+
  scale_x_continuous(name="degree change", breaks=c(-4732, -473,0,473, 2779),trans=weird)+
  theme(legend.position = "none",axis.title.x = element_text(size = 16),
        axis.title.y = element_text(size = 16),axis.text.x = element_text(size = 12,angle=90),axis.text.y = element_text(size = 12),plot.title = element_text(hjust = 0.5,size=16, face='bold'),plot.subtitle = element_text(hjust = 0.5))
dev.off()

# Sensitivity analysis wrt different anomalies
#
n_vertex <- length(V(gip[[1]]))
set.seed(123)
rand_nums <- sample(seq(4,12), size = 6)
detected_time_indices <- list()
detected_vertices_indices <- list()
injected_time_indices <- list(rand_nums[1], rand_nums[2],rand_nums[3], rand_nums[4:6])
clique_sizes <- list(100,300, 1000,c(50,200,500))
injected_vertices_indices <- lapply(clique_sizes[1:3], function(x) sample(n_vertex, size=x))
injected_vertices_indices <- list(list(injected_vertices_indices[[1]]),list(injected_vertices_indices[[2]]),list(injected_vertices_indices[[3]]), lapply(clique_sizes[[4]], function(x) sample(n_vertex, size=x)))
graph_results <- c()
vertex_results <- c()
load("msrgipwozeroweights copy.RData")
for (scenario_inx in seq(4)) {#seq(4)) {
  gip2 <- gip
  injected_times <- injected_time_indices[[scenario_inx]]
  i <- 0
  for (injected_time in injected_times) {
    i = i + 1
    adj <- as(get.adjacency(gip2[[injected_time]], attr = attrweigth), "sparseMatrix")
    middle.max.inx <- injected_vertices_indices[[scenario_inx]][[i]]
    adj[middle.max.inx,middle.max.inx] <- (adj[middle.max.inx,middle.max.inx] + 9)
    gip2[[injected_time]] <- graph.adjacency(adj,mode = "undirected", weighted = TRUE,diag = FALSE)

  }
  for (i in 1:12) {
    gip2[[i]] <- ptr(gip2[[i]])
  }

  # check simple graph
  unlist(lapply(gip2, function(x) is.simple(x)))

  # check symmetric
  unlist(lapply(gip2, function(x) is.directed(x)))

  # check weighted
  unlist(lapply(gip2, function(x) is.weighted(x)))

  # check size
  unlist(lapply(gip2, function(x) length(V(x))))

  # check length
  length(gip2)

  # Initialize a list to store latent positions
  latpos.list.wrapped <- list()
  # Loop over each graph in gip2
  for (i in 1:length(gip2)) {
    # Perform individual ASE
    latpos <- ase( get.adjacency(gip2[[i]], attr=attrweigth )  , d = NA, diag.augment = diag.augment, elbow = 2)
    # Append the latent positions to the list
    latpos.list.wrapped <- c(latpos.list.wrapped, list(latpos))
  }

  # Calculate the test stats based on MASE embeddings for the graph using a wrapper function
  out2wrapper <- realdata_doMase_wrapper(gip2, latpos.list=latpos.list.wrapped, fixedd=NULL, elbow_graph=2, embed_span=2, graph_attr="weight")
  out12wrapper <- realdata_doMase_wrapper(gip2, latpos.list=latpos.list.wrapped, fixedd=NULL, elbow_graph=2, embed_span=12, graph_attr="weight")
  # Calculate the test stats based on Omni embeddings for the graph using another wrapper function
  out2omniwrapper <- realdata_doOmni2_wrapper(gip2, fixedd=round(sqrt(vcount(gip2[[1]])*2)), elbow_graph=2, graph_attr="weight")

  # Get adjusted p-values based on different embeddings test stats and spans
  adj_pval_m2<-get_graph_adj_pvals(out2wrapper, t_list,latpos.list=NULL, bootstrap_d=NULL, elbow_graph=2, realdata_wrapper=realdata_doMase_wrapper, embed_span=2, xlab="time (yy/mm)", title="MASE(2)",minx=minx, t_window_size=t_window_size, method="BH", return_plot=FALSE, bootstrap_verbose=FALSE)
  adj_pval_m12<-get_graph_adj_pvals(out12wrapper, t_list,latpos.list=NULL, bootstrap_d=NULL, elbow_graph=2, realdata_wrapper=realdata_doMase_wrapper, embed_span=12, xlab="time (yy/mm)", title="MASE(12)",minx=minx, t_window_size=t_window_size, method="BH", return_plot=FALSE, bootstrap_verbose=FALSE)
  adj_pval_o2<-get_graph_adj_pvals(out2omniwrapper, t_list,latpos.list=NULL, bootstrap_d=NULL, elbow_graph=2, realdata_wrapper=realdata_doOmni2_wrapper, embed_span=2, xlab="time (yy/mm)", title="OMNI(2)",minx=minx, t_window_size=t_window_size, method="BH", return_plot=FALSE, bootstrap_verbose=FALSE)

  # Compute quality control charts based on different embeddings test stats and spans
  c1m2 <- get_graph_qcc(gip2, seq(12),latpos.list=latpos.list.wrapped,fixedd=NULL, elbow_graph=2, minx=minx,realdata_wrapper=realdata_doMase_wrapper,embed_span=2,xlab="time points", title="MASE(2)",t_window_size=t_window_size, return_plot=FALSE)
  c1m12 <- get_graph_qcc(gip2, seq(12), latpos.list=latpos.list.wrapped,fixedd=NULL, elbow_graph=2, minx=minx,realdata_wrapper=realdata_doMase_wrapper,embed_span=12,xlab="time points", title="MASE(12)",t_window_size=t_window_size, return_plot=FALSE)
  c1o2 <- get_graph_qcc(gip2, seq(12), fixedd=round(sqrt(vcount(gip2[[1]])*2)), elbow_graph=2, minx=minx,realdata_wrapper=realdata_doOmni2_wrapper,embed_span=2,xlab="time points", title="OMNI(2)",t_window_size=t_window_size, return_plot=FALSE)

  # Attempt to run SCAN method, capturing any errors
  err_message <- try({

    outscanwrapper <- realdata_doScan_wrapper(gip2, fixedd=NULL, graph_attr="weight")
    adj_pval_s2<-get_graph_adj_pvals(outscanwrapper, t_list,latpos.list=NULL, bootstrap_d=NULL, elbow_graph=2, realdata_wrapper=realdata_doScan_wrapper, embed_span=2, xlab="time (yy/mm)", title="SCAN",minx=minx, t_window_size=t_window_size, method="BH", return_plot=FALSE, bootstrap_verbose=FALSE)
    c1s2 <- get_graph_qcc(gip2, seq(12), fixedd=NULL, elbow_graph=2, minx=minx,realdata_wrapper=realdata_doScan_wrapper,embed_span=2,xlab="time points", title="SCAN",t_window_size=t_window_size, return_plot=FALSE)
    'no error'
    # return('no error') # Return this message if no errors occur
  }, silent=FALSE)

  # Handle errors from the try block
  if (inherits(err_message, "try-error")) {
    err_message <- "Error Occurred"
    c1s2 <- list()
    c1s2$violations <- list()
    c1s2$violations$beyond.limits <- NA
  }

  ## Update results table with violation timings for different methods
  ## Table 2

  if (length(c1m2$violations$beyond.limits)>0) {graph_results <- rbind(graph_results, data.frame(time=c1m2$violations$beyond.limits+s,scenario=scenario_inx, type="MASE(2)", approach="cc"))}
  if (length(c1m12$violations$beyond.limits)>0) {graph_results <- rbind(graph_results, data.frame(time=c1m12$violations$beyond.limits+s,scenario=scenario_inx, type="MASE(12)", approach="cc"))}
  if (length(c1o2$violations$beyond.limits)>0) {graph_results <- rbind(graph_results, data.frame(time=c1o2$violations$beyond.limits+s,scenario=scenario_inx, type="OMNI(2)", approach="cc"))}
  if (length(c1s2$violations$beyond.limits)>0) {graph_results <- rbind(graph_results, data.frame(time=c1s2$violations$beyond.limits+s,scenario=scenario_inx, type="SCAN", approach="cc"))}

  # pvalues

  # Checking the error message to decide on which p-value results to consider
  if (err_message != "Error Occurred") {
    # If no error occurred, then all methods including SCAN are considered
    pval_results <- list(adj_pval_m2, adj_pval_m12, adj_pval_o2, adj_pval_s2)
    method_names <- c("MASE(2)","MASE(12)","OMNI(2)","SCAN")
  } else {
    # If error occurred, the SCAN method results are excluded
    pval_results <- list(adj_pval_m2, adj_pval_m12, adj_pval_o2)
    method_names <- c("MASE(2)","MASE(12)","OMNI(2)")
  }

  # Loop through the p-value results
  for (j in 1:(length(pval_results))) {
    # Get the name of the current method being processed
    current_method <- method_names[j]
    # Adjust the p-values (Note: Benjamini-Hochberg)
    p.adj.pvalues <- pval_results[[j]]

    # Table 1: Update results table with significant p-values
    # For each significant p-value (less than 0.05), append its details to the results table
    graph_results <- rbind(graph_results, data.frame(time=which(p.adj.pvalues<0.05)+s,scenario=scenario_inx, type=method_names[j], approach="adj-pval"))
  }
  
  # Save the graph_results object to a file for future use or in case of errors in further steps
  save(graph_results, file =paste0("graph_results_weight_9_",format(Sys.time(), "%Y_%m_%d"),".RData"))
  

}
save(graph_results, file =paste0("graph_results_weight_9_",format(Sys.time(), "%Y_%m_%d"),".RData"))
print(graph_results)
print(injected_time_indices)
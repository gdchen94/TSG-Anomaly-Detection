# R version 4.3.1
options(warn = -1)  # suppress warnings
# options(warn = 0)  # show warnings
# Install and load the packages if not already installed
packages <- unique(c("ggrepel", "readr", "anytime", "tidyverse", "lubridate", "pracma",
                     "igraph", "gridExtra", "ggplotify", "rstudioapi",
                     "gtools", "rARPACK", "metafolio", "latex2exp", "qcc", "rstudioapi",
                     "doParallel", "irlba", "pdp", "reshape2",
                     "pROC", "kernlab", "graphkernels"))

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


out <- read.table("../Data/enron.csv", quote="\"", comment.char="")
out <- tibble(out)
names(out) <- c("time","from","to")
# sort the enron email dataframe by timestamps
out <- out %>% arrange(.[[1]])

# Convert 'time' column to date format
sorted_out <- out %>%
  mutate(time = as_datetime(time))

# Extract month and year information
sorted_out <- sorted_out %>%
  mutate(week = week(time), month = month(time),
         year = year(time))

sorted_out <- sorted_out %>% arrange(time)

# Split dataframe into multiple months
df_list <- sorted_out %>%
  group_split(year)

df_list_ym <- lapply(df_list, function(x) group_split(x, month))

# Flatten the list of dataframes
flat_list_of_df <- unlist(df_list_ym, recursive = FALSE)

# skip the first month because it is not even time gap
flat_list_of_df <- flat_list_of_df[2:length(flat_list_of_df)]

# Convert each dataframe in the flattened list to an igraph object
g_list <- lapply(flat_list_of_df, function(df) graph_from_data_frame(df%>% select(from, to), directed = FALSE))
t_list <- lapply(flat_list_of_df, function(df) strcat(c(as.character(df$year[1]),"/", as.character(df$month[1]))))
A_list <- lapply(g_list, function(x) get.adjacency(x))

for (i in 1:length(g_list)) {
  E(g_list[[i]])$weight <- 1
}
g_simplified_list <- lapply(g_list, function(x)  simplify(x, remove.multiple = TRUE, edge.attr.comb = list(weight="sum", "ignore")))
g_simplified_list <- lapply(g_simplified_list, function(x) ptr(x))

# Get all unique vertices from all graphs
all_vertices <- unique(unlist(lapply(g_simplified_list, function(g) V(g)$name)))
all_vertices <- as.character(sort(as.numeric(all_vertices)))
# Add missing vertices to each graph
g_simplified_list <- lapply(g_simplified_list, add_missing_vertices, all_vertices = all_vertices)

# match the vertices
gip <- lapply(g_simplified_list, function(g) igraph::permute.vertices(g, match(V(g)$name, all_vertices)))
gip <- lapply(gip, function(g) simplify(g, remove.multiple = TRUE, remove.loops = TRUE))

aip <- lapply(gip, function(g) as_adjacency_matrix(g, attr = "weight", sparse = FALSE))

# check simple graph
unlist(lapply(gip, function(x) is.simple(x)))

# Specify the directory where you want to save the file
output_dir <- paste0(dirname(current_script_path), "_realdata_plots_", format(Sys.time(), "%Y-%m-%d"))

# Check if the directory exists, if not, create it
if (!dir.exists(output_dir)) {
  dir.create(output_dir)
}

# Wrapper function to retrieve the latent position list based on the given graphs and attributes
latpos.list.wrapped <- realdata_latpos_list_wrapper(gip, fixedd=NULL, elbow_graph=1, graph_attr="weight")

# Iterate over different time window sizes for analysis
for (t_window_size in seq(10, 12)) {
  
  # Apply MASE, OMNI, and SCAN wrappers to process and prepare the graph to get test stats
  out2wrapper <- realdata_doMase_wrapper(gip, fixedd=NULL, elbow_graph=1, embed_span=2, graph_attr="weight")
  out12wrapper <- realdata_doMase_wrapper(gip, fixedd=NULL, elbow_graph=1, embed_span=12, graph_attr="weight")
  out2omniwrapper <- realdata_doOmni2_wrapper(gip, fixedd=NULL, elbow_graph=1, graph_attr="weight")
  out12omniwrapper <- realdata_doOmni12_wrapper(gip, fixedd=NULL, elbow_graph=1, graph_attr="weight")
  outscanwrapper <- realdata_doScan_wrapper(gip, fixedd=NULL, graph_attr="weight")
  
  # Compute the adjusted p-values for various methods
  adj_pval_pm2 <- get_graph_adj_pvals(out2wrapper, t_list, latpos.list=NULL, bootstrap_d=NULL, elbow_graph=1, realdata_wrapper=realdata_doMase_wrapper, embed_span=2, xlab="time (yy/mm)", title="MASE(2)", minx=NULL, t_window_size=t_window_size, method="BH", return_plot=TRUE, bootstrap_verbose=FALSE)
  adj_pval_pm12 <- get_graph_adj_pvals(out12wrapper, t_list, latpos.list=NULL, bootstrap_d=NULL, elbow_graph=1, realdata_wrapper=realdata_doMase_wrapper, embed_span=12, xlab="time (yy/mm)", title="MASE(12)", minx=NULL, t_window_size=t_window_size, method="BH", return_plot=TRUE, bootstrap_verbose=FALSE)
  adj_pval_po2 <- get_graph_adj_pvals(out2omniwrapper, t_list, latpos.list=NULL, bootstrap_d=NULL, elbow_graph=1, realdata_wrapper=realdata_doOmni2_wrapper, embed_span=2, xlab="time (yy/mm)", title="OMNI(2)", minx=NULL, t_window_size=t_window_size, method="BH", return_plot=TRUE, bootstrap_verbose=FALSE)
  adj_pval_po12 <- get_graph_adj_pvals(out12omniwrapper, t_list, latpos.list=NULL, bootstrap_d=NULL, elbow_graph=1, realdata_wrapper=realdata_doOmni12_wrapper, embed_span=12, xlab="time (yy/mm)", title="OMNI(12)", minx=NULL, t_window_size=t_window_size, method="BH", return_plot=TRUE, bootstrap_verbose=FALSE)
  adj_pval_ps2 <- get_graph_adj_pvals(outscanwrapper, t_list, latpos.list=NULL, bootstrap_d=NULL, elbow_graph=1, realdata_wrapper=realdata_doScan_wrapper, embed_span=2, xlab="time (yy/mm)", title="SCAN", minx=NULL, t_window_size=t_window_size, method="BH", return_plot=TRUE, bootstrap_verbose=FALSE)
  
  # Generate quality control charts for visualization using different methods
  enron_pm2 <- get_graph_qcc(gip, t_list, latpos.list=latpos.list.wrapped, fixedd=NULL, elbow_graph=1, realdata_wrapper=realdata_doMase_wrapper, embed_span=2, xlab="time points", title="MASE(2)", t_window_size=t_window_size, return_plot=TRUE)
  enron_pm12 <- get_graph_qcc(gip, t_list, latpos.list=latpos.list.wrapped, fixedd=NULL, elbow_graph=1, realdata_wrapper=realdata_doMase_wrapper, embed_span=12, xlab="time points", title="MASE(12)", t_window_size=t_window_size, return_plot=TRUE)
  enron_po2 <- get_graph_qcc(gip, t_list, fixedd=NULL, elbow_graph=1, realdata_wrapper=realdata_doOmni2_wrapper, embed_span=2, xlab="time points", title="OMNI(2)", t_window_size=t_window_size, return_plot=TRUE)
  enron_po12 <- get_graph_qcc(gip, t_list, fixedd=NULL, elbow_graph=1, realdata_wrapper=realdata_doOmni12_wrapper, embed_span=12, xlab="time points", title="OMNI(12)", t_window_size=t_window_size, return_plot=TRUE)
  enron_ps <- get_graph_qcc(gip, t_list,fixedd=NULL, realdata_wrapper=realdata_doScan_wrapper,embed_span=2,xlab="time points", title="SCAN",t_window_size=t_window_size, return_plot=TRUE)
  
  # Figures 11, 12, 13
  plots <- list(adj_pval_pm2, adj_pval_pm12, adj_pval_po2, adj_pval_po12, adj_pval_ps2)
  
  for (plot in plots) {
    print(plot)
  }
  png(paste0(output_dir, paste0("/enron_g_pval_t_window_size",as.character(t_window_size),"_",format(Sys.time(), "%Y_%m_%d")),".png"), width = 10, height = 9, units="in", res=400)
  grid.arrange(adj_pval_pm2,adj_pval_pm12,adj_pval_po2,adj_pval_po12,adj_pval_ps2, nrow=3, heights=c(4,4,4))
  dev.off()
  
  # Figures 8, 9, 10
  plots <- list(enron_pm2,enron_pm12,enron_po2,enron_po12,enron_ps)
  
  for (plot in plots) {
    print(plot)
  }
  png(paste0(output_dir, paste0("/enron_g_cc_t_window_size",as.character(t_window_size),"_",format(Sys.time(), "%Y_%m_%d")),".png"), width = 10, height = 9, units="in", res=400)
  grid.arrange(enron_pm2,enron_pm12,enron_po2,enron_po12,enron_ps, nrow=3, heights=c(4,4,4))
  dev.off()
}


# Specify a time index for analysis
temp.anomalous.time.inx <- 5

# Define a window size for analysis
t_window_size=11

# Get control chart details for given window size and time index
c1m2 <- get_graph_qcc(gip, t_list,latpos.list=latpos.list.wrapped,fixedd=NULL, elbow_graph=1, realdata_wrapper=realdata_doMase_wrapper,embed_span=2,xlab="time points", title="MASE(2)",t_window_size=t_window_size, return_plot=FALSE)
c1m2v <- get_vertex_qcc(gip, t_list, middle.max.inx=NULL,latpos.list=NULL, fixedd=NULL, elbow_graph=1, realdata_wrapper=realdata_doMase_wrapper,embed_span=2,xlab="time points", title="MASE(2)", minx=NULL,t_window_size=t_window_size, return_plot=FALSE)
n.vertex <- vcount(gip[[1]])
tempinx <- (c1m2$violations$beyond.limits+t_window_size)[temp.anomalous.time.inx]
temp.latpos.list <- jrdpg.latent(aip[tempinx:(tempinx+1)], latpos.list = latpos.list.wrapped[tempinx:(tempinx+1)], d1=NA,d2=NA, d3=NA, center=FALSE, python=FALSE, SVD=2)$Xhat

temp <- rep(FALSE, n.vertex)
temp[c1m2v[[c1m2$violations$beyond.limits[temp.anomalous.time.inx]]]$violations$beyond.limits] = "anomalous"
df <- data.frame(x=temp.latpos.list[[1]][,1], y=temp.latpos.list[[1]][,2],class=temp, label=seq(n.vertex))

# Create ggplot visualizations for the latent positions before and after the anomaly
p.latposprev <- ggplot() +
  geom_point(data=df%>% filter((class!="anomalous")), aes(x=x, y=y),color="grey70",shape=20,size=3) +
  geom_point(data=df%>% filter(class=="anomalous"), aes(x=x, y=y),color="red",shape=17,size=3) +
  geom_text_repel(data=df%>% filter(class=="anomalous"),aes(x=x, y=y,label=label, size=16))+#, hjust=-.05) +
  scale_color_manual(values=c("black")) +
  theme_void() + theme(legend.position="none",plot.title = element_text(hjust = 0.5, face = "italic", size=16))+ggtitle(strcat("t = ", t_list[[tempinx]]))

tempinx <- (c1m2$violations$beyond.limits+t_window_size)[temp.anomalous.time.inx] + 1
temp <- rep(FALSE, n.vertex)
temp[c1m2v[[c1m2$violations$beyond.limits[temp.anomalous.time.inx]]]$violations$beyond.limits] = "anomalous"
df <- data.frame(x=temp.latpos.list[[2]][,1], y=temp.latpos.list[[2]][,2],class=temp,label=seq(n.vertex))#
p.latposafter <- ggplot() +
  geom_point(data=df%>% filter((class!="anomalous")), aes(x=x, y=y),color="grey70",shape=20,size=3) +
  geom_point(data=df%>% filter(class=="anomalous"), aes(x=x, y=y),color="red",shape=17,size=3) +
  geom_text_repel(data=df%>% filter(class=="anomalous"),aes(x=x, y=y,label=label, size=16))+#, hjust=-.05) +
  scale_color_manual(values=c("black")) +
  theme_void() + theme(legend.position="none",plot.title = element_text(hjust = 0.5, face = "italic", size=16))+ggtitle(strcat("t = ", t_list[[tempinx]]))


# Figure 6 b
grid.arrange(p.latposprev, p.latposafter, nrow=2, heights=c(2,2))
png(paste0(output_dir, paste0("/enron_vertex_illus_cc_t_window_size",as.character(t_window_size),"_",format(Sys.time(), "%Y_%m_%d")),".png"), width = 5, height = 4, units="in", res=400)
grid.arrange(p.latposprev, p.latposafter, nrow=2, heights=c(2,2))
dev.off()

# Set the size of the window for the time series
t_window_size <- 11

# Calculate quality control charts (qcc) for different embedding spans and methods
# Get qcc for MASE(2)
enron_c_m2 <- get_graph_qcc(gip, t_list,latpos.list=latpos.list.wrapped,fixedd=NULL, elbow_graph=1, realdata_wrapper=realdata_doMase_wrapper,embed_span=2,xlab="time points", title="MASE(2)",t_window_size=t_window_size, return_plot=FALSE)

# Get qcc for MASE(12)
enron_c_m12 <- get_graph_qcc(gip, t_list,latpos.list=latpos.list.wrapped,fixedd=NULL, elbow_graph=1, realdata_wrapper=realdata_doMase_wrapper,embed_span=12,xlab="time points", title="MASE(12)",t_window_size=t_window_size, return_plot=FALSE)

# Get qcc for OMNI(2)
enron_c_o2 <- get_graph_qcc(gip, t_list,fixedd=NULL, elbow_graph=1, realdata_wrapper=realdata_doOmni2_wrapper,embed_span=2,xlab="time points", title="OMNI(2)",t_window_size=t_window_size, return_plot=FALSE)

# Get qcc for OMNI(12)
enron_c_o12 <- get_graph_qcc(gip, t_list,fixedd=NULL, elbow_graph=1, realdata_wrapper=realdata_doOmni12_wrapper,embed_span=12,xlab="time points", title="OMNI(12)",t_window_size=t_window_size, return_plot=FALSE)

# Get qcc for SCAN
enron_c_s <- get_graph_qcc(gip, t_list,fixedd=NULL, realdata_wrapper=realdata_doScan_wrapper,embed_span=2,xlab="time points", title="SCAN",t_window_size=t_window_size, return_plot=FALSE)

# Construct a dataframe by combining statistics from different methods
# Start with MASE(2)
df <- data.frame(time=seq(length(enron_c_m2$statistics)), y= enron_c_m2$statistics,ucl=enron_c_m2$limits[,2],center=enron_c_m2$center,lim=replace(rep(0,length(enron_c_m2$statistics)),enron_c_m2$violations$beyond.limits,1),run=replace(rep(0,length(enron_c_m2$statistics)),enron_c_m2$violations$violating.runs,1), Method="MASE")

# Add OMNI(2) data
df <- rbind(df, data.frame(time=seq(length(enron_c_o2$statistics)), y= enron_c_o2$statistics,ucl=enron_c_o2$limits[,2],center=enron_c_o2$center,lim=replace(rep(0,length(enron_c_o2$statistics)),enron_c_o2$violations$beyond.limits,1),run=replace(rep(0,length(enron_c_o2$statistics)),enron_c_o2$violations$violating.runs,1), Method="OMNI"))

# Add SCAN data
df <- rbind(df, data.frame(time=seq(length(enron_c_s$statistics)), y= enron_c_s$statistics,ucl=enron_c_s$limits[,2],center=enron_c_s$center,lim=replace(rep(0,length(enron_c_s$statistics)),enron_c_s$violations$beyond.limits,1),run=replace(rep(0,length(enron_c_s$statistics)),enron_c_s$violations$violating.runs,1), Method="SCAN"))

# Re-order the levels in the Method factor
df$Method <- factor(df$Method, levels = c("MASE","OMNI","SCAN"))

# Generate labels for the x-axis using time periods
minx <- lapply(t_list, function(s) substr(s, 3, nchar(s)))
minx <- paste(minx, lead(minx), sep = ":")
minx <- minx[(t_window_size+1):(length(gip)-1)]

# Plot the data in the dataframe using ggplot2
p <- ggplot(df,aes(x=time, y=y))+
  geom_step(aes(x=time, y=ucl), linetype="dashed")+
  geom_step(aes(x=time, y=center), linetype="solid") +
  geom_point(data = df, alpha=1, color="grey70",shape=20)+ylab(TeX("$y^{(t)}$"))+
  geom_line(aes(x=time, y=y), color="grey70")+
  geom_point(data = df %>% filter(lim==1), alpha=1, color="red",shape=17) + 
  scale_x_discrete(name ="time (yy/mm)", 
                   limits=c(minx))+theme_bw()+
  theme_classic(base_size = 18)+
  facet_wrap(~Method,nrow = 3,scales = "free_y", strip.position = "left")+
  theme(axis.text.x = element_text(size = 14,angle=90),legend.position = "none",plot.title = element_text(hjust = 0.5,size=20, face='bold'),plot.subtitle = element_text(hjust = 0.5)) #+ggtitle("Control Chart")

# Display the plot
p

# Save the plot as a PNG file with a timestamp in the filename
png(paste0(output_dir, paste0("/enron_comb_cc_t_window_size_",as.character(t_window_size),"_",format(Sys.time(), "%Y_%m_%d")),".png"), width =7, height = 4, units="in", res=400)
p
dev.off()



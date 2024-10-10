# R version 4.3.1
# Install and load the packages if not already installed
packages <- c("gtools", "rARPACK", "metafolio", "latex2exp", "qcc", "rstudioapi",
              "doParallel", "irlba", "pdp", "reshape2", "tidyverse", "invgamma",
              "igraph", "pracma", "gridExtra", "pROC", "kernlab", "graphkernels")

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

# registerDoParallel(detectCores()-1) ## uncomment this to use multicore

# Running this code without GP will take about 12 hrs.
# Running GP results will take about 60 hrs.
#power analysis at all timepoints
d_true <- 4
#p-values
method3 = 2
nullnmc <- 50
pvalnmc <- 100
tmax <- 12
tvec <- 1:tmax
m2 <- names(running(tvec, width=2))
n_values <- seq(80, 400, by=80) # Increase the number of vertices from 100 to 1000 in steps of 100
nperturb <- 100
rmin <- 0.2
rmax <- 0.8
center <- FALSE
cperturb <- .04

approx <- FALSE
center <- FALSE

d.max <- "sqrt"
plot <- FALSE
dASE <- 4
dSVD <- 4
dSVD12 <- 4
elbow_graph <- 2

# Initialize the final data frame to store results
df <- data.frame()

# df_xx represent the cts of rejecting the null using thres .05 under MC
df_tnorm_combined <- c()
df_pdist_combined <- c()
df_rr_combined <- data.frame()

# df_xx represent the cts of rejecting the null using thres .05 under MC
df_tnorm <- c()
df_pdist <- c()
for(n in n_values){
  alpha = 0 # Since we are not changing alpha
  
  # Initialize output data frames
  # out_xx represent the empirical null distributions
  out_tnorm <- c()
  out_pdist <- c()
  
  # Initialize vectors to hold the count of significant p-values
  sig_count_tnorm <- matrix(0, tmax-1, 2) # 2 for MASE(2) and MASE(12)
  sig_count_pdist <- array(0, dim = c(n, tmax - 1, 2)) # 2 for MASE(2) and MASE(12)
  
  # Combined Monte Carlo Simulation
  for(i in 1:(nullnmc + pvalnmc)) {
    set.seed(123 + i - 1)
    glist <- genTSGonepar(n, nperturb, cperturb * as.numeric(i > nullnmc), rmin, rmax, tmax, d_true, alpha)
    
    # Do MASE for different dimensions
    out2 <- doMase(glist, latpos.list=NULL, nmase=2, dSVD=dSVD, dASE=dASE, d.max=d.max, center=center, approx=approx, elbow_graph=elbow_graph, plot=plot, latent.form=method3)
    out12 <- doMase(glist, latpos.list=NULL, nmase=12, dSVD=dSVD12, dASE=dASE, d.max=d.max, center=center, approx=approx, elbow_graph=elbow_graph, plot=plot, latent.form=method3)
    
    # Process for tnorm and pdist if within the first nullnmc simulations
    if (i <= nullnmc) {
      # Process tnorm data
      out.p2_tnorm <- data.frame(as.vector(t(out2$tnorm[-c(5,6,7)])), Method="MASE(2)")
      out.p12_tnorm <- data.frame(as.vector(t(out12$tnorm[-c(5,6,7)])), Method="MASE(12)")
      
      colnames(out.p2_tnorm)[1] <- "tnorm"
      colnames(out.p12_tnorm)[1] <- "tnorm"
      out_tnorm <- rbind(out_tnorm, rbind(out.p2_tnorm, out.p12_tnorm))
      
      # Process pdist data
      out.p2_pdist <- data.frame(t(out2$pdist[,-c(5,6,7)]), Method="MASE(2)")
      out.p12_pdist <- data.frame(t(out12$pdist[,-c(5,6,7)]), Method="MASE(12)")
      
      out_pdist <- rbind(out_pdist, rbind(out.p2_pdist, out.p12_pdist))
    }
    
    # Process for p-values if beyond the first nullnmc simulations
    if (i > nullnmc && i <= (nullnmc + pvalnmc)) {
      
      rr2 <- apply(out2$pdist, 2, function(x) rank(1 / x))
      rr12 <- apply(out12$pdist, 2, function(x) rank(1 / x))
      
      # Process tnorm p-values
      pvalues2_tnorm <- rep(0, (tmax-1))
      pvalues12_tnorm <- rep(0, (tmax-1))
      for (j in 1:(tmax-1)) {
        pvalues2_tnorm[j] <- sum(out2$tnorm[j] < (out_tnorm[which(out_tnorm$Method == "MASE(2)"), 1])) / (nullnmc*8)
        pvalues12_tnorm[j] <- sum(out12$tnorm[j] < (out_tnorm[which(out_tnorm$Method == "MASE(12)"), 1])) / (nullnmc*8)
      }
      # Adjust p-values for tnorm
      pvalues2_tnorm <- p.adjust(pvalues2_tnorm, "BH")
      pvalues12_tnorm <- p.adjust(pvalues12_tnorm, "BH")
      
      # Count significant p-values for tnorm
      sig_count_tnorm[,1] <- sig_count_tnorm[,1] + (pvalues2_tnorm < 0.05)
      sig_count_tnorm[,2] <- sig_count_tnorm[,2] + (pvalues12_tnorm < 0.05)
      
      # Process pdist p-values
      # Assuming similar logic for pdist p-values
      # Calculate p-values for pdist
      pvalues2_pdist <- matrix(0, nrow = n, ncol = tmax - 1)
      pvalues12_pdist <- matrix(0, nrow = n, ncol = tmax - 1)
      
      for (j in 1:(tmax - 1)) {
        for (w in 1:n) {
          pvalues2_pdist[w, j] <- sum(out2$pdist[w, j] < (out_pdist[which(out_pdist$Method == "MASE(2)"), w])) / (nullnmc*8)
          pvalues12_pdist[w, j] <- sum(out12$pdist[w, j] < (out_pdist[which(out_pdist$Method == "MASE(12)"), w])) / (nullnmc*8)
        }
      }
      
      # Adjust p-values for pdist
      pvalues2_pdist_adjusted <- matrix(p.adjust(pvalues2_pdist, "BH"), n, tmax-1)
      pvalues12_pdist_adjusted <- matrix(p.adjust(pvalues12_pdist, "BH"), n, tmax-1)
      
      # Count significant p-values at each vertex for each time point
      sig_count_pdist[, , 1] <- sig_count_pdist[, , 1] + (pvalues2_pdist_adjusted < 0.05)
      sig_count_pdist[, , 2] <- sig_count_pdist[, , 2] + (pvalues12_pdist_adjusted < 0.05)
      
      df_rr2 <- melt(1 / rr2)
      names(df_rr2) <- c("vertex", "time", "rr")
      df_rr2$time <- rep(factor(m2, levels = m2), each = n)
      df_rr2 <- cbind(df_rr2, mc = i, Method = "MASE(2)", n = factor(n))
      
      df_rr12 <- melt(1 / rr12)
      names(df_rr12) <- c("vertex", "time", "rr")
      df_rr12$time <- rep(factor(m2, levels = m2), each = n)
      df_rr12 <- cbind(df_rr12, mc = i, Method = "MASE(12)", n = factor(n))
      
      # Combine the results    
      df_rr_combined <- rbind(df_rr_combined, rbind(df_rr2, df_rr12))
    }
  }
  
  df2 <- data.frame(time=factor(m2, levels=m2), power=sig_count_tnorm[,1]/pvalnmc, se = sqrt( (sig_count_tnorm[,1]/pvalnmc) * (1 - sig_count_tnorm[,1]/pvalnmc)/pvalnmc ), Method="MASE(2)", n=factor(n))
  df12 <- data.frame(time=factor(m2, levels=m2), power=sig_count_tnorm[,2]/pvalnmc, se = sqrt( (sig_count_tnorm[,2]/pvalnmc) * (1 - sig_count_tnorm[,2]/pvalnmc)/pvalnmc ), Method="MASE(12)", n=factor(n))
  
  df_tnorm <- rbind(df_tnorm, rbind(df2, df12))
  
  df2_pdist <- melt(sig_count_pdist[, , 1])
  names(df2_pdist) <- c("vertex", "time", "cts")
  df2_pdist$time <- rep(factor(m2, levels=m2), each = n)
  df2_pdist <- cbind(df2_pdist, Method="MASE(2)", n=factor(n))
  
  df12_pdist <- melt(sig_count_pdist[, , 2])
  names(df12_pdist) <- c("vertex", "time", "cts")
  df12_pdist$time <- rep(factor(m2, levels=m2), each = n)
  df12_pdist <- cbind(df12_pdist, Method="MASE(12)", n=factor(n))
  
  df_pdist <- rbind(df_pdist, rbind(df2_pdist, df12_pdist))
}

df_tnorm_null <- df_tnorm %>% filter(time!="5:6"&time!="6:7" & time!="7:8")
df_tnorm_type_one_err <- df_tnorm_null[,c(2,4,5)]
df_tnorm_type_one_err$power <- df_tnorm_type_one_err$power * pvalnmc
# Summarizing the data
df_tnorm_type_one_err <- df_tnorm_type_one_err %>%
  group_by(Method, n) %>%
  summarize(
    across(
      everything(),
      list(
        fpr = ~mean(.) / pvalnmc,
        se = ~sqrt((mean(.) / pvalnmc) * (1 - mean(.) / pvalnmc) / pvalnmc)
      )
    )
  )

# Renaming the first column
df_tnorm_type_one_err <- df_tnorm_type_one_err %>%
  rename(fpr=power_fpr, se=power_se)

# plots without GP
df_tnorm_power <- df_tnorm %>% filter(time=="6:7")
colnames(df_tnorm_power)[4] <- "Method"
colnames(df_tnorm_power)[2] <- "mean"
colnames(df_tnorm_type_one_err)[3] <- "mean"
df_tnorm_power$metric <- "power"
df_tnorm_type_one_err$metric <- "FPR"

df_tnorm_power <- select(df_tnorm_power, colnames(df_tnorm_type_one_err))

cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

method_colors <- c("MASE(2)" = cbbPalette[1], "MASE(12)" = cbbPalette[2])
method_shapes <- c("MASE(2)" = 17, "MASE(12)" = 19)

# Merging and renaming for final graph
df.graph <- bind_rows(df_tnorm_power, df_tnorm_type_one_err)
df.graph$metric <- factor(df.graph$metric, levels = c("power", "FPR"))

pg <- ggplot(df.graph, aes(x=n, y=mean, color=Method,group=Method)) +
  geom_errorbar(aes(x=n,ymin=mean-se, ymax=mean+se, width=.1,linetype=Method, color=Method)) +
  geom_line(aes(x=n,y=mean,linetype=Method)) +theme_bw()+xlab("Number of Vertices")+
  geom_point(aes(x=n,y=mean, shape=Method)) +
  scale_color_manual(values = method_colors) + 
  scale_shape_manual(values = method_shapes) +
  facet_wrap(~metric,nrow = 2,scales = "free_y",strip.position="left")+
  theme(plot.title = element_text(size=10, face = 'bold'),legend.text = element_text(size=8, face='bold'),axis.title = element_text(face="bold"),axis.text = element_text(face="bold"))  + ylab("")

print(pg)

png(paste0("pg_aymp","_",format(Sys.time(), "%Y_%m_%d.png")), width = 4, height = 4, units="in", res=400)
pg
dev.off()
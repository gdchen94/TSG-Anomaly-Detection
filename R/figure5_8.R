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
nullnmc <- 50#25#50
pvalnmc <- 100#200#12#00#200
tmax <- 12
tvec <- 1:tmax
m2 <- names(running(tvec, width=2))
n <- 400
nperturb <- 100
rmin <- 0.2
rmax <- 0.8
center <- FALSE
cperturb <- .12

approx <- FALSE
center <- FALSE

d.max <- "sqrt"
plot <- FALSE
dASE <- 4
dSVD <- 4
dSVD12 <- 4
elbow_graph <- 2
#pairwise mase
# Initialize the final data frame to store results
df <- data.frame()

# df_xx represent the cts of rejecting the null using thres .05 under MC
df_tnorm_combined <- c()
df_pdist_combined <- c()
df_rr_combined <- data.frame()

# df_xx represent the cts of rejecting the null using thres .05 under MC
df_tnorm <- c()
df_pdist <- c()
alpha_values = c(0, seq(8)/8)
for(a in 1:length(alpha_values)){
  alpha = alpha_values[a]
  
  
  # Initialize output data frames
  # out_xx represent the empirical null distributions
  out_tnorm <- c()
  out_pdist <- c()
  
  # Initialize vectors to hold the count of significant p-values
  sig_count_tnorm <- matrix(0, tmax-1, 6)
  sig_count_pdist <- array(0, dim = c(n, tmax - 1, 6)) # 6 for MASE(2), MASE(12), OMNI(2), OMNI(12), SCAN and DIST
  
  
  # Combined Monte Carlo Simulation
  for(i in 1:(nullnmc + pvalnmc)) {
    set.seed(123 + i - 1)
    glist <- genTSGonepar(n, nperturb, cperturb * as.numeric(i > nullnmc), rmin, rmax, tmax, d_true, alpha)
    
    # Do MASE for different dimensions
    out2 <- doMase(glist, latpos.list=NULL, nmase=2, dSVD=dSVD, dASE=dASE, d.max=d.max, center=center, approx=approx, elbow_graph=elbow_graph, plot=plot, latent.form=method3)
    out12 <- doMase(glist, latpos.list=NULL, nmase=12, dSVD=dSVD12, dASE=dASE, d.max=d.max, center=center, approx=approx, elbow_graph=elbow_graph, plot=plot, latent.form=method3)
    outo2 <- doOmni(glist, omniase=NULL, nomni=2, dASE=dASE, d.max=d.max, center=center, approx=approx, elbow_graph=elbow_graph, plot=plot)
    outo12 <- doOmni(glist, omniase=NULL, nomni=12, dASE=dASE, d.max=d.max, center=center, approx=approx, elbow_graph=elbow_graph, plot=plot)
    outs <- doScan(glist)
    outadj <- doAdj(glist)
    
    # Process for tnorm and pdist if within the first nullnmc simulations
    if (i <= nullnmc) {
      # Process tnorm data
      out.p2_tnorm <- data.frame(as.vector(t(out2$tnorm[-c(5,6,7)])), Method="MASE(2)")
      out.p12_tnorm <- data.frame(as.vector(t(out12$tnorm[-c(5,6,7)])), Method="MASE(12)")
      out.po2_tnorm <- data.frame(as.vector(t(outo2$tnorm[-c(5,6,7)])), Method="OMNI(2)")
      out.po12_tnorm <- data.frame(as.vector(t(outo12$tnorm[-c(5,6,7)])), Method="OMNI(12)")
      out.ps_tnorm <- data.frame(as.vector(t(outs$tnorm[-c(5,6,7)])), Method="SCAN")
      out.padj_tnorm <- data.frame(as.vector(t(outadj$tnorm[-c(5,6,7)])), Method="DIST")
      
      colnames(out.p2_tnorm)[1] <- "tnorm"
      colnames(out.p12_tnorm)[1] <- "tnorm"
      colnames(out.po2_tnorm)[1] <- "tnorm"
      colnames(out.po12_tnorm)[1] <- "tnorm"
      colnames(out.ps_tnorm)[1] <- "tnorm"
      colnames(out.padj_tnorm)[1] <- "tnorm"
      out_tnorm <- rbind(out_tnorm, rbind(out.p2_tnorm, out.p12_tnorm, out.po2_tnorm, out.po12_tnorm, out.ps_tnorm, out.padj_tnorm))
      
      # Process pdist data
      out.p2_pdist <- data.frame(t(out2$pdist[,-c(5,6,7)]), Method="MASE(2)")
      out.p12_pdist <- data.frame(t(out12$pdist[,-c(5,6,7)]), Method="MASE(12)")
      out.po2_pdist <- data.frame(t(outo2$pdist[,-c(5,6,7)]), Method="OMNI(2)")
      out.po12_pdist <- data.frame(t(outo12$pdist[,-c(5,6,7)]), Method="OMNI(12)")
      out.ps_pdist <- data.frame(t(outs$pdist[,-c(5,6,7)]), Method="SCAN")
      out.padj_pdist <- data.frame(t(outadj$pdist[,-c(5,6,7)]), Method="DIST")
      
      out_pdist <- rbind(out_pdist, rbind(out.p2_pdist, out.p12_pdist, out.po2_pdist, out.po12_pdist, out.ps_pdist, out.padj_pdist))
    }
    
    # Process for p-values if beyond the first nullnmc simulations
    if (i > nullnmc && i <= (nullnmc + pvalnmc)) {
      
      rr2 <- apply(out2$pdist, 2, function(x) rank(1 / x))
      rr12 <- apply(out12$pdist, 2, function(x) rank(1 / x))
      rro2 <- apply(outo2$pdist, 2, function(x) rank(1 / x))
      rro12 <- apply(outo12$pdist, 2, function(x) rank(1 / x))
      rrs <- apply(outs$pdist, 2, function(x) rank(1 / x))
      rradj <- apply(outadj$pdist, 2, function(x) rank(1 / x))
      
      # Process tnorm p-values
      pvalues2_tnorm <- rep(0, (tmax-1))
      pvalues12_tnorm <- rep(0, (tmax-1))
      pvalueso2_tnorm <- rep(0, (tmax-1))
      pvalueso12_tnorm <- rep(0, (tmax-1))
      pvaluess_tnorm <- rep(0, (tmax-1))
      pvaluesadj_tnorm <- rep(0, (tmax-1))
      for (j in 1:(tmax-1)) {
        pvalues2_tnorm[j] <- sum(out2$tnorm[j] < (out_tnorm[which(out_tnorm$Method == "MASE(2)"), 1])) / (nullnmc*8)
        pvalues12_tnorm[j] <- sum(out12$tnorm[j] < (out_tnorm[which(out_tnorm$Method == "MASE(12)"), 1])) / (nullnmc*8)
        pvalueso2_tnorm[j] <- sum(outo2$tnorm[j] < (out_tnorm[which(out_tnorm$Method == "OMNI(2)"), 1])) / (nullnmc*8)
        pvalueso12_tnorm[j] <- sum(outo12$tnorm[j] < (out_tnorm[which(out_tnorm$Method == "OMNI(12)"), 1])) / (nullnmc*8)
        pvaluess_tnorm[j] <- sum(outs$tnorm[j] < (out_tnorm[which(out_tnorm$Method == "SCAN"), 1])) / (nullnmc*8)
        pvaluesadj_tnorm[j] <- sum(outadj$tnorm[j] < (out_tnorm[which(out_tnorm$Method == "DIST"), 1])) / (nullnmc*8)
      }
      # Adjust p-values for tnorm
      pvalues2_tnorm <- p.adjust(pvalues2_tnorm, "BH")
      pvalues12_tnorm <- p.adjust(pvalues12_tnorm, "BH")
      pvalueso2_tnorm <- p.adjust(pvalueso2_tnorm, "BH")
      pvalueso12_tnorm <- p.adjust(pvalueso12_tnorm, "BH")
      pvaluess_tnorm <- p.adjust(pvaluess_tnorm, "BH")
      pvaluesadj_tnorm <- p.adjust(pvaluesadj_tnorm, "BH")
      
      # Count significant p-values for tnorm
      sig_count_tnorm[,1] <- sig_count_tnorm[,1] + (pvalues2_tnorm < 0.05)
      sig_count_tnorm[,2] <- sig_count_tnorm[,2] + (pvalues12_tnorm < 0.05)
      sig_count_tnorm[,3] <- sig_count_tnorm[,3] + (pvalueso2_tnorm < 0.05)
      sig_count_tnorm[,4] <- sig_count_tnorm[,4] + (pvalueso12_tnorm < 0.05)
      sig_count_tnorm[,5] <- sig_count_tnorm[,5] + (pvaluess_tnorm < 0.05)
      sig_count_tnorm[,6] <- sig_count_tnorm[,6] + (pvaluesadj_tnorm < 0.05)
      
      
      # Process pdist p-values
      # Assuming similar logic for pdist p-values
      # Calculate p-values for pdist
      pvalues2_pdist <- matrix(0, nrow = n, ncol = tmax - 1)
      pvalues12_pdist <- matrix(0, nrow = n, ncol = tmax - 1)
      pvalueso2_pdist <- matrix(0, nrow = n, ncol = tmax - 1)
      pvalueso12_pdist <- matrix(0, nrow = n, ncol = tmax - 1)
      pvaluess_pdist <- matrix(0, nrow = n, ncol = tmax - 1)
      pvaluesadj_pdist <- matrix(0, nrow = n, ncol = tmax - 1)
      
      for (j in 1:(tmax - 1)) {
        for (w in 1:n) {
          pvalues2_pdist[w, j] <- sum(out2$pdist[w, j] < (out_pdist[which(out_pdist$Method == "MASE(2)"), w])) / (nullnmc*8)
          pvalues12_pdist[w, j] <- sum(out12$pdist[w, j] < (out_pdist[which(out_pdist$Method == "MASE(12)"), w])) / (nullnmc*8)
          pvalueso2_pdist[w, j] <- sum(outo2$pdist[w, j] < (out_pdist[which(out_pdist$Method == "OMNI(2)"), w])) / (nullnmc*8)
          pvalueso12_pdist[w, j] <- sum(outo12$pdist[w, j] < (out_pdist[which(out_pdist$Method == "OMNI(12)"), w])) / (nullnmc*8)
          pvaluess_pdist[w, j] <- sum(outs$pdist[w, j] < (out_pdist[which(out_pdist$Method == "SCAN"), w])) / (nullnmc*8)
          pvaluesadj_pdist[w, j] <- sum(outadj$pdist[w, j] < (out_pdist[which(out_pdist$Method == "DIST"), w])) / (nullnmc*8)
        }
      }
      
      # Adjust p-values for pdist
      pvalues2_pdist_adjusted <- matrix(p.adjust(pvalues2_pdist, "BH"),n, tmax-1)
      pvalues12_pdist_adjusted <- matrix(p.adjust(pvalues12_pdist, "BH"),n, tmax-1)
      pvalueso2_pdist_adjusted <- matrix(p.adjust(pvalueso2_pdist, "BH"),n, tmax-1)
      pvalueso12_pdist_adjusted <- matrix(p.adjust(pvalueso12_pdist, "BH"),n, tmax-1)
      pvaluess_pdist_adjusted <- matrix(p.adjust(pvaluess_pdist, "BH"),n, tmax-1)
      pvaluesadj_pdist_adjusted <- matrix(p.adjust(pvaluesadj_pdist, "BH"),n, tmax-1)
      
      # Count significant p-values at each vertex for each time point
      sig_count_pdist[, , 1] <- sig_count_pdist[, , 1] + (pvalues2_pdist_adjusted < 0.05)
      sig_count_pdist[, , 2] <- sig_count_pdist[, , 2] + (pvalues12_pdist_adjusted < 0.05)
      sig_count_pdist[, , 3] <- sig_count_pdist[, , 3] + (pvalueso2_pdist_adjusted < 0.05)
      sig_count_pdist[, , 4] <- sig_count_pdist[, , 4] + (pvalueso12_pdist_adjusted < 0.05)
      sig_count_pdist[, , 5] <- sig_count_pdist[, , 5] + (pvaluess_pdist_adjusted < 0.05)
      sig_count_pdist[, , 6] <- sig_count_pdist[, , 6] + (pvaluesadj_pdist_adjusted < 0.05)
      
      df_rr2 <- melt(1 / rr2)
      names(df_rr2) <- c("vertex", "time", "rr")
      df_rr2$time <- rep(factor(m2, levels = m2), each = n)
      df_rr2 <- cbind(df_rr2, mc = i, Method = "MASE(2)", alpha = factor(alpha))
      
      df_rr12 <- melt(1 / rr12)
      names(df_rr12) <- c("vertex", "time", "rr")
      df_rr12$time <- rep(factor(m2, levels = m2), each = n)
      df_rr12 <- cbind(df_rr12, mc = i, Method = "MASE(12)", alpha = factor(alpha))
      
      df_rro2 <- melt(1 / rro2)
      names(df_rro2) <- c("vertex", "time", "rr")
      df_rro2$time <- rep(factor(m2, levels = m2), each = n)
      df_rro2 <- cbind(df_rro2, mc = i, Method = "OMNI(2)", alpha = factor(alpha))
      
      df_rro12 <- melt(1 / rro12)
      names(df_rro12) <- c("vertex", "time", "rr")
      df_rro12$time <- rep(factor(m2, levels = m2), each = n)
      df_rro12 <- cbind(df_rro12, mc = i, Method = "OMNI(12)", alpha = factor(alpha))
      
      df_rrs <- melt(1 / rrs)
      names(df_rrs) <- c("vertex", "time", "rr")
      df_rrs$time <- rep(factor(m2, levels = m2), each = n)
      df_rrs <- cbind(df_rrs, mc = i, Method = "SCAN", alpha = factor(alpha))
      
      df_rradj <- melt(1 / rradj)
      names(df_rradj) <- c("vertex", "time", "rr")
      df_rradj$time <- rep(factor(m2, levels = m2), each = n)
      df_rradj <- cbind(df_rradj, mc = i, Method = "DIST", alpha = factor(alpha))
      
      # Combine the results    
      df_rr_combined <- rbind(df_rr_combined, rbind(df_rr2, df_rr12, df_rro2, df_rro12, df_rrs, df_rradj))
    }
    
    
  }
  df2 <- data.frame(time=factor(m2, levels=m2), power=sig_count_tnorm[,1]/pvalnmc, se = sqrt( (sig_count_tnorm[,1]/pvalnmc) * (1 - sig_count_tnorm[,1]/pvalnmc)/pvalnmc ), Method="MASE(2)", alpha=factor(alpha))
  df12 <- data.frame(time=factor(m2, levels=m2), power=sig_count_tnorm[,2]/pvalnmc, se = sqrt( (sig_count_tnorm[,2]/pvalnmc) * (1 - sig_count_tnorm[,2]/pvalnmc)/pvalnmc ), Method="MASE(12)", alpha=factor(alpha))
  dfo2 <- data.frame(time=factor(m2, levels=m2), power=sig_count_tnorm[,3]/pvalnmc, se = sqrt( (sig_count_tnorm[,3]/pvalnmc) * (1 - sig_count_tnorm[,3]/pvalnmc)/pvalnmc ), Method="OMNI(2)", alpha=factor(alpha))
  dfo12 <- data.frame(time=factor(m2, levels=m2), power=sig_count_tnorm[,4]/pvalnmc, se = sqrt( (sig_count_tnorm[,4]/pvalnmc) * (1 - sig_count_tnorm[,4]/pvalnmc)/pvalnmc ), Method="OMNI(12)", alpha=factor(alpha))
  dfs <- data.frame(time=factor(m2, levels=m2), power=sig_count_tnorm[,5]/pvalnmc, se = sqrt( (sig_count_tnorm[,5]/pvalnmc) * (1 - sig_count_tnorm[,5]/pvalnmc)/pvalnmc ), Method="SCAN", alpha=factor(alpha))
  dfadj <- data.frame(time=factor(m2, levels=m2), power=sig_count_tnorm[,6]/pvalnmc, se = sqrt( (sig_count_tnorm[,6]/pvalnmc) * (1 - sig_count_tnorm[,6]/pvalnmc)/pvalnmc ), Method="DIST", alpha=factor(alpha))
  
  df_tnorm <- rbind(df_tnorm, rbind(df2, df12, dfo2, dfo12, dfs, dfadj))
  
  df2_pdist <- melt(sig_count_pdist[, , 1])
  names(df2_pdist) <- c("vertex", "time", "cts")
  df2_pdist$time <- rep(factor(m2, levels=m2), each = n)
  df2_pdist <- cbind(df2_pdist, Method="MASE(2)", alpha=factor(alpha))
  
  df12_pdist <- melt(sig_count_pdist[, , 2])
  names(df12_pdist) <- c("vertex", "time", "cts")
  df12_pdist$time <- rep(factor(m2, levels=m2), each = n)
  df12_pdist <- cbind(df12_pdist, Method="MASE(12)", alpha=factor(alpha))
  
  dfo2_pdist <- melt(sig_count_pdist[, , 3])
  names(dfo2_pdist) <- c("vertex", "time", "cts")
  dfo2_pdist$time <- rep(factor(m2, levels=m2), each = n)
  dfo2_pdist <- cbind(dfo2_pdist, Method="OMNI(2)", alpha=factor(alpha))
  
  dfo12_pdist <- melt(sig_count_pdist[, , 4])
  names(dfo12_pdist) <- c("vertex", "time", "cts")
  dfo12_pdist$time <- rep(factor(m2, levels=m2), each = n)
  dfo12_pdist <- cbind(dfo12_pdist, Method="OMNI(12)", alpha=factor(alpha))
  
  dfs_pdist <- melt(sig_count_pdist[, , 5])
  names(dfs_pdist) <- c("vertex", "time", "cts")
  dfs_pdist$time <- rep(factor(m2, levels=m2), each = n)
  dfs_pdist <- cbind(dfs_pdist, Method="SCAN", alpha=factor(alpha))
  
  dfadj_pdist <- melt(sig_count_pdist[, , 6])
  names(dfadj_pdist) <- c("vertex", "time", "cts")
  dfadj_pdist$time <- rep(factor(m2, levels=m2), each = n)
  dfadj_pdist <- cbind(dfadj_pdist, Method="DIST", alpha=factor(alpha))
  
  df_pdist <- rbind(df_pdist, rbind(df2_pdist, df12_pdist, dfo2_pdist, dfo12_pdist, dfs_pdist, dfadj_pdist))
  
}


df_tnorm_null <- df_tnorm %>% filter(time!="5:6"&time!="6:7" & time!="7:8")
df_tnorm_type_one_err <- df_tnorm_null[,c(2,4,5)]
df_tnorm_type_one_err$power <- df_tnorm_type_one_err$power * pvalnmc
# Summarizing the data
df_tnorm_type_one_err <- df_tnorm_type_one_err %>%
  group_by(Method, alpha) %>%
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

method_colors <- c("MASE(12)" = cbbPalette[1], "MASE(2)" = cbbPalette[2], "OMNI(12)" = cbbPalette[3], 
                   "OMNI(2)" = cbbPalette[4], "SCAN" = cbbPalette[5], "DIST" = cbbPalette[6], "GP" = cbbPalette[7])
method_shapes <- c("MASE(12)" = 19, "MASE(2)" = 17, "OMNI(12)" = 15, 
                   "OMNI(2)" = 13, "SCAN" = 12, "DIST" = 8, "GP" = 25)


# Merging and renaming for final graph
df.graph <- bind_rows(df_tnorm_power, df_tnorm_type_one_err)
df.graph$metric <- factor(df.graph$metric, levels = c("power", "FPR"))

cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

method_colors <- c("MASE(12)" = cbbPalette[1], "MASE(2)" = cbbPalette[2], "OMNI(12)" = cbbPalette[3], 
                   "OMNI(2)" = cbbPalette[4], "SCAN" = cbbPalette[5], "DIST" = cbbPalette[6], "GP" = cbbPalette[7])
method_shapes <- c("MASE(12)" = 19, "MASE(2)" = 17, "OMNI(12)" = 15, 
                   "OMNI(2)" = 13, "SCAN" = 12, "DIST" = 8, "GP" = 25)


pg <- ggplot(df.graph, aes(x=alpha, y=mean, color=Method,group=Method)) +
  geom_errorbar(aes(x=alpha,ymin=mean-se, ymax=mean+se, width=.1,linetype=Method, color=Method)) +
  geom_line(aes(x=alpha,y=mean,linetype=Method)) +theme_bw()+xlab(TeX("${\\theta}$"))+
  geom_point(aes(x=alpha,y=mean, shape=Method)) +
  scale_color_manual(values = method_colors) + 
  scale_shape_manual(values = method_shapes) +
  facet_wrap(~metric,nrow = 2,scales = "free_y",strip.position="left")+
  # theme(plot.title = element_text(size=10, face = 'bold'),legend.position = "none",axis.title = element_text(face="bold"),axis.text = element_text(face="bold"))  + ylab("")
  theme(plot.title = element_text(size=10, face = 'bold'),legend.text = element_text(size=8, face='bold'),axis.title = element_text(face="bold"),axis.text = element_text(face="bold"))  + ylab("") #+ 


dfi <- df_rr_combined
dfi$Method <- factor(dfi$Method)
df.auc <- c()
for (checked_mc in seq(nullnmc+1,nullnmc + pvalnmc)) {
  for (checked_alpha in c(0,seq(8)/8)) {
    for (checked_method in c("MASE(2)","OMNI(2)","MASE(12)","OMNI(12)","SCAN", "DIST")) {#"UASE(2)", "UASE(12)")) {#)) {
      tempdf <- dfi %>% filter(time=="6:7"&mc==checked_mc&alpha==checked_alpha&Method==checked_method)
      temp <- data.frame(auc=pROC::roc(tempdf$vertex<=100, tempdf$rr, quiet=TRUE)$auc, mc=checked_mc, Method=checked_method, alpha=factor(checked_alpha))
      df.auc <- rbind(df.auc, temp)
    }
  }
}

dfommasi.auc <- df.auc %>%
  group_by(Method, alpha) %>%
  select(-2) %>%
  summarize(
    across(
      everything(),
      list(mean = mean, sd = sd, se = ~ sd(.)/sqrt(n()))
    )
  ) %>%
  rename(AUC = auc_mean, Method = Method, sd=auc_sd, se=auc_se)

mycol <- gg_color_hue(2)[2]


dfpert <- dfi %>% filter(vertex<=100&time=="6:7")
dfnonpert <- dfi %>% filter(vertex>100&time=="6:7")
dfpert <- cbind(dfpert, rep("pert",dim(dfpert)[1]))
colnames(dfpert)[6+1] <- "vertex.class"
dfnonpert <- cbind(dfnonpert, rep("nonpert",dim(dfnonpert)[1]))
colnames(dfnonpert)[6+1] <- "vertex.class"
dfall <- rbind(dfpert,dfnonpert)
dfall$vertex.class <- factor(dfall$vertex.class, c("nonpert", "pert"),ordered = TRUE)
a = dfall[, -(1:2)] %>%
  group_by(alpha, mc, Method, vertex.class) %>%
  summarize(across(everything(), mean))

dfa = a %>%
  group_by(alpha, mc, Method) %>%
  mutate(rr_Diff = c(NA, diff(rr)))

dfa = dfa %>%
  filter(vertex.class == "pert")

dfommai = dfa %>%
  group_by(Method, alpha) %>%
  select(-c(2,4,5)) %>%
  summarize(
    across(
      everything(),
      list(mean = mean, sd = sd, se = ~ sd(.) / sqrt(n()))
    )
  ) %>%
  rename(mrr = rr_Diff_mean, Method = Method, sd=rr_Diff_sd, se=rr_Diff_se)

mycol <- gg_color_hue(2)[2]

colnames(dfommai)[3] <- "mean"
colnames(dfommasi.auc)[3] <- "mean"
dfommai$metric <- "reciprocal rank difference"
dfommasi.auc$metric <- "AUC"
df.vertex <- rbind(dfommai, dfommasi.auc)
df.vertex$Method <- factor(df.vertex$Method, levels = c("MASE(12)","MASE(2)","OMNI(12)","OMNI(2)","SCAN","DIST","GP"))



pv <- ggplot(df.vertex, aes(x=alpha, y=mean, color=Method,group=Method)) +
  geom_errorbar(aes(x=alpha,ymin=mean-se, ymax=mean+se, width=.1,linetype=Method, color=Method)) +
  facet_wrap(~metric,nrow = 2,scales = "free_y",strip.position="left")+
  geom_line(aes(x=alpha,y=mean,linetype=Method)) +theme_bw()+xlab(TeX("$\\theta$"))+
  geom_point(aes(x=alpha,y=mean, shape=Method)) + 
  scale_color_manual(values = method_colors) + 
  scale_shape_manual(values = method_shapes) +
  # theme(plot.title = element_text(size=10, face = 'bold'), legend.position = "none",axis.title = element_text(face="bold"),axis.text = element_text(face="bold"))  + ylab("")#+
theme(plot.title = element_text(size=10, face = 'bold'), legend.text = element_text(size=8, face='bold'),axis.title = element_text(face="bold"),axis.text = element_text(face="bold"))  + ylab("")#+



# Filter the data to exclude certain time periods from df_pdist_combined
df_pdist_null <- df_pdist %>% 
  filter(!(time %in% c("5:6", "6:7", "7:8")))

# Calculate type one error rate
df_pdist_type_one_err <- df_pdist_null %>% 
  transmute(time, Method, alpha, fpr = cts / pvalnmc) %>% 
  group_by(Method, alpha) %>% 
  summarize(fpr = mean(fpr), se = sqrt(fpr * (1 - fpr) / pvalnmc), .groups = 'drop')

# Calculating power for a specific time period in df_pdist_combined
df_pdist_power <- df_pdist %>% 
  filter((vertex<=100)&(time=="6:7")) %>% 
  transmute(time, Method, alpha, power = cts / pvalnmc) %>% 
  group_by(Method, alpha) %>% 
  summarize(power = mean(power), se = sqrt(power * (1 - power) / pvalnmc), .groups = 'drop')


colnames(df_pdist_type_one_err)[3] <- "mean"
colnames(df_pdist_power)[3] <- "mean"
df_pdist_power$metric <- "power"
df_pdist_type_one_err$metric <- "FPR"
df_vertex_test <- rbind(df_pdist_power, df_pdist_type_one_err)
df_vertex_test$metric <- factor(df_vertex_test$metric, levels = c("power", "FPR"))

pv_test_errbar <- ggplot(df_vertex_test, aes(x=alpha, y=mean, color=Method,group=Method)) +
  geom_errorbar(aes(x=alpha,ymin=mean-se, ymax=mean+se, width=.1,linetype=Method, color=Method)) +
  facet_wrap(~metric,nrow = 2,scales = "free_y",strip.position="left")+
  geom_line(aes(x=alpha,y=mean,linetype=Method)) +theme_bw()+xlab(TeX("$\\theta$"))+
  geom_point(aes(x=alpha,y=mean, shape=Method)) +
  scale_color_manual(values = method_colors) +
  scale_shape_manual(values = method_shapes) +
  theme(plot.title = element_text(size=10, face = 'bold'), legend.position = "none",axis.title = element_text(face="bold"),axis.text = element_text(face="bold"))  + ylab("")#+
# theme(plot.title = element_text(size=10, face = 'bold'), legend.text = element_text(size=8, face='bold'),axis.title = element_text(face="bold"),axis.text = element_text(face="bold"))  + ylab("")#+

df_v_fpr <- df_pdist %>% 
  filter(!(time %in% c("5:6", "6:7", "7:8"))) %>%
  mutate(metric = "FPR",
         value = cts/pvalnmc) # This should be adjusted to how you calculate FPR

# Calculate Power (for an experimental scenario)
df_v_power <- df_pdist %>% 
  filter((vertex<=100)&(time=="6:7")) %>%
  mutate(metric = "power",
         value = cts/pvalnmc)

df_vertex_test_violin <- rbind(df_v_fpr, df_v_power)
df_vertex_test_violin$Method <- factor(df_vertex_test_violin$Method, levels = c("MASE(12)","MASE(2)","OMNI(12)","OMNI(2)","SCAN","DIST","GP"))

pv_test_violin <- ggplot(df_vertex_test_violin, aes(x=alpha, y=value, fill=Method, color=Method)) +
  geom_violin() + # Add violin plot
  facet_wrap(~metric, nrow = 2, scales = "free_y", strip.position = "left") +
  theme_bw() +
  xlab(TeX("$\\theta$")) +
  scale_fill_manual(values = method_colors) + # Use fill for violin color
  scale_color_manual(values = method_colors) + # Adjust color for consistency
  # scale_shape_manual(values = method_shapes) +
  theme(plot.title = element_text(size=10, face = 'bold'), 
        # legend.position = "none",
        legend.text = element_text(size=8, face='bold'),
        axis.title = element_text(face="bold"),
        axis.text = element_text(face="bold")) +
  ylab("")


dfgp_larger <- c()
dfgp_smaller <- c()
df_lists <- list(
  mu_larger = list(), mu_smaller = list(),
  sigma_larger = list(), sigma_smaller = list(),
  pi_larger = list(), pi_smaller = list(),
  H_larger = list(), H_smaller = list()
)

pvalnmc <- 20
for(alpha in c(0,seq(8)/8)) {#c(0)) {#,seq(8)/8)) {
  print(alpha)
  # print(dfgp_smaller)
  trainG <- list()
  for(i in 1:nullnmc) {
    set.seed(123+i-1)
    glist <- genTSGonepar(n, nperturb, cperturb=0, rmin, rmax, tmax, d_true,alpha )
    trainG <- c(trainG, lapply(glist[-c(6,7)], function(x) as_adj(x)))
  }
  testG <- list()
  for(i in 1:pvalnmc){
    set.seed(123+nullnmc+i-1)
    glist <- genTSGonepar(n, nperturb, cperturb, rmin, rmax, tmax, d_true,alpha)
    testG <- c(testG, lapply(glist, function(x) as_adj(x)))
  }
  
  Y <- c(rep(replace(rep(-1, 12), c(6,7),1),pvalnmc), rep(-1, length(trainG)))
  # Unbalanced split ----
  split_index <- c(rep(rep("validate", 12),pvalnmc), rep("train", length(trainG))) # OCC
  
  G <- c(testG,trainG)
  
  # GP-F ----
  # GP with squared-exponential kernel and Frobenius distance
  D <- dist.frobenius(list(G))
  a <- c(10, 10)
  b <- c(50, 50)
  
  ns <- 2500
  nchains <- 1
  m.train <- sum(split_index == "train")
  p <- dim(D)[3]
  params <- list(m.train + p + 2)
  for (i in 1:m.train) {
    params[i] <- paste0("f", i)
  }
  params[m.train + 1] <- "sigma"
  for (pp in 1:p) {
    params[m.train + 1 + pp]  <- paste0("l", pp)
  }
  params[m.train + p + 2] <- "lp.f"
  sims <- mcmc_array(ns, nchains = nchains, params)
  
  for (ic in 1:nchains) {
    temp_sims <- t(gp.class(dist = D, Y = Y, split_index = split_index
                            , a = a, b = b, ns = ns, monitor = FALSE))
    sims[, ic, ] <- temp_sims
  }
  pred.GP_F <- gp.class_pred(sims[, 1 ,], D, split_index, p
                             , burn = .2, thin = 10, avg = TRUE)
  occ.GP_F <- gp.occ(sims[,1,], D, burn = .2, thin = 10)
  occ.scores <- occ.GP_F
  gp_mu <- (occ.scores$mu)
  mu_thres <- mean(sort(occ.scores$mu)[c(which.max(diff(sort(occ.scores$mu)))
                                      , which.max(diff(sort(occ.scores$mu))) + 1)])
  
  gp_pi <- (occ.scores$pi)
  pi_thres <- mean(sort(occ.scores$pi)[c(which.max(diff(sort(occ.scores$pi)))
                                            , which.max(diff(sort(occ.scores$pi))) + 1)])
  gp_sigma <- (occ.scores$sigma)
  sigma_thres <- mean(sort(occ.scores$sigma)[c(which.max(diff(sort(occ.scores$sigma)))
                                , which.max(diff(sort(occ.scores$sigma))) + 1)])
  
  gp_H <- occ.scores$H
  H_thres <- mean(sort(occ.scores$H)[c(which.max(diff(sort(occ.scores$H)))
                            , which.max(diff(sort(occ.scores$H))) + 1)])
  thresholds <- list(
    mu = mean(sort(gp_mu)[c(which.max(diff(sort(gp_mu))), which.max(diff(sort(gp_mu))) + 1)]),
    sigma = mean(sort(gp_sigma)[c(which.max(diff(sort(gp_sigma))), which.max(diff(sort(gp_sigma))) + 1)]),
    pi = mean(sort(gp_pi)[c(which.max(diff(sort(gp_pi))), which.max(diff(sort(gp_pi))) + 1)]),
    H = mean(sort(gp_H)[c(which.max(diff(sort(gp_H))), which.max(diff(sort(gp_H))) + 1)])
  )
  
  count_list <- list(
    mu_larger = rep(0, 12),
    sigma_larger = rep(0, 12),
    pi_larger = rep(0, 12),
    H_larger = rep(0, 12),
    mu_smaller = rep(0, 12),
    sigma_smaller = rep(0, 12),
    pi_smaller = rep(0, 12),
    H_smaller = rep(0, 12)
  )
  
  # Define the number of elements to pick in each iteration
  chunk_size <- 12
  
  # Loop through Y in chunks of size 12
  c2_larger <- rep(0,12)
  c2_smaller <- rep(0,12)
  # dfgp_larger <- c()
  for(i in seq(1, length(Y[split_index == "validate"]), by = chunk_size)) {
    chunk_end <- min(i + chunk_size - 1, length(Y))  # Ensure the chunk does not exceed the length of Y

    c2_larger <- c2_larger + (gp_mu[i:chunk_end]>mu_thres|gp_pi[i:chunk_end]>pi_thres|gp_sigma[i:chunk_end]>sigma_thres|gp_H[i:chunk_end]>H_thres)
    c2_smaller <- c2_smaller + (gp_mu[i:chunk_end]<mu_thres|gp_pi[i:chunk_end]<pi_thres|gp_sigma[i:chunk_end]<sigma_thres|gp_H[i:chunk_end]<H_thres)
    
    # For mu
    count_list$mu_larger <- count_list$mu_larger + (gp_mu[i:chunk_end] > thresholds$mu)
    count_list$mu_smaller <- count_list$mu_smaller + (gp_mu[i:chunk_end] < thresholds$mu)
    
    # For sigma
    count_list$sigma_larger <- count_list$sigma_larger + (gp_sigma[i:chunk_end] > thresholds$sigma)
    count_list$sigma_smaller <- count_list$sigma_smaller + (gp_sigma[i:chunk_end] < thresholds$sigma)
    
    # For pi
    count_list$pi_larger <- count_list$pi_larger + (gp_pi[i:chunk_end] > thresholds$pi)
    count_list$pi_smaller <- count_list$pi_smaller + (gp_pi[i:chunk_end] < thresholds$pi)
    
    # For H
    count_list$H_larger <- count_list$H_larger + (gp_H[i:chunk_end] > thresholds$H)
    count_list$H_smaller <- count_list$H_smaller + (gp_H[i:chunk_end] < thresholds$H)
    
  }
  temp <- data.frame(time=factor(tvec, levels=tvec),power=c2_larger/pvalnmc, se = sqrt( (c2_larger/pvalnmc) * (1 - c2_larger/pvalnmc)/pvalnmc ) ,type="GP", alpha=factor(alpha))
  dfgp_larger <- rbind(dfgp_larger, temp)
  temp <- data.frame(time=factor(tvec, levels=tvec),power=c2_smaller/pvalnmc, se = sqrt( (c2_smaller/pvalnmc) * (1 - c2_smaller/pvalnmc)/pvalnmc ) ,type="GP", alpha=factor(alpha))
  dfgp_smaller <- rbind(dfgp_smaller, temp)
  
  for(param in c("mu", "sigma", "pi", "H")) {
    for(thres_type in c("larger", "smaller")) {
      param_thres <- paste0(param, "_", thres_type)
      count_val <- count_list[[param_thres]]
      temp <- data.frame(
        time = factor(tvec, levels = tvec),
        power = count_val / pvalnmc,
        se = sqrt((count_val / pvalnmc) * (1 - count_val / pvalnmc) / pvalnmc),
        type = "GP",
        alpha = factor(alpha)
      )
      df_lists[[param_thres]] <- append(df_lists[[param_thres]], list(temp))
    }
  }
  
}

# Combine results for mu_larger
df_mu_larger <- do.call(rbind, df_lists$mu_larger)

# Combine results for mu_smaller
df_mu_smaller <- do.call(rbind, df_lists$mu_smaller)

# Combine results for sigma_larger
df_sigma_larger <- do.call(rbind, df_lists$sigma_larger)

# Combine results for sigma_smaller
df_sigma_smaller <- do.call(rbind, df_lists$sigma_smaller)

# Combine results for pi_larger
df_pi_larger <- do.call(rbind, df_lists$pi_larger)

# Combine results for pi_smaller
df_pi_smaller <- do.call(rbind, df_lists$pi_smaller)

# Combine results for H_larger
df_H_larger <- do.call(rbind, df_lists$H_larger)

# Combine results for H_smaller
df_H_smaller <- do.call(rbind, df_lists$H_smaller)

# List all the data frames for easy looping
df_list <- list(
  mu_larger = df_mu_larger,
  mu_smaller = df_mu_smaller,
  sigma_larger = df_sigma_larger,
  sigma_smaller = df_sigma_smaller,
  pi_larger = df_pi_larger,
  pi_smaller = df_pi_smaller,
  H_larger = df_H_larger,
  H_smaller = df_H_smaller
)

dfgp <- c()
# Loop through each alpha
for (theta in c(0, seq(8)/8)) {
  
  max_power_for_alpha <- -Inf  # Initialize with a very low value
  best_statistic_for_alpha <- NA
  
  # Loop through each dataframe in df_list
  for (name in names(df_list)) {
    df <- df_list[[name]]
    
    temp <- df %>% filter(alpha == theta)
    
    if ((temp%>%filter(time==1))$power>0.05) {
          temp$power <- 1-temp$power
        }
    
    # Calculate average power at times 6 and 7
    avg_power_6_7 <- mean((temp %>% filter(time %in% c(6, 7)))$power)
    
    # Update max power and best statistic for the current alpha, if necessary
    if (avg_power_6_7 > max_power_for_alpha) {
      max_power_for_alpha <- avg_power_6_7
      best_statistic_for_alpha <- temp
    }
  }
  
  # Append the results for current alpha to the results dataframe
  dfgp <- rbind(dfgp, best_statistic_for_alpha)
}

# Summarizing data for type one error
dfgp_type_one_err <- dfgp %>%
  filter(time != "6" & time != "7") %>%
  select(2, 4, 5) %>%
  group_by(type, alpha) %>%
  summarize(
    across(
      everything(),
      list(
        fpr = ~ mean(.),
        se = ~ sqrt(mean(.) * (1 - mean(.)) / pvalnmc)
      )
    )
  ) %>%
  rename(Method = type, mean = power_fpr, se=power_se)

dfgp_type_one_err$metric <- "FPR"
dfgp_type_one_err <- select(dfgp_type_one_err, colnames(df_tnorm_type_one_err))

# Summarizing data for power
dfgp_power <- dfgp %>%
  filter(time == "6" | time == "7") %>%
  select(2, 3, 4, 5) %>%
  group_by(type, alpha) %>%
  summarize(
    across(
      everything(),
      list(power = ~sum(.))
    )
  ) %>%
  rename(Method = type, mean = power_power, se = se_power)

dfgp_power$time <- "6:7"
dfgp_power$metric <- "power"
dfgp_power <- select(dfgp_power, colnames(df_tnorm_type_one_err))

# Merging and renaming for final graph
df.graph <- bind_rows(df_tnorm_power, df_tnorm_type_one_err, dfgp_power, dfgp_type_one_err)
df.graph$metric <- factor(df.graph$metric, levels = c("power", "FPR"))
df.graph$Method <- str_replace_all(df.graph$Method, c("MASE12" = "MASE(12)", "MASE2" = "MASE(2)", "OMNI12" = "OMNI(12)", "OMNI2" = "OMNI(2)"))
df.graph$Method <- factor(df.graph$Method, levels = c("MASE(12)", "MASE(2)", "OMNI(12)", "OMNI(2)", "SCAN", "DIST", "GP"))


cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

method_colors <- c("MASE(12)" = cbbPalette[1], "MASE(2)" = cbbPalette[2], "OMNI(12)" = cbbPalette[3], 
                   "OMNI(2)" = cbbPalette[4], "SCAN" = cbbPalette[5], "DIST" = cbbPalette[6], "GP" = cbbPalette[7])
method_shapes <- c("MASE(12)" = 19, "MASE(2)" = 17, "OMNI(12)" = 15, 
                   "OMNI(2)" = 13, "SCAN" = 12, "DIST" = 8, "GP" = 25)


pg <- ggplot(df.graph, aes(x=alpha, y=mean, color=Method,group=Method)) +
  geom_errorbar(aes(x=alpha,ymin=mean-se, ymax=mean+se, width=.1,linetype=Method, color=Method)) +
  geom_line(aes(x=alpha,y=mean,linetype=Method)) +theme_bw()+xlab(TeX("${\\theta}$"))+
  geom_point(aes(x=alpha,y=mean, shape=Method)) +
  scale_color_manual(values = method_colors) + 
  scale_shape_manual(values = method_shapes) +
  facet_wrap(~metric,nrow = 2,scales = "free_y",strip.position="left")+
  # theme(plot.title = element_text(size=10, face = 'bold'),legend.position = "none",axis.title = element_text(face="bold"),axis.text = element_text(face="bold"))  + ylab("")
  theme(plot.title = element_text(size=10, face = 'bold'),legend.text = element_text(size=8, face='bold'),axis.title = element_text(face="bold"),axis.text = element_text(face="bold"))  + ylab("") #+ 
  


# Figure 5

grid.arrange(pg, pv_test_errbar, nrow=1, widths=c(5,4))
png(paste0("comb","_",format(Sys.time(), "%Y_%m_%d.png")), width = 9, height = 4, units="in", res=400)
grid.arrange(pg, pv_test_errbar, nrow=1, widths=c(5,4))
dev.off()

# Figure 18
pv
png(paste0("vertex_auc_mrr","_",format(Sys.time(), "%Y_%m_%d.png")), width = 5, height = 4, units="in", res=400)
pv
dev.off()

# Figure 19
pv_test_violin
png(paste0("pv_test_violin","_",format(Sys.time(), "%Y_%m_%d.png")), width = 8, height = 4, units="in", res=400)
pv_test_violin
dev.off()

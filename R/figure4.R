# R version 4.3.1
# Install and load the packages if not already installed

packages <- c("igraph", "gtools", "tidyverse", "irlba", "doParallel", "rstudioapi",
              "qcc", "latex2exp", "gridExtra", "reshape2","metafolio")

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


cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
colors =  cbbPalette[c(7, 8, 4, 6)]


#same dimension case
d_true <- 4
#p-values
nullnmc <- 50
pvalnmc <- 200
tmax <- 12
tvec <- 1:tmax
m2 <- names(running(tvec, width=2))
n <- 400
nperturb <- 100
rmin <- 0.2
rmax <- 0.8
center <- FALSE
cperturb <- .22

approx <- FALSE
center <- FALSE

dSVD=2
dASE=NA#
dSVD12=NA#9#NA
method3 = 2
alpha = .875

# Initialize output data frames
# out_xx represent the empirical null distributions
out_tnorm <- c()
out_pdist <- c()

# df_xx represent the empirical p-val distributions
df_tnorm <- c()
df_pdist <- c()
attr_weight <-  NULL
norm_choice <- "F"
d.max <- "sqrt"
plot <- FALSE
elbow_graph <- 1
# Combined Monte Carlo Simulation
for(i in 1:(nullnmc + pvalnmc)) {
  set.seed(123 + i - 1)
  glist <- genTSGdd(n, nperturb, cperturb * as.numeric(i > nullnmc), rmin, rmax, tmax, d_true, alpha)
  
  # Do MASE for different dimensions
  # doMase <- function(glist, latpos.list=NULL, nmase=2, dSVD=NA, dASE=NA, d.max="full", center=FALSE, approx=TRUE, elbow_graph=1, plot=TRUE,latent.form=3, attr_weight=NULL, norm_choice="F")
  
  out2 <- doMase(glist, latpos.list=NULL, nmase=2, dSVD=dSVD, dASE=dASE, d.max=d.max, center=center, approx=approx, elbow_graph=elbow_graph, plot=plot, latent.form=method3, attr_weight=attr_weight, norm_choice=norm_choice)
  out12 <- doMase(glist, latpos.list=NULL, nmase=12, dSVD=dSVD12, dASE=dASE, d.max=d.max, center=center, approx=approx, elbow_graph=elbow_graph, plot=plot, latent.form=method3, attr_weight=attr_weight, norm_choice=norm_choice)
  outs <- doScan(glist)
  
  # Process for tnorm and pdist if within the first nullnmc simulations
  if (i <= nullnmc) {
    # Process tnorm data
    out.p2_tnorm <- data.frame(as.vector(t(out2$tnorm[-c(5,6,7)])), Method="MASE(2)")
    out.p12_tnorm <- data.frame(as.vector(t(out12$tnorm[-c(5,6,7)])), Method="MASE(12)")
    out.ps_tnorm <- data.frame(as.vector(t(outs$tnorm[-c(5,6,7)])), Method="SCAN")
    colnames(out.p2_tnorm)[1] <- "tnorm"
    colnames(out.p12_tnorm)[1] <- "tnorm"
    colnames(out.ps_tnorm)[1] <- "tnorm"
    out_tnorm <- rbind(out_tnorm, rbind(out.p2_tnorm, out.p12_tnorm, out.ps_tnorm))
    
    # Process pdist data
    out.p2_pdist <- data.frame(t(out2$pdist[,-c(5,6,7)]), Method="MASE(2)")
    out.p12_pdist <- data.frame(t(out12$pdist[,-c(5,6,7)]), Method="MASE(12)")
    out.ps_pdist <- data.frame(t(outs$pdist[,-c(5,6,7)]), Method="SCAN")
    
    out_pdist <- rbind(out_pdist, rbind(out.p2_pdist, out.p12_pdist, out.ps_pdist))
  }
  
  # Process for p-values if beyond the first nullnmc simulations
  if (i > nullnmc && i <= (nullnmc + pvalnmc)) {
    # Process tnorm p-values
    pvalues2_tnorm <- rep(0, (tmax-1))
    pvalues12_tnorm <- rep(0, (tmax-1))
    pvaluess_tnorm <- rep(0, (tmax-1))
    for (j in 1:(tmax-1)) {
      pvalues2_tnorm[j] <- sum(out2$tnorm[j] < (out_tnorm[which(out_tnorm$Method == "MASE(2)"), 1])) / (nullnmc*8)
      pvalues12_tnorm[j] <- sum(out12$tnorm[j] < (out_tnorm[which(out_tnorm$Method == "MASE(12)"), 1])) / (nullnmc*8)
      pvaluess_tnorm[j] <- sum(outs$tnorm[j] < (out_tnorm[which(out_tnorm$Method == "SCAN"), 1])) / (nullnmc*8)
    }
    # Adjust p-values for tnorm
    pvalues2_tnorm <- p.adjust(pvalues2_tnorm, "BH")
    pvalues12_tnorm <- p.adjust(pvalues12_tnorm, "BH")
    pvaluess_tnorm <- p.adjust(pvaluess_tnorm, "BH")
    
    # Add to df_tnorm
    df2_tnorm <- data.frame(time=factor(m2, levels=m2), pvalues=pvalues2_tnorm, mc=i, Method="MASE(2)")
    df12_tnorm <- data.frame(time=factor(m2, levels=m2), pvalues=pvalues12_tnorm, mc=i, Method="MASE(12)")
    dfs_tnorm <- data.frame(time=factor(m2, levels=m2), pvalues=pvaluess_tnorm, mc=i, Method="SCAN")
    df_tnorm <- rbind(df_tnorm, rbind(df2_tnorm, df12_tnorm, dfs_tnorm))
    
    # Process pdist p-values
    # Assuming similar logic for pdist p-values
    # Calculate p-values for pdist
    pvalues2_pdist <- matrix(0, nrow = n, ncol = tmax - 1)
    pvalues12_pdist <- matrix(0, nrow = n, ncol = tmax - 1)
    pvaluess_pdist <- matrix(0, nrow = n, ncol = tmax - 1)
    
    for (j in 1:(tmax - 1)) {
      for (w in 1:n) {
        pvalues2_pdist[w, j] <- sum(out2$pdist[w, j] < (out_pdist[which(out_pdist$Method == "MASE(2)"), w])) / (nullnmc*8)
        pvalues12_pdist[w, j] <- sum(out12$pdist[w, j] < (out_pdist[which(out_pdist$Method == "MASE(12)"), w])) / (nullnmc*8)
        pvaluess_pdist[w, j] <- sum(outs$pdist[w, j] < (out_pdist[which(out_pdist$Method == "SCAN"), w])) / (nullnmc*8)
      }
    }
    
    # Adjust p-values for pdist
    pvalues2_pdist_adjusted <- matrix(p.adjust(pvalues2_pdist, "BH"),n, tmax-1)
    pvalues12_pdist_adjusted <- matrix(p.adjust(pvalues12_pdist, "BH"),n, tmax-1)
    pvaluess_pdist_adjusted <- matrix(p.adjust(pvaluess_pdist, "BH"),n, tmax-1)
    
    df2_pdist <- melt(pvalues2_pdist_adjusted)
    names(df2_pdist) <- c("vertex", "time", "pvalues")
    df2_pdist$time <- rep(factor(m2, levels=m2), each = n)
    df2_pdist <- cbind(df2_pdist,mc=i, Method="MASE(2)")
    
    df12_pdist <- melt(pvalues12_pdist_adjusted)
    names(df12_pdist) <- c("vertex", "time", "pvalues")
    df12_pdist$time <- rep(factor(m2, levels=m2), each = n)
    df12_pdist <- cbind(df12_pdist,mc=i, Method="MASE(12)")
    
    dfs_pdist <- melt(pvaluess_pdist_adjusted)
    names(dfs_pdist) <- c("vertex", "time", "pvalues")
    dfs_pdist$time <- rep(factor(m2, levels=m2), each = n)
    dfs_pdist <- cbind(dfs_pdist,mc=i, Method="SCAN")
    
    df_pdist <- rbind(df_pdist, rbind(df2_pdist, df12_pdist, dfs_pdist))
  }
}
df <- df_tnorm[,-3]
df <- df %>%
  group_by(time, Method) %>%
  summarize(
    across(
      everything(),
      list(
        median = ~median(.),
        lci = ~sort(.)[round(n()/2 + qnorm(0.025/(tmax-1)) * sqrt(n()/4))],
        uci = ~sort(.)[round(n()/2 - qnorm(0.025/(tmax-1)) * sqrt(n()/4))]
      )
    )
  )
dfmdd <- df %>% rename(pvalues = pvalues_median, lci = pvalues_lci, uci = pvalues_uci)

mycol <- gg_color_hue(2)[2]

df <- df_pdist[,-4]
df <- df %>%
  group_by(vertex, time, Method) %>%
  summarize(
    across(
      everything(),
      list(
        median = ~median(.),
        lci = ~sort(.)[round(n()/2 + qnorm(0.025/(tmax-1)) * sqrt(n()/4))],
        uci = ~sort(.)[round(n()/2 - qnorm(0.025/(tmax-1)) * sqrt(n()/4))]
      )
    )
  )
dfmidd <- df %>% rename(pvalues = pvalues_median, lci = pvalues_lci, uci = pvalues_uci)
mycol <- gg_color_hue(2)[2]

colnames(dfmdd)[2] <- "type"
dfmdd <- ungroup(dfmdd)
dfommadd <- dfmdd
dfommawoom12dd <- dfommadd 
dfommawoom12dd <- dfommawoom12dd  #%>% filter(type!="SCAN")
pddo <- ggplot(dfommawoom12dd, aes(x=time, y=pvalues, group=type)) + 
  geom_errorbar(aes(ymin=lci, ymax=uci), width=.1, color=mycol) +
  facet_wrap(~type,nrow = 3,strip.position ="left", scales = "free") +
  geom_line(color="grey70") +ylab("Adjusted p-value")+
  geom_point(data=dfommawoom12dd%>%filter(time=="6:7"&pvalues<0.05), color="red",shape=17) +
  geom_point(data=dfommawoom12dd%>%filter(!(time=="6:7"&pvalues<0.05)), color="grey",shape=20) +
  theme_bw()+ theme(text = element_text(size=12, face='bold'), strip.text = element_text(size = 12)) + geom_hline(aes(yintercept=.05))



# fdr control on the pvalues in vertices control charts diff dimension
dfmidd <- dfmidd
#only mase
colnames(dfmidd)[3] <- "type"
dfmidd <- ungroup(dfmidd)
dfmaseip <- dfmidd
dfmasewoom12i <- dfmaseip %>% filter(time=="2:3"|time=="6:7") #%>% filter(type!="SCAN")
pddi <- ggplot() +   geom_errorbar(data=dfmasewoom12i, aes(ymin=lci, ymax=uci,x=vertex), width=.5, color=mycol) +
  geom_point(data=dfmasewoom12i %>%filter(!vertex %in% c(1:(nperturb))),aes(x=vertex,y=pvalues),alpha=1, size=.5,color="grey70",shape=20)+ 
  geom_hline(data=dfmasewoom12i,aes(yintercept=.05))+geom_vline(data=dfmasewoom12i,aes(xintercept=100), linetype="dashed",color="red")+
  geom_point(data=dfmasewoom12i %>%filter(vertex %in% c(1:(nperturb))),aes(x=vertex,y=pvalues), color="red", alpha=1,size=.5,shape=17)+
  facet_grid(type~time,switch = "both",scales = "free")+theme_bw() +ylab("Adjusted p-value")+
  theme(legend.position = "none",plot.title = element_text(size=12),plot.subtitle = element_text(hjust = 0.5), text = element_text(size=12, face='bold'),strip.text = element_text(size = 12))



alpha = 1-.875

# Initialize output data frames
# out_xx represent the empirical null distributions
out_tnorm <- c()
out_pdist <- c()

# df_xx represent the empirical p-val distributions
df_tnorm <- c()
df_pdist <- c()
# Combined Monte Carlo Simulation
for(i in 1:(nullnmc + pvalnmc)) {
  set.seed(123 + i - 1)
  glist <- genTSGdd(n, nperturb, cperturb * as.numeric(i > nullnmc), rmin, rmax, tmax, d_true, alpha)
  
  # Do MASE for different dimensions
  # doMase <- function(glist, latpos.list=NULL, nmase=2, dSVD=NA, dASE=NA, d.max="full", center=FALSE, approx=TRUE, elbow_graph=1, plot=TRUE,latent.form=3, attr_weight=NULL, norm_choice="F")
    
  out2 <- doMase(glist, latpos.list=NULL, nmase=2, dSVD=dSVD, dASE=dASE, d.max=d.max, center=center, approx=approx, elbow_graph=elbow_graph, plot=plot, latent.form=method3, attr_weight=attr_weight, norm_choice=norm_choice)
  out12 <- doMase(glist, latpos.list=NULL, nmase=12, dSVD=dSVD12, dASE=dASE, d.max=d.max, center=center, approx=approx, elbow_graph=elbow_graph, plot=plot, latent.form=method3, attr_weight=attr_weight, norm_choice=norm_choice)
  outs <- doScan(glist)
  
  # Process for tnorm and pdist if within the first nullnmc simulations
  if (i <= nullnmc) {
    # Process tnorm data
    out.p2_tnorm <- data.frame(as.vector(t(out2$tnorm[-c(5,6,7)])), Method="MASE(2)")
    out.p12_tnorm <- data.frame(as.vector(t(out12$tnorm[-c(5,6,7)])), Method="MASE(12)")
    out.ps_tnorm <- data.frame(as.vector(t(outs$tnorm[-c(5,6,7)])), Method="SCAN")
    colnames(out.p2_tnorm)[1] <- "tnorm"
    colnames(out.p12_tnorm)[1] <- "tnorm"
    colnames(out.ps_tnorm)[1] <- "tnorm"
    out_tnorm <- rbind(out_tnorm, rbind(out.p2_tnorm, out.p12_tnorm, out.ps_tnorm))
    
    # Process pdist data
    out.p2_pdist <- data.frame(t(out2$pdist[,-c(5,6,7)]), Method="MASE(2)")
    out.p12_pdist <- data.frame(t(out12$pdist[,-c(5,6,7)]), Method="MASE(12)")
    out.ps_pdist <- data.frame(t(outs$pdist[,-c(5,6,7)]), Method="SCAN")
    
    out_pdist <- rbind(out_pdist, rbind(out.p2_pdist, out.p12_pdist, out.ps_pdist))
  }
  
  # Process for p-values if beyond the first nullnmc simulations
  if (i > nullnmc && i <= (nullnmc + pvalnmc)) {
    # Process tnorm p-values
    pvalues2_tnorm <- rep(0, (tmax-1))
    pvalues12_tnorm <- rep(0, (tmax-1))
    pvaluess_tnorm <- rep(0, (tmax-1))
    for (j in 1:(tmax-1)) {
      pvalues2_tnorm[j] <- sum(out2$tnorm[j] < (out_tnorm[which(out_tnorm$Method == "MASE(2)"), 1])) / (nullnmc*8)
      pvalues12_tnorm[j] <- sum(out12$tnorm[j] < (out_tnorm[which(out_tnorm$Method == "MASE(12)"), 1])) / (nullnmc*8)
      pvaluess_tnorm[j] <- sum(outs$tnorm[j] < (out_tnorm[which(out_tnorm$Method == "SCAN"), 1])) / (nullnmc*8)
    }
    # Adjust p-values for tnorm
    pvalues2_tnorm <- p.adjust(pvalues2_tnorm, "BH")
    pvalues12_tnorm <- p.adjust(pvalues12_tnorm, "BH")
    pvaluess_tnorm <- p.adjust(pvaluess_tnorm, "BH")
    
    # Add to df_tnorm
    df2_tnorm <- data.frame(time=factor(m2, levels=m2), pvalues=pvalues2_tnorm, mc=i, Method="MASE(2)")
    df12_tnorm <- data.frame(time=factor(m2, levels=m2), pvalues=pvalues12_tnorm, mc=i, Method="MASE(12)")
    dfs_tnorm <- data.frame(time=factor(m2, levels=m2), pvalues=pvaluess_tnorm, mc=i, Method="SCAN")
    df_tnorm <- rbind(df_tnorm, rbind(df2_tnorm, df12_tnorm, dfs_tnorm))
    
    # Process pdist p-values
    # Assuming similar logic for pdist p-values
    # Calculate p-values for pdist
    pvalues2_pdist <- matrix(0, nrow = n, ncol = tmax - 1)
    pvalues12_pdist <- matrix(0, nrow = n, ncol = tmax - 1)
    pvaluess_pdist <- matrix(0, nrow = n, ncol = tmax - 1)
    
    for (j in 1:(tmax - 1)) {
      for (w in 1:n) {
        pvalues2_pdist[w, j] <- sum(out2$pdist[w, j] < (out_pdist[which(out_pdist$Method == "MASE(2)"), w])) / (nullnmc*8)
        pvalues12_pdist[w, j] <- sum(out12$pdist[w, j] < (out_pdist[which(out_pdist$Method == "MASE(12)"), w])) / (nullnmc*8)
        pvaluess_pdist[w, j] <- sum(outs$pdist[w, j] < (out_pdist[which(out_pdist$Method == "SCAN"), w])) / (nullnmc*8)
      }
    }
    
    # Adjust p-values for pdist
    pvalues2_pdist_adjusted <- matrix(p.adjust(pvalues2_pdist, "BH"),n, tmax-1)
    pvalues12_pdist_adjusted <- matrix(p.adjust(pvalues12_pdist, "BH"),n, tmax-1)
    pvaluess_pdist_adjusted <- matrix(p.adjust(pvaluess_pdist, "BH"),n, tmax-1)
    
    df2_pdist <- melt(pvalues2_pdist_adjusted)
    names(df2_pdist) <- c("vertex", "time", "pvalues")
    df2_pdist$time <- rep(factor(m2, levels=m2), each = n)
    df2_pdist <- cbind(df2_pdist,mc=i, Method="MASE(2)")
    
    df12_pdist <- melt(pvalues12_pdist_adjusted)
    names(df12_pdist) <- c("vertex", "time", "pvalues")
    df12_pdist$time <- rep(factor(m2, levels=m2), each = n)
    df12_pdist <- cbind(df12_pdist,mc=i, Method="MASE(12)")
    
    dfs_pdist <- melt(pvaluess_pdist_adjusted)
    names(dfs_pdist) <- c("vertex", "time", "pvalues")
    dfs_pdist$time <- rep(factor(m2, levels=m2), each = n)
    dfs_pdist <- cbind(dfs_pdist,mc=i, Method="SCAN")
    
    df_pdist <- rbind(df_pdist, rbind(df2_pdist, df12_pdist, dfs_pdist))
  }
}
df <- df_tnorm[,-3]
df <- df %>%
  group_by(time, Method) %>%
  summarize(
    across(
      everything(),
      list(
        median = ~median(.),
        lci = ~sort(.)[round(n()/2 + qnorm(0.025/(tmax-1)) * sqrt(n()/4))],
        uci = ~sort(.)[round(n()/2 - qnorm(0.025/(tmax-1)) * sqrt(n()/4))]
      )
    )
  )
dfmdd <- df %>% rename(pvalues = pvalues_median, lci = pvalues_lci, uci = pvalues_uci)

mycol <- gg_color_hue(2)[2]

df <- df_pdist[,-4]
df <- df %>%
  group_by(vertex, time, Method) %>%
  summarize(
    across(
      everything(),
      list(
        median = ~median(.),
        lci = ~sort(.)[round(n()/2 + qnorm(0.025/(tmax-1)) * sqrt(n()/4))],
        uci = ~sort(.)[round(n()/2 - qnorm(0.025/(tmax-1)) * sqrt(n()/4))]
      )
    )
  )
dfmidd <- df %>% rename(pvalues = pvalues_median, lci = pvalues_lci, uci = pvalues_uci)
mycol <- gg_color_hue(2)[2]

colnames(dfmdd)[2] <- "type"
dfmdd <- ungroup(dfmdd)
dfommadd <- dfmdd
dfommawoom12dd <- dfommadd 
dfommawoom12dd <- dfommawoom12dd  #%>% filter(type!="SCAN")
psdo <- ggplot(dfommawoom12dd, aes(x=time, y=pvalues, group=type)) + 
  geom_errorbar(aes(ymin=lci, ymax=uci), width=.1, color=mycol) +
  facet_wrap(~type,nrow = 3,strip.position ="left", scales = "free") +
  geom_line(color="grey70") +ylab("Adjusted p-value")+
  geom_point(data=dfommawoom12dd%>%filter(time=="6:7"&pvalues<0.05), color="red",shape=17) +
  geom_point(data=dfommawoom12dd%>%filter(!(time=="6:7"&pvalues<0.05)), color="grey",shape=20) +
  theme_bw()+ theme(text = element_text(size=12, face='bold'), strip.text = element_text(size = 12)) + geom_hline(aes(yintercept=.05))



# fdr control on the pvalues in vertices control charts diff dimension
dfmidd <- dfmidd
#only mase
colnames(dfmidd)[3] <- "type"
dfmidd <- ungroup(dfmidd)
dfmaseip <- dfmidd
dfmasewoom12i <- dfmaseip %>% filter(time=="2:3"|time=="6:7") #%>% filter(type!="SCAN")
psdi <- ggplot() +   geom_errorbar(data=dfmasewoom12i, aes(ymin=lci, ymax=uci,x=vertex), width=.5, color=mycol) +
  geom_point(data=dfmasewoom12i %>%filter(!vertex %in% c(1:(nperturb))),aes(x=vertex,y=pvalues),alpha=1, size=.5,color="grey70",shape=20)+ 
  geom_hline(data=dfmasewoom12i,aes(yintercept=.05))+geom_vline(data=dfmasewoom12i,aes(xintercept=100), linetype="dashed",color="red")+
  geom_point(data=dfmasewoom12i %>%filter(vertex %in% c(1:(nperturb))),aes(x=vertex,y=pvalues), color="red", alpha=1,size=.5,shape=17)+
  facet_grid(type~time,switch = "both",scales = "free")+theme_bw() +ylab("Adjusted p-value")+
  theme(legend.position = "none",plot.title = element_text(size=12),plot.subtitle = element_text(hjust = 0.5), text = element_text(size=12, face='bold'),strip.text = element_text(size = 12))




# Figure 4 a
grid.arrange(psdo,psdi, ncol=1, heights=c(4,4))
png("sdmcomultidn400FDR.png", width = 6, height = 7, units="in", res=400)
grid.arrange(psdo,psdi, ncol=1, heights=c(4,4))
dev.off()

# Figure 4 b
grid.arrange(pddo,pddi, ncol=1, heights=c(4,4))
png("ddmcomultidn400FDR.png", width = 6, height = 7, units="in", res=400)
grid.arrange(pddo,pddi, ncol=1, heights=c(4,4))
dev.off()
# 
# grid.arrange(pddo,pddi, ncol=1, heights=c(4,4))
# png("ddmcomultid4june18n400maseRevisionFDR.png", width = 6, height = 6, units="in", res=400)
# grid.arrange(pddo,pddi, ncol=1, heights=c(4,4))
# dev.off()



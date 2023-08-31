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
alpha = .125
d_true <- 4
#p-values
py = FALSE
method3 = 2
dSVD=5
dASE=4
dSVD12=5
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

approx <- TRUE
center <- FALSE

#pairwise mase
out <- foreach(i = 1:nullnmc, .combine='rbind') %dopar% {
  set.seed(123+i-1)
  glist <- genTSGdd(n, nperturb, cperturb=0, rmin, rmax, tmax, d_true, alpha)
  out2 <- doMase(glist, 2, dSVD,dASE,center, approx, py, method3)
  out12 <- doMase(glist, 12, dSVD12,dASE,center, approx, py, method3)
  out.p2 <- data.frame(as.vector(t(out2$tnorm[-c(5,6,7)])), mase="MASE(2)")
  colnames(out.p2)[1] <- "tnorm"
  out.p12 <- data.frame(as.vector(t(out12$tnorm[-c(5,6,7)])), mase="MASE(12)")
  colnames(out.p12)[1] <- "tnorm"
  out.p <- rbind(out.p2, out.p12)
}
df <- foreach(i = 1:pvalnmc, .combine='rbind') %dopar% {
  set.seed(123+nullnmc+i-1)
  glist <- genTSGdd(n, nperturb, cperturb, rmin, rmax, tmax, d_true, alpha)
  out2 <- doMase(glist, 2, dSVD,dASE,center, approx, py, method3)
  out12 <- doMase(glist, 12, dSVD12,dASE, center, approx, py, method3)
  pvalues2 <- rep(0, (tmax-1))
  pvalues12 <- rep(0, (tmax-1))
  for (j in 1:11) {
    pvalues2[j] <- sum(out2$tnorm[j]<sort(out[which(out$mase=="MASE(2)"),1]))/(nullnmc*8)
    pvalues12[j] <- sum(out12$tnorm[j]<sort(out[which(out$mase=="MASE(12)"),1]))/(nullnmc*8)
  }
  pvalues2 <- p.adjust(pvalues2, "BH")
  pvalues12 <- p.adjust(pvalues12, "BH")
  df2 <- data.frame(time=factor(m2, levels=m2), pvalues=pvalues2, mc=i,mase="MASE(2)")
  df12 <- data.frame(time=factor(m2, levels=m2), pvalues=pvalues12, mc=i,mase="MASE(12)")
  df <- rbind(df2, df12)
  df
}
df <- df[,-3]
df <- df %>%
  group_by(time, mase) %>%
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

dfmsd <- df %>%
  rename(pvalues = pvalues_median, lci = pvalues_lci, uci = pvalues_uci)

mycol <- gg_color_hue(2)[2]


dSVD=5
dASE=4
dSVD12=NA


#pairwise mase
out <- foreach(i = 1:nullnmc, .combine='rbind') %dopar% {
  set.seed(123+i-1)
  glist <- genTSGdd(n, nperturb, cperturb=0, rmin, rmax, tmax, d_true, alpha)
  out2 <- doMase(glist, 2, dSVD,dASE,center, approx, py, method3)
  out12 <- doMase(glist, 12, dSVD12,dASE,center, approx, py, method3)
  #out2$pdist[,1]
  out.p2 <- data.frame(t(out2$pdist[,1]), mase="MASE(2)")
  out.p12 <- data.frame(t(out12$pdist[,1]), mase="MASE(12)")
  out.p <- rbind(out.p2, out.p12)
}

df <- foreach(i = 1:pvalnmc, .combine='rbind') %dopar% {
  set.seed(123+nullnmc+i-1)
  glist <- genTSGdd(n, nperturb, cperturb, rmin, rmax, tmax, d_true, alpha)
  out2 <- doMase(glist, 2, dSVD,dASE,center, approx, py, method3)
  out12 <- doMase(glist, 12, dSVD12, dASE, center, approx, py, method3)
  pvalues2 <- matrix(0, n, tmax-1)#rep(0, 11)
  pvalues12 <- matrix(0, n, tmax-1)
  for (j in 1:(tmax-1)) {
    for (w in 1:n) {
      pvalues2[w,j] <- sum(out2$pdist[w,j]<sort(out[which(out$mase=="MASE(2)"),w]))/(nullnmc)
      pvalues12[w,j] <- sum(out12$pdist[w,j]<sort(out[which(out$mase=="MASE(12)"),w]))/(nullnmc)
    }
    
  }
  pvalues2.adjust <- matrix(p.adjust(pvalues2, "BH"),n, tmax-1)
  pvalues12.adjust <- matrix(p.adjust(pvalues12, "BH"),n, tmax-1)
  df.p2 <- melt(pvalues2.adjust)
  names(df.p2) <- c("vertex", "time", "pvalues")
  df.p2$time <- rep(factor(m2, levels=m2), each = n)
  df.p2 <- cbind(df.p2,mc=i, mase="MASE(2)")
  df.p12 <- melt(pvalues12.adjust)
  names(df.p12) <- c("vertex", "time", "pvalues")
  df.p12$time <- rep(factor(m2, levels=m2), each = n)
  df.p12 <- cbind(df.p12,mc=i, mase="MASE(12)")
  df.p <- rbind(df.p2, df.p12)
  df.p

}
df <- df[,-4]
df <- df %>%
  group_by(vertex, time, mase) %>%
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

dfmisd <- df %>%  rename(pvalues = pvalues_median, lci = pvalues_lci, uci = pvalues_uci)
mycol <- gg_color_hue(2)[2]


colnames(dfmsd)[2] <- "type"
dfmsd <- ungroup(dfmsd)
dfommasd <- dfmsd
dfommawoom12sd <- dfommasd 
dfommawoom12sd <- dfommawoom12sd %>% filter(type!="MASE(12)")
# psdo <- #ggplot(dfommawoom12sd, aes(x=time, y=pvalues, group=type)) + 
  # geom_errorbar(aes(ymin=lci, ymax=uci), width=.1, color=mycol) +
  # facet_wrap(~type,nrow = 2,switch  = "y", scales = "free") +
  # geom_line(color=mycol) +ylab("Adjusted p-value")+
  # geom_point(color=mycol) +theme_bw()+ theme(legend.position = "none",plot.title = element_text(size=25),plot.subtitle = element_text(hjust = 0.5), text = element_text(size=(17)))+
  # geom_hline(aes(yintercept=.05))
psdo <- ggplot(dfommawoom12sd, aes(x=time, y=pvalues, group=type)) + 
  geom_errorbar(aes(ymin=lci, ymax=uci), width=.1, color=mycol) +
  facet_wrap(~type,nrow = 2,strip.position ="left", scales = "free") +
  geom_line(color="grey70") +ylab("Adjusted p-value")+
  geom_point(data=dfommawoom12sd%>%filter(time=="6:7"&pvalues<0.05), color="red",shape=17) +
  geom_point(data=dfommawoom12sd%>%filter(!(time=="6:7"&pvalues<0.05)), color="grey",shape=20) +
  theme_bw()+ theme(text = element_text(size=12, face='bold'), strip.text = element_text(size = 12)) + geom_hline(aes(yintercept=.05))

# fdr control on the pvalues in vertices control charts same dimension
#only mase
colnames(dfmisd)[3] <- "type"
dfmisd <- ungroup(dfmisd)
dfmaseip <- dfmisd
dfmasewoom12i <- dfmaseip %>% filter(time=="2:3"|time=="6:7") %>% filter(type!="MASE(12)")
# psdi <- ggplot() +   geom_errorbar(data=dfmasewoom12i, aes(ymin=lci, ymax=uci,x=vertex), width=.5, color=mycol) +
#   geom_point(data=dfmasewoom12i,aes(x=vertex,y=pvalues),alpha=1, size=.5)+ 
#   geom_hline(data=dfmasewoom12i,aes(yintercept=.05))+geom_vline(data=dfmasewoom12i,aes(xintercept=100), linetype="dashed",color="red")+
#   geom_point(data=dfmasewoom12i %>%filter(vertex %in% c(1:(nperturb))),aes(x=vertex,y=pvalues), color="red", alpha=1,size=.5)+
#   facet_grid(type~time,switch = "both",scales = "free")+theme_bw() +ylab("Adjusted p-value")+
#   theme(legend.position = "none",plot.title = element_text(size=25),plot.subtitle = element_text(hjust = 0.5), text = element_text(size=(17)))
psdi <- ggplot() +   geom_errorbar(data=dfmasewoom12i, aes(ymin=lci, ymax=uci,x=vertex), width=.5, color=mycol) +
  geom_point(data=dfmasewoom12i %>%filter(!vertex %in% c(1:(nperturb))),aes(x=vertex,y=pvalues),alpha=1, size=.5,color="grey70",shape=20)+ 
  geom_hline(data=dfmasewoom12i,aes(yintercept=.05))+geom_vline(data=dfmasewoom12i,aes(xintercept=100), linetype="dashed",color="red")+
  geom_point(data=dfmasewoom12i %>%filter(vertex %in% c(1:(nperturb))),aes(x=vertex,y=pvalues), color="red", alpha=1,size=.5,shape=17)+
  facet_grid(type~time,switch = "both",scales = "free")+theme_bw() +ylab("Adjusted p-value")+
  theme(legend.position = "none",plot.title = element_text(size=12),plot.subtitle = element_text(hjust = 0.5), text = element_text(size=12, face='bold'),strip.text = element_text(size = 12))

#-------------------------------------------
#diff dimension case
alpha = 1-.125
dSVD=5
dASE=4
dSVD12=8

#pairwise mase
out <- foreach(i = 1:nullnmc, .combine='rbind') %dopar% {
  set.seed(123+i-1)
  glist <- genTSGdd(n, nperturb, cperturb=0, rmin, rmax, tmax, d_true, alpha)
  out2 <- doMase(glist, 2, dSVD,dASE,center, approx, py, method3)
  out12 <- doMase(glist, 12, dSVD12,dASE,center, approx, py, method3)
  out.p2 <- data.frame(as.vector(t(out2$tnorm[-c(5,6,7)])), mase="MASE(2)")
  colnames(out.p2)[1] <- "tnorm"
  out.p12 <- data.frame(as.vector(t(out12$tnorm[-c(5,6,7)])), mase="MASE(12)")
  colnames(out.p12)[1] <- "tnorm"
  out.p <- rbind(out.p2, out.p12)
}
df <- foreach(i = 1:pvalnmc, .combine='rbind') %dopar% {
  set.seed(123+nullnmc+i-1)
  glist <- genTSGdd(n, nperturb, cperturb, rmin, rmax, tmax, d_true, alpha)
  out2 <- doMase(glist, 2, dSVD,dASE,center, approx, py, method3)
  out12 <- doMase(glist, 12, dSVD12,dASE, center, approx, py, method3)
  pvalues2 <- rep(0, (tmax-1))
  pvalues12 <- rep(0, (tmax-1))
  for (j in 1:11) {
    pvalues2[j] <- sum(out2$tnorm[j]<sort(out[which(out$mase=="MASE(2)"),1]))/(nullnmc*8)
    pvalues12[j] <- sum(out12$tnorm[j]<sort(out[which(out$mase=="MASE(12)"),1]))/(nullnmc*8)
  }
  pvalues2 <- p.adjust(pvalues2, "BH")
  pvalues12 <- p.adjust(pvalues12, "BH")
  df2 <- data.frame(time=factor(m2, levels=m2), pvalues=pvalues2, mc=i,mase="MASE(2)")
  df12 <- data.frame(time=factor(m2, levels=m2), pvalues=pvalues12, mc=i,mase="MASE(12)")
  df <- rbind(df2, df12)
  df
}
df <- df[,-3]
df <- df %>%
  group_by(time, mase) %>%
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



dSVD=5
dASE=4
dSVD12=10

#pairwise mase
out <- foreach(i = 1:nullnmc, .combine='rbind') %dopar% {
  set.seed(123+i-1)
  glist <- genTSGdd(n, nperturb, cperturb=0, rmin, rmax, tmax, d_true, alpha)
  out2 <- doMase(glist, 2, dSVD,dASE,center, approx, py, method3)
  out12 <- doMase(glist, 12, dSVD12,dASE,center, approx, py, method3)
  out.p2 <- data.frame(t(out2$pdist[,1]), mase="MASE(2)")
  out.p12 <- data.frame(t(out12$pdist[,1]), mase="MASE(12)")
  out.p <- rbind(out.p2, out.p12)
}

df <- foreach(i = 1:pvalnmc, .combine='rbind') %dopar% {
  set.seed(123+nullnmc+i-1)
  glist <- genTSGdd(n, nperturb, cperturb, rmin, rmax, tmax, d_true, alpha)
  out2 <- doMase(glist, 2, dSVD,dASE,center, approx, py, method3)
  out12 <- doMase(glist, 12, dSVD12, dASE, center, approx, py, method3)
  pvalues2 <- matrix(0, n, tmax-1)#rep(0, 11)
  pvalues12 <- matrix(0, n, tmax-1)
  for (j in 1:(tmax-1)) {
    for (w in 1:n) {
      pvalues2[w,j] <- sum(out2$pdist[w,j]<sort(out[which(out$mase=="MASE(2)"),w]))/(nullnmc)
      pvalues12[w,j] <- sum(out12$pdist[w,j]<sort(out[which(out$mase=="MASE(12)"),w]))/(nullnmc)
    }
  }

  pvalues2.adjust <- matrix(p.adjust(pvalues2, "BH"),n, tmax-1)
  pvalues12.adjust <- matrix(p.adjust(pvalues12, "BH"),n, tmax-1)
  df.p2 <- melt(pvalues2.adjust)
  names(df.p2) <- c("vertex", "time", "pvalues")
  df.p2$time <- rep(factor(m2, levels=m2), each = n)
  df.p2 <- cbind(df.p2,mc=i, mase="MASE(2)")
  df.p12 <- melt(pvalues12.adjust)
  names(df.p12) <- c("vertex", "time", "pvalues")
  df.p12$time <- rep(factor(m2, levels=m2), each = n)
  df.p12 <- cbind(df.p12,mc=i, mase="MASE(12)")
  df.p <- rbind(df.p2, df.p12)
  df.p
  
}
df <- df[,-4]
df <- df %>%
  group_by(vertex, time, mase) %>%
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
dfommawoom12dd <- dfommawoom12dd 
# pddo <- ggplot(dfommawoom12dd, aes(x=time, y=pvalues, group=type)) + 
#   geom_errorbar(aes(ymin=lci, ymax=uci), width=.1, color=mycol) +
#   facet_wrap(~type,nrow = 2,switch  = "y", scales = "free") +
#   geom_line(color=mycol) +ylab("Adjusted p-value")+
#   
#   geom_point(color=mycol) +theme_bw()+ theme(legend.position = "none",plot.title = element_text(size=25),plot.subtitle = element_text(hjust = 0.5), text = element_text(size=(17)))+
#   geom_hline(aes(yintercept=.05))
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
dfmasewoom12i <- dfmaseip %>% filter(time=="2:3"|time=="6:7") #%>% filter(type!="MASE(12)")
# pddi <- ggplot() +   geom_errorbar(data=dfmasewoom12i, aes(ymin=lci, ymax=uci,x=vertex), width=.5, color=mycol) +
#   geom_point(data=dfmasewoom12i,aes(x=vertex,y=pvalues),alpha=1, size=.5)+ 
#   geom_hline(data=dfmasewoom12i,aes(yintercept=.05))+geom_vline(data=dfmasewoom12i,aes(xintercept=100), linetype="dashed",color="red")+
#   geom_point(data=dfmasewoom12i %>%filter(vertex %in% c(1:(nperturb))),aes(x=vertex,y=pvalues), color="red", alpha=1,size=.5)+
#   facet_grid(type~time,switch = "both",scales = "free")+theme_bw() +ylab("Adjusted p-value")+
#   theme(legend.position = "none",plot.title = element_text(size=25),plot.subtitle = element_text(hjust = 0.5), text = element_text(size=(17)))# +ggtitle(TeX(paste0("center=",center,", n=",n, ", $ \\delta_x=$",cperturb)),subtitle = "time" ) 
pddi <- ggplot() +   geom_errorbar(data=dfmasewoom12i, aes(ymin=lci, ymax=uci,x=vertex), width=.5, color=mycol) +
  geom_point(data=dfmasewoom12i %>%filter(!vertex %in% c(1:(nperturb))),aes(x=vertex,y=pvalues),alpha=1, size=.5,color="grey70",shape=20)+ 
  geom_hline(data=dfmasewoom12i,aes(yintercept=.05))+geom_vline(data=dfmasewoom12i,aes(xintercept=100), linetype="dashed",color="red")+
  geom_point(data=dfmasewoom12i %>%filter(vertex %in% c(1:(nperturb))),aes(x=vertex,y=pvalues), color="red", alpha=1,size=.5,shape=17)+
  facet_grid(type~time,switch = "both",scales = "free")+theme_bw() +ylab("Adjusted p-value")+
  theme(legend.position = "none",plot.title = element_text(size=12),plot.subtitle = element_text(hjust = 0.5), text = element_text(size=12, face='bold'),strip.text = element_text(size = 12))


# Figure 4 a
grid.arrange(psdo,psdi, ncol=1, heights=c(4,4))
png("sdmcomultid4june19n400maseRevisionFDR.png", width = 6, height = 6, units="in", res=400)
grid.arrange(psdo,psdi, ncol=1, heights=c(4,4))
dev.off()

# Figure 4 b
grid.arrange(pddo,pddi, ncol=1, heights=c(4,4))
png("ddmcomultid4june18n400maseRevisionFDR.png", width = 6, height = 6, units="in", res=400)
grid.arrange(pddo,pddi, ncol=1, heights=c(4,4))
dev.off()



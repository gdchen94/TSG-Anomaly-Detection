# Install and load the packages if not already installed
packages <- c("igraph", "gtools", "tidyverse", "irlba", "doParallel", 
              "qcc", "latex2exp", "gridExtra", "reshape2","metafolio", "rstudioapi")

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
# registerDoParallel(detectCores()-1) ## uncomment this to use multicore

cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
colors =  cbbPalette[c(7, 8, 4, 6)]
py = FALSE
method3 = 2
dSVD=2
dASE=1
nullnmc <- 100
pvalnmc <- 100
tmax <- 12
tvec <- 1:tmax
m2 <- names(running(tvec, width=2))
n <- 100
d <- 1
nperturb <- 20
rmin <- 0.2
rmax <- 0.8
center <- FALSE
cperturb <- .15
approx <- TRUE
center <- FALSE

#pairwise mase
out <- foreach(i = 1:nullnmc, .combine='rbind') %dopar% {
  set.seed(123+i-1)
  glist <- genTSG(n, nperturb, cperturb=0, rmin, rmax, tmax)
  out2 <- doMase(glist, 2, dSVD,dASE,center, approx, py, method3)
  out12 <- doMase(glist, 12, 3,dASE,center, approx, py, method3)
  out.p2 <- data.frame(as.vector(t(out2$tnorm[-c(5,6,7)])), mase="MASE(2)")
  colnames(out.p2)[1] <- "tnorm"
  out.p12 <- data.frame(as.vector(t(out12$tnorm[-c(5,6,7)])), mase="MASE(12)")
  colnames(out.p12)[1] <- "tnorm"
  out.p <- rbind(out.p2, out.p12)
}
df <- foreach(i = 1:pvalnmc, .combine='rbind') %dopar% {
  set.seed(123+nullnmc+i-1)
  glist <- genTSG(n, nperturb, cperturb, rmin, rmax, tmax)
  out2 <- doMase(glist, 2, dSVD,dASE,center, approx, py, method3)
  out12 <- doMase(glist, 12, 3,dASE, center, approx, py, method3)
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

dfmo <- df %>%
  rename(pvalues = pvalues_median, lci = pvalues_lci, uci = pvalues_uci)

mycol <- gg_color_hue(2)[2]

#only mase
colnames(dfmo)[2] <- "type"
dfmo <- ungroup(dfmo)
dfommao <- dfmo
dfommawoom12o <- dfommao #%>% filter(type!="MASE(12)")
p <-ggplot(dfommawoom12o, aes(x=time, y=pvalues, group=type)) + 
  geom_errorbar(aes(ymin=lci, ymax=uci), width=.1, color=mycol) +
  facet_wrap(~type,nrow = 3,strip.position ="left", scales = "fixed") +
  geom_line(color="grey70") +ylab("Adjusted p-value")+
  geom_point(data=dfommawoom12o%>%filter(time=="6:7"&pvalues<0.05), color="red",shape=17) +
  geom_point(data=dfommawoom12o%>%filter(!(time=="6:7"&pvalues<0.05)), color="grey",shape=20) +
  theme_bw()+ theme(text = element_text(size=(16), face = 'bold')) + geom_hline(aes(yintercept=.05))

#------------------------vertex pvalues

#pairwise mase
out <- foreach(i = 1:nullnmc, .combine='rbind') %dopar% {
  set.seed(123+i-1)
  glist <- genTSG(n, nperturb, cperturb=0, rmin, rmax, tmax)
  out2 <- doMase(glist, 2, dSVD,dASE,center, approx, py, method3)
  #out2$pdist[,1]
  out12 <- doMase(glist, 12, 3,dASE,center, approx, py, method3)
  #out2$pdist[,1]
  out.p2 <- data.frame(t(out2$pdist[,-c(5,6,7)]), mase="MASE(2)")
  out.p12 <- data.frame(t(out12$pdist[,-c(5,6,7)]), mase="MASE(12)")
  out.p <- rbind(out.p2, out.p12)
}

df <- foreach(i = 1:pvalnmc, .combine='rbind') %dopar% {
  set.seed(123+nullnmc+i-1)
  glist <- genTSG(n, nperturb, cperturb, rmin, rmax, tmax)
  out2 <- doMase(glist, 2, dSVD,dASE,center, approx, py, method3)
  out12 <- doMase(glist, 12, 3,dASE, center, approx, py, method3)
  pvalues2 <- matrix(0, n, tmax-1)#rep(0, 11)
  pvalues12 <- matrix(0, n, tmax-1)
  for (j in 1:11) {
    for (w in 1:n) {
      pvalues2[w,j] <-sum(out2$pdist[w,j]<sort(out[which(out$mase=="MASE(2)"),w]))/(8*nullnmc)
      pvalues12[w,j] <- sum(out12$pdist[w,j]<sort(out[which(out$mase=="MASE(12)"),w]))/(8*nullnmc)
    }
  }
  pvalues2.adjust <- matrix(p.adjust(pvalues2, "BH"),n, tmax-1)
  pvalues12.adjust <- matrix(p.adjust(pvalues12, "BH"),n, tmax-1)

  df.p2 <- reshape2::melt(pvalues2.adjust)
  names(df.p2) <- c("vertex", "time", "pvalues")
  df.p2$time <- rep(factor(m2, levels=m2), each = n)
  df.p2 <- cbind(df.p2,mc=i, mase="MASE(2)")
  df.p12 <- reshape2::melt(pvalues12.adjust)
  names(df.p12) <- c("vertex", "time", "pvalues")
  df.p12$time <- rep(factor(m2, levels=m2), each = n)
  df.p12 <- cbind(df.p12,mc=i, mase="MASE(12)")
  df.p <- rbind(df.p2, df.p12)
  df.p
}
dfmi <- df[,-4]
# corrected for pvalues, change n() to n in z-norm
dfmi <- dfmi %>%
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

dfmi <- dfmi %>%
  rename(pvalues = pvalues_median, lci = pvalues_lci, uci = pvalues_uci)
# ppp1mi <-  ggplot() + 
#   geom_errorbar(data=dfmi, aes(ymin=pvalues-lci, ymax=pvalues+uci,x=vertex), width=.05, color=mycol) +
#   geom_point(data=dfmi %>%filter(vertex %in% c((nperturb+1):n)),aes(x=vertex,y=pvalues),alpha=1, size=.1)+
#   geom_point(data=dfmi %>%filter(vertex %in% c(1:nperturb)),aes(x=vertex,y=pvalues), color="red", alpha=0.3)+
#   geom_hline(data=dfmi,aes(yintercept=.05))+facet_grid(mase~time,scales = "free") +theme_bw()+ylab("p-value")+
#   theme(legend.position = "none",plot.title = element_text(size=10, face='bold'),plot.subtitle = element_text(hjust = 0.5)) +ggtitle(TeX(paste0("center=",center,", n=",n, ", $ \\delta_x=$",cperturb)),subtitle = "time" ) 




#only mase
colnames(dfmi)[3] <- "type"
dfmi <- ungroup(dfmi)
dfommaip <- dfmi
dfommawoom12i <- dfommaip %>% filter(time=="2:3"|time=="6:7")
pi <- ggplot() +   geom_errorbar(data=dfommawoom12i, aes(ymin=lci, ymax=uci,x=vertex), width=.5, color=mycol) +
  geom_point(data=dfommawoom12i %>%filter(!vertex %in% c(1:(nperturb))),aes(x=vertex,y=pvalues),alpha=1, size=.5,color="grey70",shape=20)+ 
  geom_hline(data=dfommawoom12i,aes(yintercept=.05))+geom_vline(data=dfommawoom12i,aes(xintercept=20), linetype="dashed",color="red")+
  geom_point(data=dfommawoom12i %>%filter(vertex %in% c(1:(nperturb))),aes(x=vertex,y=pvalues), color="red", alpha=1,size=.5,shape=17)+
  facet_grid(type~time,switch = "both",scales = "free")+theme_bw() +ylab("Adjusted p-value")+
  theme(legend.position = "none",plot.title = element_text(size=14, face='bold'),plot.subtitle = element_text(hjust = 0.5), text = element_text(size=(14), face='bold'),strip.text = element_text(size = 14))


# Figure 3
grid.arrange(p,pi, ncol=2, heights=c(3))
png("pvaluesexample2methodMASEFDR.png", width = 12, height = 3, units="in", res=400)
grid.arrange(p,pi, ncol=2, heights=c(3))
dev.off()
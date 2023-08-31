# Install and load the packages if not already installed
packages <- c("igraph", "gtools", "tidyverse", "irlba", "doParallel", 
              "qcc", "latex2exp", "gridExtra", "reshape2", "rstudioapi")

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

#mase
py = FALSE
method3 = 2
dSVD=2
dASE=1
nmc <- 1
approx <- TRUE
center <- FALSE

set.seed(59)
tmax <- 12
tvec <- 1:tmax
m2 <- names(running(tvec, width=2))
n <- 100
d <- 1
nperturb <- 20
rmin <- 0.2
rmax <- 0.8

cperturb <- .15
mean2 <- rep(1,(tmax-1))
std2 <- rep(1,(tmax-1))
mean12 <- rep(1, (tmax-1))
std12 <- rep(1,(tmax-1))
for (w in 1:(tmax-1)) {
  out <- foreach(i = 1:nmc, .combine='rbind') %dopar% {
    set.seed(123+i-1)
    tempglist <- conchartgenTSG(n, nperturb, cperturb, rmin, rmax, tmax+10)
    glist <- tempglist[(w:(w+10))]
    out2 <- doMase(glist, 2, dSVD,dASE,center, approx, py, method3)
    out12 <- doMase(glist, 12, 6,dASE, center, approx, py, method3)
    df.norm2 <- data.frame(time=factor(m2[1:(tmax-2)], levels=m2[1:(tmax-2)]), norm=out2$tnorm, mc=i, mase="mase2")
    df.norm12 <- data.frame(time=factor(m2[1:(tmax-2)], levels=m2[1:(tmax-2)]), norm=out12$tnorm, mc=i, mase="mase12")
    df.norm <- rbind(df.norm2, df.norm12)
    df.norm
  }
  df <- out  %>% filter(mase=="mase2")
  df <- df[,c(-1,-3,-4)]
  mase2 <- matrix(df, tmax-2, nmc)
  mean2[w] <- mean(mase2)
  std2[w] <- sd.xbar.one(mase2, rep(nmc,tmax-2), "MR")#sd.xbar(mase2, rep(nmc,tmax-2), "MR")
  df <- out  %>% filter(mase=="mase12")
  df <- df[,c(-1,-3,-4)]
  mase12 <- matrix(df, tmax-2, nmc)
  mean12[w] <- mean(mase12)
  std12[w] <- sd.xbar.one(mase12, rep(nmc,tmax-2), "MR")#sd.xbar(mase12, rep(nmc,tmax-2), "MR")
}

out <- foreach(i = 1:nmc, .combine='rbind') %dopar% {
  set.seed(123+i-1)
  glist <- genTSG(n, nperturb, cperturb, rmin, rmax, tmax)
  out2 <- doMase(glist, 2, dSVD,dASE,center, approx, py, method3)
  out12 <- doMase(glist, 12, 6,dASE, center, approx, py, method3)
  df.norm2 <- data.frame(time=factor(m2, levels=m2), norm=out2$tnorm, mc=i, mase="mase2")
  df.norm12 <- data.frame(time=factor(m2, levels=m2), norm=out12$tnorm, mc=i, mase="mase12")
  df.norm <- rbind(df.norm2, df.norm12)
  df.norm
}
df <- out  %>% filter(mase=="mase2")
df <- df[,c(-1,-3,-4)]
mase2 <- matrix(df, tmax-1, nmc)
c1m2 <- qcc(mase2,type = "xbar.one",center = mean2 ,std.dev = std2, nsigmas=3, plot=FALSE)
df <- out  %>% filter(mase=="mase12")
df <- df[,c(-1,-3,-4)]
mase12 <- matrix(df, tmax-1, nmc)
c1m12 <- qcc(mase12,type = "xbar.one",center = mean12,std.dev = std12, nsigmas=3, plot=FALSE)


df <- data.frame(time=seq(length(c1m2$statistics)), y= c1m2$statistics,ucl=c1m2$limits[,2],center=c1m2$center,lim=replace(rep(0,length(c1m2$statistics)),c1m2$violations$beyond.limits,1),run=replace(rep(0,length(c1m2$statistics)),c1m2$violations$violating.runs,1), Method="MASE(2)")
df <- rbind(df, data.frame(time=seq(length(c1m12$statistics)), y= c1m12$statistics,ucl=c1m12$limits[,2],center=c1m12$center,lim=replace(rep(0,length(c1m12$statistics)),c1m12$violations$beyond.limits,1),run=replace(rep(0,length(c1m12$statistics)),c1m12$violations$violating.runs,1), Method="MASE(12)"))
df$Method <-  factor(df$Method, levels = c("MASE(2)","MASE(12)"))
# Figure 2 left
pg <- ggplot(df,aes(x=time, y=y))+
  geom_step(aes(x=time, y=ucl), linetype="dashed")+
  geom_step(aes(x=time, y=center), linetype="solid") +
  geom_point(data = df%>% filter(lim!=1), alpha=1, color="grey70")+ylab(TeX("$y^{(t)}$"))+
  geom_line(aes(x=time, y=y), color="grey70", shape=20)+
  # geom_point(data = df %>% filter(run==1), alpha=1, color="yellow")+
  geom_point(data = df %>% filter(lim==1), alpha=1, color="red", shape=17) + 
  scale_x_discrete(name ="time points", 
                   limits=c(m2))+theme_bw()+
  theme_classic(base_size = 18)+
  facet_wrap(~Method,nrow = 2,scales = "free_y", strip.position = "left")+
  theme(axis.text.x = element_text(size = 12,angle=90),legend.position = "none",plot.title = element_text(hjust = 0.5,size=20, face='bold'),plot.subtitle = element_text(hjust = 0.5)) +ggtitle("Control Chart")

#mase
method3 = 3
set.seed(123)
mean2 <- rep(1,(tmax-1))
std2 <- rep(1,(tmax-1))
mean12 <- rep(1, (tmax-1))
std12 <- rep(1,(tmax-1))
for (w in 1:(tmax-1)) {
  out <- foreach(i = 1:nmc, .combine='rbind') %dopar% {
    set.seed(123+i-1)
    tempglist <- conchartgenTSG(n, nperturb, cperturb, rmin, rmax, tmax+10)
    glist <- tempglist[(w:(w+10))]
    out2 <- doMase(glist, 2, dSVD,dASE,center, approx, py, method3)
    out12 <- doMase(glist, tmax, 6,dASE,center, approx, py, method3)
    df.dist2 <- melt(out2$pdist)
    names(df.dist2) <- c("vertex", "time", "pdist")
    df.dist2$time <- rep(factor(m2[1:(tmax-2)], levels=m2[1:(tmax-2)]), each = n)
    df.dist2 <- cbind(df.dist2,mc=i, mase="mase2")
    df.dist12 <- melt(out12$pdist)
    names(df.dist12) <- c("vertex", "time", "pdist")
    df.dist12$time <- rep(factor(m2[1:(tmax-2)], levels=m2[1:(tmax-2)]), each = n)
    df.dist12 <- cbind(df.dist12,mc=i, mase="mase12")
    df.dist <- rbind(df.dist2,df.dist12)
  }
  df <- out  %>% filter(mase=="mase2")
  df <- df[,c(-4,-5)]
  mase2 <- matrix(df[,3], 10, nmc*n, byrow = TRUE)
  mean2[w] <- mean(mase2)
  std2[w] <- sd.xbar(mase2, rep(nmc*n,10), "UWAVE-SD")
  df <- out  %>% filter(mase=="mase12")
  df <- df[,c(-4,-5)]
  mase12 <- matrix(df[,3], 10, nmc*n, byrow = TRUE)
  mean12[w] <- mean(mase12)
  std12[w] <- sd.xbar(mase12, rep(nmc*n,10), "UWAVE-SD")#sd.xbar(mase12, rep(nmc,10), "MR")
}

out <- foreach(i = 1:nmc, .combine='rbind') %dopar% {
  set.seed(123+i-1)
  glist <- genTSG(n, nperturb, cperturb, rmin, rmax, tmax)
  out2 <- doMase(glist, 2, dSVD,dASE,center, approx, py, method3)
  out12 <- doMase(glist, tmax, 6,dASE,center, approx, py, method3)
  df.dist2 <- melt(out2$pdist)
  names(df.dist2) <- c("vertex", "time", "pdist")
  df.dist2$time <- rep(factor(m2, levels=m2), each = n)
  df.dist2 <- cbind(df.dist2,mc=i, mase="mase2")
  df.dist12 <- melt(out12$pdist)
  names(df.dist12) <- c("vertex", "time", "pdist")
  df.dist12$time <- rep(factor(m2, levels=m2), each = n)
  df.dist12 <- cbind(df.dist12,mc=i, mase="mase12")
  df.dist <- rbind(df.dist2,df.dist12)
}
df <- out  %>% filter(mase=="mase2")
df <- df[,c(-4,-5)]
mase2 <- matrix(df[,3], 11, nmc*n, byrow = TRUE)

c1m2v <- list()
for (i in 1:(tmax-1)) {
  c1m2v[[i]] <- qcc(mase2[i,],type = "xbar.one",center = mean2[i],std.dev = std2[i], nsigmas=3, plot=FALSE)
}
df <- out  %>% filter(mase=="mase12")
df <- df[,c(-4,-5)]
mase12 <- matrix(df[,3], 11, nmc*n, byrow = TRUE)
c1m12v <- list()
for (i in 1:(tmax-1)) {
  c1m12v[[i]] <- qcc(mase12[i,],type = "xbar.one",center = mean12[i],std.dev = std12[i], nsigmas=3, plot=FALSE)
}


p1v= plot.qcc.vertex(c1m2v, add.stats = FALSE, chart.all = TRUE, 
                     label.limits = c("LCL ", "UCL"), title="Control Chart MASE(2)", xlab="time points",m2=seq(n), ylab="y", 
                     axes.las = 0, digits = getOption("digits"),
                     restore.par = FALSE)

p2v=plot.qcc.vertex(c1m12v, add.stats = FALSE, chart.all = TRUE, 
                    label.limits = c("LCL ", "UCL"), title="Control Chart MASE(12)", xlab="time points",m2=seq(n), ylab="y", 
                    axes.las = 0, digits = getOption("digits"),
                    restore.par = FALSE)

pv <- grid.arrange(p1v,p2v, ncol=1)

# Figure 2 
grid.arrange(pg,pv, ncol=2, widths=c(4,4))
png("qccIllusMase.png", width = 8, height = 4, units="in", res=400)
grid.arrange(pg,pv, ncol=2, widths=c(4,4))#, heights=c(4,4),widths=c(6,6))
dev.off()

packages <- c("gtools", "tidyverse", "doParallel", "qcc", 
              "latex2exp", "gridExtra", "ggalt", "rstudioapi")

for(pkg in packages) {
  if(!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg)
  }
  library(pkg, character.only = TRUE)
}

# registerDoParallel(detectCores()-1) ## uncomment this to use multicore

# Get the path of the currently active document in RStudio
current_script_path <- getSourceEditorContext()$path

# Set the working directory to the directory containing the current script
setwd(dirname(current_script_path))
source("utils.R")

seed=155
set.seed(seed)
tmax <- 12
d_true <- 4
n <- 400#200#200#500
alpha <- .125
rmin <- 0.2
rmax <- 0.8
cperturb <- .12
nperturb <- 100

k <- d_true
Z <- matrix(0, n, d_true)
initialZ <- matrix(0, n, d_true)
tempinitialZ <- matrix(0, n, d_true)
if(alpha<.1^4){
  for (i in 1:k) {
    Z[((n/k)*(i-1)+1):((n/k)*(i)), i] = 1
  }
}else{
  tempinitialZ <- na.omit(rdirichlet(n, alpha * rep(1, d_true) ))
  while(dim(tempinitialZ)[1]!=n){
    i <- n - dim(tempinitialZ)[1]
    tempinitialZ <- rbind(tempinitialZ, na.omit(rdirichlet(i, alpha * rep(1, d_true) )))
  }
  initialZ <- tempinitialZ
  Z <- initialZ
}

q <- .5 * runif(1, rmin, rmax)
B <- (  .5 * runif(1, rmin, rmax) - q + .5 )* diag(rep(1,d_true) ) +  q* rep(1, d_true) %*% t(rep(1, d_true))
temp <-  svd(Z%*%B %*% t(Z))
X1 <- temp$u[,1:d_true] %*%  diag( sqrt(temp$d[1:d_true]) )
order_knn <- order(apply( X1, MARGIN = 1, function(x) sum((x-X1[1,])^2) ))
temp <- X1
temp[(1:(n/k)),] <- X1[order_knn[(1:(n/k))],]
temp[-(1:(n/k)),] <- X1[-order_knn[(1:(n/k))],]
X1 <- temp
Xlist <- rep(list(X1), tmax)
Z_gamma <- matrix(0, n, k)
Z_pert <- t(matrix(rep(.6*rdirichlet(1, rep(1,  d_true) )+.2, times = n/k), d_true,n/k) )
Z_gamma <- rbind(Z_pert, matrix(0, n-n/k, d_true))


Xlist[[6]] <- Xlist[[1]] + cperturb * Z_gamma
Xlist[[7]] <- Xlist[[1]] - cperturb * Z_gamma
xlims = c(-1, 1)
xlims = c(-0.2, 0.2)
ylims = c(-1, 1)
ylims = c(-0.2, 0.2)
V <- rbind(X1,Xlist[[6+1]])
vec <- data.frame(xstart=X1[1:nperturb,1],ystart=X1[1:nperturb,2],xe=Xlist[[6+1]][1:nperturb,1],ye=Xlist[[6+1]][1:nperturb,2])
cols <- c(rep("grey90", n),rep("red", n))
sp <- factor(rep(c(rep("anomalous",nperturb),rep("normal",n-nperturb)),times=2))
df <- data.frame(x = V[,1], y = V[,2], col=cols, class=sp)
p0 <- ggplot() +
        geom_point(data=df%>% filter((class!="anomalous")|(col!="red")), aes(x=x, y=y),color="grey70",shape=20,size=3) +
        geom_encircle( data=df%>% filter(col=="grey90"),aes(x=x,y=y),linetype = 2,color="black")+
        ylab("") +
        scale_color_manual(values=c("black"))+
  theme_void() + theme(legend.position="none",plot.title = element_text(hjust = 0.5, face = "italic"))+ggtitle("t = 5")
V <- rbind(X1,Xlist[[6]])
vec <- data.frame(xstart=X1[1:nperturb,1],ystart=X1[1:nperturb,2],xe=Xlist[[6]][1:nperturb,1],ye=Xlist[[6]][1:nperturb,2])
cols <- c(rep("grey90", n),rep("red", n))
sp <- factor(rep(c(rep("anomalous",nperturb),rep("normal",n-nperturb)),times=2))
df <- data.frame(x = V[,1], y = V[,2], col=cols, class=sp)
p1 <- ggplot() +
  geom_point(data=df%>% filter(class=="anomalous"&col=="red"), aes(x=x, y=y),color="red",shape=17,size=3) +
  geom_point(data=df%>% filter((class!="anomalous")&(col=="red")), aes(x=x, y=y),color="grey70",shape=20,size=3) +
  geom_encircle( data=df%>%  filter(class=="anomalous"&col=="red"),aes(x=x,y=y),linetype = 2,color="red")+
  geom_encircle( data=df%>% filter((class!="anomalous")&(col=="red")),aes(x=x,y=y),linetype = 2,color="black")+
  geom_segment(data=vec[5,],aes(x = xstart, y = ystart, xend = xe, yend =ye),
               lineend = , linejoin = "round",col="blue",
               size =1, arrow = arrow(length = unit(.2, "inches"))
  )+
  ylab("") +
  scale_color_manual(values=c("black"))+
  theme_void() + theme(legend.position="none",plot.title = element_text(hjust = 0.5, face = "italic"))+ggtitle("t = 6")
V <- rbind(X1,Xlist[[7]])
vec <- data.frame(xstart=X1[1:nperturb,1],ystart=X1[1:nperturb,2],xe=Xlist[[7]][1:nperturb,1],ye=Xlist[[7]][1:nperturb,2])
cols <- c(rep("grey90", n),rep("red", n))
sp <- factor(rep(c(rep("anomalous",nperturb),rep("normal",n-nperturb)),times=2))
df <- data.frame(x = V[,1], y = V[,2], col=cols, class=sp)
p2 <- ggplot() +
  geom_point(data=df%>% filter(class=="anomalous"&col=="red"), aes(x=x, y=y),color="red",shape=17,size=3) +
  geom_point(data=df%>% filter((class!="anomalous")&(col=="red")), aes(x=x, y=y),color="grey70",shape=20,size=3) +
  geom_encircle( data=df%>%  filter(class=="anomalous"&col=="red"),aes(x=x,y=y),linetype = 2,color="red")+
  geom_encircle( data=df%>% filter((class!="anomalous")&(col=="red")),aes(x=x,y=y),linetype = 2,color="black")+
  geom_segment(data=vec[5,],aes(x = xstart, y = ystart, xend = xe, yend =ye),
               lineend = , linejoin = "round",col="blue",
               size =1, arrow = arrow(length = unit(.2, "inches"))
  )+
  ylab("") +
  scale_color_manual(values=c("black"))+
  theme_void() + theme(legend.position="none",plot.title = element_text(hjust = 0.5, face = "italic"))+ggtitle("t = 7")
p3 <- p0 +ggtitle("t = 8")

p_latpos_illustration<-grid.arrange(p0,p1,p2,p3, ncol=2, heights=c(4,4))


# independent of device size
nmc <- 100
tmax <- 12
tvec <- 1:tmax
m2 <- names(running(tvec, width=2))
n <- 100#100
d <- 1#
nperturb <- 20 # was 2
rmin <- 0.2
rmax <- 0.8
center <- FALSE

cperturb <- .15
meant <- rep(1,(tmax-1))
stdt <- rep(1,(tmax-1))
for (w in 1:(tmax-1)) {
  out <- foreach(i = 1:nmc, .combine='rbind') %dopar% {
    set.seed(123+i-1)
    out2 <- (truenorm.aux(n, nperturb, cperturb, rmin, rmax, tmax+10)$tnorm)[(w:(w+9))]
    df.norm <- data.frame(time=factor(m2[1:(tmax-2)], levels=m2[1:(tmax-2)]), norm=out2, mc=i, mase="true")
    df.norm
  }
  df <- out  
  df <- df[,c(-1,-3,-4)]
  true <- matrix(df, tmax-2, nmc)
  meant[w] <- mean(true)
  stdt[w] <- sd.xbar.one(true, rep(nmc,tmax-2), "MR")
  df <- out  
  df <- df[,c(-1,-3,-4)]
}

out <- foreach(i = 1:nmc, .combine='rbind') %dopar% {
  set.seed(123+i-1)
  out2 <- truenorm(n, nperturb, cperturb, rmin, rmax, tmax)
  df.norm2 <- data.frame(time=factor(m2, levels=m2), norm=out2$tnorm, mc=i, mase="true")
  df.norm2
}

out <- out[,-3]
df <- out %>% 
  group_by(time, mase) %>% 
  summarize(across(everything(), list(mean = ~mean(.), sd = ~sd(.), se = ~sd(.)/sqrt(n()))))
true <- df$norm_mean
c1t <- qcc(true, type = "xbar.one",center = meant,std.dev = stdt, nsigmas=3, plot=FALSE)


# Figure 1a
plot.qcc.true(c1t,m2=m2)
png("illuswolclRevision.png", width = 6, height = 2.5, units="in", res=400)
plot.qcc.true(c1t,m2=m2)
dev.off()

# Figure 1b
grid.arrange(p0,p1,p2,p3, ncol=2, heights=c(4,4))
png("alpha0125illus.png", width = 4, height = 2.5, units="in", res=400)
grid.arrange(p0,p1,p2,p3, ncol=2, heights=c(4,4))
dev.off()


# Figure 1 combined
grid.arrange(plot.qcc.true(c1t,m2=m2),p_latpos_illustration, ncol=2, widths=c(4,3))
png("illuswolclRevision_alpha0125illus.png", width = 10, height = 2.5, units="in", res=400)
grid.arrange(plot.qcc.true(c1t,m2=m2),p_latpos_illustration, ncol=2, widths=c(4,3))
dev.off()

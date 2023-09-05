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
py = FALSE
method3 = 2
d1=4
d2=4
d3=4
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
cperturb <- .12

approx <- TRUE
center <- FALSE

#pairwise mase
# Initialize the final data frame to store results
df <- data.frame()
alpha_values = c(0, seq(8)/8)
for(a in 1:length(alpha_values)){
  alpha = alpha_values[a]
  
  # Initialize the intermediate 'out' data frame to collect results for each alpha
  out = data.frame()
  
  for(i in 1:nullnmc){
    # Setting the seed for reproducibility
    set.seed(123 + i - 1)
    
    # Generating the time series data
    glist <- genTSGonepar(n, nperturb, cperturb=0, rmin, rmax, tmax, d_true, alpha)
    
    # Calculating MASE for two different periods (2 and 12)
    out2 <- doMase(glist, 2, d1, d3, center, approx, py, method3)
    out12 <- doMase(glist, 12, d1, d3, center, approx, py, method3)
    
    # Storing the results in data frames and setting column names
    out.p2 <- data.frame(as.vector(t(out2$tnorm[-c(5,6,7)])), mase="MASE2", alpha=factor(alpha))
    colnames(out.p2)[1] <- "tnorm"
    out.p12 <- data.frame(as.vector(t(out12$tnorm[-c(5,6,7)])), mase="MASE12", alpha=factor(alpha))
    colnames(out.p12)[1] <- "tnorm"
    
    # Combining the results for MASE2 and MASE12
    out.p <- rbind(out.p2, out.p12)
    
    # Adding the combined results to the 'out' data frame
    out <- rbind(out, out.p)
  }
  
  # Initialize vectors to hold p-values and counts
  c2 <- rep(0, 11)
  c12 <- rep(0, 11)
  
  # Calculating p-values and counts for the MASE
  for(i in 1:pvalnmc){
    set.seed(123+nullnmc+i-1)
    glist <- genTSGonepar(n, nperturb, cperturb, rmin, rmax, tmax, d_true, alpha)
    out2 <- doMase(glist, 2, d1, d3, center, approx, py, method3)
    out12 <- doMase(glist, 12, d1, d3, center, approx, py, method3)
    
    pvalues2 <- rep(0, (tmax-1))
    pvalues12 <- rep(0, (tmax-1))
    for (j in 1:11) {
      pvalues2[j] <- sum(out2$tnorm[j]<sort(out[which(out$mase=="MASE2"&out$alpha==alpha),1]))/(nullnmc*8)
      pvalues12[j] <- sum(out12$tnorm[j]<sort(out[which(out$mase=="MASE12"&out$alpha==alpha),1]))/(nullnmc*8)
    }
    
    c2 <- c2 + (p.adjust(pvalues2, "BH")<.05)
    c12 <- c12 + (p.adjust(pvalues12, "BH")<.05)
  }
  
  # Create data frames to store power and standard error for MASE2 and MASE12
  df2 <- data.frame(time=factor(m2, levels=m2), power=c2/pvalnmc, se = sqrt( (c2/pvalnmc) * (1 - c2/pvalnmc)/pvalnmc ), mase="MASE2", alpha=factor(alpha))
  df12 <- data.frame(time=factor(m2, levels=m2), power=c12/pvalnmc, se = sqrt( (c12/pvalnmc) * (1 - c12/pvalnmc)/pvalnmc ), mase="MASE12", alpha=factor(alpha))
  
  # Combine the data frames for MASE2 and MASE12
  df_temp <- rbind(df2, df12)
  
  # Add the combined data frame to the final result 'df'
  df <- rbind(df, df_temp)
}
dfm <- df
dfm <- dfm  %>% group_by(mase)

#pairwise omni
dmax <- 4#4
# Initialize the final data frame to store results
df <- data.frame()
for(a in 1:length(alpha_values)){
  alpha = alpha_values[a]
  
  # Initialize an 'out' data frame to collect results for each alpha
  out = data.frame()
  
  for(i in 1:nullnmc){
    # Setting the seed for reproducibility
    set.seed(123 + i - 1)
    
    # Generating the time series data
    glist <- genTSGonepar(n, nperturb, cperturb=0, rmin, rmax, tmax, d_true, alpha)
    
    # Calculating Omni for two different periods (2 and 12)
    out2 <- doOmni(glist, nomni=2, dmax, center=center, approx=approx)
    out12 <- doOmni(glist, nomni=12, dmax, center=center, approx=approx)
    
    # Storing the results in data frames and setting column names
    out.p2 <- data.frame(as.vector(t(out2$tnorm[-c(5,6,7)])), omni="OMNI2", alpha=factor(alpha))
    colnames(out.p2)[1] <- "tnorm"
    out.p12 <- data.frame(as.vector(t(out12$tnorm[-c(5,6,7)])), omni="OMNI12", alpha=factor(alpha))
    colnames(out.p12)[1] <- "tnorm"
    
    # Combining the results for OMNI2 and OMNI12
    out.p <- rbind(out.p2, out.p12)
    
    # Adding the combined results to the 'out' data frame
    out <- rbind(out, out.p)
  }
  
  # Initialize vectors to hold p-values and counts
  c2 <- rep(0, 11)
  c12 <- rep(0, 11)
  
  # Calculating p-values and counts for the Omni
  for(i in 1:pvalnmc){
    set.seed(123+nullnmc+i-1)
    glist <- genTSGonepar(n, nperturb, cperturb, rmin, rmax, tmax, d_true, alpha)
    out2 <- doOmni(glist, nomni=2, dmax, center=center, approx=approx)
    out12 <- doOmni(glist, nomni=12, dmax, center=center, approx=approx)
    
    pvalues2 <- rep(0, (tmax-1))
    pvalues12 <- rep(0, (tmax-1))
    for (j in 1:11) {
      pvalues2[j] <- sum(out2$tnorm[j]<sort(out[which(out$omni=="OMNI2"&out$alpha==alpha),1]))/(nullnmc*8)
      pvalues12[j] <- sum(out12$tnorm[j]<sort(out[which(out$omni=="OMNI12"&out$alpha==alpha),1]))/(nullnmc*8)
    }
    
    c2 <- c2 + (p.adjust(pvalues2, "BH")<.05)
    c12 <- c12 + (p.adjust(pvalues12, "BH")<.05)
  }
  
  # Create data frames to store power and standard error for OMNI2 and OMNI12
  df2 <- data.frame(time=factor(m2, levels=m2), power=c2/pvalnmc, se = sqrt((c2/pvalnmc) * (1 - c2/pvalnmc)/pvalnmc), omni="OMNI2", alpha=factor(alpha))
  df12 <- data.frame(time=factor(m2, levels=m2), power=c12/pvalnmc, se = sqrt((c12/pvalnmc) * (1 - c12/pvalnmc)/pvalnmc), omni="OMNI12", alpha=factor(alpha))
  
  # Combine the data frames for OMNI2 and OMNI12
  df_temp <- rbind(df2, df12)
  
  # Add the combined data frame to the final result 'df'
  df <- rbind(df, df_temp)
}
dfo <- df
dfo <- dfo  %>% group_by(omni)

#scan stat
#pairwise scan stat
# Initialize the final data frame to store results
df <- data.frame()
for(a in 1:length(alpha_values)){
  alpha = alpha_values[a]
  
  # Initialize an 'out' data frame to collect results for each alpha
  out = data.frame()
  
  for(i in 1:nullnmc){
    # Setting the seed for reproducibility
    set.seed(123 + i - 1)
    
    # Generating the time series data
    glist <- genTSGonepar(n, nperturb, cperturb=0, rmin, rmax, tmax, d_true, alpha)
    
    # Applying the 'doScan' function on the generated data
    out2 <- doScan(glist)
    
    # Creating a data frame to hold the results and setting column names
    out.p <- data.frame(as.vector(t(out2$tnorm[-c(5,6,7)])), type="SCAN", alpha=factor(alpha))
    colnames(out.p)[1] <- "tnorm"
    
    # Adding the results to the 'out' data frame
    out <- rbind(out, out.p)
  }
  
  # Initialize a vector to hold p-values and counts
  c2 <- rep(0, 11)
  
  # Calculating p-values and counts for the SCAN metric
  for(i in 1:pvalnmc){
    set.seed(123 + nullnmc + i - 1)
    
    # Generating the time series data
    glist <- genTSGonepar(n, nperturb, cperturb, rmin, rmax, tmax, d_true, alpha)
    
    # Applying the 'doScan' function on the generated data
    out2 <- doScan(glist)
    
    pvalues2 <- rep(0, (tmax - 1))
    for (j in 1:11) {
      pvalues2[j] <- sum(out2$tnorm[j] < sort(out[which(out$type == "SCAN" & out$alpha == alpha), 1])) / (nullnmc * 8)
    }
    
    c2 <- c2 + (p.adjust(pvalues2, "BH") < 0.05)
  }
  
  # Create a data frame to store power and standard error for SCAN
  df2 <- data.frame(time = factor(m2, levels = m2), power = c2 / pvalnmc, se = sqrt((c2 / pvalnmc) * (1 - c2 / pvalnmc) / pvalnmc), type = "SCAN", alpha = factor(alpha))
  
  # Add the results to the final data frame 'df'
  df <- rbind(df, df2)
}
dfs <- df
dfs <- dfs  %>% group_by(type)
mycol <- gg_color_hue(2)[2]

#adjacency matrices
#pairwise adjacency matrices
# Initialize an empty data frame to store the final results
df <- data.frame()
# Loop through different alpha values;
for(a in 1:length(alpha_values)){
  alpha = alpha_values[a]
  
  # Initialize a data frame to collect inner loop results for each alpha value
  out <- data.frame()
  
  for(i in 1:nullnmc){
    # Set the seed for reproducibility
    set.seed(123 + i - 1)
    
    # Generate time series data using the 'genTSGonepar' function
    glist <- genTSGonepar(n, nperturb, cperturb=0, rmin, rmax, tmax, d_true, alpha)
    
    # Run the 'doAdj' function on the generated data
    out2 <- doAdj(glist)
    
    # Create a temporary data frame to store the results and set its column names
    out.p <- data.frame(as.vector(t(out2$tnorm[-c(5,6,7)])), Method="DIST", alpha=factor(alpha))
    colnames(out.p)[1] <- "tnorm"
    
    # Append the results to the 'out' data frame
    out <- rbind(out, out.p)
  }
  
  # Initialize a vector to keep track of counts for each time point
  c2 <- rep(0, 11)
  
  # Loop to calculate p-values for each time point
  for(i in 1:pvalnmc){
    # Set the seed for reproducibility
    set.seed(123 + nullnmc + i - 1)
    
    # Generate time series data
    glist <- genTSGonepar(n, nperturb, cperturb, rmin, rmax, tmax, d_true, alpha)
    
    # Run the 'doAdj' function on the generated data
    out2 <- doAdj(glist)
    
    pvalues2 <- rep(0, (tmax-1))
    
    # Calculate p-values
    for(j in 1:11){
      pvalues2[j] <- sum(out2$tnorm[j] < sort(out[which(out$Method == "DIST" & out$alpha == alpha), 1])) / (nullnmc * 8)
    }
    
    # Adjust p-values and count the number of significant results
    c2 <- c2 + (p.adjust(pvalues2, "BH") < 0.05)
  }
  
  # Create a data frame to store the results for this alpha value
  df2 <- data.frame(time=factor(m2, levels=m2), power=c2/pvalnmc, se=sqrt((c2/pvalnmc)*(1 - c2/pvalnmc)/pvalnmc), Method="DIST", alpha=factor(alpha))
  
  # Append the results to the final 'df' data frame
  df <- rbind(df, df2)
}
dfadj <- df
dfadj <- dfadj  %>% group_by(Method)
mycol <- gg_color_hue(2)[2]
colnames(dfadj)[4] <- "type"
colnames(dfm)[4] <- "type"
dfm <- ungroup(dfm)
colnames(dfo)[4] <- "type"
dfo <- ungroup(dfo)


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


# plots without GP
dfmnull <- dfm %>% filter(time!="5:6"&time!="6:7" & time!="7:8")
dfonull <- dfo %>% filter(time!="5:6"&time!="6:7" & time!="7:8")
dfsnull <- dfs %>% filter(time!="5:6"&time!="6:7" & time!="7:8")
dfadjnull <- dfadj %>% filter(time!="5:6"&time!="6:7" & time!="7:8")

dfommasnull <- rbind(dfmnull, dfonull, dfsnull, dfadjnull)
dfommas_type_one_err <- dfommasnull[,c(2,4,5)]
dfommas_type_one_err$power <- dfommas_type_one_err$power * pvalnmc
# Summarizing the data
dfommas_type_one_err <- dfommas_type_one_err %>%
  group_by(type, alpha) %>%
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
dfommas_type_one_err <- dfommas_type_one_err %>%
  rename(Method = type, fpr=power_fpr, se=power_se)

# Creating the ggplot
ppp1_type_one_err <- ggplot(dfommas_type_one_err, aes(x=alpha, y=fpr, color=Method, group=Method)) +
  geom_errorbar(aes(ymin = fpr - se, ymax = fpr + se, width = .1, linetype = Method)) +
  geom_line(aes(linetype = Method)) +
  geom_point(aes(shape = Method)) +
  theme_bw() +
  xlab(TeX("$\\theta$")) +
  theme(
    plot.title = element_text(size = 10, face = 'bold'),
    legend.position = c(.8, 0.8),
    legend.text = element_text(size = 8, face = 'bold'),
    axis.title = element_text(face = "bold"),
    axis.text = element_text(face = "bold")
  )
dfs <- ungroup(dfs)

dfommas <- rbind(dfm, dfo, dfs,dfadj)
dfommas <- dfommas
dfommas <- dfommas %>% filter(time=="6:7")
colnames(dfommas)[4] <- "Method"

ppp1spower <- ggplot(dfommas, aes(x=alpha, y=power, color=Method,group=Method)) +
  geom_errorbar(aes(x=alpha,ymin=power-se, ymax=power+se, width=.1,linetype=Method, color=Method)) +
  geom_line(aes(x=alpha,y=power,linetype=Method)) +theme_bw()+xlab(TeX("${\\theta}$"))+
  geom_point(aes(x=alpha,y=power, shape=Method)) +
  theme(plot.title = element_text(size=10, face = 'bold'),legend.position = "none",axis.title = element_text(face="bold"),axis.text = element_text(face="bold"))  #+
# theme(plot.title = element_text(size=10, face = 'bold'),axis.title = element_text(face="bold"),axis.text = element_text(face="bold"))  #+

# plots with GP
dfommas <- rbind(dfm, dfo, dfs,dfadj)
dfommas <- dfommas %>% filter(time=="6:7")
colnames(dfommas)[4] <- "Method"
colnames(dfommas)[2] <- "mean"
colnames(dfommas_type_one_err)[3] <- "mean"
dfommas$metric <- "power"
dfommas_type_one_err$metric <- "FPR"

dfommas <- select(dfommas, colnames(dfommas_type_one_err))

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
dfgp_type_one_err <- select(dfgp_type_one_err, colnames(dfommas_type_one_err))

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
dfgp_power <- select(dfgp_power, colnames(dfommas_type_one_err))

# Merging and renaming for final graph
df.graph <- bind_rows(dfommas, dfommas_type_one_err, dfgp_power, dfgp_type_one_err)
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
  


# mrr
#pairwise mase
# Run the nested foreach loop using %:% for nesting and .combine to combine the results
df <- foreach(alpha = c(0, seq(8) / 8), .combine = 'rbind') %:% 
  foreach(i = 1:nullnmc, .combine = 'rbind') %dopar% {
    
    set.seed(123 + i - 1)
    
    # Generate time-series data
    glist <- genTSGonepar(n, nperturb, cperturb, rmin, rmax, tmax, d_true, alpha)
    
    # Apply MASE methods
    out2 <- doMase(glist, 2, d1, d3, center, approx, py, method3)
    out12 <- doMase(glist, 12, d1, d3, center, approx, py, method3)
    
    # Rank transformation
    rr2 <- apply(out2$pdist, 2, function(x) rank(1 / x))
    rr12 <- apply(out12$pdist, 2, function(x) rank(1 / x))
    
    # Melt and rename for MASE(2)
    df.rr2 <- melt(1 / (rr2))
    names(df.rr2) <- c("vertex", "time", "rr")
    df.rr2$time <- rep(factor(m2, levels = m2), each = n)
    df.rr2 <- cbind(df.rr2, mc = i, type = "MASE(2)", alpha = factor(alpha))
    
    # Melt and rename for MASE(12)
    df.rr12 <- melt(1 / (rr12))
    names(df.rr12) <- c("vertex", "time", "rr")
    df.rr12$time <- rep(factor(m2, levels = m2), each = n)
    df.rr12 <- cbind(df.rr12, mc = i, type = "MASE(12)", alpha = factor(alpha))
    
    # Combine the results
    df.p <- rbind(df.rr2, df.rr12)
    df.p
  }
dfmi <- df

#pairwise omni
# Run the nested foreach loop using %:% for nesting and .combine to combine the results
df <- foreach(alpha = c(0, seq(8) / 8), .combine = 'rbind') %:%
  foreach(i = 1:nullnmc, .combine = 'rbind') %dopar% {
    
    set.seed(123 + i - 1)
    
    # Generate time-series data
    glist <- genTSGonepar(n, nperturb, cperturb, rmin, rmax, tmax, d_true, alpha)
    
    # Apply Omni methods
    out2 <- doOmni(glist, nomni=2, dmax, center=center, approx=approx)
    out12 <- doOmni(glist, 12, dmax, center, approx)
    
    # Rank transformation
    rr2 <- apply(out2$pdist, 2, function(x) rank(1 / x))
    rr12 <- apply(out12$pdist, 2, function(x) rank(1 / x))
    
    # Melt and rename for OMNI(2)
    df.rr2 <- melt(1 / rr2)
    names(df.rr2) <- c("vertex", "time", "rr")
    df.rr2$time <- rep(factor(m2, levels = m2), each = n)
    df.rr2 <- cbind(df.rr2, mc = i, type = "OMNI(2)", alpha = factor(alpha))
    
    # Melt and rename for OMNI(12)
    df.rr12 <- melt(1 / rr12)
    names(df.rr12) <- c("vertex", "time", "rr")
    df.rr12$time <- rep(factor(m2, levels = m2), each = n)
    df.rr12 <- cbind(df.rr12, mc = i, type = "OMNI(12)", alpha = factor(alpha))
    
    # Combine the results
    df.p <- rbind(df.rr2, df.rr12)
    df.p
  }
# Assign the result to dfoi
dfoi <- df

#pairwise scan
# Run the nested foreach loop using %:% for nesting and .combine to combine the results
df <- foreach(alpha = c(0, seq(8) / 8), .combine = 'rbind') %:%
  foreach(i = 1:nullnmc, .combine = 'rbind') %dopar% {
    
    # Set seed
    set.seed(123 + i - 1)
    
    # Generate time-series data
    glist <- genTSGonepar(n, nperturb, cperturb, rmin, rmax, tmax, d_true, alpha)
    
    # Apply Scan method
    out2 <- doScan(glist)
    
    # Rank transformation
    rr2 <- apply(out2$pdist, 2, function(x) rank(1 / x))
    
    # Melt and rename
    df.rr2 <- melt(1 / rr2)
    names(df.rr2) <- c("vertex", "time", "rr")
    df.rr2$time <- rep(factor(m2, levels = m2), each = n)
    
    # Combine the results
    df.p <- cbind(df.rr2, mc = i, type = "SCAN", alpha = factor(alpha))
    df.p
  }
dfsi <- df

#pairwise Adj
df <- foreach(alpha = c(0, seq(8) / 8), .combine = 'rbind') %:% 
  foreach(i = 1:nullnmc, .combine = 'rbind') %dopar% {
    
    # Set seed for reproducibility
    set.seed(123 + i - 1)
    
    # Generate time-series graph
    glist <- genTSGonepar(n, nperturb, cperturb, rmin, rmax, tmax, d_true, alpha)
    
    # Apply Adj method
    out2 <- doAdj(glist)
    
    # Rank transformation
    rr2 <- apply(out2$pdist, 2, function(x) rank(1 / x))
    
    # Melt and rename
    df.rr2 <- melt(1 / rr2)
    names(df.rr2) <- c("vertex", "time", "rr")
    df.rr2$time <- rep(factor(m2, levels = m2), each = n)
    
    # Combine the results
    df.p <- cbind(df.rr2, mc = i, type = "DIST", alpha = factor(alpha))
    df.p
  }
dfadji <- df



dfi <- rbind(dfmi,dfoi, dfsi, dfadji)
dfi$type <- factor(dfi$type)
df.auc <- c()
for (checked_mc in seq(nullnmc)) {
  for (checked_alpha in c(0,seq(8)/8)) {
    for (checked_type in c("MASE(2)","OMNI(2)","MASE(12)","OMNI(12)","SCAN", "DIST")) {#"UASE(2)", "UASE(12)")) {#)) {
      tempdf <- dfi %>% filter(time=="6:7"&mc==checked_mc&alpha==checked_alpha&type==checked_type)
      temp <- data.frame(auc=pROC::roc(tempdf$vertex<=100, tempdf$rr, quiet=TRUE)$auc, mc=checked_mc, type=checked_type, alpha=factor(checked_alpha))
      df.auc <- rbind(df.auc, temp)
    }
  }
}

dfommasi.auc <- df.auc %>%
  group_by(type, alpha) %>%
  select(-2) %>%
  summarize(
    across(
      everything(),
      list(mean = mean, sd = sd, se = ~ sd(.)/sqrt(n()))
    )
  ) %>%
  rename(AUC = auc_mean, Method = type, sd=auc_sd, se=auc_se)

mycol <- gg_color_hue(2)[2]
ppp1auc <- ggplot(dfommasi.auc, aes(x=alpha, y=AUC, color=Method,group=Method)) +
  geom_errorbar(aes(x=alpha,ymin=AUC-se, ymax=AUC+se, width=.1,linetype=Method, color=Method)) +
  geom_line(aes(x=alpha,y=AUC,linetype=Method)) +theme_bw()+xlab(TeX("$\\theta$"))+
  # geom_point(aes(x=alpha,y=AUC, shape=Method)) + theme(plot.title = element_text(size=10, face = 'bold'),legend.position = c(0.8, 0.85), legend.text = element_text(size=8, face='bold'),axis.title = element_text(face="bold"),axis.text = element_text(face="bold"))  #+
  geom_point(aes(x=alpha,y=AUC, shape=Method)) + theme(plot.title = element_text(size=10, face = 'bold'), legend.text = element_text(size=8, face='bold'),axis.title = element_text(face="bold"),axis.text = element_text(face="bold"))  #+


dfpert <- dfi %>% filter(vertex<=100&time=="6:7")
dfnonpert <- dfi %>% filter(vertex>100&time=="6:7")
dfpert <- cbind(dfpert, rep("pert",dim(dfpert)[1]))
colnames(dfpert)[6+1] <- "vertex.class"
dfnonpert <- cbind(dfnonpert, rep("nonpert",dim(dfnonpert)[1]))
colnames(dfnonpert)[6+1] <- "vertex.class"
dfall <- rbind(dfpert,dfnonpert)
dfall$vertex.class <- factor(dfall$vertex.class, c("nonpert", "pert"),ordered = TRUE)
a = dfall[, -(1:2)] %>%
  group_by(alpha, mc, type, vertex.class) %>%
  summarize(across(everything(), mean))

dfa = a %>%
  group_by(alpha, mc, type) %>%
  mutate(rr_Diff = c(NA, diff(rr)))

dfa = dfa %>%
  filter(vertex.class == "pert")

dfommai = dfa %>%
  group_by(type, alpha) %>%
  select(-c(2,4,5)) %>%
  summarize(
    across(
      everything(),
      list(mean = mean, sd = sd, se = ~ sd(.) / sqrt(n()))
    )
  ) %>%
  rename(mrr = rr_Diff_mean, Method = type, sd=rr_Diff_sd, se=rr_Diff_se)
# colnames(dfommai)[1] <- "Method"
mycol <- gg_color_hue(2)[2]
ppp1mrr <- ggplot(dfommai, aes(x=alpha, y=mrr, color=Method,group=Method)) +
  geom_errorbar(aes(x=alpha,ymin=mrr-se, ymax=mrr+se, width=.1,linetype=Method, color=Method)) +
  geom_line(aes(x=alpha,y=mrr,linetype=Method)) +theme_bw()+xlab(TeX("$\\theta$"))+
  # geom_point(aes(x=alpha,y=mrr, shape=Method)) + theme(plot.title = element_text(size=10),legend.position = c(0.8, 0.6),legend.text = element_text(size=8)) +ylab("reciprocal rank difference")
  # geom_point(aes(x=alpha,y=mrr, shape=Method)) + theme(plot.title = element_text(size=10, face='bold'),legend.position = c(0.1, 0.25),legend.text = element_text(size=8, face='bold'),axis.title = element_text(face="bold"),axis.text = element_text(face="bold")) +ylab("reciprocal rank difference")
  geom_point(aes(x=alpha,y=mrr, shape=Method)) + theme(plot.title = element_text(size=10, face='bold'),legend.position = "none",legend.text = element_text(size=8, face='bold'),axis.title = element_text(face="bold"),axis.text = element_text(face="bold")) +ylab("reciprocal rank difference")

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
  theme(plot.title = element_text(size=10, face = 'bold'), legend.position = "none",axis.title = element_text(face="bold"),axis.text = element_text(face="bold"))  + ylab("")#+
# theme(plot.title = element_text(size=10, face = 'bold'), legend.text = element_text(size=8, face='bold'),axis.title = element_text(face="bold"),axis.text = element_text(face="bold"))  + ylab("")#+


# Figure 5
grid.arrange(pg, pv, nrow=1, widths=c(5,4))
png(paste0("comb","_",format(Sys.time(), "%Y_%m_%d.png")), width = 9, height = 4, units="in", res=400)
grid.arrange(pg, pv, nrow=1, widths=c(5,4))
dev.off()

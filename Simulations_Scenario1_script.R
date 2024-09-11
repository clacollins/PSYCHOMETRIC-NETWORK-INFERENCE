##################################
### Mental Health Graphs Inference
### Comparing different methods for graph inference in the context of Mental Health Networks
### Simulations - SCENARIO 1 
##################################
rm(list=ls())
set.seed(1234)

# install.packages("optparse")
library("optparse")
library("igraph")
library('BDgraph')
library("qgraph")
library('bootnet')
library('GGMnonreg')
library('GeneNet')
library("PAsso")
library("psychonetrics")
library("BGGM")
library('e1071') 




##########
# Parameters 
nsims <- 100 # number of repetitions - for each simulation setting - Should be fixed to 100

#######
# Parameters to modulate through the Rscript command 

option_list <- list(
make_option( c("--nbind","-n"), type = "integer", default = 100,
            help = "number of individuals"),
make_option(c("--symptom","-p"), type = "integer", default = 6,
            help = "number of ordinal symptoms"),
make_option(c("--Levels","-L"), type = "integer", default = 4,
            help = "number of levels per symptom"),
make_option(c("--wishart","-w"), type = "numeric", default = 1,
            help = "density of the Wishart graph")
)


opt <- parse_args(OptionParser(option_list = option_list))

n <- opt$n
p <- opt$symptom
L <- opt$L
prob <- opt$w

# Usage :
# Rscript Simulations_Scenario1_script.R --n 100 --symptom 6 --L 4 --wishart 0.5
##############


# defining names for the variables below
names<-c()
for (i in 1:p){
  names <- c(names,paste0("V",i))
}

#########
# Conventions 
# Sigma are covariance matrices ; Rho are correlation matrices =-cov2cor(Sigma)
# Omega are precision matrices = inverse of Sigma ; Theta are partial correlation matrices =-cov2cor(Omega) =-cov2cor(solve(Rho))


##########
# Initialize values for performance measures

# First simulation setting - Gaussian setting + binning
# MSE for all methods 
MSE.poly.mle <- rep(NA, nsims)
MSE.poly.wls <- rep(NA, nsims)
MSE.pearson <- rep(NA, nsims)
MSE.Glasso.poly <- rep(NA, nsims)
MSE.ggmModSel.poly <- rep(NA, nsims)
MSE.Glasso.pearson <- rep(NA, nsims)
MSE.ggmModSel.pearson <- rep(NA, nsims)
MSE.ggmSS <- rep(NA, nsims)
MSE.GGMnonreg.neighsel <- rep(NA, nsims)
MSE.GGMnonreg.boot.pearson <- rep(NA, nsims)
MSE.GGMnonreg.boot.poly <- rep(NA, nsims)
MSE.BGGM.explore <- rep(NA, nsims)
MSE.BGGM.estimate <- rep(NA, nsims)
MSE.passo <- rep(NA, nsims)

# True positives 
TP.poly.mle <- rep(NA, nsims)
TP.poly.wls <- rep(NA, nsims)
TP.pearson <- rep(NA, nsims)
TP.Glasso.poly <- rep(NA, nsims)
TP.ggmModSel.poly <- rep(NA, nsims)
TP.Glasso.pearson <- rep(NA, nsims)
TP.ggmModSel.pearson <- rep(NA, nsims)
TP.ggmSS <- rep(NA, nsims)
TP.GGMnonreg.neighsel <- rep(NA, nsims)
TP.GGMnonreg.boot.pearson <- rep(NA, nsims)
TP.GGMnonreg.boot.poly <- rep(NA, nsims)
TP.BGGM.explore <- rep(NA, nsims)
TP.BGGM.estimate <- rep(NA, nsims)
TP.passo <- rep(NA, nsims)

# True negatives 
TN.poly.mle <- rep(NA, nsims)
TN.poly.wls <- rep(NA, nsims)
TN.pearson <- rep(NA, nsims)
TN.Glasso.poly <- rep(NA, nsims)
TN.ggmModSel.poly <- rep(NA, nsims)
TN.Glasso.pearson <- rep(NA, nsims)
TN.ggmModSel.pearson <- rep(NA, nsims)
TN.ggmSS <- rep(NA, nsims)
TN.GGMnonreg.neighsel <- rep(NA, nsims)
TN.GGMnonreg.boot.pearson <- rep(NA, nsims)
TN.GGMnonreg.boot.poly <- rep(NA, nsims)
TN.BGGM.explore <- rep(NA, nsims)
TN.BGGM.estimate <- rep(NA, nsims)
TN.passo <- rep(NA, nsims)

# False positives 
FP.poly.mle <- rep(NA, nsims)
FP.poly.wls <- rep(NA, nsims)
FP.pearson <- rep(NA, nsims)
FP.Glasso.poly <- rep(NA, nsims)
FP.ggmModSel.poly <- rep(NA, nsims)
FP.Glasso.pearson <- rep(NA, nsims)
FP.ggmModSel.pearson <- rep(NA, nsims)
FP.ggmSS <- rep(NA, nsims)
FP.GGMnonreg.neighsel <- rep(NA, nsims)
FP.GGMnonreg.boot.pearson <- rep(NA, nsims)
FP.GGMnonreg.boot.poly <- rep(NA, nsims)
FP.BGGM.explore <- rep(NA, nsims)
FP.BGGM.estimate <- rep(NA, nsims)
FP.passo <- rep(NA, nsims)

# sensitivity
sensitivity.poly.mle <- rep(NA, nsims)
sensitivity.poly.wls <- rep(NA, nsims)
sensitivity.pearson <- rep(NA, nsims)
sensitivity.Glasso.poly <- rep(NA, nsims)
sensitivity.ggmModSel.poly <- rep(NA, nsims)
sensitivity.Glasso.pearson <- rep(NA, nsims)
sensitivity.ggmModSel.pearson <- rep(NA, nsims)
sensitivity.ggmSS <- rep(NA, nsims)
sensitivity.GGMnonreg.neighsel <- rep(NA, nsims)
sensitivity.GGMnonreg.boot.pearson <- rep(NA, nsims)
sensitivity.GGMnonreg.boot.poly <- rep(NA, nsims)
sensitivity.BGGM.explore <- rep(NA, nsims)
sensitivity.BGGM.estimate <- rep(NA, nsims)
sensitivity.passo <- rep(NA, nsims)

# specificity
specificity.poly.mle <- rep(NA, nsims)
specificity.poly.wls <- rep(NA, nsims)
specificity.pearson <- rep(NA, nsims)
specificity.Glasso.poly <- rep(NA, nsims)
specificity.ggmModSel.poly <- rep(NA, nsims)
specificity.Glasso.pearson <- rep(NA, nsims)
specificity.ggmModSel.pearson <- rep(NA, nsims)
specificity.ggmSS <- rep(NA, nsims)
specificity.GGMnonreg.neighsel <- rep(NA, nsims)
specificity.GGMnonreg.boot.pearson <- rep(NA, nsims)
specificity.GGMnonreg.boot.poly <- rep(NA, nsims)
specificity.BGGM.explore <- rep(NA, nsims)
specificity.BGGM.estimate <- rep(NA, nsims)
specificity.passo <- rep(NA, nsims)

# precision
precision.poly.mle <- rep(NA, nsims)
precision.poly.wls <- rep(NA, nsims)
precision.pearson <- rep(NA, nsims)
precision.Glasso.poly <- rep(NA, nsims)
precision.ggmModSel.poly <- rep(NA, nsims)
precision.Glasso.pearson <- rep(NA, nsims)
precision.ggmModSel.pearson <- rep(NA, nsims)
precision.ggmSS <- rep(NA, nsims)
precision.GGMnonreg.neighsel <- rep(NA, nsims)
precision.GGMnonreg.boot.pearson <- rep(NA, nsims)
precision.GGMnonreg.boot.poly <- rep(NA, nsims)
precision.BGGM.explore <- rep(NA, nsims)
precision.BGGM.estimate <- rep(NA, nsims)
precision.passo <- rep(NA, nsims)


#########
# Initialize skewness of the dataset
skew <- data.frame(matrix(NA,nrow=nsims,ncol=p))

# initialize errors - in case a method does not work
poly.mle.errors <- 0
poly.wls.errors <- 0
pearson.errors <- 0
Glasso.poly.errors <- 0
ggmModSel.poly.errors <- 0
Glasso.pearson.errors <- 0
ggmModSel.pearson.errors <- 0
ggmSS.errors <- 0
GGMnonreg.neighsel.errors <- 0
GGMnonreg.boot.pearson.errors <- 0
GGMnonreg.boot.poly.errors <- 0
BGGM.explore.errors <- 0 
BGGM.estimate.errors <- 0
passo.errors <- 0

####################
# Start the simulations 
for (repet in 1:nsims){
  
  #######################
  ### SCENARIO 1: based on G-Wishart distributions for polychoric partial correlations
  #######################
  ### 1. Simulate the data 
  # Parameters of Wishart and density of the graph 
  b <- 3 # degree of freedom of Wishart distr. - default value 
  D <- diag(p) # scale matrix - Identity - default value 
  # graph skeleton 
  skeleton <- as.matrix(as_adjacency_matrix(erdos.renyi.game(n=p,p=prob)))
  # The TRUE precision matrix is set as  
  Omega <- round(rgwish(n=1,adj=skeleton,b=b,D=D),2)# it's needed to threshold the values to obtain something sparse
  # And corresponding TRUE partial polychoric correlation matrix is 
  Theta <- -cov2cor(Omega)
  # Simulate the data 
  X <- bdgraph.sim(p=p,graph='fixed',n=n,type="categorical", cut= L,K=Omega)
  
  # But we consider only the upper diagonal part - to remove redundant information   
  tru.par.cor <- Theta[upper.tri(Theta,diag=FALSE)]
  # We record the sparsity of true graph  
  tru.nb.edges <- sum(tru.par.cor !=0)
  # and skewness of distribution
  skew[repet,] <-lapply(as.data.frame(X$data),skewness)
  
  ##########
  ### 2. Estimate the partial correlations 
  
  ###
  ## 2.1 First family of methods: Directly estimate polychoric correlations
  # 2.1.a First variant : Maximum likelihood estimation of polychoric correlations
  # Below variables are transformed to ordinal - necessary to run cor_auto correctly
  Rho.estim.poly.mle <- cor_auto(data=as.data.frame(lapply(as.data.frame(X$data),ordered))) # estimated polychoric correlations 
  res <- try(-cov2cor(solve(Rho.estim.poly.mle)))  # estimated partial polychoric correlations
  if(inherits(res, "try-error")){
    #error handling code, count number of errors
    mle.errors <- mle.errors+1
    Theta.estim.poly.mle <- matrix(NA,p,p)
  } else {
    Theta.estim.poly.mle <- res
  }
  
  # 2.1.b Second variant: Weighted least squares estimation of polychoric correlations 
  mod <- try(ggm(data=X$data, ordered = TRUE, missing = "pairwise"))
  if(inherits(mod, "try-error")){
    #error handling code, count number of errors
    poly.wls.errors <- poly.wls.errors+1
    Theta.estim.poly.wls <- matrix(NA,p,p)
  } else {
    # Estimate model:
    mod <- try(mod %>% prune %>% modelsearch) # If L is large, outputs an error here
    if(inherits(mod, "try-error")){
      #error handling code, count number of errors
      poly.wls.errors <- poly.wls.errors+1
      Theta.estim.poly.wls <- matrix(NA,p,p)
    } else {
      # Obtain network:
      Theta.estim.poly.wls <- getmatrix(mod, "omega")
    }
  }
  
  
  # 2.1.c Instead of polychoric, look at Pearson's correlation (then treat the variables as continuous)
  Rho.estim.pearson <- cor(X$data)
  res <- try(-cov2cor(solve(Rho.estim.pearson)))
  if(inherits(res, "try-error")){
    #error handling code, count number of errors
    pearson.errors <- pearson.errors+1
    Theta.estim.pearson <- matrix(NA,p,p)
  } else {
    Theta.estim.pearson <- res
  }  
  
  ###
  ## 2.2 Second family of methods: Regularized GGM relying on polychoric correlations
  # 2.2.a First variant : EBICglasso 
  # the input is the polychoric correlation computed above - we choose MLE   
  res <- try(EBICglasso(Rho.estim.poly.mle, n=n, 
                        gamma=0.5,  # EBIC tuning parameter
                        nlambda = 100, # Number of lambda values to test
                        lambda.min.ratio = 0.01, # Ratio of lowest lambda value compared to maximal lambda
  ))
  if(inherits(res, "try-error")){
    #error handling code, count number of errors
    Glasso.poly.errors <- Glasso.poly.errors+1
    Theta.estim.Glasso.poly <- matrix(NA,p,p)
  } else {
    Theta.estim.Glasso.poly <- res
  }  
  
  
  # 2.2.b Second variant: ggmModSelect 
  # the input is the polychoric correlation computed above  
  res <- try(ggmModSelect(S=Rho.estim.poly.mle,n=n)$graph)
  if(inherits(res, "try-error")){
    #error handling code, count number of errors
    ggmModSel.poly.errors <- ggmModSel.poly.errors+1
    Theta.estim.ggmModSel.poly <- matrix(NA,p,p)
  } else {
    Theta.estim.ggmModSel.poly <- res
  }  
  
  
  ###
  ## 2.3 Third family of methods: Regularized GGM relying on Pearson's correlations
  # 2.3.a First variant : EBICglasso 
  # the input is the Pearson correlation computed above  
  res <- try(EBICglasso(Rho.estim.pearson, n=n, 
                        gamma=0.5,  # EBIC tuning parameter
                        nlambda = 100, # Number of lambda values to test
                        lambda.min.ratio = 0.01, # Ratio of lowest lambda value compared to maximal lambda
  ))
  if(inherits(res, "try-error")){
    #error handling code, count number of errors
    Glasso.pearson.errors <- Glasso.pearson.errors+1
    Theta.estim.Glasso.pearson <- matrix(NA,p,p)
  } else {
    Theta.estim.Glasso.pearson <- res
  }
  
  
  
  # 2.3.b Second variant: ggmModSelect 
  # the input is the Pearson correlation computed above  
  res <- try(ggmModSelect(S=Rho.estim.pearson,n=n)$graph)
  if(inherits(res, "try-error")){
    #error handling code, count number of errors
    ggmModSel.pearson.errors <- ggmModSel.pearson.errors+1
    Theta.estim.ggmModSel.pearson <- matrix(NA,p,p)
  } else {
    Theta.estim.ggmModSel.pearson <- res
  }
  
  
  # 2.3.c Third variant: method of Schafer and Strimmer
  # Input is original data - so we cannot choose to input polychoric correlations. It relies on Pearson's correlations
  res <- try(ggm.estimate.pcor(X$data))
  if(inherits(res, "try-error")){
    #error handling code, count number of errors
    ggmSS.errors <- ggmSS.errors+1
    Theta.estim.ggmSS <- matrix(NA,p,p)
  } else {
    Theta.estim.ggmSS <- res
  }
  
  ###
  ## 2.4 Fourth family of methods: Non Regularized GGM   
  
  # 2.4.a First approach - do p regressions and combine them with neighborhood selection - relying on Pearson correlation
  res <- try(ggm_inference(X$data,boot=FALSE, method="pearson")$pcors)
  if(inherits(res, "try-error")){
    #error handling code, count number of errors
    GGMnonreg.neighsel.errors <- GGMnonreg.neighsel.errors+1
    Theta.estim.GGMnonreg.neighsel <- matrix(NA,p,p)
  } else {
    Theta.estim.GGMnonreg.neighsel <- res
  }
  
  # 2.4.b Second approach - relying on bootstrap and Pearson corr
  res <- ggm_inference(X$data,boot=T, method="pearson")$pcors
  if(inherits(res, "try-error")){
    #error handling code, count number of errors
    GGMnonreg.boot.pearson.errors <- GGMnonreg.boot.pearson.errors+1
    Theta.estim.GGMnonreg.boot.pearson <- matrix(NA,p,p)
  } else {
    Theta.estim.GGMnonreg.boot.pearson <- res
  }
  
  # 2.4.c Third approach - relying on bootstrap and polychoric corr
  # WARNING THIS IS VERY LONG !!!!
  res <- try(ggm_inference(X$data,boot=T, method="polychoric")$pcors)
  if(inherits(res, "try-error"))
  {
    #error handling code, count number of errors
    GGMnonreg.boot.poly.errors <- GGMnonreg.boot.poly.errors+1
    Theta.estim.GGMnonreg.boot.poly <- matrix(NA,p,p)
  } else {
    Theta.estim.GGMnonreg.boot.poly <- res
  }
  
  # 2.4.d Fourth approach - package BGGM - same methods as in Isvoranu & Esk
  # With this method ordinal data should be categories not including 0 values - so we add 1 to the data matrix (does not change anything) 
  fit1 <- try(BGGM::explore(Y=X$data+1,type="ordinal",impute = FALSE))
  if(inherits(fit1, "try-error"))
  {
    #error handling code, count number of errors
    BGGM.explore.errors <- BGGM.explore.errors+1
    Theta.estim.BGGM.explore <- matrix(NA,p,p)
  } else {
    fit1 <- BGGM::select(fit1)
    Theta.estim.BGGM.explore <- fit1$pcor_mat_zero
  }
  
  
  fit2 <- try(BGGM::estimate(Y=X$data+1,type="ordinal",impute = FALSE)) 
  if(inherits(fit2, "try-error"))
  {
    #error handling code, count number of errors
    BGGM.estimate.errors <- BGGM.estimate.errors+1
    Theta.estim.BGGM.estimate <- matrix(NA,p,p)
  } else {
    fit2 <- BGGM::select(fit2)
    Theta.estim.BGGM.estimate <- fit2$pcor_adj
  }
  
  
  ###
  ## .2.5 Fifth method: Passo - Surrogate residuals method by Liu et al 
  X$data <- as.data.frame(X$data) # data needs to be structured as data.frame
  Theta.estim.passo <- diag(p)/2 # when symetrizing, diagonal doubles
  # partial correlations need to be computed pairwisely
  for (i in 1:(p-1)){
    for (j in (i+1):p){
      resp <-c(names[i],names[j])
      adj <-names[-c(i,j)] 
      res <- try(PAsso(responses =resp ,adjustments =adj, data = X$data))
      if(inherits(res, "try-error"))
      {
        #error handling code, count number of errors
        #warning these errors are counted for every (i,j)
        passo.errors <- passo.errors+1
        Theta.estim.passo[i,j] <- NA
      } else {Theta.estim.passo[i,j] <- res$corr[1,2]}
    }
  }
  # symetrize the matrix 
  Theta.estim.passo <- Theta.estim.passo+t(Theta.estim.passo)
  
  
  ### .3. Performance evaluation
  # Compute MSE between true vector and estimated ones, for each method 
  MSE.poly.mle[repet] <- sum(Theta.estim.poly.mle[upper.tri(Theta.estim.poly.mle,diag=FALSE)] - tru.par.cor)^2
  MSE.poly.wls[repet] <- sum(Theta.estim.poly.wls[upper.tri(Theta.estim.poly.wls,diag=FALSE)] - tru.par.cor)^2
  MSE.pearson[repet] <- sum(Theta.estim.pearson[upper.tri(Theta.estim.pearson,diag=FALSE)] - tru.par.cor)^2
  MSE.Glasso.poly[repet] <- sum(Theta.estim.Glasso.poly[upper.tri(Theta.estim.Glasso.poly,diag=FALSE)] - tru.par.cor)^2
  MSE.ggmModSel.poly[repet] <- sum(Theta.estim.ggmModSel.poly[upper.tri(Theta.estim.ggmModSel.poly,diag=FALSE)] - tru.par.cor)^2
  MSE.Glasso.pearson[repet] <- sum(Theta.estim.Glasso.pearson[upper.tri(Theta.estim.Glasso.pearson,diag=FALSE)] - tru.par.cor)^2
  MSE.ggmModSel.pearson[repet] <- sum(Theta.estim.ggmModSel.pearson[upper.tri(Theta.estim.ggmModSel.pearson,diag=FALSE)] - tru.par.cor)^2
  MSE.ggmSS[repet] <- sum(Theta.estim.ggmSS[upper.tri(Theta.estim.ggmSS,diag=FALSE)] - tru.par.cor)^2
  MSE.GGMnonreg.neighsel[repet] <- sum(Theta.estim.GGMnonreg.neighsel[upper.tri(Theta.estim.GGMnonreg.neighsel,diag=FALSE)] - tru.par.cor)^2
  MSE.GGMnonreg.boot.pearson[repet] <- sum(Theta.estim.GGMnonreg.boot.pearson[upper.tri(Theta.estim.GGMnonreg.boot.pearson,diag=FALSE)] - tru.par.cor)^2
  MSE.GGMnonreg.boot.poly[repet] <- sum(Theta.estim.GGMnonreg.boot.poly[upper.tri(Theta.estim.GGMnonreg.boot.poly,diag=FALSE)] - tru.par.cor)^2
  MSE.BGGM.explore[repet] <- sum(Theta.estim.BGGM.explore[upper.tri(Theta.estim.BGGM.explore,diag=FALSE)] - tru.par.cor)^2
  MSE.BGGM.estimate[repet] <- sum(Theta.estim.BGGM.estimate[upper.tri(Theta.estim.BGGM.estimate,diag=FALSE)] - tru.par.cor)^2
  MSE.passo[repet] <- sum(Theta.estim.passo[upper.tri(Theta.estim.passo,diag=FALSE)] - tru.par.cor)^2
  # Compute TP and TN and FP - note that FN=n-TP-TN-FP
  TP.poly.mle[repet] <- sum((Theta.estim.poly.mle[upper.tri(Theta.estim.poly.mle,diag=FALSE)]!=0)*(tru.par.cor!=0))
  TP.poly.wls[repet] <- sum((Theta.estim.poly.wls[upper.tri(Theta.estim.poly.wls,diag=FALSE)]!=0)*(tru.par.cor!=0))
  TP.pearson[repet] <- sum((Theta.estim.pearson[upper.tri(Theta.estim.pearson,diag=FALSE)]!=0)*(tru.par.cor!=0))
  TP.Glasso.poly[repet] <- sum((Theta.estim.Glasso.poly[upper.tri(Theta.estim.Glasso.poly,diag=FALSE)]!=0)*(tru.par.cor!=0))
  TP.ggmModSel.poly[repet] <- sum((Theta.estim.ggmModSel.poly[upper.tri(Theta.estim.ggmModSel.poly,diag=FALSE)]!=0)*(tru.par.cor!=0))
  TP.Glasso.pearson[repet] <- sum((Theta.estim.Glasso.pearson[upper.tri(Theta.estim.Glasso.pearson,diag=FALSE)]!=0)*(tru.par.cor!=0))
  TP.ggmModSel.pearson[repet] <- sum((Theta.estim.ggmModSel.pearson[upper.tri(Theta.estim.ggmModSel.pearson,diag=FALSE)]!=0)*(tru.par.cor!=0))
  TP.ggmSS[repet] <- sum((Theta.estim.ggmSS[upper.tri(Theta.estim.ggmSS,diag=FALSE)]!=0)*(tru.par.cor!=0))
  TP.GGMnonreg.neighsel[repet] <- sum((Theta.estim.GGMnonreg.neighsel[upper.tri(Theta.estim.GGMnonreg.neighsel,diag=FALSE)]!=0)*(tru.par.cor!=0))
  TP.GGMnonreg.boot.pearson[repet] <- sum((Theta.estim.GGMnonreg.boot.pearson[upper.tri(Theta.estim.GGMnonreg.boot.pearson,diag=FALSE)]!=0)*(tru.par.cor!=0))
  TP.GGMnonreg.boot.poly[repet] <- sum((Theta.estim.GGMnonreg.boot.poly[upper.tri(Theta.estim.GGMnonreg.boot.poly,diag=FALSE)]!=0)*(tru.par.cor!=0))
  TP.BGGM.explore[repet] <- sum((Theta.estim.BGGM.explore[upper.tri(Theta.estim.BGGM.explore,diag=FALSE)]!=0)*(tru.par.cor!=0))
  TP.BGGM.estimate[repet] <- sum((Theta.estim.BGGM.estimate[upper.tri(Theta.estim.BGGM.estimate,diag=FALSE)]!=0)*(tru.par.cor!=0))
  TP.passo[repet] <- sum((Theta.estim.passo[upper.tri(Theta.estim.passo,diag=FALSE)]!=0)*(tru.par.cor!=0))
  ## TN
  TN.poly.mle[repet] <- sum((Theta.estim.poly.mle[upper.tri(Theta.estim.poly.mle,diag=FALSE)]==0)*(tru.par.cor==0))
  TN.poly.wls[repet] <- sum((Theta.estim.poly.wls[upper.tri(Theta.estim.poly.wls,diag=FALSE)]==0)*(tru.par.cor==0))
  TN.pearson[repet] <- sum((Theta.estim.pearson[upper.tri(Theta.estim.pearson,diag=FALSE)]==0)*(tru.par.cor==0))
  TN.Glasso.poly[repet] <- sum((Theta.estim.Glasso.poly[upper.tri(Theta.estim.Glasso.poly,diag=FALSE)]==0)*(tru.par.cor==0))
  TN.ggmModSel.poly[repet] <- sum((Theta.estim.ggmModSel.poly[upper.tri(Theta.estim.ggmModSel.poly,diag=FALSE)]==0)*(tru.par.cor==0))
  TN.Glasso.pearson[repet] <- sum((Theta.estim.Glasso.pearson[upper.tri(Theta.estim.Glasso.pearson,diag=FALSE)]==0)*(tru.par.cor==0))
  TN.ggmModSel.pearson[repet] <- sum((Theta.estim.ggmModSel.pearson[upper.tri(Theta.estim.ggmModSel.pearson,diag=FALSE)]==0)*(tru.par.cor==0))
  TN.ggmSS[repet] <- sum((Theta.estim.ggmSS[upper.tri(Theta.estim.ggmSS,diag=FALSE)]==0)*(tru.par.cor==0))
  TN.GGMnonreg.neighsel[repet] <- sum((Theta.estim.GGMnonreg.neighsel[upper.tri(Theta.estim.GGMnonreg.neighsel,diag=FALSE)]==0)*(tru.par.cor==0))
  TN.GGMnonreg.boot.pearson[repet] <- sum((Theta.estim.GGMnonreg.boot.pearson[upper.tri(Theta.estim.GGMnonreg.boot.pearson,diag=FALSE)]==0)*(tru.par.cor==0))
  TN.GGMnonreg.boot.poly[repet] <- sum((Theta.estim.GGMnonreg.boot.poly[upper.tri(Theta.estim.GGMnonreg.boot.poly,diag=FALSE)]==0)*(tru.par.cor==0))
  TN.BGGM.explore[repet] <- sum((Theta.estim.BGGM.explore[upper.tri(Theta.estim.BGGM.explore,diag=FALSE)]==0)*(tru.par.cor==0))
  TN.BGGM.estimate[repet] <- sum((Theta.estim.BGGM.estimate[upper.tri(Theta.estim.passo,diag=FALSE)]==0)*(tru.par.cor==0))
  TN.passo[repet] <- sum((Theta.estim.passo[upper.tri(Theta.estim.passo,diag=FALSE)]==0)*(tru.par.cor==0))
  # FP
  FP.poly.mle[repet] <- sum((Theta.estim.poly.mle[upper.tri(Theta.estim.poly.mle,diag=FALSE)]!=0)*(tru.par.cor==0))
  FP.poly.wls[repet] <- sum((Theta.estim.poly.wls[upper.tri(Theta.estim.poly.wls,diag=FALSE)]!=0)*(tru.par.cor==0))
  FP.pearson[repet] <- sum((Theta.estim.pearson[upper.tri(Theta.estim.pearson,diag=FALSE)]!=0)*(tru.par.cor==0))
  FP.Glasso.poly[repet] <- sum((Theta.estim.Glasso.poly[upper.tri(Theta.estim.Glasso.poly,diag=FALSE)]!=0)*(tru.par.cor==0))
  FP.ggmModSel.poly[repet] <- sum((Theta.estim.ggmModSel.poly[upper.tri(Theta.estim.ggmModSel.poly,diag=FALSE)]!=0)*(tru.par.cor==0))
  FP.Glasso.pearson[repet] <- sum((Theta.estim.Glasso.pearson[upper.tri(Theta.estim.Glasso.pearson,diag=FALSE)]!=0)*(tru.par.cor==0))
  FP.ggmModSel.pearson[repet] <- sum((Theta.estim.ggmModSel.pearson[upper.tri(Theta.estim.ggmModSel.pearson,diag=FALSE)]!=0)*(tru.par.cor==0))
  FP.ggmSS[repet] <- sum((Theta.estim.ggmSS[upper.tri(Theta.estim.ggmSS,diag=FALSE)]!=0)*(tru.par.cor==0))
  FP.GGMnonreg.neighsel[repet] <- sum((Theta.estim.GGMnonreg.neighsel[upper.tri(Theta.estim.GGMnonreg.neighsel,diag=FALSE)]!=0)*(tru.par.cor==0))
  FP.GGMnonreg.boot.pearson[repet] <- sum((Theta.estim.GGMnonreg.boot.pearson[upper.tri(Theta.estim.GGMnonreg.boot.pearson,diag=FALSE)]!=0)*(tru.par.cor==0))
  FP.GGMnonreg.boot.poly[repet] <- sum((Theta.estim.GGMnonreg.boot.poly[upper.tri(Theta.estim.GGMnonreg.boot.poly,diag=FALSE)]!=0)*(tru.par.cor==0))
  FP.BGGM.explore[repet] <- sum((Theta.estim.BGGM.explore[upper.tri(Theta.estim.BGGM.explore,diag=FALSE)]!=0)*(tru.par.cor==0))
  FP.BGGM.estimate[repet] <- sum((Theta.estim.BGGM.estimate[upper.tri(Theta.estim.BGGM.estimate,diag=FALSE)]!=0)*(tru.par.cor==0))
  FP.passo[repet] <- sum((Theta.estim.passo[upper.tri(Theta.estim.passo,diag=FALSE)]!=0)*(tru.par.cor==0))
  
  sensitivity.poly.mle[repet] <- TP.poly.mle[repet]/tru.nb.edges 
  sensitivity.poly.wls[repet] <- TP.poly.wls[repet]/tru.nb.edges
  sensitivity.pearson[repet] <- TP.pearson[repet]/tru.nb.edges
  sensitivity.Glasso.poly[repet] <- TP.Glasso.poly[repet]/tru.nb.edges
  sensitivity.ggmModSel.poly[repet] <- TP.ggmModSel.poly[repet]/tru.nb.edges
  sensitivity.Glasso.pearson[repet] <- TP.Glasso.pearson[repet]/tru.nb.edges
  sensitivity.ggmModSel.pearson[repet] <- TP.ggmModSel.pearson[repet]/tru.nb.edges
  sensitivity.ggmSS[repet] <- TP.ggmSS[repet]/tru.nb.edges
  sensitivity.GGMnonreg.neighsel[repet] <- TP.GGMnonreg.neighsel[repet]/tru.nb.edges
  sensitivity.GGMnonreg.boot.pearson[repet] <- TP.GGMnonreg.boot.pearson[repet]/tru.nb.edges
  sensitivity.GGMnonreg.boot.poly[repet] <- TP.GGMnonreg.boot.poly[repet]/tru.nb.edges
  sensitivity.BGGM.explore[repet] <- TP.BGGM.explore[repet]/tru.nb.edges
  sensitivity.BGGM.estimate[repet] <- TP.BGGM.estimate[repet]/tru.nb.edges
  sensitivity.passo[repet] <- TP.passo[repet]/tru.nb.edges
  
  specificity.poly.mle[repet] <- TN.poly.mle[repet]/(p*(p-1)/2-tru.nb.edges) 
  specificity.poly.wls[repet] <- TN.poly.wls[repet]/(p*(p-1)/2-tru.nb.edges)
  specificity.pearson[repet] <- TN.pearson[repet]/(p*(p-1)/2-tru.nb.edges)
  specificity.Glasso.poly[repet] <- TN.Glasso.poly[repet]/(p*(p-1)/2-tru.nb.edges)
  specificity.ggmModSel.poly[repet] <- TN.ggmModSel.poly[repet]/(p*(p-1)/2-tru.nb.edges)
  specificity.Glasso.pearson[repet] <- TN.Glasso.pearson[repet]/(p*(p-1)/2-tru.nb.edges)
  specificity.ggmModSel.pearson[repet] <- TN.ggmModSel.pearson[repet]/(p*(p-1)/2-tru.nb.edges)
  specificity.ggmSS[repet] <- TN.ggmSS[repet]/(p*(p-1)/2-tru.nb.edges)
  specificity.GGMnonreg.neighsel[repet] <- TN.GGMnonreg.neighsel[repet]/(p*(p-1)/2-tru.nb.edges)
  specificity.GGMnonreg.boot.pearson[repet] <- TN.GGMnonreg.boot.pearson[repet]/(p*(p-1)/2-tru.nb.edges)
  specificity.GGMnonreg.boot.poly[repet] <- TN.GGMnonreg.boot.poly[repet]/(p*(p-1)/2-tru.nb.edges)
  specificity.BGGM.explore[repet] <- TN.BGGM.explore[repet]/(p*(p-1)/2-tru.nb.edges)
  specificity.BGGM.estimate[repet] <- TN.BGGM.estimate[repet]/(p*(p-1)/2-tru.nb.edges)
  specificity.passo[repet] <- TN.passo[repet]/(p*(p-1)/2-tru.nb.edges)
  
  precision.poly.mle[repet] <- TP.poly.mle[repet]/sum(Theta.estim.poly.mle[upper.tri(Theta.estim.poly.mle,diag=FALSE)]!=0)
  precision.poly.wls[repet] <- TP.poly.wls[repet]/sum(Theta.estim.poly.wls[upper.tri(Theta.estim.poly.wls,diag=FALSE)]!=0)
  precision.pearson[repet] <- TP.pearson[repet]/sum(Theta.estim.pearson[upper.tri(Theta.estim.pearson,diag=FALSE)]!=0)
  precision.Glasso.poly[repet] <- TP.Glasso.poly[repet]/sum(Theta.estim.Glasso.poly[upper.tri(Theta.estim.Glasso.poly,diag=FALSE)]!=0)
  precision.ggmModSel.poly[repet] <- TP.ggmModSel.poly[repet]/sum(Theta.estim.ggmModSel.poly[upper.tri(Theta.estim.ggmModSel.poly,diag=FALSE)]!=0)
  precision.Glasso.pearson[repet] <- TP.Glasso.pearson[repet]/sum(Theta.estim.Glasso.pearson[upper.tri(Theta.estim.Glasso.pearson,diag=FALSE)]!=0)
  precision.ggmModSel.pearson[repet] <- TP.ggmModSel.pearson[repet]/sum(Theta.estim.ggmModSel.pearson[upper.tri(Theta.estim.ggmModSel.pearson,diag=FALSE)]!=0)
  precision.ggmSS[repet] <- TP.ggmSS[repet]/sum(Theta.estim.ggmSS[upper.tri(Theta.estim.ggmSS,diag=FALSE)]!=0)
  precision.GGMnonreg.neighsel[repet] <- TP.GGMnonreg.neighsel[repet]/sum(Theta.estim.GGMnonreg.neighsel[upper.tri(Theta.estim.GGMnonreg.neighsel,diag=FALSE)]!=0)
  precision.GGMnonreg.boot.pearson[repet] <- TP.GGMnonreg.boot.pearson[repet]/sum(Theta.estim.GGMnonreg.boot.pearson[upper.tri(Theta.estim.GGMnonreg.boot.pearson,diag=FALSE)]!=0)
  precision.GGMnonreg.boot.poly[repet] <- TP.GGMnonreg.boot.poly[repet]/sum(Theta.estim.GGMnonreg.boot.poly[upper.tri(Theta.estim.GGMnonreg.boot.poly,diag=FALSE)]!=0)
  precision.BGGM.explore[repet] <- TP.BGGM.explore[repet]/sum(Theta.estim.BGGM.explore[upper.tri(Theta.estim.BGGM.explore,diag=FALSE)]!=0)
  precision.BGGM.estimate[repet] <- TP.BGGM.estimate[repet]/sum(Theta.estim.BGGM.estimate[upper.tri(Theta.estim.BGGM.estimate,diag=FALSE)]!=0)
  precision.passo[repet] <- TP.passo[repet]/sum(Theta.estim.passo[upper.tri(Theta.estim.passo,diag=FALSE)]!=0)
  
  ########
  ########
  
#save(list=ls(),file=paste0("Simu_repet",repet,".Rdata"))
}

filename<- paste0("Simu_Scenario1_n=",n,"_p=",p,"_L=",L,"_prob=",prob,".Rdata")
#filename
save(list=ls(),file=filename)


#######################
#######################
# Simulation protocol:  
#######################
#######################

# General parameters to vary: 
# n in 100, 250, 500, 1000, 2500
# p in 6:10 ; 20 
# L in 3:5

# Parameters to vary - First scenario
# prob in 1; 0.5 ; 0.4 ; 0.3 ; 0.2; 0.1

######

# Plots to compare the methods: 
# - compare all the methods in one plot (for each performance measure and each simulation setting)
# - compare groups of methods : 
#   (poly.mle ; poly.wls ; pearson)  
#   (Glasso.poly ; Glasso.pearson) 
#   (ggmModSel.poly ; ggmModSel.pearson)
#   (Glasso.poly ; Glasso.pearson ; ggmModSel.poly ; ggmModSel.pearson ; ggmSS)
#   (GGMnonreg.neighsel ; GGMnonreg.boot.pearson ; GGMnonreg.boot.poly)
#   (BGGM.explore, BGGM.fit)
# Depending on these results, we might discard some methods and compare other groups 

#######
# Other things to look at :
# - look at the skewness of the distributions for each variable 

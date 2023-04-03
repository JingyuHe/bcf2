library(dbarts)
library(XBART) ## !! this package needs to be installed from https://github.com/JingyuHe/XBART
library(bcf2)
library(rnn)
library(nnet)
library(foreach)
library(doParallel)

# main code
options(width = 10000)
n_cores <- 8
start.at <- 1
mc <- 100 # number of mc iterations

dgp <- 3 # 3 = hetero, nonlinear, large; 8 = hetero, linear, small

for(dgp.number in dgp)
{
  
  dgp <- as.numeric(int2bin(dgp.number-1,length=3))
  
  small.sample <- dgp[1]
  hetero <- dgp[2]
  linear <- dgp[3]
  
  # initialize storages
  if(start.at == 1){
    results.warm <- matrix(0,mc,7)
    results.xbcf <- matrix(0,mc,7)
  } 
  
  colnames(results.warm)<- c('rmse.ate','rmse.cate','cover.ate','coverage.cate', 'IL.ate','AIL.cate','rt.call')  
  colnames(results.xbcf)<- c('rmse.ate','rmse.cate','cover.ate','coverage.cate', 'IL.ate','AIL.cate','rt.call')
  
  for(iter in start.at:mc)
  {
    # generate the dataframe
    n <- 500 # sample size
    
    x1 <- rnorm(n)
    x2 <- sample(1:2,n,replace=TRUE)
    x3 <- sample(1:3,n,replace=TRUE,prob = c(0.3,0.4,0.3))
    x4 <- rnorm(n)
    x5 <- rnorm(n)
    
    df <- data.frame(x1,x4,x5,x2,as.factor(x3))
    
    # small/large sample
    if (small.sample == 1){ n <- 250
    } else { n <- 500 }
    x <- df[1:n,]
    
    # heterogeneous/homogeneous effects
    if (hetero == 1){alpha <- 1 + 2*x[,2]*x[,4]}else{alpha <- rep(3,n)}
    ate.true <- mean(alpha)
    
    # define mu, tau, pi for each dgp
    # linear/non-linear prognostic function
    mu <- function(x){
      lev <- c(2,-1,-4)
      if (linear == 1){
        result <- 1 + lev[x[,5]] + x[,1]*x[,3]
      } else {
          result <- -6 + lev[x[,5]] + 6*abs(x[,3] - 1)
      }
      return(result)
    }
    
    # prognostic scores
    pi <- 0.8*pnorm(3*mu(x)/sd(mu(x))-0.5*x[,1],0,1) + 0.05 + 0.1*runif(n)
    
    # generate the treatment variable z
    z <- rbinom(n,1,pi)

    # generate the outcome variable y
    mu_true <- mu(x)
    Ey <- mu(x) + alpha * z
    sig <- 0.5 * sd(Ey)
    y <- Ey + sig * rnorm(n)
    
    x.mod <- makeModelMatrixFromDataFrame(x)
    for(i in c(1:7)) {x.mod[,i] <- as.numeric(x.mod[,i])}
    
    # compute pihat here 
    # (for our simulation study we assume that we don't know the prognostic scores, so we compute them)
    compute_pi = TRUE
    if(compute_pi){
      { sink("/dev/null"); fitz <- nnet(z~.,data = x.mod, size = 3,rang = 0.1, maxit = 1000, abstol = 1.0e-8, decay = 5e-2); sink();}
      pihat <- fitz$fitted.values
    } else { pihat = pi}
    
    
    ### method 1: XBCF ###
    # XBCF parameters
    burnin <- 20
    sweeps <- 60
    treesmu <- 30
    treestau <- 10
    p_cat <- ncol(x)-3
    
    tau1 <- 0.6*var(y)/treesmu
    tau2 <- 0.1*var(y)/treestau
    
    t.x = proc.time()
    xbcf_fit = XBART::XBCF.discrete(y, z, x.mod, x.mod, pihat = pihat, num_sweeps = sweeps, burnin = burnin,
                                    max_depth = 250, Nmin = 1, num_cutpoints = 50, no_split_penality = "Auto",
                                    mtry_con = ncol(x.mod)+1, mtry_mod = ncol(x.mod),
                                    pcat_con = pcat,  pcat_mod = pcat,
                                    n_trees_con = treesmu, alpha_con = 0.95, beta_con = 1.25, tau_con = tau1, kap_con = 1, s_con = 1,
                                    n_trees_mod = treestau, alpha_mod = 0.25, beta_mod = 3, tau_mod = tau2, kap_mod = 1, s_mod = 1,
                                    update_tau = TRUE, verbose = FALSE, a_scaling = TRUE, b_scaling = TRUE)
    t.x <- proc.time() - t.x
    pred <- predict.XBCFdiscrete(xbcf_fit, X_con = x.mod, X_mod = x.mod, Z = z, pihat = pihat, burnin = burnin)
    
    # compute rmse, coverage and average interval length for CATE and ATE
    th_xbcf = pred$tau.adj[,(burnin+1):sweeps]
    tauhats_xbcf <- rowMeans(th_xbcf)

    rmse.cate <- sqrt(mean((alpha - tauhats_xbcf)^2))
    rmse.ate <- sqrt(mean((mean(alpha)-mean(tauhats_xbcf))^2))

    lb <- as.numeric(quantile(colMeans(th_xbcf),0.025))
    ub <- as.numeric(quantile(colMeans(th_xbcf),0.975))

    cover.ate <- as.numeric(ub > ate.true & lb < ate.true)

    lbs <- as.numeric(apply(th_xbcf,1,quantile,0.025))
    ubs <- as.numeric(apply(th_xbcf,1,quantile,0.975))

    coverage.cate = mean(ubs > alpha & lbs < alpha)

    IL.ate <- ub-lb
    AIL.cate <- mean(ubs-lbs)

    #total = proc.time() - start
    results.xbcf[iter,] <- c(rmse.ate,rmse.cate,cover.ate,coverage.cate, IL.ate,AIL.cate, as.list(t.x)$elapsed)

    ### END OF XBCF ###

    ### method 2: WARM-START BCF###
    
    ## FIT wsBCF
    
    ################  warmstart parameters

    n_draw_warmstart <- 100
    burnin_warmstart <- 10
    mod_tree_scaling <- 1

    t.ws <- proc.time() # begin time tracking for ws-BCF

    fit_warmstart <- bcf2::bcf.warmstart(y = y, z = z, x_control = x.mod, x_moderate = x.mod, pihat = pihat,
                                         nburn=burnin_warmstart, nsim=n_draw_warmstart, include_pi = 'control',
                                         use_tauscale = TRUE, ntree_control = treesmu, ntree_moderate = treestau, ini_bcf = FALSE,
                                         verbose = FALSE, n_cores = n_cores, warm_start_fit = xbcf_fit
    )

    t.ws <- proc.time() - t.ws # end time tracking for warmstart

    # compute rmse, coverage and average interval length for CATE and ATE
    taus <- rowMeans(fit_warmstart$tauhat)

    rmse.cate <- sqrt(mean((alpha - taus)^2))
    rmse.ate <- sqrt(mean((mean(alpha)-mean(taus))^2))

    lb <- as.numeric(quantile(colMeans(fit_warmstart$tauhat),0.025))
    ub <- as.numeric(quantile(colMeans(fit_warmstart$tauhat),0.975))

    cover.ate <- as.numeric(ub > ate.true & lb < ate.true)

    lbs <- as.numeric(apply(fit_warmstart$tauhat,1,quantile,0.025))
    ubs <- as.numeric(apply(fit_warmstart$tauhat,1,quantile,0.975))

    coverage.cate <- mean(ubs > alpha & lbs < alpha)

    IL.ate <- ub-lb
    AIL.cate <- mean(ubs-lbs)

    # store results
    results.warm[iter,] <- c(rmse.ate,rmse.cate,cover.ate,coverage.cate, IL.ate,AIL.cate,as.list(t.ws)$elapsed)

    ### END WARM-START BCF ###
    

    # save to file
    # save(file = paste('~/test_sim_study_dgp_', dgp.number, sep=''),
    #      results.xbcf, results.warm, results.bcf, results.crf,
    #      results.bart, results.bartf0f1, results.psbart, results.lm,
    #      n, small.sample, hetero, linear)
  
  }
  
}

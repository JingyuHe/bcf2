bcf_ini <- function(
                treedraws_con, treedraws_mod, muscale_ini, bscale0_ini, bscale1_ini, sigma_ini, pi_con_tau, pi_con_sigma,
                pi_mod_tau, pi_mod_sigma, mod_tree_scaling,
                y, z, x_control, x_moderate=x_control, pihat, w = NULL, 
                random_seed = sample.int(.Machine$integer.max, 1),
                n_chains = 4,
                n_cores  = n_chains,
                n_threads = max((RcppParallel::defaultNumThreads()-2)/n_cores,1), #max number of threads, minus a arbitrary holdback, over the number of cores
                nburn, nsim, nthin = 1, update_interval = 100,
                ntree_control = 250,
                sd_control = 2*sd(y),
                base_control = 0.95,
                power_control = 2,
                ntree_moderate = 50,
                sd_moderate = sd(y),
                base_moderate = 0.25,
                power_moderate = 3,
                save_tree_directory = '..',
                nu = 3, lambda = NULL, sigq = .9, sighat = NULL,
                include_pi = "control", use_muscale=TRUE, use_tauscale=TRUE, verbose=FALSE,
                ini_bcf = FALSE, update_mu_loading_tree = FALSE,
                x_c = NULL, x_m = NULL, cutpoint_list_c = NULL, cutpoint_list_m = NULL
) {

  
  if(is.null(w)){
    w <- matrix(1, ncol = 1, nrow = length(y))
    }

  pihat = as.matrix(pihat)
  # if(!.ident(length(y),
  #            length(z),
  #            length(w),
  #            nrow(x_control),
  #            nrow(x_moderate),
  #            nrow(pihat))
  #   ) {
  #   stop("Data size mismatch. The following should all be equal:
  #        length(y): ", length(y), "\n",
  #        "length(z): ", length(z), "\n",
  #        "length(w): ", length(w), "\n",
  #        "nrow(x_control): ", nrow(x_control), "\n",
  #        "nrow(x_moderate): ", nrow(x_moderate), "\n",
  #        "nrow(pihat): ", nrow(pihat),"\n"
  #   )
  # }

  # if(any(is.na(y))) stop("Missing values in y")
  # if(any(is.na(z))) stop("Missing values in z")
  # if(any(is.na(w))) stop("Missing values in w")
  # if(any(is.na(x_control))) stop("Missing values in x_control")
  # if(any(is.na(x_moderate))) stop("Missing values in x_moderate")
  # if(any(is.na(pihat))) stop("Missing values in pihat")
  # if(any(!is.finite(y))) stop("Non-numeric values in y")
  # if(any(!is.finite(z))) stop("Non-numeric values in z")
  # if(any(!is.finite(w))) stop("Non-numeric values in w")
  # if(any(!is.finite(x_control))) stop("Non-numeric values in x_control")
  # if(any(!is.finite(x_moderate))) stop("Non-numeric values in x_moderate")
  # if(any(!is.finite(pihat))) stop("Non-numeric values in pihat")
  # if(!all(sort(unique(z)) == c(0,1))) stop("z must be a vector of 0's and 1's, with at least one of each")

  # if(length(unique(y))<5) warning("y appears to be discrete")

  # if(nburn<0) stop("nburn must be positive")
  # if(nsim<0) stop("nsim must be positive")
  # if(nthin<0) stop("nthin must be positive")
  # if(nthin>nsim+1) stop("nthin must be < nsim")
  # if(nburn<1000) warning("A low (<1000) value for nburn was supplied")
  # if(nsim*nburn<1000) warning("A low (<1000) value for total iterations after burn-in was supplied")

  ### TODO range check on parameters

  ###
  if(is.null(x_c)){
    return_sort_index = TRUE

    x_c = matrix(x_control, ncol=ncol(x_control))
    if(include_pi=="both" | include_pi=="control") {
      x_c = cbind(pihat, x_control)
    }
  }else{
    return_sort_index = FALSE
  }

  if(is.null(x_m)){
    x_m = matrix(x_moderate, ncol=ncol(x_moderate))
    if(include_pi=="both" | include_pi=="moderate") {
      x_m = cbind(pihat, x_moderate)
    }
  }

  if(is.null(cutpoint_list_c)){
    cutpoint_list_c = lapply(1:ncol(x_c), function(i) .cp_quantile(x_c[,i]))
  }

  if(is.null(cutpoint_list_m)){
    cutpoint_list_m = lapply(1:ncol(x_m), function(i) .cp_quantile(x_m[,i]))
  }

  sdy = sqrt(Hmisc::wtd.var(y, w))
  muy = stats::weighted.mean(y, w)
  yscale = (y-muy)/sdy


  if(is.null(lambda)) {
    if(is.null(sighat)) {
      lmf = lm(yscale~z+as.matrix(x_c), weights = w)
      sighat = summary(lmf)$sigma #sd(y) #summary(lmf)$sigma
    }
    qchi = qchisq(1.0-sigq,nu)
    lambda = (sighat*sighat*qchi)/nu
  }

  dir = tempdir()

  perm = order(z, decreasing=TRUE)

  con_sd = ifelse(abs(2*sdy - sd_control)<1e-6, 2, sd_control/sdy)
  mod_sd = ifelse(abs(sdy - sd_moderate)<1e-6, 1, sd_moderate/sdy)/ifelse(use_tauscale,0.674,1) # if HN make sd_moderate the prior median

  RcppParallel::setThreadOptions(numThreads=n_threads)
  
  do_type_config <- .get_do_type(n_cores)
  `%doType%` <- do_type_config$doType
  
  chain_out <- foreach::foreach(iChain=1:n_chains) %doType% {
    
    this_seed = random_seed + iChain - 1
    
    cat("Calling bcfoverparRcppClean From R\n")
    set.seed(this_seed)
    
    tree_files = .get_chain_tree_files(save_tree_directory, iChain)
    
    print(tree_files)

    fitbcf = bcfoverparRcppClean_ini(
      ini_bcf, treedraws_con, treedraws_mod, muscale_ini, bscale0_ini, bscale1_ini, sigma_ini, pi_con_tau, pi_con_sigma,
                        pi_mod_tau, pi_mod_sigma, mod_tree_scaling,
                                 y_ = yscale[perm], z_ = z[perm], w_ = w[perm],
                                 x_con_ = t(x_c[perm,,drop=FALSE]), x_mod_ = t(x_m[perm,,drop=FALSE]), 
                                 x_con_info_list = cutpoint_list_c, 
                                 x_mod_info_list = cutpoint_list_m,
                                 random_des = matrix(1),
                                 random_var = matrix(1),
                                 random_var_ix = matrix(1),
                                 random_var_df = 3,
                                 burn = nburn, nd = nsim, thin = nthin,
                                 ntree_mod = ntree_moderate, ntree_con = ntree_control, 
                                 lambda = lambda, nu = nu,
                                 con_sd = con_sd,
                                 mod_sd = mod_sd, # if HN make sd_moderate the prior median
                                 mod_alpha = base_moderate, 
                                 mod_beta = power_moderate, 
                                 con_alpha = base_control, 
                                 con_beta = power_control,
                                 treef_con_name_ = tree_files$con_trees, 
                                 treef_mod_name_ = tree_files$mod_trees, 
                                 status_interval = update_interval,
                                 use_mscale = use_muscale, use_bscale = use_tauscale, 
                                 b_half_normal = TRUE, verbose_sigma=verbose, update_mu_loading_tree = update_mu_loading_tree, trt_init = 1.0)
    
    cat("bcfoverparRcppClean returned to R\n")

    ac = fitbcf$m_post[,order(perm)]

    Tm = fitbcf$b_post[,order(perm)] * (1.0/ (fitbcf$b1 - fitbcf$b0))

    Tc = ac * (1.0/fitbcf$msd) 

    tau_post = sdy*fitbcf$b_post[,order(perm)]

    mu_post  = muy + sdy*(Tc*fitbcf$msd + Tm*fitbcf$b0)
    
    list(sigma = sdy*fitbcf$sigma,
         yhat = muy + sdy*fitbcf$yhat_post[,order(perm)],
         sdy = sdy,
         con_sd = con_sd,
         mod_sd = mod_sd,
         muy = muy,
         mu  = mu_post,
         tau = tau_post,
         mu_scale = fitbcf$msd,
         tau_scale = fitbcf$bsd,
         b0 = fitbcf$b0,
         b1 = fitbcf$b1,
         perm = perm,
         include_pi = include_pi,
         random_seed=this_seed
    )
    
  }


  all_sigma = c()
  all_mu_scale = c()
  all_tau_scale = c()

  all_b0 = c()
  all_b1 = c()
  
  all_yhat = c()
  all_mu   = c()
  all_tau  = c()
  
  chain_list=list()

  n_iter = length(chain_out[[1]]$sigma)
  
  
  for (iChain in 1:n_chains){
    sigma        <- chain_out[[iChain]]$sigma
    mu_scale     <- chain_out[[iChain]]$mu_scale
    tau_scale    <- chain_out[[iChain]]$tau_scale
    
    b0          <- chain_out[[iChain]]$b0
    b1          <- chain_out[[iChain]]$b1

    yhat         <- chain_out[[iChain]]$yhat
    tau          <- chain_out[[iChain]]$tau
    mu           <- chain_out[[iChain]]$mu

    # -----------------------------    
    # Support Old Output
    # -----------------------------
    all_sigma       = c(all_sigma,     sigma)
    all_mu_scale    = c(all_mu_scale,  mu_scale)
    all_tau_scale   = c(all_tau_scale, tau_scale)
    all_b0 = c(all_b0, b0)
    all_b1 = c(all_b1, b1)

    all_yhat = rbind(all_yhat, yhat)
    all_mu   = rbind(all_mu,   mu)
    all_tau  = rbind(all_tau,  tau)

    # -----------------------------    
    # Make the MCMC Object
    # -----------------------------

    scalar_df <- data.frame("sigma"     = sigma,
                            "tau_bar"   = matrixStats::rowWeightedMeans(tau, w),
                            "mu_bar"    = matrixStats::rowWeightedMeans(mu, w),
                            "yhat_bar"  = matrixStats::rowWeightedMeans(yhat, w),
                            "mu_scale"  = mu_scale, 
                            # "tau_scale" = tau_scale,
                            "b0"  = b0, 
                            "b1"  = b1)
    
    # y_df <- as.data.frame(chain$yhat)
    # colnames(y_df) <- paste0('y',1:ncol(y_df))
    # 
    # mu_df <- as.data.frame(chain$mu)
    # colnames(mu_df) <- paste0('mu',1:ncol(mu_df))
    # 
    # tau_df <- as.data.frame(chain$tau)
    # colnames(tau_df) <- paste0('tau',1:ncol(tau_df))
    
    chain_list[[iChain]] <- coda::as.mcmc(scalar_df)
    # -----------------------------    
    # Sanity Check Constants Accross Chains
    # -----------------------------
    if(chain_out[[iChain]]$sdy        != chain_out[[1]]$sdy)        stop("sdy not consistent between chains for no reason")
    if(chain_out[[iChain]]$con_sd     != chain_out[[1]]$con_sd)     stop("con_sd not consistent between chains for no reason")
    if(chain_out[[iChain]]$mod_sd     != chain_out[[1]]$mod_sd)     stop("mod_sd not consistent between chains for no reason")
    if(chain_out[[iChain]]$muy        != chain_out[[1]]$muy)        stop("muy not consistent between chains for no reason")
    if(chain_out[[iChain]]$include_pi != chain_out[[1]]$include_pi) stop("include_pi not consistent between chains for no reason")
    if(any(chain_out[[iChain]]$perm   != chain_out[[1]]$perm))      stop("perm not consistent between chains for no reason")
  }

  fitObj <- list(sigma = all_sigma,
                 yhat = all_yhat,
                 sdy = chain_out[[1]]$sdy,
                 muy = chain_out[[1]]$muy,
                 mu  = all_mu,
                 tau = all_tau,
                 mu_scale = all_mu_scale,
                 tau_scale = all_tau_scale,
                 b0 = all_b0,
                 b1 = all_b1,
                 perm = perm,
                 include_pi = chain_out[[1]]$include_pi,
                 random_seed = chain_out[[1]]$random_seed,
                 coda_chains = coda::as.mcmc.list(chain_list),
                 raw_chains = chain_out)
  
  attr(fitObj, "class") <- "bcf"
  
  .cleanup_after_par(do_type_config)
  
  return(fitObj)
}

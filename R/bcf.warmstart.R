# .ident <- function(...){
# # courtesy https://stackoverflow.com/questions/19966515/how-do-i-test-if-three-variables-are-equal-r
#   args <- c(...)
#   if( length( args ) > 2L ){
#     #  recursively call ident()
#     out <- c( identical( args[1] , args[2] ) , .ident(args[-1]))
#   }else{
#     out <- identical( args[1] , args[2] )
#   }
#   return( all( out ) )
# }
#

.comb <- function(x, ...) {
  lapply(seq_along(x),
         function(i) c(x[[i]], lapply(list(...), function(y) y[[i]])))
}

.cp_list= function(x){
  if(length(unique(x)) == 1) {
    ret = x[1]
    warning("A supplied covariate contains a single distinct value.")
  } else if(identical(sort(unique(x)),c(0,1))) {
    ret <- 0.5
  } else {
    ret <- sort(unique(x))[-1]
    }
  return(ret)
}

#' Fit warm-start Bayesian Causal Forests
#'
#' @references Hahn, Murray, and Carvalho (2020). Bayesian Regression Tree Models for Causal Inference: Regularization, Confounding, and Heterogeneous Effects (with Discussion).
#' Bayesian Anal. 15(3): 965-1056 (September 2020). DOI: 10.1214/19-BA1195
#'
#' @details Fits multiple Bayesian Causal Forest models (Hahn et. al., 2020),
#' that are independently initialized with forests extracted from a fitted Accelerated Bayesian Causal Forests heuristic (Krantsevich et al., 2023).
#' For a response variable y, binary treatment z, and covariates x,
#' \deqn{y_i = a \mu(x_i, \pi_i) + b_{z_i}\tau(x_i, \pi_i) + \epsilon_i}
#' where \eqn{\pi_i} is an estimate of the propensity score \eqn{\Pr(Z_i=1 | X_i=x_i)};
#' \eqn{a, b_0, b_1} are scaling parameters; and \eqn{\epsilon_i \sim N(0,\sigma^2)}.
#'
#' Some notes:
#' \itemize{
#'    \item x_control and x_moderate must be numeric matrices. See e.g. the makeModelMatrix function in the
#'    dbarts package for appropriately constructing a design matrix from a data.frame.
#'    For proper initialization, these matrices must be identical to the ones used when fitting XBCF model.
#'    \item sd_control and sd_moderate are the prior SD(mu(x)) and SD(tau(x)) at a given value of x (respectively). If
#'    use_muscale = FALSE, then this is the parameter \eqn{\sigma_\mu} from the original BART paper, where the leaf parameters
#'    have prior distribution \eqn{N(0, \sigma_\mu/m)}, where m is the number of trees.
#'    If use_muscale=TRUE then sd_control is the prior median of a half Cauchy prior for SD(mu(x)). If use_tauscale = TRUE,
#'    then sd_moderate is the prior median of a half Normal prior for SD(tau(x)).
#'    \item By default the prior on \eqn{\sigma^2} is calibrated as in Chipman, George and McCulloch (2008).
#'
#' }
#' @param y Response variable
#' @param z Treatment variable
#' @param x_control Design matrix for the "prognostic" function mu(x)
#' @param x_moderate Design matrix for the covariate-dependent treatment effects tau(x)
#' @param pihat Length n estimates of
#' @param w An optional vector of weights. When present, BCF fits a model \eqn{y | x ~ N(f(x), \sigma^2 / w)}, where \eqn{f(x)} is the unknown function.
#' @param nburn Number of burn-in MCMC iterations
#' @param nsim Number of MCMC iterations to save after burn-in
#' @param nthin Save every nthin'th MCMC iterate. The total number of MCMC iterations will be nsim*nthin + nburn.
#' @param update_interval Print status every update_interval MCMC iterations
#' @param ntree_control Number of trees in mu(x)
#' @param sd_control SD(mu(x)) marginally at any covariate value (or its prior median if use_muscale=TRUE)
#' @param base_control Base for tree prior on mu(x) trees (see details)
#' @param power_control Power for the tree prior on mu(x) trees
#' @param ntree_moderate Number of trees in tau(x)
#' @param sd_moderate SD(tau(x)) marginally at any covariate value (or its prior median if use_tauscale=TRUE)
#' @param base_moderate Base for tree prior on tau(x) trees (see details)
#' @param power_moderate Power for the tree prior on tau(x) trees (see details)
#' @param nu Degrees of freedom in the chisq prior on \eqn{sigma^2}
#' @param lambda Scale parameter in the chisq prior on \eqn{sigma^2}
#' @param sigq Calibration quantile for the chisq prior on \eqn{sigma^2}
#' @param sighat Calibration estimate for the chisq prior on \eqn{sigma^2}
#' @param include_pi Takes values "control", "moderate", "both" or "none". Whether to
#' include pihat in mu(x) ("control"), tau(x) ("moderate"), both or none. Values of "control"
#' or "both" are HIGHLY recommended with observational data.
#' @param use_muscale Use a half-Cauchy hyperprior on the scale of mu.
#' @param use_tauscale Use a half-Normal prior on the scale of tau.
#' @param n_cores An optional integer of the number of cores to run your warm start-BCF.
#' @param warm_start_fit An XBCF fit object, that will be used to initialize warm start-BCF.
#' @param cutpoint_grid The scheme for cutpoint grid generation:
#' "unifrom" [default] option generates a unifromly distributed cutpoint grid for each continuous variable (size is 10*number of observations);
#' "exact" option uses sets of unique values at each categorical or continuous variable as availible cutpoints.
#' @return A list with elements
#' \item{tau}{\code{nsim} by \code{n} matrix of posterior samples of individual treatment effects}
#' \item{yhat}{\code{nsim} by \code{n} matrix of posterior samples of individual outcome}
#'
#' @useDynLib bcf2
#' @import Rcpp
#' @importFrom stats approxfun lm qchisq quantile sd
#' @export
bcf.warmstart <- function(
                y, z, x_control, x_moderate=x_control, pihat, w = NULL,
                nburn, nsim, nthin = 1, update_interval = 100,
                ntree_control = 200,
                sd_control = 2*sd(y),
                base_control = 0.95,
                power_control = 2,
                ntree_moderate = 50,
                sd_moderate = sd(y),
                base_moderate = 0.25,
                power_moderate = 3,
                nu = 3, lambda = NULL, sigq = .9, sighat = NULL, randeff = FALSE,
                include_pi = "control", use_muscale=TRUE, use_tauscale=TRUE, ini_bcf = FALSE, update_mu_loading_tree = FALSE,
                verbose = FALSE,n_cores = NULL, warm_start_fit = NULL, cutpoint_grid = 'uniform'
) {
  if (is.null(warm_start_fit)) {
    stop("bcf.warmstart requires an XBCF fit object for initialization. Stopping.")
  }
  if (class(warm_start_fit) != "XBCFdiscrete") {
    stop("bcf.warmstart can only be initialized with an object of class XBCF. Stopping.")
  }
  if (!(cutpoint_grid == 'uniform' || cutpoint_grid == 'exact')) {
    stop("Invalid input: cutpoint_grid should be either 'unifrom' or 'exact'. Stopping.")
  }
  if (is.null(n_cores)) {
    n_cores <- 2
    warning("n_cores is not provided. Using 2 as default for parallelization.")
  }

  if(is.null(w)){
    w <- matrix(1, ncol = 1, nrow = length(y))
    }

  pihat = as.matrix(pihat)


  ### TODO range check on parameters

  ###
  n_sweeps = warm_start_fit$model_params$n_sweeps
  n_burnin = warm_start_fit$model_params$n_burnin
  ntree_control = warm_start_fit$model_params$n_trees_con
  ntree_moderate = warm_start_fit$model_params$n_trees_mod
  include_pi = warm_start_fit$model_params$pihat_status

  x_c = matrix(x_control, ncol=ncol(x_control))
  if(include_pi=="both" | include_pi=="control") {
    x_c = cbind(pihat, x_control)
  }

  x_m = matrix(x_moderate, ncol=ncol(x_moderate))
  if(include_pi=="both" | include_pi=="moderate") {
    x_m = cbind(pihat, x_moderate)
  }

  if(cutpoint_grid == "exact") {
    cutpoint_list_c = lapply(1:ncol(x_c), function(i) .cp_list(x_c[,i]))
    cutpoint_list_m = lapply(1:ncol(x_m), function(i) .cp_list(x_m[,i]))
  } else {
    cutpoint_list_c = lapply(1:ncol(x_c), function(i) .cp_quantile(x_c[,i], num = 10*nrow(x_c)))
    cutpoint_list_m = lapply(1:ncol(x_m), function(i) .cp_quantile(x_m[,i], num = 10*nrow(x_m)))
  }

  sdy = sqrt(Hmisc::wtd.var(y, w))
  muy = stats::weighted.mean(y, w)
  yscale = (y-muy)/sdy

  if(is.null(sighat)) {
    lmf = lm(yscale~z+as.matrix(x_c))
    sighat = summary(lmf)$sigma #sd(y) #summary(lmf)$sigma
  }
  qchi = qchisq(1.0-sigq,nu)
  lambda = (sighat*sighat*qchi)/nu

  dir = tempdir()

  perm = order(z, decreasing=TRUE)

  ## parameter matching


  ## set up parallel
  cl <- parallel::makeCluster(n_cores)
  doParallel::registerDoParallel(cl)

  if(verbose){
    cat(paste0("Calling bcfoverparRcppClean_ini from R in parallel on ", n_cores, " cores.\n"))
  }
  n_chains <- n_sweeps - n_burnin

  out <- foreach::foreach(chain_num=1:n_chains,
                                .combine='.comb', .multicombine=TRUE, .init=list(list(), list())) %dopar% {

    i <- n_burnin + chain_num

    bscale0_ini = warm_start_fit$b[i, 1]
    bscale1_ini = warm_start_fit$b[i, 2]
    #sigma_ini = warm_start_fit$sigma0_draws[(ntree_control+ntree_moderate),i]

    treedraws_con = as.vector(warm_start_fit$tree_string_con[i])
    treedraws_mod = as.vector(warm_start_fit$tree_string_mod[i])
    muscale_ini = warm_start_fit$a[i, 1]

    sigma_ini = warm_start_fit$sigma0[ntree_control,i]
    pi_con_sigma_ini = warm_start_fit$sigma0[ntree_control,i] / warm_start_fit$a[i, 1]
    pi_mod_sigma_ini = warm_start_fit$sigma0[ntree_control,i]
    pi_con_tau = sqrt(warm_start_fit$tau_con[i, 1])
    pi_mod_tau = sqrt(warm_start_fit$tau_mod[i, 1])
    mod_tree_scaling = 1
    ini_bcf = FALSE

    fitbcf = bcfoverparRcppClean_ini( ini_bcf, treedraws_con, treedraws_mod, muscale_ini, bscale0_ini, bscale1_ini,
                                      sigma_ini, pi_con_tau, pi_con_sigma_ini, pi_mod_tau, pi_mod_sigma_ini, mod_tree_scaling,
                                      yscale[perm], z[perm], w[perm],
                                      t(x_c[perm,]), t(x_m[perm,,drop=FALSE]), t(x_m[1,,drop=FALSE]),
                                      cutpoint_list_c, cutpoint_list_m,
                                      random_des = matrix(1),
                                      random_var = matrix(1),
                                      random_var_ix = matrix(1),
                                      random_var_df = 3,
                                      nburn, nsim, nthin,
                                      ntree_moderate, ntree_control, lambda, nu,
                                      con_sd = ifelse(abs(2*sdy - sd_control)<1e-6, 2, sd_control/sdy),
                                      mod_sd = ifelse(abs(sdy - sd_moderate)<1e-6, 1, sd_moderate/sdy)/ifelse(use_tauscale,0.674,1), # if HN make sd_moderate the prior median
                                      base_control, power_control, base_moderate, power_moderate,
                                      "tmp", status_interval = update_interval, randeff = randeff,
                                      use_mscale = use_muscale, use_bscale = use_tauscale, b_half_normal = TRUE, update_mu_loading_tree = update_mu_loading_tree, trt_init = 1.0, verbose = verbose)


    return(list(t(sdy*fitbcf$b_post[,order(perm)]),
                t(muy + sdy*fitbcf$yhat_post[,order(perm)])
                     )
          )
                                }

  if(verbose){
    cat("Loop complete, back to R.\n")
  }

  ## parallel clean-up
  parallel::stopCluster(cl)

  return(list(tauhat = matrix(unlist(out[[1]]), nrow = length(y), byrow = FALSE),
              yhat = matrix(unlist(out[[2]]), nrow = length(y), byrow = FALSE)
  ))
}

#' @export
verify_install <- function() {
    cat("BCF2 Installed Correctly\n")
}
#include <RcppArmadillo.h>

#include <iostream>
#include <fstream>
#include <vector>
#include <ctime>

#include "rng.h"
#include "tree.h"
#include "info.h"
#include "funs.h"
#include "bd.h"
#include "logging.h"

using namespace Rcpp;

#pragma GCC diagnostic ignored "-Wunused-parameter"
#pragma GCC diagnostic ignored "-Wcomment"
#pragma GCC diagnostic ignored "-Wformat"
#pragma GCC diagnostic ignored "-Wsign-compare"

// y = m(x) + b(x)z + e, e~N(0, sigma^2_y

//x_con is the design matrix for m. It should have n = rows
//x_mod is the design matrix for b. It should have n = rows
//data should come in sorted with all trt first, then control cases

// [[Rcpp::export]]
List bcfoverparRcppClean_ini(bool ini_bcf, SEXP treedraws_con, SEXP treedraws_mod, double muscale_ini, double bscale0_ini, double bscale1_ini, double sigma_ini, double pi_con_tau, double pi_con_sigma, double pi_mod_tau, double pi_mod_sigma,
                             NumericVector y_, NumericVector z_, NumericVector w_,
                             NumericVector x_con_, NumericVector x_mod_, NumericVector x_mod_est_,
                             List x_con_info_list, List x_mod_info_list,
                             arma::mat random_des,                          //needs to come in with n rows no matter what(?)
                             arma::mat random_var, arma::mat random_var_ix, //random_var_ix*random_var = diag(Var(random effects))
                             double random_var_df,
                             int burn, int nd, int thin, //Draw nd*thin + burn samples, saving nd draws after burn-in
                             int ntree_mod, int ntree_con,
                             double lambda, double nu, //prior pars for sigma^2_y
                             double con_sd,            // Var(m(x)) = con_sd^2 marginally a priori (approx)
                             double mod_sd,            // Var(b(x)) = mod_sd^2 marginally a priori (approx)
                             double con_alpha, double con_beta,
                             double mod_alpha, double mod_beta,
                             CharacterVector treef_name_,
                             int status_interval = 100,
                             bool RJ = false, bool randeff = false, bool use_mscale = true, bool use_bscale = true, bool b_half_normal = true,
                             double trt_init = 1.0)
{

    // bool randeff = true;
    // if (random_var_ix.n_elem == 1)
    // {
    //     randeff = false;
    // }

    if (randeff){
        Rcout << "Using random effects." << std::endl;
    }else{
        Rcout << "NOT using random effects." << std::endl;
    }

    std::string treef_name = as<std::string>(treef_name_);
    std::ofstream treef(treef_name.c_str());

    RNGScope scope;
    RNG gen; //this one random number generator is used in all draws

    //double lambda = 1.0; //this one really needs to be set
    //double nu = 3.0;
    //double kfac=2.0; //original is 2.0

    Logger logger = Logger();
    char logBuff[100];

    bool log_level = false;

    logger.setLevel(log_level);
    logger.log("============================================================");
    logger.log(" Starting up BCF: ");
    logger.log("============================================================");
    if (log_level)
    {
        logger.getVectorHead(y_, logBuff);
        Rcout << "y: " << logBuff << "\n";
        logger.getVectorHead(z_, logBuff);
        Rcout << "z: " << logBuff << "\n";
        logger.getVectorHead(w_, logBuff);
        Rcout << "w: " << logBuff << "\n";
    }

    logger.log("BCF is Weighted");

    // sprintf(logBuff, "Updating Moderate Tree: %d of %d");
    // logger.log(logBuff);
    logger.log("");

    /*****************************************************************************
  /* Read, format y
  *****************************************************************************/
    std::vector<double> y; //storage for y
    double miny = INFINITY, maxy = -INFINITY;
    sinfo allys; //sufficient stats for all of y, use to initialize the bart trees.

    for (NumericVector::iterator it = y_.begin(); it != y_.end(); ++it)
    {
        y.push_back(*it);
        if (*it < miny)
            miny = *it;
        if (*it > maxy)
            maxy = *it;
        allys.sy += *it;            // sum of y
        allys.sy2 += (*it) * (*it); // sum of y^2
    }
    size_t n = y.size();
    allys.n = n;

    double ybar = allys.sy / n;                                  //sample mean
    double shat = sqrt((allys.sy2 - n * ybar * ybar) / (n - 1)); //sample standard deviation
    /*****************************************************************************
  /* Read, format  weights 
  *****************************************************************************/
    double *w = new double[n]; //y-(allfit-ftemp) = y-allfit+ftemp

    for (int j = 0; j < n; j++)
    {
        w[j] = w_[j];
    }

    /*****************************************************************************
  /* Read, format X_con
  *****************************************************************************/
    //the n*p numbers for x are stored as the p for first obs, then p for second, and so on.
    std::vector<double> x_con;
    for (NumericVector::iterator it = x_con_.begin(); it != x_con_.end(); ++it)
    {
        x_con.push_back(*it);
    }
    size_t p_con = x_con.size() / n;

    Rcout << "Using " << p_con << " control variables." << std::endl;

    //x cutpoints
    xinfo xi_con;

    xi_con.resize(p_con);
    for (int i = 0; i < p_con; ++i)
    {
        NumericVector tmp = x_con_info_list[i];
        std::vector<double> tmp2;
        for (size_t j = 0; j < tmp.size(); ++j)
        {
            tmp2.push_back(tmp[j]);
        }
        xi_con[i] = tmp2;
    }

    /*****************************************************************************
  /* Read, format X_mod
  *****************************************************************************/
    int ntrt = 0;
    for (size_t i = 0; i < n; ++i)
    {
        if (z_[i] > 0)
            ntrt += 1;
    }
    std::vector<double> x_mod;
    for (NumericVector::iterator it = x_mod_.begin(); it != x_mod_.end(); ++it)
    {
        x_mod.push_back(*it);
    }
    size_t p_mod = x_mod.size() / n;

    Rcout << "Using " << p_mod << " potential effect moderators." << std::endl;

    //x cutpoints
    xinfo xi_mod;

    // optimize, get rid of .push_back(). should declare vector size in the beginning
    xi_mod.resize(p_mod);
    for (int i = 0; i < p_mod; ++i)
    {
        NumericVector tmp = x_mod_info_list[i];
        std::vector<double> tmp2;
        for (size_t j = 0; j < tmp.size(); ++j)
        {
            tmp2.push_back(tmp[j]);
        }
        xi_mod[i] = tmp2;
    }

    //  Rcout <<"\nburn,nd,number of trees: " << burn << ", " << nd << ", " << m << endl;
    //  Rcout <<"\nlambda,nu,kfac: " << lambda << ", " << nu << ", " << kfac << endl;

    // cout << "before loading trees " << endl;
    // load trees
    Rcpp::CharacterVector itrees_con(Rcpp::wrap(treedraws_con));
    Rcpp::CharacterVector itrees_mod(Rcpp::wrap(treedraws_mod));

    std::string itv_con(itrees_con[0]);
    std::stringstream ttss_con(itv_con);
    size_t mm_con, pp_con;
    ttss_con >> mm_con >> pp_con;
    std::vector<tree> t_con(mm_con);
    std::vector<tree::tree_p> temp_node;
    size_t index_in_full;
    for (size_t j = 0; j < mm_con; j++)
    {
        // load from input trees, cutpoint are raw values
        ttss_con >> t_con[j];

        // search xinfo matrix in bcf package, find corresponding index of cutpoints
        temp_node.clear();
        t_con[j].getnodes(temp_node);
        cout << " print some nodes " << endl;
        cout << "variable index " << temp_node[0]->getv() << " cutpoint index " << temp_node[0]->getc() << " cutpoint value " << temp_node[0]->getc_value() << endl;

        if (!ini_bcf)
        {
            // if not initialize at BCF
            // initialize at XBCF
            // need to figure out indicies of cutpoint in xi matrix

            cout << "call initialization of XBCF " << endl;
            for (size_t kk = 0; kk < temp_node.size(); kk++)
            {
                if (temp_node[kk]->getl())
                {
                    // if left / right child exist, it is not a leaf node

                    if (xi_con[temp_node[kk]->getv()].size() == 1)
                    {
                        // if there is only one cutpoint candidate (binary variable)
                        // you have to split at 0.5
                        // index is 0 and value is 0.5
                        temp_node[kk]->setc(0);
                        temp_node[kk]->setc_value(0.5);
                    }
                    else
                    {
                        index_in_full = 0;
                        while (xi_con[temp_node[kk]->getv()][index_in_full] < temp_node[kk]->getc_value() && index_in_full <= xi_con[temp_node[kk]->getv()].size())
                        {
                            index_in_full++;
                        }

                        if (index_in_full > xi_con[temp_node[kk]->getv()].size())
                        {
                            index_in_full = xi_con[temp_node[kk]->getv()].size() - 1;
                        }
                        temp_node[kk]->setc(index_in_full - 1);
                    }
                }
            }
        }

        cout << " print some nodes, after " << endl;
        cout << "variable index " << temp_node[0]->getv() << " cutpoint index " << temp_node[0]->getc() << " cutpoint value " << temp_node[0]->getc_value() << endl;
        cout << " last cutpoint " << xi_con[temp_node[0]->getv()][temp_node[0]->getc() - 1] << " current cutpoint " << xi_con[temp_node[0]->getv()][temp_node[0]->getc()] << " next cutpoint " << xi_con[temp_node[0]->getv()][temp_node[0]->getc() + 1] << endl;
    }

    // cout << "load mod trees " << endl;

    std::string itv_mod(itrees_mod[0]);
    std::stringstream ttss_mod(itv_mod);
    size_t mm_mod, pp_mod;
    ttss_mod >> mm_mod >> pp_mod;
    std::vector<tree> t_mod(mm_mod);
    for (size_t j = 0; j < mm_mod; j++)
    {
        ttss_mod >> t_mod[j];

        // search xinfo matrix in bcf package, find corresponding index of cutpoints
        temp_node.clear();
        t_mod[j].getnodes(temp_node);
        cout << " print some nodes " << endl;
        cout << "variable index " << temp_node[0]->getv() << " cutpoint index " << temp_node[0]->getc() << " cutpoint value " << temp_node[0]->getc_value() << endl;

        if (!ini_bcf)
        {
            cout << "call initialization of XBCF " << endl;

            for (size_t kk = 0; kk < temp_node.size(); kk++)
            {
                if (temp_node[kk]->getl())
                {
                    // if left / right child exist, it is not a leaf node

                    if (xi_mod[temp_node[kk]->getv()].size() == 1)
                    {
                        // if there is only one cutpoint candidate (binary variable)
                        // you have to split at 0.5
                        // index is 0 and value is 0.5
                        temp_node[kk]->setc(0);
                        temp_node[kk]->setc_value(0.5);
                    }
                    else
                    {
                        index_in_full = 0;
                        while (xi_mod[temp_node[kk]->getv()][index_in_full] < temp_node[kk]->getc_value() && index_in_full <= xi_mod[temp_node[kk]->getv()].size())
                        {
                            index_in_full++;
                        }

                        if (index_in_full > xi_mod[temp_node[kk]->getv()].size())
                        {
                            index_in_full = xi_mod[temp_node[kk]->getv()].size() - 1;
                        }
                        temp_node[kk]->setc(index_in_full - 1);
                    }
                }
            }
        }
        cout << " print some nodes, after " << endl;
        cout << "variable index " << temp_node[0]->getv() << " cutpoint index " << temp_node[0]->getc() << " cutpoint value " << temp_node[0]->getc_value() << endl;
        cout << " last cutpoint " << xi_con[temp_node[0]->getv()][temp_node[0]->getc() - 1] << " current cutpoint " << xi_con[temp_node[0]->getv()][temp_node[0]->getc()] << " next cutpoint " << xi_con[temp_node[0]->getv()][temp_node[0]->getc() + 1] << endl;
    }

    // cout << "load all trees " << endl;

    cout << "print con trees " << endl;
    for (size_t tt = 0; tt < mm_con; tt++)
    {
        cout << "index " << tt << endl;
        cout << t_con[tt] << endl;
    }

    cout << "print mod trees " << endl;
    for (size_t tt = 0; tt < mm_mod; tt++)
    {
        cout << "index " << tt << endl;
        cout << t_mod[tt] << endl;
    }
    /*****************************************************************************
  /* Setup the model
  *****************************************************************************/
    //--------------------------------------------------
    //trees
    //make trt_init a parameter later
    // std::vector<tree> t_mod(ntree_mod);
    // for(size_t i=0;i<ntree_mod;i++) t_mod[i].setm(trt_init/(double)ntree_mod);

    // std::vector<tree> t_con(ntree_con);
    // for(size_t i=0;i<ntree_con;i++) t_con[i].setm(ybar/(double)ntree_con);

    //--------------------------------------------------
    //prior parameters
    // PX scale parameter for b: b(x) = bscale*b_0(x), bscale ~ N(0,1/bscale_prec)
    double bscale_prec = 2;
    // double bscale0 = -0.5;
    // double bscale1 = 0.5;
    double bscale0 = bscale0_ini;
    double bscale1 = bscale1_ini;

    double mscale_prec = 1.0;
    // double mscale = 1.0;
    double mscale = muscale_ini;
    double delta_con = 1.0;
    double delta_mod = 1.0;

    pinfo pi_mod;
    pi_mod.pbd = 1.0; //prob of birth/death move
    pi_mod.pb = .5;   //prob of birth given  birth/death

    pi_mod.alpha = mod_alpha; //prior prob a bot node splits is alpha/(1+d)^beta, d is depth of node
    pi_mod.beta = mod_beta;   //2 for bart means it is harder to build big trees.
    // pi_mod.tau = con_sd / (sqrt(delta_mod) * sqrt((double)ntree_mod)); //sigma_mu, variance on leaf parameters

    cout << "pi mod tau " << con_sd / (sqrt(delta_mod) * sqrt((double)ntree_mod)) << " " << pi_mod_tau << endl;
    pi_mod.tau = pi_mod_tau;

    shat = sigma_ini;

    // pi_mod.sigma = shat; //resid variance is \sigma^2_y/bscale^2 in the backfitting update
    cout << "pi mod sigma " << shat << " " << pi_mod_sigma << endl;
    pi_mod.sigma = pi_mod_sigma;

    pinfo pi_con;
    pi_con.pbd = 1.0; //prob of birth/death move
    pi_con.pb = .5;   //prob of birth given  birth/death

    pi_con.alpha = con_alpha;
    pi_con.beta = con_beta;
    // pi_con.tau   = con_sd/(sqrt(delta_con)*sqrt((double) ntree_con)); //sigma_mu, variance on leaf parameters

    // pi_con.sigma = shat/fabs(mscale); //resid variance in backfitting is \sigma^2_y/mscale^2

    pi_con.tau = pi_con_tau;
    pi_con.sigma = pi_con_sigma;

    cout << "pi con tau " << con_sd/(sqrt(delta_con)*sqrt((double) ntree_con)) << " " << pi_con_tau << endl;
    cout << "pi con sigma " << shat/fabs(mscale) << " " << pi_con_sigma << endl;

    // double sigma = shat;
    double sigma = sigma_ini;

    //   // @Peter This is where dinfo is initialized

    //   //--------------------------------------------------
    //   //dinfo for control function m(x)
    // //  Rcout << "ybar " << ybar << endl;
    //   double* allfit_con = new double[n]; //sum of fit of all trees
    //   for(size_t i=0;i<n;i++) allfit_con[i] = ybar;
    //   double* r_con = new double[n]; //y-(allfit-ftemp) = y-allfit+ftemp
    //   dinfo di_con;
    //   di_con.n=n;
    //   di_con.p = p_con;
    //   di_con.x = &x_con[0];
    //   di_con.y = r_con; //the y for each draw will be the residual
    //   di_con.w = w;

    //   //--------------------------------------------------
    //   //dinfo for trt effect function b(x)
    //   double* allfit_mod = new double[n]; //sum of fit of all trees
    //   for(size_t i=0;i<n;i++) allfit_mod[i] = (z_[i]*bscale1 + (1-z_[i])*bscale0)*trt_init;
    //   double* r_mod = new double[n]; //y-(allfit-ftemp) = y-allfit+ftemp
    //   dinfo di_mod;
    //   di_mod.n=n;
    //   di_mod.p=p_mod;
    //   di_mod.x = &x_mod[0];
    //   di_mod.y = r_mod; //the y for each draw will be the residual
    //   di_mod.w = w;

    // allocate storage
    double *ftemp = new double[n]; //fit of current tree

    // never used???
    // double* precs = new double[n]; // temp storage for conditional ''precisions'' in heteroskedastic updates

    double *allfit_con = new double[n]; //sum of fit of all trees
    double *allfit_mod = new double[n]; //sum of fit of all trees
    double *allfit = new double[n];     //yhat
    double *r_con = new double[n];      //y-(allfit-ftemp) = y-allfit+ftemp
    double *r_mod = new double[n];      //y-(allfit-ftemp) = y-allfit+ftemp
    dinfo di_con;
    di_con.n = n;
    di_con.p = p_con;
    di_con.x = &x_con[0]; //the y for each draw will be the residual
    dinfo di_mod;
    di_mod.n = n;
    di_mod.p = p_mod;
    di_mod.x = &x_mod[0]; //the y for each draw will be the residual
    di_con.w = w;
    di_mod.w = w;

    // initialize allfit, residuals from trees read in
    for (size_t kk = 0; kk < n; kk++)
    {
        allfit[kk] = 0;
        allfit_con[kk] = 0;
        allfit_mod[kk] = 0;
    }

    for (size_t j = 0; j < ntree_con; j++)
    {
        fit(t_con[j], xi_con, di_con, ftemp);
        for (size_t k = 0; k < n; k++)
        {
            allfit[k] += mscale * ftemp[k];
            allfit_con[k] += mscale * ftemp[k];
            r_con[k] = (y[k] - allfit[k]) / mscale;
        }
    }

    for (size_t j = 0; j < ntree_mod; j++)
    {
        fit(t_mod[j], xi_mod, di_mod, ftemp);
        for (size_t k = 0; k < ntrt; k++)
        {
            allfit[k] += bscale1 * ftemp[k];
            allfit_mod[k] += bscale1 * ftemp[k];
            r_mod[k] = (y[k] - allfit[k]) / bscale1;
        }
        for (size_t k = ntrt; k < n; k++)
        {
            allfit[k] += bscale0 * ftemp[k];
            allfit_mod[k] += bscale0 * ftemp[k];
            r_mod[k] = (y[k] - allfit[k]) / bscale0;
        }
    }

    di_con.y = r_con;
    di_mod.y = r_mod;

    //--------------------------------------------------
    //dinfo and design for trt effect function out of sample
    //x for predictions
    dinfo di_mod_est;              //data information for prediction
    std::vector<double> x_mod_est; //stored like x
    size_t n_mod_est;
    //  if(x_mod_est_.size()) {
    for (NumericVector::iterator it = x_mod_est_.begin(); it != x_mod_est_.end(); ++it)
    {
        x_mod_est.push_back(*it);
    }
    n_mod_est = x_mod_est.size() / p_mod;
    //  Rcout << "n_mod_est " << n_mod_est << std::endl;
    if (x_mod_est.size() != n_mod_est * p_mod)
        stop("error, wrong number of elements in effect estimate data set\n");
    //if(n_mod_est)
    di_mod_est.n = n_mod_est;
    di_mod_est.p = p_mod;
    di_mod_est.x = &x_mod_est[0];
    di_mod_est.y = 0; //there are no y's!
    di_mod_est.w = w;

    //  }
    //--------------------------------------------------
    //storage for ouput
    //in sample fit
    /*
  //out of sample fit
  double* ppredmean=0; //posterior mean for prediction
  double* fpredtemp=0; //temporary fit vector to compute prediction
  if(dip.n) {
  ppredmean = new double[dip.n];
  fpredtemp = new double[dip.n];
  for(size_t i=0;i<dip.n;i++) ppredmean[i]=0.0;
  }
  //for sigma draw
  double rss, restemp;
  */

    //--------------------------------------------------
    //setup for random effects
    size_t random_dim = random_des.n_cols;
    int nr = 1;
    if (randeff)
        nr = n;

    arma::vec r(nr);           //working residuals
    arma::vec Wtr(random_dim); // W'r

    arma::mat WtW = random_des.t() * random_des; //W'W
    arma::mat Sigma_inv_random = diagmat(1 / (random_var_ix * random_var));

    // PX parameters
    arma::vec eta(random_var_ix.n_cols); //random_var_ix is num random effects by num variance components
    eta.fill(1.0);

    for (size_t k = 0; k < nr; ++k)
    {
        r(k) = y[k] - allfit_con[k] - allfit_mod[k];
    }

    Wtr = random_des.t() * r;
    arma::vec gamma = solve(WtW / (sigma * sigma) + Sigma_inv_random, Wtr / (sigma * sigma));
    arma::vec allfit_random = random_des * gamma;
    if (!randeff)
        allfit_random.fill(0);

    //--------------------------------------------------
    //storage for the fits
    // double* allfit = new double[n]; //yhat
    for (size_t i = 0; i < n; i++)
    {
        allfit[i] = allfit_mod[i] + allfit_con[i];
        if (randeff)
            allfit[i] += allfit_random[i];
    }
    // double* ftemp  = new double[n]; //fit of current tree

    NumericVector sigma_post(nd);
    NumericVector msd_post(nd);
    NumericVector bsd_post(nd);
    NumericMatrix m_post(nd, n);
    NumericMatrix yhat_post(nd, n);
    NumericMatrix b_post(nd, n);
    NumericMatrix b_est_post(nd, n_mod_est);
    arma::mat gamma_post(nd, gamma.n_elem);
    arma::mat random_var_post(nd, random_var.n_elem);
    NumericVector b0_vec(nd);
    NumericVector b1_vec(nd);
    NumericVector m_vec(nd);
    //  NumericMatrix spred2(nd,dip.n);

    /*
  //save stuff to tree file
  treef << xi << endl; //cutpoints
  treef << m << endl;  //number of trees
  treef << p << endl;  //dimension of x's
  treef << (int)(nd/thin) << endl;
  */

    //*****************************************************************************
    /* MCMC
   * note: the allfit objects are all carrying the appropriate scales
   */
    //*****************************************************************************
    Rcout << "\n============================================================\nBeginning MCMC:\n============================================================\n";
    time_t tp;
    int time1 = time(&tp);

    size_t save_ctr = 0;
    bool verbose_itr = false;

    // double* allfit_mod = new double[n]; //sum of fit of all trees

    double *weight = new double[n];
    double *weight_het = new double[n];

    logger.setLevel(0);

    // cout << "check fitted values " << endl;
    // for(size_t k=0; k < 100; k ++ ){
    //   // cout << ntrt << " " << k << " " << allfit[k] << " " << allfit_mod[k] + allfit_con[k] << " " << allfit_mod[k] << " " << allfit_con[k] << endl;
    //   cout << di_con.y[k] << " " << di_mod.y[k] << endl;
    // }

    for (size_t iIter = 0; iIter < (nd * thin + burn); iIter++)
    {
        // verbose_itr = iIter>=burn;
        verbose_itr = false;

        if (iIter % status_interval == 0)
        {
            Rcout << "iteration: " << iIter << " sigma: " << sigma << endl;
        }

        logger.setLevel(verbose_itr);

        logger.log("==============================================");
        sprintf(logBuff, "MCMC iteration: %d of %d, sigma %f", iIter + 1, nd * thin + burn, sigma);
        logger.log(logBuff);
        logger.log("==============================================");
        if (verbose_itr)
        {
            logger.getVectorHead(y, logBuff);
            Rcout << "           y: " << logBuff << "\n";

            logger.getVectorHead(allfit, logBuff);
            Rcout << "Current Fit : " << logBuff << "\n";
        }

        logger.log("=====================================");
        logger.log("- Tree Processing");
        logger.log("=====================================");

        //draw trees for m(x)
        for (size_t iTreeCon = 0; iTreeCon < ntree_con; iTreeCon++)
        {

            logger.log("==================================");
            sprintf(logBuff, "Updating Control Tree: %d of %d", iTreeCon + 1, ntree_con);
            logger.log(logBuff);
            logger.log("==================================");
            logger.startContext();

            logger.log("Attempting to Print Tree Pre Update \n");
            if (verbose_itr)
            {
                t_con[iTreeCon].pr(xi_con);
                Rcout << "\n\n";
            }

            // This seems to populate fTemp with the appopriate
            // mu_i for every observation in di_con
            // based on the tree defined by t and xi_con
            // it seems like fTemp is the only guy that is modified
            fit(t_con[iTreeCon], // tree& t
                xi_con,          // xinfo& xi
                di_con,          // dinfo& di
                ftemp);          // std::vector<double>& fv
                                 // if(iIter == 0 && iIter == 0){
                                 //   for(size_t ll = 0; ll < 10; ll ++ ){
                                 //     cout << "ftemp t_con " << ll << " " << ftemp[ll] << endl;
                                 //   }
                                 // }

            logger.log("Attempting to Print Tree Post first call to fit \n");
            if (verbose_itr)
            {
                t_con[iTreeCon].pr(xi_con);
                Rcout << "\n\n";
            }
            // I'd expect a residual calculation somewhere here.
            // so that's what this probably is?
            for (size_t k = 0; k < n; k++)
            {
                if (ftemp[k] != ftemp[k])
                {
                    Rcout << "control tree " << iTreeCon << " obs " << k << " " << endl;
                    Rcout << t_con[iTreeCon] << endl;
                    stop("nan in ftemp");
                }
                // Whatever these allfits used to be, we're subtracting off
                // ftempt from them. so assumping allfits started at the real y_hat
                // would be y_hat without this tree
                allfit[k] = allfit[k] - mscale * ftemp[k];
                allfit_con[k] = allfit_con[k] - mscale * ftemp[k];

                // If the above is correct, then r_con would be tree_y or equavalently Rj
                // I suspect (but have not confirmed) that because of pointers
                // updating r_con updates di_con.y, looking at how it's intialized
                // The name r_con is never used again till it's deleted, so that would make sense
                r_con[k] = (y[k] - allfit[k]) / mscale;

                if (r_con[k] != r_con[k])
                {
                    Rcout << (y[k] - allfit[k]) << endl;
                    Rcout << mscale << endl;
                    Rcout << r_con[k] << endl;
                    stop("NaN in resid");
                }
            }

            for (int k = 0; k < n; ++k)
            {
                weight[k] = w[k] * mscale * mscale / (sigma * sigma); // for non-het case, weights need to be divided by sigma square to make it similar to phi
            }
            if (verbose_itr)
            {
                logger.getVectorHead(weight, logBuff);
                Rcout << "\n weight: " << logBuff << "\n\n";
            }
            // It seems like with 100% probability, either a
            // birth or death will be attempted
            // but the change can still be rejected or accepted
            logger.log("Starting Birth / Death Processing");
            logger.startContext();
            bd(t_con[iTreeCon], // tree& x
               xi_con,          // xinfo& xi
               di_con,          // dinfo& di
               weight,          // phi
               pi_con,          // pinfo& pi
               gen,
               logger); // RNG& gen
            logger.stopContext();

            logger.log("Attempting to Print Tree Post db \n");
            if (verbose_itr)
            {
                t_con[iTreeCon].pr(xi_con);
                Rcout << "\n";
            }

            if (verbose_itr)
            {
                logger.log("Printing Current Status of Fit");

                logger.getVectorHead(z_, logBuff);
                // logger.log(logBuff);
                Rcout << "\n          z : " << logBuff << "\n";

                logger.getVectorHead(y, logBuff);
                Rcout << "          y : " << logBuff << "\n";

                logger.getVectorHead(allfit, logBuff);
                Rcout << "Fit - Tree  : " << logBuff << "\n";

                logger.getVectorHead(r_con, logBuff);
                Rcout << "     r_con  : " << logBuff << "\n\n";

                Rcout << " MScale: " << mscale << "\n";

                Rcout << " bscale0 : " << bscale0 << "\n";

                Rcout << " bscale1 : " << bscale1 << "\n\n";
            }
            logger.log("Starting To Draw Mu");
            logger.startContext();

            drmu(t_con[iTreeCon], // tree& x
                 xi_con,          // xinfo& xi
                 di_con,          // dinfo& di
                 pi_con,          // pinfo& pi,
                 weight,
                 gen); // RNG& gen

            logger.stopContext();

            logger.log("Attempting to Print Tree Post drmu \n");
            if (verbose_itr)
            {
                t_con[iTreeCon].pr(xi_con);
                Rcout << "\n";
            }

            fit(t_con[iTreeCon],
                xi_con,
                di_con,
                ftemp);
            // if(iIter == 0 && iIter == 0){
            //   for(size_t ll = 0; ll < 10; ll ++ ){
            //     cout << "ftemp t_con 2 " << ll << " " << ftemp[ll] << endl;
            //   }
            // }
            for (size_t k = 0; k < n; k++)
            {
                allfit[k] += mscale * ftemp[k];
                allfit_con[k] += mscale * ftemp[k];
            }

            logger.log("Attempting to Print tree Post second call to fit \n");

            if (verbose_itr)
            {
                t_con[iTreeCon].pr(xi_con);
                Rcout << "\n";
            }
            logger.stopContext();
        }

        //draw trees for b(x)
        for (size_t k = 0; k < ntrt; ++k)
        {
            weight_het[k] = w[k] * bscale1 * bscale1 / (sigma * sigma);
        }
        for (size_t k = ntrt; k < n; ++k)
        {
            weight_het[k] = w[k] * bscale0 * bscale0 / (sigma * sigma);
        }

        for (size_t iTreeMod = 0; iTreeMod < ntree_mod; iTreeMod++)
        {
            logger.log("==================================");
            sprintf(logBuff, "Updating Moderate Tree: %d of %d", iTreeMod + 1, ntree_mod);
            logger.log(logBuff);
            logger.log("==================================");
            logger.startContext();

            logger.log("Attempting to Print Tree Pre Update \n");
            if (verbose_itr)
            {
                t_mod[iTreeMod].pr(xi_mod);
                Rcout << "\n";
            }

            fit(t_mod[iTreeMod],
                xi_mod,
                di_mod,
                ftemp);
            // if(iIter == 0 && iIter == 0){
            //   for(size_t ll = 0; ll < 10; ll ++ ){
            //     cout << "ftemp t_mod " << ll << " " << ftemp[ll] << endl;
            //   }
            // }
            logger.log("Attempting to Print Tree Post first call to fit");
            if (verbose_itr)
            {
                t_mod[iTreeMod].pr(xi_mod);
                Rcout << "\n";
            }

            for (size_t k = 0; k < n; k++)
            {
                if (ftemp[k] != ftemp[k])
                {
                    Rcout << "moderator tree " << iTreeMod << " obs " << k << " " << endl;
                    Rcout << t_mod[iTreeMod] << endl;
                    stop("nan in ftemp");
                }
                double bscale = (k < ntrt) ? bscale1 : bscale0;
                allfit[k] = allfit[k] - bscale * ftemp[k];
                allfit_mod[k] = allfit_mod[k] - bscale * ftemp[k];
                r_mod[k] = (y[k] - allfit[k]) / bscale;
            }
            logger.log("Starting Birth / Death Processing");
            logger.startContext();
            bd(t_mod[iTreeMod],
               xi_mod,
               di_mod,
               weight_het,
               pi_mod,
               gen,
               logger);
            logger.stopContext();

            logger.log("Attempting to Print Tree  Post bd \n");
            if (verbose_itr)
            {
                t_mod[iTreeMod].pr(xi_mod);
                Rcout << "\n";
            }

            if (verbose_itr)
            {
                logger.log("Printing Status of Fit");

                logger.getVectorHead(z_, logBuff);
                Rcout << "\n          z : " << logBuff << "\n";

                logger.getVectorHead(y, logBuff);
                Rcout << "          y : " << logBuff << "\n";

                logger.getVectorHead(allfit, logBuff);
                Rcout << "Fit - Tree  : " << logBuff << "\n";

                logger.getVectorHead(r_mod, logBuff);
                Rcout << "     r_mod  : " << logBuff << "\n\n";

                Rcout << " MScale: " << mscale << "\n";

                Rcout << " bscale0 : " << bscale0 << "\n";

                Rcout << " bscale1 : " << bscale1 << "\n\n";
            }
            logger.log("Starting To Draw Mu");
            logger.startContext();
            drmu(t_mod[iTreeMod],
                 xi_mod,
                 di_mod,
                 pi_mod,
                 weight_het,
                 gen);
            logger.stopContext();

            logger.log("Attempting to Print Tree Post drmuhet \n");
            if (verbose_itr)
            {
                t_mod[iTreeMod].pr(xi_mod);
                Rcout << "\n";
            }

            fit(t_mod[iTreeMod],
                xi_mod,
                di_mod,
                ftemp);

            for (size_t k = 0; k < ntrt; k++)
            {
                allfit[k] += bscale1 * ftemp[k];
                allfit_mod[k] += bscale1 * ftemp[k];
            }
            for (size_t k = ntrt; k < n; k++)
            {
                allfit[k] += bscale0 * ftemp[k];
                allfit_mod[k] += bscale0 * ftemp[k];
            }

            logger.log("Attempting to Print Tree Post second call to fit");

            if (verbose_itr)
            {
                t_mod[iTreeMod].pr(xi_mod);
                Rcout << "\n";
            }
            logger.stopContext();

        } // end tree lop
        logger.log("=====================================");
        logger.log("- MCMC iteration Cleanup");
        logger.log("=====================================");

        if (use_bscale)
        {
            double ww0 = 0.0, ww1 = 0.;
            double rw0 = 0.0, rw1 = 0.;
            double s2 = sigma * sigma;
            for (size_t k = 0; k < n; ++k)
            {
                double bscale = (k < ntrt) ? bscale1 : bscale0;
                double w = s2 * bscale * bscale / (allfit_mod[k] * allfit_mod[k]);

                if (w != w)
                {
                    Rcout << " w " << w << endl;
                    stop("");
                }

                double randeff_contrib = randeff ? allfit_random[k] : 0.0;

                double r = (y[k] - allfit_con[k] - randeff_contrib) * bscale / allfit_mod[k];

                if (r != r)
                {
                    Rcout << "bscale " << k << " r " << r << " mscale " << mscale << " b*z " << allfit_mod[k] * z_[k] << " bscale " << bscale0 << " " << bscale1 << endl;
                    stop("");
                }
                if (k < ntrt)
                {
                    ww1 += 1 / w;
                    rw1 += r / w;
                }
                else
                {
                    ww0 += 1 / w;
                    rw0 += r / w;
                }
            }
            logger.log("Drawing bscale 1");
            logger.startContext();
            double bscale1_old = bscale1;
            double bscale_fc_var = 1 / (ww1 + bscale_prec);
            bscale1 = bscale_fc_var * rw1 + gen.normal(0., 1.) * sqrt(bscale_fc_var);
            if (verbose_itr)
            {

                Rcout << "Original bscale1 : " << bscale1_old << "\n";
                Rcout << "bscale_prec : " << bscale_prec << ", ww1 : " << ww1 << ", rw1 : " << rw1 << "\n";
                Rcout << "New  bscale1 : " << bscale1 << "\n\n";
            }
            logger.stopContext();

            logger.log("Drawing bscale 0");
            logger.startContext();
            double bscale0_old = bscale0;
            bscale_fc_var = 1 / (ww0 + bscale_prec);
            bscale0 = bscale_fc_var * rw0 + gen.normal(0., 1.) * sqrt(bscale_fc_var);
            if (verbose_itr)
            {
                Rcout << "Original bscale0 : " << bscale0_old << "\n";
                Rcout << "bscale_prec : " << bscale_prec << ", ww0 : " << ww0 << ", rw0 : " << rw0 << "\n";
                Rcout << "New  bscale0 : " << bscale0 << "\n\n";
            }
            logger.stopContext();

            for (size_t k = 0; k < ntrt; ++k)
            {
                allfit_mod[k] = allfit_mod[k] * bscale1 / bscale1_old;
            }
            for (size_t k = ntrt; k < n; ++k)
            {
                allfit_mod[k] = allfit_mod[k] * bscale0 / bscale0_old;
            }

            // update delta_mod
            if (b_half_normal)
            {
                logger.log("Updating delta_mod because b_half_normal");
                double ssq = 0.0;
                tree::npv bnv;
                typedef tree::npv::size_type bvsz;
                double endnode_count = 0.0;

                for (size_t iTreeMod = 0; iTreeMod < ntree_mod; iTreeMod++)
                {
                    bnv.clear();
                    t_mod[iTreeMod].getbots(bnv);
                    bvsz nb = bnv.size();
                    for (bvsz ii = 0; ii < nb; ++ii)
                    {
                        double mm = bnv[ii]->getm(); //node parameter
                        ssq += mm * mm / (pi_mod.tau * pi_mod.tau);
                        endnode_count += 1.0;
                    }
                }
                delta_mod = gen.gamma(0.5 * (1. + endnode_count), 1.0) / (0.5 * (1 + ssq));
            }
            if (verbose_itr)
            {
                Rcout << "Original pi_mod.tau : " << pi_mod.tau << "\n";
                pi_mod.tau = con_sd / (sqrt(delta_mod) * sqrt((double)ntree_mod));
                Rcout << "New pi_mod.tau : " << pi_mod.tau << "\n\n";
            }
        }
        else
        {
            bscale0 = -0.5;
            bscale1 = 0.5;
        }
        pi_mod.sigma = sigma;

        if (use_mscale)
        {
            double ww = 0.;
            double rw = 0.;
            double s2 = sigma * sigma;
            for (size_t k = 0; k < n; ++k)
            {
                double w = s2 * mscale * mscale / (allfit_con[k] * allfit_con[k]);
                if (w != w)
                {
                    Rcout << " w " << w << endl;
                    stop("");
                }

                double randeff_contrib = randeff ? allfit_random[k] : 0.0;

                double r = (y[k] - allfit_mod[k] - randeff_contrib) * mscale / allfit_con[k];
                if (r != r)
                {
                    Rcout << "mscale " << k << " r " << r << " mscale " << mscale << " b*z " << allfit_mod[k] * z_[k] << " bscale " << bscale0 << " " << bscale1 << endl;
                    stop("");
                }
                ww += 1 / w;
                rw += r / w;
            }

            logger.log("Drawing mscale");

            double mscale_old = mscale;
            double mscale_fc_var = 1 / (ww + mscale_prec);
            mscale = mscale_fc_var * rw + gen.normal(0., 1.) * sqrt(mscale_fc_var);
            if (verbose_itr)
            {
                Rcout << "Original mscale : " << mscale_old << "\n";
                Rcout << "mscale_prec : " << mscale_prec << ", ww : " << ww << ", rw : " << rw << "\n";
                Rcout << "New  mscale : " << mscale << "\n\n";
            }

            //Rcout<< mscale_fc_var << " " << rw <<" " << mscale << endl;

            for (size_t k = 0; k < n; ++k)
            {
                allfit_con[k] = allfit_con[k] * mscale / mscale_old;
            }

            // update delta_con

            double ssq = 0.0;
            tree::npv bnv;
            typedef tree::npv::size_type bvsz;
            double endnode_count = 0.0;

            for (size_t iTreeCon = 0; iTreeCon < ntree_con; iTreeCon++)
            {
                bnv.clear();
                t_con[iTreeCon].getbots(bnv);
                bvsz nb = bnv.size();
                for (bvsz ii = 0; ii < nb; ++ii)
                {
                    double mm = bnv[ii]->getm(); //node parameter
                    ssq += mm * mm / (pi_con.tau * pi_con.tau);
                    endnode_count += 1.0;
                }
            }

            delta_con = gen.gamma(0.5 * (1. + endnode_count), 1.0) / (0.5 * (1 + ssq));
            if (verbose_itr)
            {
                logger.log("Updating pi_con.tau");
                Rcout << "Original pi_con.tau : " << pi_con.tau << "\n";
                pi_con.tau = con_sd / (sqrt(delta_con) * sqrt((double)ntree_con));
                Rcout << "New pi_con.tau : " << pi_con.tau << "\n\n";
            }
        }
        else
        {
            mscale = 1.0;
        }
        pi_con.sigma = sigma / fabs(mscale); //should be sigma/abs(mscale) for backfitting

        //sync allfits after scale updates, if necessary. Could do smarter backfitting updates inline
        if (use_mscale || use_bscale)
        {
            logger.log("Sync allfits after scale updates");

            for (size_t k = 0; k < n; ++k)
            {
                double randeff_contrib = randeff ? allfit_random[k] : 0.0;
                allfit[k] = allfit_con[k] + allfit_mod[k] + randeff_contrib;
            }
        }

        if (randeff)
        {
            Rcout << "==================================\n";
            Rcout << "- Random Effects \n";
            Rcout << "==================================\n";

            //update random effects
            for (size_t k = 0; k < n; ++k)
            {
                r(k) = y[k] - allfit_con[k] - allfit_mod[k];
                allfit[k] -= allfit_random[k];
            }

            Wtr = random_des.t() * r;

            arma::mat adj = diagmat(random_var_ix * eta);
            //    Rcout << adj << endl << endl;
            arma::mat Phi = adj * WtW * adj / (sigma * sigma) + Sigma_inv_random;
            arma::vec m = adj * Wtr / (sigma * sigma);
            //Rcout << m << Phi << endl << Sigma_inv_random;
            gamma = rmvnorm_post(m, Phi);

            //Rcout << "updated gamma";

            // Update px parameters eta

            arma::mat adj2 = diagmat(gamma) * random_var_ix;
            arma::mat Phi2 = adj2.t() * WtW * adj2 / (sigma * sigma) + arma::eye(eta.size(), eta.size());
            arma::vec m2 = adj2.t() * Wtr / (sigma * sigma);
            //Rcout << m << Phi << endl << Sigma_inv_random;
            eta = rmvnorm_post(m2, Phi2);

            //Rcout << "updated eta";

            // Update variance parameters

            arma::vec ssqs = random_var_ix.t() * (gamma % gamma);
            //Rcout << "A";
            arma::rowvec counts = sum(random_var_ix, 0);
            //Rcout << "B";
            for (size_t ii = 0; ii < random_var_ix.n_cols; ++ii)
            {
                random_var(ii) = 1.0 / gen.gamma(0.5 * (random_var_df + counts(ii)), 1.0) * 2.0 / (random_var_df + ssqs(ii));
            }
            //Rcout << "updated vars" << endl;
            Sigma_inv_random = diagmat(1 / (random_var_ix * random_var));

            allfit_random = random_des * diagmat(random_var_ix * eta) * gamma;

            //Rcout << "recom allfit vars" << endl;

            for (size_t k = 0; k < n; ++k)
            {
                allfit[k] = allfit_con[k] + allfit_mod[k] + allfit_random(k); //+= allfit_random[k];
            }
        }

        logger.log("Draw Sigma");
        //draw sigma
        double rss = 0.0;
        double restemp = 0.0;
        for (size_t k = 0; k < n; k++)
        {
            restemp = y[k] - allfit[k];
            rss += w[k] * restemp * restemp;
        }
        sigma = sqrt((nu * lambda + rss) / gen.chi_square(nu + n));
        pi_con.sigma = sigma / fabs(mscale);
        pi_mod.sigma = sigma;

        if (((iIter >= burn) & (iIter % thin == 0)))
        {

            msd_post(save_ctr) = fabs(mscale) * con_sd;
            bsd_post(save_ctr) = fabs(bscale1 - bscale0) * mod_sd;

            m_vec(save_ctr) = mscale;
            b0_vec(save_ctr) = bscale0;
            b1_vec(save_ctr) = bscale1;

            gamma_post.row(save_ctr) = (diagmat(random_var_ix * eta) * gamma).t();
            random_var_post.row(save_ctr) = (sqrt(eta % eta % random_var)).t();

            sigma_post(save_ctr) = sigma;
            for (size_t k = 0; k < n; k++)
            {
                m_post(save_ctr, k) = allfit_con[k];
                yhat_post(save_ctr, k) = allfit[k];
            }
            for (size_t k = 0; k < n; k++)
            {
                double bscale = (k < ntrt) ? bscale1 : bscale0;
                b_post(save_ctr, k) = (bscale1 - bscale0) * allfit_mod[k] / bscale;
            }
            //if(di_mod_est.n) {
            for (size_t k = 0; k < di_mod_est.n; k++)
            {
                b_est_post(save_ctr, k) = (bscale1 - bscale0) * fit_i(k, t_mod, xi_mod, di_mod_est);
            }
            //}
            save_ctr += 1;
        }
    } // end MCMC Loop

    int time2 = time(&tp);
    Rcout << "\n============================================================\n MCMC Complete \n============================================================\n";

    Rcout << "time for loop: " << time2 - time1 << endl;

    t_mod.clear();
    t_con.clear();
    delete[] allfit;
    delete[] allfit_mod;
    delete[] allfit_con;
    delete[] r_mod;
    delete[] r_con;
    delete[] ftemp;

    treef.close();

    return (List::create(_["yhat_post"] = yhat_post, _["b_post"] = b_post, _["b_est_post"] = b_est_post,
                         _["sigma"] = sigma_post, _["msd"] = msd_post, _["bsd"] = bsd_post,
                         _["gamma"] = gamma_post, _["random_var_post"] = random_var_post,
                         _["mscale_vec"] = m_vec, _["bscale0_vec"] = b0_vec, _["bscale1_vec"] = b1_vec));
}

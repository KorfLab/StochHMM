/*
 *  distributions.cpp
 *
 *  Created by Paul Lott on 4/22/09.
 *  Copyright 2009 University of California, Davis. All rights reserved.
 *
 */

#include "distributions.h"

namespace StochHMM{
    int maximum=1000000;
    
    
    
    //////////////////////// DISTRIBUTIONS ///////////////////////////////////
    
    
    ////// DISCRETE & FINITE //////
    
    
    //!Binomial Distribution Survival Function
    //!\param [out] dist
    //!\param trials Number of trials in experiment
    //!\param prob Probability of success
    void sBinomial (std::vector<double>&dist, int trials, double prob){
        //http://en.wikipedia.org/wiki/Binomial_distribution
        double cdf=0;
        for (int i=1; i<maximum;i++){
            double pmf=bin_coef(trials,i)*pow(prob,i)*pow(1-prob,trials-i);
            cdf+=pmf;
            double val=1-cdf;
            if (val<=0){
                break;
            }
            else{
                dist.push_back(val);
            }
        }
        return;
    }
    
    
    //!Binomial Cumulative Distribution
    //!\param [out] dist
    //!\param trials Number of trials in experiment
    //!\param prob Probability of success
    void cBinomial (std::vector<double>&dist, int trials, double prob){
        //http://en.wikipedia.org/wiki/Binomial_distribution
        double cdf=0;
        for (int i=1; i<maximum;i++){
            double pmf=bin_coef(trials,i)*pow(prob,i)*pow(1-prob,trials-i);
            cdf+=pmf;
            double val=cdf;
            if (val>=1.0){
                dist.push_back(1.0);
                break;
            }
            else{
                dist.push_back(val);
            }
        }
        return;
    }
    
    
    
    //!BetaBinomial Survival Function
    //!<a href = "http://en.wikipedia.org/wiki/Beta-binomial_model">
    //!\param [out] dist
    //!\param trials Number of trials
    //!\param a alpha
    //!\param b beta
    void sBetaBinomial (std::vector<double>& dist, int trials, double a, double b){
        double cdf = 0;
        for (int i=1; i<maximum;i++){
            double newAlpha = (double) i + a;
            double newBeta  = (double)(trials-i) + b;
            double pmf      = bin_coef(trials,i) * (beta(newAlpha,newBeta)/ beta(a,b));
            cdf+=pmf;
            double val=1-cdf;
            
            if (val<=0){
                dist.push_back(1.0);
                break;
            }
            else{
                dist.push_back(val);
            }
        }
        return;
        
    }
    
    //!BetaBinomial Cumulative Distribution
    //!<a href = "http://en.wikipedia.org/wiki/Beta-binomial_model">
    //!\param [out] dist
    //!\param trials Number of trials
    //!\param a alpha
    //!\param b beta
    void cBetaBinomial (std::vector<double>& dist, int trials, double a, double b){
        double cdf = 0;
        for (int i=1; i<maximum;i++){
            double newAlpha = (double) i + a;
            double newBeta  = (double)(trials-i) + b;
            double pmf      = bin_coef(trials,i) * (beta(newAlpha,newBeta)/ beta(a,b));
            cdf+=pmf;
            double val=cdf;
            
            if (val>=1.0){
                dist.push_back(1.0);
                break;
            }
            else{
                dist.push_back(val);
            }
        }
        return;
    }
    
    
    
    //!Degenerate Survival Function
    //!Support consists of only one value
    //!\param [out] dist
    //!\param value
    void sDegenerate(std::vector<double>&dist, double value){
        for (int i=0;i<value;i++){
            dist.push_back(1);
        }
        dist.push_back(0);
        return;
    }
    
    
    //!Degenerate Cumulative Distribution
    //!Support consists of only one value
    //!\param [out] dist
    //!\param value
    //!<a href = "http://en.wikipedia.org/wiki/Degenerate_distribution">
    void cDegenerate(std::vector<double>&dist, double value){
        for (int i=0;i<value;i++){
            dist.push_back(0);
        }
        dist.push_back(1);
        return;
    }
    
    
    //!Discrete Uniform Survival Function
    //!<a href = "http://en.wikipedia.org/wiki/Uniform_distribution_(discrete)">
    void sUniform (std::vector<double> &dist,int lower,int upper){
        for(int k=0;k<=upper;k++){
            if (k<lower){
                dist.push_back(1);
            }
            else if (k>=lower && k<=upper){
                double val=(k-lower+1)/(upper-lower+1);
                dist.push_back(1-val);
            }
        }
        return;
    }
    
    //!Discrete Uniform Cumulative Distribution
    //!<a href = "http://en.wikipedia.org/wiki/Uniform_distribution_(discrete)">
    void cUniform (std::vector<double> &dist,int lower,int upper){
        for(int k=0;k<=upper;k++){
            if (k<lower){
                dist.push_back(0);
            }
            else if (k>=lower && k<=upper){
                double val=(k-lower+1)/(upper-lower+1);
                dist.push_back(val);
            }
        }
        return;
    }
    
    
    
    
    //!Complement CDF Hypergeometric Distribution
    //!\param [out] dist
    //!\param N Population Size
    //!\param m Number of successes
    //!\param n Number of draws
    void sHypergeometric (std::vector<double> &dist,int N, int m, int n){
        //http://en.wikipedia.org/wiki/Hypergeometric_distribution
        double cdf=0;
        for( int k=0;k<maximum;k++){ //k = number of successes
            double pmf=((double)bin_coef(m,k) * (double)bin_coef(N-m,n-k))/(double) bin_coef(N,n);
            cdf+=pmf;
            double val=1-cdf;
            if (val<=0){
                break;
            }
            else{
                dist.push_back(val);
            }
        }
        return;
    }
    
    
    
    //! CDF Hypergeometric Distribution
    //!\param [out] dist
    //!\param N Population Size
    //!\param m Number of successes
    //!\param n Number of draws
    void cHypergeometric (std::vector<double> &dist,int N, int m, int n){
        //http://en.wikipedia.org/wiki/Hypergeometric_distribution
        double cdf=0;
        for( int k=0;k<maximum;k++){ //k = number of successes
            double pmf=((double)bin_coef(m,k) * (double)bin_coef(N-m,n-k))/(double) bin_coef(N,n);
            cdf+=pmf;
            double val=cdf;
            if (val>=1){
                dist.push_back(1.0);
                break;
            }
            else{
                dist.push_back(val);
            }
        }
        return;
    }
    
    
    
    
    //Cauchy Distribution
    void sCauchy(std::vector<double> &dist,double location, double scale){
        //http://en.wikipedia.org/wiki/Cauchy_distribution
        for (int i=0; i<maximum;i++){
            double val=1-((1/M_PI)*atan((i-location)/scale)+0.5);
            if (val==0){
                break;
            }
            else{
                dist.push_back(val);
            }
        }
        return;
    }
    
    //! CDF of Chi-Squared Distribution
    //! \param x Chi-Square Value
    //! \param df Degrees of freedom
    //! \return double Value of CDF at x
    double cChiSquared(double x,double df){
        double val = 1.0/tgamma(df/2.0);
        double div = igamma_lower(df/2.0,x/2.0);
        return val*div;
    }
    
    
    //! Complement CDF/Survival of Chi-Squared Distribution
    //! \param [out] dist Vector to store distribution in
    //! \param df Degrees of freedom
    //ChiSquared Distribution
    void sChiSquared(std::vector<double> &dist, double df){
        //http://en.wikipedia.org/wiki/Chi-square_distribution
        for (double i=1; i<maximum;i++){
            //CDF Value
            double val=1.0/tgamma(df/2.0);
            double div=igamma_lower(df/2.0, i/2.0);
            val*=div;
            
            //CCDF value
            val = 1-val;
            if (fabs(val)<0.00001){
                break;
            }
            else{
                dist.push_back(val);
            }
        }
        return;
    }
    
    
    
    
    
    
    //!Complement CDF Exponential Distribution
    //!\param [out] dist
    //!\param lambda Value of lambda to use
    //!<a href = "http://en.wikipedia.org/wiki/Exponential_distribution">
    void sExponential (std::vector<double> &dist, double lambda){
        for (int i=0; i<maximum;i++){
            double val=exp(-1*lambda*i);
            if (val==0){
                break;
            }
            else{
                dist.push_back(val);
            }
        }
        return;
    }
    
    //!CDF Exponential Distribution
    //!\param [out] dist
    //!\param lambda Value of lambda to use
    //!<a href = "http://en.wikipedia.org/wiki/Exponential_distribution">
    void cExponential (std::vector<double> &dist, double lambda){
        for (int i=0; i<maximum;i++){
            double val=1-exp(-1*lambda*i);
            if (val==0){
                break;
            }
            else{
                dist.push_back(val);
            }
        }
        return;
    }
    
    
    
    
    
    //
    //void sExtremeValue(vector<double> &dist, double location, double shape, double scale){
    //    if (scale<=0){
    //        std::cerr << "Extreme value doesn't support scale <=0" <<std::endl;
    //        return;
    //    }
    //	//http://en.wikipedia.org/wiki/Extreme_value_distribution
    //}
    //
    //void cExtremeValue(vector<double> &dist, double location, double shape, double scale){
    //    if (scale<=0){
    //        std::cerr << "Extreme value doesn't support scale <=0" <<std::endl;
    //        return;
    //    }
    //
    //    if (shape == 0){
    //
    //    }
    //    else{
    //        for(int i=0;i<maximum;i++){
    //            double tx = exp((i-location)/shape
    //            double val = 1-exp(-1*tx);
    //        }
    //    }
    //	//http://en.wikipedia.org/wiki/Extreme_value_distribution
    //}
    
    
    
    //F-Distribution Distribution
    void sFDistribution(std::vector<double> &dist, double dOne, double dTwo){
        //http://en.wikipedia.org/wiki/F_distribution
        for (int i=0;i<maximum;i++){
            double val=1-ribeta((dOne*i)/((dOne*i)+dTwo), dOne/2, dTwo/2);
            if (val==0){
                break;
            }
            else{
                dist.push_back(val);
            }
        }
        return;
    }
    
    
    //F-Distribution Distribution
    void cFDistribution(std::vector<double> &dist, double dOne, double dTwo){
        //http://en.wikipedia.org/wiki/F_distribution
        for (int i=0;i<maximum;i++){
            double val=ribeta((dOne*i)/((dOne*i)+dTwo), dOne/2, dTwo/2);
            if (val==1){
                dist.push_back(1.0);
                break;
            }
            else{
                dist.push_back(val);
            }
        }
        return;
    }
    
    
    
    
    
    //Gamma distribution
    void sGamma(std::vector<double> &dist, double shape, double scale){
        //http://en.wikipedia.org/wiki/Gamma_distribution
        for (int i=0;i<maximum;i++){
            double val=1-igamma_lower(shape, i/scale)/tgamma(shape);
            if (val==0){
                break;
            }
            else{
                dist.push_back(val);
            }
        }
        return;
    }
    
    //Gamma distribution
    void cGamma(std::vector<double> &dist, double shape, double scale){
        //http://en.wikipedia.org/wiki/Gamma_distribution
        for (int i=0;i<maximum;i++){
            double val=igamma_lower(shape, i/scale)/tgamma(shape);
            if (val==1){
                dist.push_back(1.0);
                break;
            }
            else{
                dist.push_back(val);
            }
        }
        return;
    }
    
    
    //Geometric Distribution
    void sGeometric (std::vector<double> &dist, double p){
        //http://en.wikipedia.org/wiki/Geometric_distribution
        for (int i=0;i<maximum;i++){
            double val=pow((1-p),i+1);
            if (val==0){
                break;
            }
            else{
                dist.push_back(val);
            }
        }
        return;
    }
    
    void cGeometric (std::vector<double> &dist, double p){
        //http://en.wikipedia.org/wiki/Geometric_distribution
        for (int i=0;i<maximum;i++){
            double val=1-pow((1-p),i+1);
            if (val==1){
                dist.push_back(1.0);
                break;
            }
            else{
                dist.push_back(val);
            }
        }
        return;
    }
    
    
    
    
    //Laplace Distribution
    void sLaplace(std::vector<double> &dist , double location, double scale){
        //http://en.wikipedia.org/wiki/Laplace_distribution
        for( int i=0;i<maximum;i++){
            double val=(1/(2*scale))*exp(-1*(fabs(i-location)/scale));
            if (val==0){
                break;
            }
            else{
                dist.push_back(val);
            }
        }
        return;
    }
    
    //Log-Normal Distribution
    void sLogNormal(std::vector<double> &dist,double location, double scale){
        //http://en.wikipedia.org/wiki/Lognormal_distribution
        for( int i=0;i<maximum;i++){
            double val=0.5+0.5*(erf((log(i-location)/sqrt(2*pow(scale,2)))));
            if (val==0){
                break;
            }
            else{
                dist.push_back(val);
            }
        }
        return;
    }
    
    
    //Logarithmic Distribution
    void sLogarithmic(std::vector<double> &dist, double prob){
        //http://en.wikipedia.org/wiki/Logarithmic_distribution
        double cdf=0;
        for( int i=1;i<maximum;i++){
            double pmf=(-1/log(1-prob))*((pow(prob,i)/i));
            cdf+=pmf;
            double val=1-cdf;
            if (val==0){
                break;
            }
            else{
                dist.push_back(val);
            }
        }
        return;
    }
    
    //Logistic Distribution
    void sLogistic(std::vector<double> &dist, double location, double scale){
        //http://en.wikipedia.org/wiki/Logistic_distribution
        for( int i=0;i<maximum;i++){
            double val=1-(1/(1+exp(-1*(i-location)/scale)));
            if (val==0){
                break;
            }
            else{
                dist.push_back(val);
            }
        }
        return;
    }
    
    
    
    void sMaxwellBoltzman(std::vector<double> &dist, double scale){
        //http://en.wikipedia.org/wiki/Maxwell–Boltzmann_distribution
        for( int i=0;i<maximum;i++){
            double error=erf(i/(sqrt(2)*scale));
            double i_squared=-1*pow((double)i,2);
            double val=sqrt(2*(1/M_PI))*((double)i*exp(i_squared/(2*pow(scale,2)))/scale);
            double boltz=1-(error-val);
            if (boltz==0){
                break;
            }
            else{
                dist.push_back(boltz);
            }
        }
        return;
    }
    
    //Negative Binomial Distribution
    void sNegativeBinomial(std::vector<double> &dist,int r, double p){
        //http://en.wikipedia.org/wiki/Negative_binomial_distribution
        for( int i=0;i<maximum;i++){
            double val=ribeta(p, i+1, r);
            if (val==0){
                break;
            }
            else{
                dist.push_back(val);
            }
        }
        return;
    }
    
    
    
    /* 
     void sNonCentralChiSquared(vector<double> &dist, double v, double lmda){
     //http://en.wikipedia.org/wiki/Noncentral_chi-square_distribution
     double cdf=0;
     for( int i=1;i<maximum;i++){
     double pmf=(-1/log(1-prob))*((pow(prob,i)/i));
     cdf+=pmf;
     double val=1-cdf;
     if (val==0){
     break;
     }
     else{
     dist.push_back(val);
     }
     }
     return;
     
     }
     
     
     
     //Non-Central F distribution
     void sNonCentralF(vector<double> &dist,double vone, double vtwo, double lmda){
     //http://en.wikipedia.org/wiki/Noncentral_F-distribution
     }
     
     
     //Non-Central T distribution
     void sNonCentralT(vector<double> &dist, double v, double lambda){
     //http://en.wikipedia.org/wiki/Noncentral_t-distribution
     }
     */
    
    //Normal distribution
    void sNormal(std::vector<double> &dist, double mean, double stdev){
        //http://en.wikipedia.org/wiki/Normal_distribution
        for( int i=0;i<maximum;i++){
            double val=1-(0.5*(1+erf((i-mean)/sqrt(2*pow(stdev,2)))));
            if (val==0){
                break;
            }
            else{
                dist.push_back(val);
            }
        }
        return;
        
    }
    
    //Pareto Distribution
    void sPareto(std::vector<double> &dist,double shape, double scale){
        for( int i=0;i<maximum;i++){
            double val=pow((scale/(double) i),shape);
            if (val==0){
                break;
            }
            else{
                dist.push_back(val);
            }
        }
        return;
    }
    
    //Poisson Distribution
    void sPoisson(std::vector<double> &dist, double lambda){
        //http://en.wikipedia.org/wiki/Poisson_distribution
        for( int i=0;i<maximum;i++){
            double val=igamma_upper(i+1,lambda)/factorial(i);
            if (val==0){
                break;
            }
            else{
                dist.push_back(val);
            }
        }
        return;
    }
    
    
    //Rayleigh Distribution
    void sRayleigh(std::vector<double> &dist, double sigma){
        //http://en.wikipedia.org/wiki/Rayleigh_distribution
        for( int i=0;i<maximum;i++){
            double val=exp((-1*(pow((double) i, 2)))/(2*pow(sigma,2)));;
            if (val==0){
                break;
            }
            else{
                dist.push_back(val);
            }
        }
        return;
    }
    
    //Triangular distribution
    void sTriangular(std::vector<double> &dist,int a, int b, int c){
        //http://en.wikipedia.org/wiki/Triangular_distribution
        double val=0;
        for(int x=0;x<=b;x++){
            if (x<a){
                dist.push_back(1);
            }
            else if (x>=a && x<=c){
                val=pow((double)(x-a),2)/((b-a)*(c-a));
                dist.push_back(1-val);
            }
            else{
                val=pow((double)(b-x),2)/((b-a)*(b-c));
                dist.push_back(val);
            }
        }
        return;
    }
    
    
    
    
    //User defined distribution - Converts pdf to survival distrib.
    void sUser (std::vector<double>&dist, std::vector<double> &prob_dist){
        double c_df=0;
        //calculate CDF
        for(int i=0;i<prob_dist.size();i++){
            c_df+=prob_dist[i];
            dist.push_back(c_df);
        }
        //Convert CDF to survival function
        for(int i=0; i<dist.size();i++){
            dist[i]=1-dist[i];
        }
        return;
    }
    
    //Weibull-Survival distribution 
    void sWeibull(std::vector<double> &dist,double shape, double scale){
        //http://en.wikipedia.org/wiki/Weibull_distribution
        for (int i=0; i<maximum;i++){
            double val=exp(-1*pow(i/scale,shape));
            if (val==0){
                break;
            }
            else{
                dist.push_back(val);
            }
        }
        return;
    }
};




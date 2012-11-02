//
//  PDF.cpp
//  StochHMM
//
//  Created by Paul Lott on 10/11/12.
//  Copyright (c) 2012 Korf Lab, Genome Center, UC Davis, Davis, CA. All rights reserved.
//

#include "PDF.h"

namespace StochHMM{
    
    
    /*-------------------- Discrete with finite support ------------- */
    
    //!Binomial probability mass function
    //!param k Number of successes
    //!param n Number of trials
    //!param p Probability of success
    float binomial_pdf(int k,int n, float p){
        return bin_coef(n,k)*pow(p,k)*pow(1-p,n-k);
    }
    
    
    //!Binomial probaility mass function
    //!param k Number of successes
    //!param param std::vector<float> (n=number of trials, p=probability of success)
    float binomial_pdf(int k, std::vector<float> param){
        return binomial_pdf(k, param[0],param[1]);
    }
    
    
    //!BetaBinomial probability mass function
    //!<a href = "http://en.wikipedia.org/wiki/Beta-binomial_model">
    //!\param [out] dist
    //!\param n Number of trials
    //!\param a alpha
    //!\param b beta
    float beta_binomial_pdf(int n, int k, double a, double b){
        float newAlpha = (double) k + a;
        float newBeta  = (double)(n-k) + b;
        float pmf      = bin_coef(n,k) * (beta(newAlpha,newBeta)/ beta(a,b));
        return pmf;
    }
    
    //!Degenerate probability mass function
    //!param value Current value
    //!param k Degenerate value
    float degenerate_pdf(float value, float k){
        if (value == k){
            return 1.0;
        }
        else{
            return 0.0;
        }
    }
    
    
    //!Discrete Uniform CDF
    //!param a Minimum position
    //!param b Maximum position
    //!param position Position to calculate
    float discrete_uniform_pdf(int a, int b, int position){
        if (position<a || position>b){
            return 0;
        }
        else{
            return (float)1/(b-a+1);
        }
    }
    
    
    //!Hypergeometric Cumulative Distribution Function
    //!param n  Number of draws from Population
    //!param N  Size of population
    //!param m  Number of successes in Population
    //!param k  Number of successes
    float hypergeometric_pdf(int n, int N, int m, int k){
        float prob = (bin_coef(m, k)*bin_coef(N-m,n-k))/bin_coef(N, n);
        return prob;
    }
    
    //TODO:  Need to implement poisson binomial
    //!Poisson Binomial Probability mass function
    //!param n Number for trials
    //!param k Number of successful trials
    //!param p Probability of Success for each trial
    float poisson_binomial_pdf(int n, int k, std::vector<float> p){
        return 0.0;
    }
    
    
    /*-------------------- Discrete with infinite support ------------- */

    //!Beta negative binomial pdf
    //!param n  Number of successful trials
    //!param k  Shape parameter
    //!param a  Shape parameter
    //!param b  Shape parameter
    float beta_negative_binomial_pdf(int n, int k, float a, float b){
        double pmf = bin_coef(n+k-1, k);
        pmf*=gamma(a+n)*gamma(b+k)*gamma(a+b)/(gamma(a+b+n+k)*gamma(a)*gamma(b));
        return pmf;
    }
    
    
    //!Maxwell-Boltzman Probability distribution function
    //!<a href = "http://en.wikipedia.org/wiki/Maxwell–Boltzmann_distribution">
    //!param a Shape parameter  a>0
    //!param x 
    float maxwell_boltzman_pdf(double a, float x){
        float pmf = sqrtf(2/PI)*(pow(x,2)*exp(-1*pow(x,2)/(2*pow(a, 2))))/pow(a,3);
        return pmf;
    }
    
    
    //!Geometric probability mass function
    //!param k Trial of first success
    //!param p probability of success
    float geometric_pdf(int k, double p){
        float pmf = pow(1-p,k-1) * p;
        return pmf;
    }
    
    //TODO: check previous function parameters
    
    //!Logarithmic probability mass function
    //!param k
    //!param p Probability
    float logarithmic_pdf(int k, double p){
        if (k<1 || p==1.0 || p==0.0){
            std::cerr << "Logarithmic PDF Function:  k must be >=1, 0<p<1\n";
            exit(2);
        }
        float pmf = (1/log(1-p))*(pow(p,k)/k);
        return pmf;
    }
    
    
    //!Negative binomial probability mass function
    //!param r  Number of failures until experiment stopped
    //!param p  Success probability in each experiment
    float negative_binomial_pdf(int k, int r, double p){
        if (r<=0 || p<0.0 || p > 1.0){
            std::cerr << "Negative Binomial PMF: Incorrect parameters\n";
            exit(2);
        }
        float pmf = bin_coef(k+r-1, k)*pow(1-p,r)*pow(p,k);
        return pmf;
    }
    
    
    //!Poisson probability mass function
    //!param k  Trial of first success
    //!param lambda Probability of success
    float poisson_pdf(int k, float lambda){
        if (k<0 || lambda<=0){
            std::cerr << "Poisson PMF: Incorrect parameters\n";
            exit(2);
        }
        float pmf = (pow(lambda,(double)k)* exp(-1*lambda))/factorial(k);
        return pmf;
    }
    
    
    //!Yule-Simon probability mass function
    //!param k  
    //!param p  Shape parameter
    float yule_simon_pdf(int k, double p){
        if (k<1 || p<=0.0){
            std::cerr << "Yule-Simon PMF: Incorrect parameters\n";
            exit(2);
        }
        float pmf = p*beta(k, p+1);
        return pmf;
    }
    
    
    //!Zipf probability mass function
    //!param k  The rank of element
    //!param N  Number of elements
    //!param s  Shape parameter (Exponent value)
    float zipf_pdf(int k,int N, float s){
        if (k<1 || N<1 || s<=0.0){
            std::cerr << "Zipf's PMF: Incorrect parameters\n";
            exit(2);
        }
        float pmf = 1/pow((double)k,s);
        float denom(0.0);
        for (int n = 1; n < N; n++){
            denom+=1/pow((double)n, s);
        }
        return pmf/denom;
    }
    
    
    //!Zipf-Mandelbrot probability mass function
    //!param k  The rank of element
    //!param N  Number of elements
    //!param s  Shape parameter (Exponent value)
    //!param q  Shape parameter
    float zipf_mandelbrot_pdf(int k,int N, float s, float q){
        if (k<1 || N<1 || s<=0.0 || q<0.0){
            std::cerr << "Zipf's PMF: Incorrect parameters\n";
            exit(2);
        }
        float pmf = 1/pow((double)k+q,s);
        float denom(0.0);
        for (int n = 1; n < N; n++){
            denom+=1/pow((double)n+q, s);
        }
        return pmf/denom;
    }
    
    
    
    /*-------------------- Continuous within Interval ------------------------ */
    
    //!Arcsine probability distribution function
    //!param x value 0<=X<=1
    float arcsine_pdf(float x){
        if (x > 1.0 || x < 0.0 ){
            std::cerr << "Arcsine PDF: Incorrect parameters\n";
            exit(2);
        }
        return 1/(PI*sqrtf(x*(1-x)));
    }
    
    
    //!Beta probability distribution function
    //!param x  Value 0<x<1
    //!param a  Shape parameter a>0
    //!param b  Shape parameter b>0
    float beta_pdf(float x, float a, float b){
        if (x > 1.0 || x < 0.0 || a <= 0.0 || b <= 0.0){
            std::cerr << "Beta PDF: Incorrect parameters\n";
            exit(2);
        }
        return betaPDF(x, a, b);
    }
    
    //Beta probability distribution function
    //!param x  Value 0<x<1
    //!param a  Shape parameter a>0
    //!param b  Shape parameter b>0
    float beta_pdf(float x, std::vector<float> param){
        return beta_pdf(x,param[0],param[1]);
    }
    
    //!Logit Normal probability distribution function
    //!param x  Value
    //!param mu Mean
    //!param sigma  Std. deviation
    float logit_normal_pdf(double x, double mu, double sigma){
        if (x >= 1.0 || x <= 0.0 || sigma <= 0.0){
            std::cerr << "Logit Normal PDF: Incorrect parameters\n";
            exit(2);
        }
        float pdf = (1/(sigma*sqrtf(2*PI)))*exp(-1*pow((logit(x)-mu),2)/(2*pow(sigma,2)))*(1/(x*(1-x)));
        return pdf;
    }
    
    
    //!Continuous Uniform probability distribution function
    //!param a Minimum position
    //!param b Maximum position
    //!param position Position to calculate
    float continuous_uniform_pdf(float a, float b, float position){
        if (b<=a){
            std::cerr << "Continuous uniform PDF: Incorrect parameters\n";
            exit(2);
        }
        else if (position<a || position>b){
            return 0.0f;
        }
        else{
            return 1/(b-a);
        }
    }
    
    //!Kumaraswamy probability distribution function
    //!param x
    //!param a Shape parameter
    //!param b Shape parameter
    float kumaraswamy_pdf(float x, float a, float b){
        if (a <= 0 || b <= 0 || x < 0 || x > 1.0 ){
            std::cerr << "Kumaraswamy PDF:  Incorrect parameters\n";
            exit(2);
        }
        float pdf = a * b * powf(x, a-1) * powf((1 - powf(x, a)),b-1);
        return pdf;
    }
    
    //!Raised cosine probability distribution function
    //!param mu
    //!param s   S>0
    //!param x  u-s<=x<=u+s
    float raised_cosine_pdf(float x, float mu, float s){
        if (s <= 0 || x < mu-s || x > mu+s){
            std::cerr << "Raised cosine PDF:  Incorrect parameters\n";
            exit(2);
        }
        float pdf = 1/(2*s)*(1+cosf(PI*(x-mu)/s));
        return pdf;
    }
    
    //!Triangular probability distribution function
    //!param a  Real number
    //!param b  a<b
    //!param c  a<=c<=b
    //!param x  a<=x<=b
    float triangular_pdf(float x, float a, float b, float c){
        if (b<=a || c<=a || c>=b){
            std::cerr << "Triangular PDF:  Incorrect parameters\n";
            exit(2);
        }
        
        
        if (x >= a && x <= c){
            return 2*(x-a)/((b-a)*(c-a));
        }
        else if (x > c && x <= b){
            return 2*(b-x)/((b-a)*(b-c));
        }
        else{
            return 0;
        }
    
    }
    
    
    //!Truncated Normal probability distribution function
    //!param a  minimum value
    //!param b  maximum value
    //!param mu Mean
    //!param sd Standard deviation
    //!param x  Value
    float truncated_normal_pdf(float x, float mu, float sd, float a, float b){
        if (sd<0 || x < a || x > b){
            std::cerr << "Truncated normal PDF:  Incorrect paramenters\n";
            exit(2);
        }
        
        float alpha = (a-mu)/sd;
        float beta  = (b-mu)/sd;
        float xi = (x-mu)/sd;
        
        double alpha_cdf = (0.5*(1+erf((alpha-mu)/sqrt(2*pow(sd,2)))));
        double beta_cdf = (0.5*(1+erf((beta-mu)/sqrt(2*pow(sd,2)))));
        float z = beta_cdf - alpha_cdf;
        
        double pdf = (0.5*(1+erf((xi-mu)/sqrt(2*pow(sd,2)))))/(sd * z);
        return pdf;
    }
    
    
    //!U-quadratic probability distribution function
    //!param a
    //!param b
    //!param x
    float u_quadratic_pdf(float x, float a, float b){
        if (b <= a || x < a || x > b){
            std::cerr << "U-quadratic PDF: Incorrect parameters\n";
            exit(2);
        }
        
        float alpha = 12 / powf(b-a, 3);
        float beta = (b + a)/2;
        return alpha * powf(x-beta, 2);;
    }
    
    
    //!Wigner semicircle probability distribution function
    //!param r radius R>0
    //!param x value
    float wigner_semicircle_pdf(float x, float r){
        if (r<=0){
            std::cerr << "Wigner semicircle PDF: Incorrect parameters\n";
            exit(2);
        }
        
        return (2/(PI*powf(r, 2)))*sqrtf(powf(r, 2) - powf(x, 2));
    }
    
    
    /*-------------- Continuous within Semi-infinite Interval ---------------- */
    
    //!Beta prime probability distribution function
    //!param x Value
    //!param a  Alpha Shape parameter
    //!param b  Beta Shape parameter
    float beta_prime_pdf(float x, float a, float b){
        if (a<=0 || b<=0 || x<=0){
            std::cerr << "Beta prime PDF: Incorrect parameters\n";
            exit(2);
        }
        float pdf = (powf(x, a-1)*powf(1+x, -1*a-b))/ beta(a, b);
        return pdf;
    }
    
    
    //!Chi probability distribution function
    //!param x Value x>=0
    //!param k degrees of freedom k>0
    float chi_pdf(float x, float k){
        if (k<=0){
            std::cerr << "Chi PDF: Incorrect parameters\n";
            exit(2);
        }
        float pdf  = (powf(2, 1-(k/2))*powf(x, k-1)*exp(-1*(powf(x, 2)/2)))/gamma(k/2);
        return pdf;
        
    }
    
    
    //!Chi-squared probability distribution function
    //!param k K>0
    //!param x x>0
    float chi_squared_pdf(float x, float k){
        if ( k < 0 || x < 0){
            std::cerr << "Chi-squared PDF: Incorrect parameters\n";
            exit(2);
        }
        
        float pdf = (powf(x, (k/2)-1)*exp(-1*x/2))/(powf(2, k/2)*gamma(k/2));
        return pdf;
    }
    
    
    
    //!Inverse-Chi-squared probability distribution function
    //!param x Value x>0
    //!param v v>0
    float inverse_chi_squared_pdf(float x, float v){
        if ( x <= 0 || v <= 0 ){
            std::cerr << "Inverse Chi squared PDF: Incorrect parameters\n";
            exit(2);
        }
        
        float pdf = (powf(2, -1*v/2)*powf(x, -1*(v/2)-1)*exp(-1/(2*x)))/gamma(v/2);
        return pdf;
    }
    
    //!Scaled Inverse Chi-squared probability distribution function
    //!param x Value x>0
    //!param v  v>0
    //!param sigma_sqrd  sigma_sqrd>0
    float scaled_inverse_chi_squared_pdf(float x, float v, float sigma_sqrd){
        if (x <= 0 || v<=0 || sigma_sqrd <= 0 ){
            std::cerr << "Scaled Inverse Chi-squared PDF: Incorrect paramenters\n";
            exit(2);
        }
        float numerator = powf(sigma_sqrd*(v/2),v/2);
        numerator *= exp((-1*v*sigma_sqrd)/(2*x));
        
        float denominator = gamma(v/2);
        denominator *= powf(x, 1+(v/2));
        
        return numerator/denominator;
    }
    
    
    //!Dagum probability distribution function
    //!param p Shape parameter p>0
    //!param a Shape parameter a>0
    //!param b Shape parameter b>0
    //!param x Value x>0
    float dagum_pdf(float x, float p, float a, float b){
        if ( p <= 0 || a <= 0 || b <= 0 || x <= 0){
            std::cerr << "Dagum PDF: Incorrect parameters\n";
            exit(2);
        }
        float pdf = (a*p)/x;
        pdf *= powf(x/b, a*p)/powf(powf(x/b, a)+1 , p+1);
        return pdf;
    }
    
    
    //!Exponential probability distribution function
    //!param lambda  Rate or Inverse scale lambda>0
    //!param x  Value x>=0
    float exponential_pdf(float x, float lambda){
        if ( lambda <= 0 || x < 0 ){
            std::cerr << "Exponential PDF: Incorrect parameters\n";
            exit(2);
        }
        
        float pdf = lambda * exp (-1 * lambda * x);
        return pdf;
    }
    
    
    //!F-Distribution
    //!param x Value x>=0
    //!param d1 Degree of freedom d1>0
    //!param d2 Degree of freedom d2>0
    float f_pdf(float x, float d1, float d2){
        if (x < 0 || d1 <= 0 || d2 <= 0){
            std::cerr << "F Distribution PDF: Incorrect parameters\n";
            exit(2);
        }
        float numerator = sqrtf((powf(d1*x, d1)*powf(d2, d2))/powf(d1*x+d2, d1+d2));
        float denominator = x * beta(d1/2, d2/2);
        return numerator / denominator;
    }
    
    
    //!Fisher's z-distribution
    //!param x Value
    //!param d1 Degree of freedom d1>0
    //!param d2 Degree of freedom d2>0
    float fishers_z_pdf(float x, float d1, float d2){
        if (d1 <= 0 || d2 <= 0){
            std::cerr << "Fisher's z-Distribution PDF: Incorrect parameters\n";
            exit(2);
        }
        float numerator = 2*powf(d1, d1/2)*powf(d2, d2/2)*exp(d1*x);
        float denominator = beta(d1/2, d2/2) * powf(d1*exp(2*x)+d2, (d1+d2)/2);
        return numerator/denominator;
    }
    
    
    //!Folded Normal probability distribution function
    //!param x Value x>=0
    //!param mu Mean(location)
    //!param sigma_sqrd Scale
    float folded_normal_pdf(float x, float mu, float sigma_sqrd){
        if (x < 0 || sigma_sqrd <= 0 ){
            std::cerr << "Folded Normal PDF: Incorrect parameters\n";
            exit(2);
        }
        float l_pdf = 1/(sqrtf(sigma_sqrd*2 * PI));
        float l_temp = exp(-1* powf(-1*x-mu,2)/(2*sigma_sqrd));
        float r_temp = exp(-1* powf(x-mu,2)/(2*sigma_sqrd));
        
        float pdf = (l_pdf*l_temp) + (l_pdf *r_temp);
        return pdf;
    }
    
    
    //!Frechet Probability distribution function
    //!param x Value X>m
    //!param alpha Shape parameter a>0
    //!param s  Scale parameter s>0 (default s=1)
    //!param m  Location minimum (default m=0)
    float frechet_pdf(float x, float alpha, float s, float m){
        if (alpha <= 0 || s <= 0 || x < m ){
            std::cerr << "Frechet PDF: Incorrect parameters\n";
            exit(2);
        }
        
        float pdf = (alpha/s)*powf((x-m)/s, -1-alpha) * exp(-1*powf((x-m)/s, -1*alpha));
        return pdf;
    }
    
}
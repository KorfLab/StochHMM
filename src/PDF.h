//
//  PDF.h
//  StochHMM
//
//  Created by Paul Lott on 10/11/12.
//  Copyright (c) 2012 Korf Lab, Genome Center, UC Davis, Davis, CA. All rights reserved.
//

#ifndef __StochHMM__PDF__
#define __StochHMM__PDF__

#include <iostream>
#include <math.h>
#include "stochMath.h"

namespace StochHMM {
    //Discrete with finite support
    float binomial_pdf(int k ,int n,float p);
    float binomial_pdf(int k, std::vector<float> param);
    float beta_binomial_pdf(int n, int k, double a, double b);
    float degenerate_pdf(float value, float k);
    float discrete_uniform_pdf(int a,int b,int position);
    float hypergeometric_pdf(int n,int N,int m,int k);
    float poisson_binomial_pdf(int n, int k,std::vector<float> p);
        
    //Discrete with infinite support
    float beta_negative_binomial_pdf(int n, int k, float a, float b);
    float maxwell_boltzman_pdf(float a, float x);
    float geometric_pdf(int k, double p);
    float logarithmic_pdf(int k, double p);
    float negative_binomial_pdf(int k, int r, double p);
    float poisson_pdf(int k, float lambda);
    float yule_simon_pdf(int k, double p);
    float zipf_pdf(int k,int N, float s);
    float zipf_mandelbrot_pdf(int k,int N, float s, float q);
    
    //Continuous Distribution, bounded interval
    float arcsine_pdf(float x);
    float beta_pdf(float x, float a , float b);
    float beta_pdf(float x, std::vector<float> param);
    float logit_normal_pdf(double x, double mu, double sigma);
    float continuous_uniform_pdf(float,float,float);
    float kumaraswamy_pdf(float x, float a, float b);
    float raised_cosine_pdf(float x, float mu, float s);
    float triangular_pdf(float x, float a, float b, float c);
    float truncated_normal_pdf(float x, float mu, float sd, float a, float b);
    float u_quadratic_pdf(float x, float a, float b);
    float wigner_semicircle_pdf(float x, float r);

    
    //Continuous Distributions, semi-bounded interval
    float beta_prime_pdf(float x, float a, float b);
    float chi_pdf(float x, float k);
    float chi_squared_pdf(float x, float k);
    float inverse_chi_squared_pdf(float x, float v);
    float scaled_inverse_chi_squared_pdf(float x, float v, float sigma_sqrd);
    float dagum_pdf(float x, float p, float a, float b);
    float exponential_pdf(float x, float lambda);
    float f_pdf(float x, float d1, float d2);
    float fishers_z_pdf(float x, float d1, float d2);
    float folded_normal_pdf(float x, float mu, float sigma_sqrd);
    float frechet_pdf(float x, float alpha, float s, float m);
    

    
    double gamma_pdf(double);
    double gamma_pdf(double, double, double);
    double gamma_pdf(double, std::vector<double>);
    
    double std_normal_pdf(double);
    double normal_pdf(double,double,double);
    double folded_normal_distribution_pdf(double);
    double half_normal_distribution_pdf(double);
    double log_normal_pdf(double);

    
    double fisher_z_pdf(double);
    
    double hotelling_t_square_pdf(double);
    double inverse_gaussian_pdf(double);
    double log_logistic_pdf(double);
    double pareto_pdf(double);
    
   
    
    
    double beta_pdf(double, std::vector<double>);
    
    double multinomial_pdf(std::vector<int>,int , std::vector<double>);
    
   
    
    
}

#endif /* defined(__StochHMM__PDF__) */

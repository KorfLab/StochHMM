//statistics.h
 //Copyright (c) 2007-2012 Paul C Lott 
 //University of California, Davis
 //Genome and Biomedical Sciences Facility
 //UC Davis Genome Center
 //Ian Korf Lab
 //Website: www.korflab.ucdavis.edu
 //Email: lottpaul@gmail.com
 //
 //Permission is hereby granted, free of charge, to any person obtaining a copy of
 //this software and associated documentation files (the "Software"), to deal in
 //the Software without restriction, including without limitation the rights to
 //use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of
 //the Software, and to permit persons to whom the Software is furnished to do so,
 //subject to the following conditions:
 //
 //The above copyright notice and this permission notice shall be included in all
 //copies or substantial portions of the Software.
 //
 //THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 //IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS
 //FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR
 //COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER
 //IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN
 //CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

#ifndef STATISTICS_H_
#define STATISTICS_H_

#include <vector>
#include <math.h>
#include <iostream>
namespace StochHMM{


    #define PI 3.145926535897932
    
    //Discrete with finite support
    float discrete_uniform_cdf(int,int,int);
    double binomial_cdf(double);
    double hypergeometric_cdf(double);
    
    
    double poisson_binomial_cdf(double);
    double degenerate_distribution(double);
    double beta_binomial_cdf(double);
    
    //Discrete with infinite support
    double boltzman_cdf(double);
    double extended_negative_binomial_cdf(double);
    double geometric_cdf(double);
    double logarithmic_cdf(double);
    double negative_binomial_cdf(double);
    double parabolic_fractal_cdf(double);
    double poisson_cdf(double);
    double skellam_cdf(double);
    double yule_simon_cdf(double);
    double zeta(double);
    
    //Continuous Distributions
    float continuous_uniform_cdf(float,float,float);
    double beta_prime_cdf(double);
    double chi_cdf(double);
    double noncentral_chi_cdf(double);
    double chi_square_cdf(double);
    double exponential_cdf(double);
    double F_cdf(double);
    double gamma_cdf(double);
    double fisher_z_cdf(double);
    double folded_normal_distribution(double);
    double half_normal_distribution(double);
    double hotelling_t_square_cdf(double);
    double inverse_gaussian_cdf(double);
    double log_logistic_cdf(double);
    double log_normal_cdf(double);
    double pareto_cdf(double);
    
    
    template <class T> T min(std::vector<T>&);
    template <class T> T max(std::vector<T>&);
    template <class T> T construct_histogram( std::vector<T>&, int);
    template <class T> T smooth_histogram( std::vector<T>& , int, int, int);
    
    
    float entropy(std::vector<float> &set);
    float rel_entropy(std::vector<float> &set, std::vector<float> &set2);
    
    
    double _integrate(double (*funct)(double, std::vector<double>),double, double, std::vector<double> );
    double integrate(double (*funct)(double, std::vector<double>),double, double, std::vector<double>, double, double);
    
    double simpson(double (*funct)(double,std::vector<double>),double alpha, double beta,double lower, double upper);
    double adapt_simpson(double (*funct)(double, double, double),double alpha, double beta, double lower, double upper, double max_error, double sum);
    
    double gamma_func(double);
    double gamma_pdf(double, double, double);
    double gamma_pdf(double, std::vector<double>);
    double gamma_cdf(double, double, double);
    
    double chi2_pdf(double, double);
    double chi2_cdf(double, double);
    
    double beta_pdf(double,double,double);
    double beta_pdf(double, std::vector<double>);
    double beta_cdf(double,double,double);
    
    double expon_pdf(double,double);
    double expon_cdf(double,double);
    
    double normal_pdf(double,double,double);
    
    float bin_coef(int,int);
    float summation(float(*funct)(int,std::vector<float>), int, int, std::vector<float>);
    
    float binomial_pdf(int,int,float);
    float binomial_pdf(int, std::vector<float>);
    float binomial_cdf(int, int, float);
    
    double multinomial_pdf(std::vector<int>,int , std::vector<double>);
    
    double _low_igamma(double, std::vector<double>);
    double low_igamma(double, double);
    double upper_igamma(double, double);
    double erf(double);
    double std_normal_cdf(double);

    
}
#endif /*INC_FILENAME_H*/

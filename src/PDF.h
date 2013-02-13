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
	/*! \defgroup PDF Probability Density Functions
	 *	Univariate and Multivariate probability density functions
	 *	for use in continuous emissions.
	 *	@{
	 */
	
    //Discrete with finite support
    double binomial_pdf(int k ,int n,double p);
    double binomial_pdf(int k, std::vector<double> param);
    double beta_binomial_pdf(int k, int n, double a, double b);
    double degenerate_pdf(double value, double k);
    double discrete_uniform_pdf(int position, int a,int b);
    double hypergeometric_pdf(int k, int n,int N,int m);
    double poisson_binomial_pdf(int k, std::vector<double>& p);
	
		
        
    //Discrete with infinite support
    double beta_negative_binomial_pdf(int k, int n, double a, double b);
    double maxwell_boltzman_pdf(double x, double a);
    double geometric_pdf(int k, double p);
    double logarithmic_pdf(int k, double p);
    double negative_binomial_pdf(int k, int r, double p);
    double poisson_pdf(int k, double lambda);
    double yule_simon_pdf(int k, double p);
    double zipf_pdf(int k,int N, double s);
    double zipf_mandelbrot_pdf(int k,int N, double s, double q);
    
	    
	//Continuous Distribution, bounded interval
    double arcsine_pdf(double x);
    double beta_pdf(double x, double a , double b);
    double beta_pdf(double x, std::vector<double> param);
    double logit_normal_pdf(double x, double mu, double sigma);
    double continuous_uniform_pdf(double,double,double);
    double kumaraswamy_pdf(double x, double a, double b);
    double raised_cosine_pdf(double x, double mu, double s);
    double triangular_pdf(double x, double a, double b, double c);
    double truncated_normal_pdf(double x, double mu, double sd, double a, double b);
    double u_quadratic_pdf(double x, double a, double b);
    double wigner_semicircle_pdf(double x, double r);
	
	

    
    //Continuous Distributions, on a semi-bounded interval
    double beta_prime_pdf(double x, double a, double b);
    double chi_pdf(double x, double k);
    double chi_squared_pdf(double x, double k);
    double inverse_chi_squared_pdf(double x, double v);
    double scaled_inverse_chi_squared_pdf(double x, double v, double sigma_sqrd);
    double dagum_pdf(double x, double p, double a, double b);
    double exponential_pdf(double x, double lambda);
    double f_pdf(double x, double d1, double d2);
    double fishers_z_pdf(double x, double d1, double d2);
    double folded_normal_pdf(double x, double mu, double sigma_sqrd);
    double frechet_pdf(double x, double alpha, double s, double m);
    double gamma_pdf(double x, double alpha, double beta);
	double inv_gamma_pdf(double x, double alpha, double beta);
	double half_normal_pdf(double x, double sigma);
	double inv_gaussian_pdf(double x, double mu, double lambda);
	double levy_pdf(double x,double mu, double scale);
	double log_cauchy_pdf(double x,double mu, double sigma);
	double log_laplace_pdf(double x, double mu, double b);
	double log_logistic_pdf(double x, double a, double b);
	double log_normal_pdf(double x, double mu, double sigma_sqrd);
	double pareto_pdf(double x, double alpha, double x_m);
	double nakagami_pdf(double x, double mu, double w);
    double rayleigh_pdf(double x, double sigma);
	double gumbel_type_two_pdf(double x, double a, double b);
	double weibull_distribution(double x, double lambda, double k);
	
	
	//Continuous Distributions, on unbounded interval
	double cauchy_pdf(double x, double x_o, double gamma);
	double gumbel_pdf(double x, double mu, double beta);
	double generalized_normal_pdf(double x, double mu, double alpha, double beta);
	double hyperbolic_secant_pdf(double x);
	double laplace_pdf(double x, double mu, double b);
	double logistic_pdf(double x, double mu, double s);
	double standard_normal_pdf(double x);
	double normal_pdf(double x, double mu, double sigma);
	double students_t_pdf(double x, double v);
	double gumbel_type_one_pdf(double x, double a, double b);
	double generalized_extreme_value_pdf(double x, double mu, double sigma, double xi);
	double generalized_pareto_pdf(double x, double mu, double sigma, double xi);
	
	
	//Multivariate Distributions
	double dirichlet_pdf(std::vector<double>& x, std::vector<double>& alpha);
	double multivariate_ewens_pdf(std::vector<double>& x, double theta);

	
	
	//Discrete with finite support
	inline double binomial_pdf(const double k,const std::vector<double>* param){
		if (param==NULL){
			std::cerr << __FUNCTION__ << ": Missing parameters\n";
			exit(2);
		}
		return log(binomial_pdf(k, (*param)[0], (*param)[2]));
	}
	
    inline double beta_binomial_pdf(const double k, const std::vector<double>* param){
		if (param==NULL){
			std::cerr << __FUNCTION__ << ": Missing parameters\n";
			exit(2);
		}
		return log(beta_binomial_pdf(k, (*param)[0], (*param)[1], (*param)[2]));
	}
	
	inline double degenerate_pdf(const double value, const std::vector<double>* param){
		if (param==NULL){
			std::cerr << __FUNCTION__ << ": Missing parameters\n";
			exit(2);
		}
		return log(degenerate_pdf(value, (*param)[0]));
	}
	
	inline double discrete_uniform_pdf(const double position, const std::vector<double>* param){
		if (param==NULL){
			std::cerr << __FUNCTION__ << ": Missing parameters\n";
			exit(2);
		}
		return log(discrete_uniform_pdf(position, (*param)[0], (*param)[1]));
	}
	
	inline double hypergeometric_pdf(const double k, const std::vector<double>* param){
		if (param==NULL){
			std::cerr << __FUNCTION__ << ": Missing parameters\n";
			exit(2);
		}
		return log(hypergeometric_pdf(k, (*param)[0], (*param)[1], (*param)[2]));
	}
	
	
	//Discrete with infinite support
	inline double beta_negative_binomial_pdf(const double k, const std::vector<double>* param){
		if (param==NULL){
			std::cerr << __FUNCTION__ << ": Missing parameters\n";
			exit(2);
		}
		return log(beta_negative_binomial_pdf(k, (*param)[0], (*param)[1], (*param)[2]));
	}
	
    inline double maxwell_boltzman_pdf(const double x, const std::vector<double>* param){
		if (param==NULL){
			std::cerr << __FUNCTION__ << ": Missing parameters\n";
			exit(2);
		}
		return log(maxwell_boltzman_pdf(x, (*param)[0]));
	}
	
	
    inline double geometric_pdf(const double k, const std::vector<double>* param){
		if (param==NULL){
			std::cerr << __FUNCTION__ << ": Missing parameters\n";
			exit(2);
		}
		return log(geometric_pdf(k, (*param)[0]));
	}
	
    inline double logarithmic_pdf(const double k, const std::vector<double>* param){
		if (param==NULL){
			std::cerr << __FUNCTION__ << ": Missing parameters\n";
			exit(2);
		}
		return log(logarithmic_pdf(k, (*param)[0]));
	}
	
    inline double negative_binomial_pdf(const double k, const std::vector<double>* param){
		if (param==NULL){
			std::cerr << __FUNCTION__ << ": Missing parameters\n";
			exit(2);
		}
		return log(negative_binomial_pdf(k, (*param)[0], (*param)[1]));
	}
	
    inline double poisson_pdf(const double k, const std::vector<double>* param){
		if (param==NULL){
			std::cerr << __FUNCTION__ << ": Missing parameters\n";
			exit(2);
		}
		return log(poisson_pdf(k, (*param)[0]));
	}
	
    inline double yule_simon_pdf(const double k, const std::vector<double>* param){
		if (param==NULL){
			std::cerr << __FUNCTION__ << ": Missing parameters\n";
			exit(2);
		}
		return log(yule_simon_pdf(k, (*param)[0]));
	}
	
    inline double zipf_pdf(const double k, const std::vector<double>* param){
		if (param==NULL){
			std::cerr << __FUNCTION__ << ": Missing parameters\n";
			exit(2);
		}
		return log(zipf_pdf(k, (*param)[0], (*param)[1]));
	}
    inline double zipf_mandelbrot_pdf(const double k, const std::vector<double>* param){
		if (param==NULL){
			std::cerr << __FUNCTION__ << ": Missing parameters\n";
			exit(2);
		}
		return log(zipf_mandelbrot_pdf(k, (*param)[0], (*param)[1], (*param)[2]));
	}

	
	//Continuous Distribution, bounded interval
	inline double arcsine_pdf(const double x, const std::vector<double>* param){
		return log(arcsine_pdf(x));
	}
    
	inline double beta_pdf(const double x, const std::vector<double>* param){
		if (param==NULL){
			std::cerr << __FUNCTION__ << ": Missing parameters\n";
			exit(2);
		}
		return log(beta_pdf(x, (*param)[0], (*param)[1]));
	}
	
    inline double logit_normal_pdf(const double x, const std::vector<double>* param){
		if (param==NULL){
			std::cerr << __FUNCTION__ << ": Missing parameters\n";
			exit(2);
		}
		return log(logit_normal_pdf(x, (*param)[0], (*param)[1]));
	}
    
	inline double continuous_uniform_pdf(const double x, const std::vector<double>* param){
		if (param==NULL){
			std::cerr << __FUNCTION__ << ": Missing parameters\n";
			exit(2);
		}
		return log(continuous_uniform_pdf(x, (*param)[0], (*param)[1]));
	}
    
	inline double kumaraswamy_pdf(const double x, const std::vector<double>* param){
		if (param==NULL){
			std::cerr << __FUNCTION__ << ": Missing parameters\n";
			exit(2);
		}
		return log(kumaraswamy_pdf(x, (*param)[0], (*param)[1]));
	}
	
    inline double raised_cosine_pdf(const double x, const std::vector<double>* param){
		if (param==NULL){
			std::cerr << __FUNCTION__ << ": Missing parameters\n";
			exit(2);
		}
		return log(raised_cosine_pdf(x, (*param)[0], (*param)[1]));
	}
	
    inline double triangular_pdf(const double x, const std::vector<double>* param){
		if (param==NULL){
			std::cerr << __FUNCTION__ << ": Missing parameters\n";
			exit(2);
		}
		return log(triangular_pdf(x, (*param)[0], (*param)[1], (*param)[2]));
	}
	
    inline double truncated_normal_pdf(const double x, const std::vector<double>* param){
		if (param==NULL){
			std::cerr << __FUNCTION__ << ": Missing parameters\n";
			exit(2);
		}
		return log(truncated_normal_pdf(x, (*param)[0], (*param)[1], (*param)[2], (*param)[3]));
	}
	
    inline double u_quadratic_pdf(const double x, const std::vector<double>* param){
		if (param==NULL){
			std::cerr << __FUNCTION__ << ": Missing parameters\n";
			exit(2);
		}
		return log(u_quadratic_pdf(x, (*param)[0], (*param)[1]));
	}
	
    inline double wigner_semicircle_pdf(const double x, const std::vector<double>* param){
		if (param==NULL){
			std::cerr << __FUNCTION__ << ": Missing parameters\n";
			exit(2);
		}
		return log(wigner_semicircle_pdf(x, (*param)[0]));
	}
	
	
	
	//Continuous Distributions, on unbounded interval
	inline double beta_prime_pdf(const double x, const std::vector<double>* param){
		if (param==NULL){
			std::cerr << __FUNCTION__ << ": Missing parameters\n";
			exit(2);
		}
		return log(beta_prime_pdf(x, (*param)[0], (*param)[1]));
	}
	
    inline double chi_pdf(const double x, const std::vector<double>* param){
		if (param==NULL){
			std::cerr << __FUNCTION__ << ": Missing parameters\n";
			exit(2);
		}
		return log(chi_pdf(x,(*param)[0]));
	}
    
	inline double chi_squared_pdf(const double x, const std::vector<double>* param){
		if (param==NULL){
			std::cerr << __FUNCTION__ << ": Missing parameters\n";
			exit(2);
		}
		return log(chi_squared_pdf(x, (*param)[0]));
	}
	
    inline double inverse_chi_squared_pdf(const double x, const std::vector<double>* param){
		if (param==NULL){
			std::cerr << __FUNCTION__ << ": Missing parameters\n";
			exit(2);
		}
		return log(inverse_chi_squared_pdf(x, (*param)[0]));
	}
	
    inline double scaled_inverse_chi_squared_pdf(const double x, const std::vector<double>* param){
		if (param==NULL){
			std::cerr << __FUNCTION__ << ": Missing parameters\n";
			exit(2);
		}
		return log(scaled_inverse_chi_squared_pdf(x, (*param)[0], (*param)[1]));
	}
	
    inline double dagum_pdf(const double x, const std::vector<double>* param){
		if (param==NULL){
			std::cerr << __FUNCTION__ << ": Missing parameters\n";
			exit(2);
		}
		return log(dagum_pdf(x, (*param)[0], (*param)[1], (*param)[2]));
	}
	
    inline double exponential_pdf(const double x, const std::vector<double>* param){
		if (param==NULL){
			std::cerr << __FUNCTION__ << ": Missing parameters\n";
			exit(2);
		}
		return log(exponential_pdf(x, (*param)[0]));
	}
	
    inline double f_pdf(const double x, const std::vector<double>* param){
		if (param==NULL){
			std::cerr << __FUNCTION__ << ": Missing parameters\n";
			exit(2);
		}
		return log(f_pdf(x, (*param)[0], (*param)[1]));
	}
	
    inline double fishers_z_pdf(const double x, const std::vector<double>* param){
		if (param==NULL){
			std::cerr << __FUNCTION__ << ": Missing parameters\n";
			exit(2);
		}
		return log(fishers_z_pdf(x, (*param)[0], (*param)[1]));
	}
	
    inline double folded_normal_pdf(const double x, const std::vector<double>* param){
		if (param==NULL){
			std::cerr << __FUNCTION__ << ": Missing parameters\n";
			exit(2);
		}
		return log(folded_normal_pdf(x, (*param)[0], (*param)[1]));
	}
	
    inline double frechet_pdf(const double x, const std::vector<double>* param){
		if (param==NULL){
			std::cerr << __FUNCTION__ << ": Missing parameters\n";
			exit(2);
		}
		return log(frechet_pdf(x, (*param)[0], (*param)[1], (*param)[2]));
	}
	
    inline double gamma_pdf(const double x, const std::vector<double>* param){
		if (param==NULL){
			std::cerr << __FUNCTION__ << ": Missing parameters\n";
			exit(2);
		}
		return log(gamma_pdf(x, (*param)[0], (*param)[1]));
	}
	
	inline double inv_gamma_pdf(const double x, const std::vector<double>* param){
		if (param==NULL){
			std::cerr << __FUNCTION__ << ": Missing parameters\n";
			exit(2);
		}
		return log(inv_gamma_pdf(x, (*param)[0], (*param)[1]));
	}
	
	inline double half_normal_pdf(const double x, const std::vector<double>* param){
		if (param==NULL){
			std::cerr << __FUNCTION__ << ": Missing parameters\n";
			exit(2);
		}
		return log(half_normal_pdf(x, (*param)[0]));
	}
	
	inline double inv_gaussian_pdf(const double x, const std::vector<double>* param){
		if (param==NULL){
			std::cerr << __FUNCTION__ << ": Missing parameters\n";
			exit(2);
		}
		return log(inv_gaussian_pdf(x, (*param)[0], (*param)[1]));
	}
	
	inline double levy_pdf(const double x, const std::vector<double>* param){
		if (param==NULL){
			std::cerr << __FUNCTION__ << ": Missing parameters\n";
			exit(2);
		}
		return log(levy_pdf(x, (*param)[0], (*param)[1]));
	}
	
	inline double log_cauchy_pdf(const double x, const std::vector<double>* param){
		if (param==NULL){
			std::cerr << __FUNCTION__ << ": Missing parameters\n";
			exit(2);
		}
		return log(log_cauchy_pdf(x, (*param)[0], (*param)[1]));
	}
	
	inline double log_laplace_pdf(const double x, const std::vector<double>* param){
		if (param==NULL){
			std::cerr << __FUNCTION__ << ": Missing parameters\n";
			exit(2);
		}
		return log(log_laplace_pdf(x, (*param)[0], (*param)[1]));
	}
	
	inline double log_logistic_pdf(const double x, const std::vector<double>* param){
		if (param==NULL){
			std::cerr << __FUNCTION__ << ": Missing parameters\n";
			exit(2);
		}
		return log(log_logistic_pdf(x, (*param)[0], (*param)[1]));
	}
	
	inline double log_normal_pdf(const double x, const std::vector<double>* param){
		if (param==NULL){
			std::cerr << __FUNCTION__ << ": Missing parameters\n";
			exit(2);
		}
		return log(log_normal_pdf(x, (*param)[0], (*param)[1]));
	}
	
	inline double pareto_pdf(const double x, const std::vector<double>* param){
		if (param==NULL){
			std::cerr << __FUNCTION__ << ": Missing parameters\n";
			exit(2);
		}
		return log(pareto_pdf(x, (*param)[0], (*param)[1]));
	}
	
	inline double nakagami_pdf(const double x, const std::vector<double>* param){
		if (param==NULL){
			std::cerr << __FUNCTION__ << ": Missing parameters\n";
			exit(2);
		}
		return log(nakagami_pdf(x, (*param)[0], (*param)[1]));
	}
	
    inline double rayleigh_pdf(const double x, const std::vector<double>* param){
		if (param==NULL){
			std::cerr << __FUNCTION__ << ": Missing parameters\n";
			exit(2);
		}
		return log(rayleigh_pdf(x, (*param)[0]));
	}
	
	inline double gumbel_type_two_pdf(const double x, const std::vector<double>* param){
		if (param==NULL){
			std::cerr << __FUNCTION__ << ": Missing parameters\n";
			exit(2);
		}
		return gumbel_type_two_pdf(x, (*param)[0], (*param)[1]);
	}
	
	inline double weibull_distribution(const double x, const std::vector<double>* param){
		if (param==NULL){
			std::cerr << __FUNCTION__ << ": Missing parameters\n";
			exit(2);
		}
		return log(weibull_distribution(x, (*param)[0], (*param)[1]));
	}
	
	
	
	//Continuous Distributions, on unbounded interval
	inline double cauchy_pdf(const double x, const std::vector<double>* param){
		if (param==NULL){
			std::cerr << __FUNCTION__ << ": Missing parameters\n";
			exit(2);
		}
		return log(cauchy_pdf(x, (*param)[0], (*param)[1]));
	}
	
	inline double gumbel_pdf(const double x, const  std::vector<double>* param){
		if (param==NULL){
			std::cerr << __FUNCTION__ << ": Missing parameters\n";
			exit(2);
		}
		return log(gumbel_pdf(x, (*param)[0], (*param)[1]));
	}
	
	
	inline double generalized_normal_pdf(const double x, const std::vector<double>* param){
		if (param==NULL){
			std::cerr << __FUNCTION__ << ": Missing parameters\n";
			exit(2);
		}
		return log(generalized_normal_pdf(x, (*param)[0], (*param)[1], (*param)[2]));
	}
	
	inline double hyperbolic_secant_pdf(const double x, const std::vector<double>* param){
		return log(hyperbolic_secant_pdf(x));
	}
	
	inline double laplace_pdf(const double x, const std::vector<double>* param){
		if (param==NULL){
			std::cerr << __FUNCTION__ << ": Missing parameters\n";
			exit(2);
		}
		return log(laplace_pdf(x, (*param)[0], (*param)[1]));
	}
	
	inline double logistic_pdf(const double x, const std::vector<double>* param){
		if (param==NULL){
			std::cerr << __FUNCTION__ << ": Missing parameters\n";
			exit(2);
		}
		return log(logistic_pdf(x, (*param)[0], (*param)[1]));
	}
	
	inline double standard_normal_pdf(const double x, const std::vector<double>* param){
		return standard_normal_pdf(x);
	}
	
	inline double normal_pdf(const double x, const std::vector<double>* param){
		if (param==NULL){
			std::cerr << __FUNCTION__ << ": Missing parameters\n";
			exit(2);
		}
		return log(normal_pdf(x, (*param)[0], (*param)[1]));
	}
	
	inline double students_t_pdf(const double x, const std::vector<double>* param){
		if (param==NULL){
			std::cerr << __FUNCTION__ << ": Missing parameters\n";
			exit(2);
		}
		return log(students_t_pdf(x, (*param)[0]));
	}
	
	inline double gumbel_type_one_pdf(const double x, const std::vector<double>* param){
		if (param==NULL){
			std::cerr << __FUNCTION__ << ": Missing parameters\n";
			exit(2);
		}
		return log(gumbel_type_one_pdf(x, (*param)[0], (*param)[1]));
	}
	
	inline double generalized_extreme_value_pdf(const double x, const std::vector<double>* param){
		if (param==NULL){
			std::cerr << __FUNCTION__ << ": Missing parameters\n";
			exit(2);
		}
		return log(generalized_extreme_value_pdf(x, (*param)[0], (*param)[1], (*param)[2]));
	}
	
	inline double generalized_pareto_pdf(const double x, const std::vector<double>* param){
		if (param==NULL){
			std::cerr << __FUNCTION__ << ": Missing parameters\n";
			exit(2);
		}
		return log(generalized_pareto_pdf(x, (*param)[0], (*param)[1], (*param)[2]));
	}

   
    /* @} */ // End of PDF group
    
}

#endif /* defined(__StochHMM__PDF__) */

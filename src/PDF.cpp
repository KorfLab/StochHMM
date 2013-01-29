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
    double binomial_pdf(int k,int n, double p){
        return bin_coef(n,k)*pow(p,k)*pow(1-p,n-k);
    }
    
    
    //!BetaBinomial probability mass function
    //!<a href = "http://en.wikipedia.org/wiki/Beta-binomial_model">
    //!\param [out] dist
    //!\param n Number of trials
    //!\param a alpha
    //!\param b beta
    double beta_binomial_pdf(int n, int k, double a, double b){
        double newAlpha = (double) k + a;
        double newBeta  = (double)(n-k) + b;
        return bin_coef(n,k) * (beta(newAlpha,newBeta)/ beta(a,b));
    }
    
    //!Degenerate probability mass function
    //!param value Current value
    //!param k Degenerate value
    double degenerate_pdf(double value, double k){
        if (value == k){
            return 1.0;
        }
        else{
            return 0.0;
        }
    }
    
    
    //!Discrete Uniform CDF
	//!param position Value or Position to calculate
    //!param a Minimum position
    //!param b Maximum position
    double discrete_uniform_pdf(int position, int a, int b ){
        if (position<a || position>b){
            return 0;
        }
        else{
            return (double)1/(b-a+1);
        }
    }
    
    
    //!Hypergeometric Cumulative Distribution Function
	//!param k  Value or Number of successes
    //!param n  Number of draws from Population
    //!param N  Size of population
    //!param m  Number of successes in Population
    double hypergeometric_pdf(int k,int n, int N, int m){
        return (bin_coef(m, k)*bin_coef(N-m,n-k))/bin_coef(N, n);
    }
    
	
    //!Poisson Binomial Probability mass function
	//!Iterative calculations don't make this appropriate for using with emissions
	//!param k Value or Number of successful trials
    //!param p Probability of Success for each trial of N trials
    double poisson_binomial_pdf(int k, std::vector<double>& p){
		if (k > p.size()){
			std::cerr << "Poisson Binomial PDF Function: k > n \n";
            exit(2);
		}
		
		//Calculate K=0
		double prob(1.0);
		for(size_t i=0; i < p.size(); ++i){
			prob*=(1-p[i]);
		}
		
		
		if(k > 0){
			
			//Calculate T(i)
			std::vector<double> t(0,p.size());
			for(size_t i=0; i < p.size(); ++i){
				for (size_t j = 0; j < p.size(); ++j){
					t[i] += pow(p[j]/(1-p[j]), i);
				}
			}
			
			double temp(0);
			for(size_t i=0; i < k; ++i){
				temp+= pow(-1,i-1) * prob * t[i];
				prob = temp;
			}
			
			prob *= 1/k;
		}
		
        return prob;
    }
    
    
    /*-------------------- Discrete with infinite support ------------- */

    //!Beta negative binomial pdf
	//!param k  Value or number of failures to get N successes
    //!param n  Number of successful trials
    //!param a  Shape parameter
    //!param b  Shape parameter
    double beta_negative_binomial_pdf(int k, int n, double a, double b){
        double pmf = bin_coef(n+k-1, k);
        pmf*=gamma(a+n)*gamma(b+k)*gamma(a+b)/(gamma(a+b+n+k)*gamma(a)*gamma(b));
        return pmf;
    }
    
    
    //!Maxwell-Boltzman Probability distribution function
    //!<a href = "http://en.wikipedia.org/wiki/Maxwell–Boltzmann_distribution">
    //!param a Shape parameter  a>0
    //!param x Value
    double maxwell_boltzman_pdf(double x, double a){
        return sqrtf(2/PI)*(pow(x,2)*exp(-1*pow(x,2)/(2*pow(a, 2))))/pow(a,3);
    }
    
    
    //!Geometric probability mass function
    //!param k Value or Trial of first success
    //!param p probability of success
    double geometric_pdf(int k, double p){
        return pow(1-p,k-1) * p;
    }
    
    
    //!Logarithmic probability mass function
    //!param k	Value
    //!param p	Probability
    double logarithmic_pdf(int k, double p){
        if (k<1 || p==1.0 || p==0.0){
            std::cerr << "Logarithmic PDF Function:  k must be >=1, 0<p<1\n";
            exit(2);
        }
        return (1/log(1-p))*(pow(p,k)/k);
    }
    
    
    //!Negative binomial probability mass function
    //!param r  Value or Number of failures until experiment stopped
    //!param p  Success probability in each experiment
    double negative_binomial_pdf(int k, int r, double p){
        if (r<=0 || p<0.0 || p > 1.0){
            std::cerr << "Negative Binomial PMF: Incorrect parameters\n";
            exit(2);
        }
        return bin_coef(k+r-1, k)*pow(1-p,r)*pow(p,k);
    }
    
    
    //!Poisson probability mass function
    //!param k  Trial of first success
    //!param lambda Probability of success
    double poisson_pdf(int k, double lambda){
        if (k<0 || lambda<=0){
            std::cerr << "Poisson PMF: Incorrect parameters\n";
            exit(2);
        }
        return (pow(lambda,(double)k)* exp(-1*lambda))/factorial(k);
    }
    
    
    //!Yule-Simon probability mass function
    //!param k  
    //!param p  Shape parameter
    double yule_simon_pdf(int k, double p){
        if (k<1 || p<=0.0){
            std::cerr << "Yule-Simon PMF: Incorrect parameters\n";
            exit(2);
        }
        return p*beta(k, p+1);
    }
    
    
    //!Zipf probability mass function
    //!param k  The rank of element
    //!param N  Number of elements
    //!param s  Shape parameter (Exponent value)
    double zipf_pdf(int k,int N, double s){
        if (k<1 || N<1 || s<=0.0){
            std::cerr << "Zipf's PMF: Incorrect parameters\n";
            exit(2);
        }
        double pmf = 1/pow((double)k,s);
        double denom(0.0);
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
    double zipf_mandelbrot_pdf(int k,int N, double s, double q){
        if (k<1 || N<1 || s<=0.0 || q<0.0){
            std::cerr << "Zipf's PMF: Incorrect parameters\n";
            exit(2);
        }
        double pmf = 1/pow((double)k+q,s);
        double denom(0.0);
        for (int n = 1; n < N; n++){
            denom+=1/pow((double)n+q, s);
        }
        return pmf/denom;
    }
    
    
    
    /*-------------------- Continuous within Interval ------------------------ */
    
    //!Arcsine probability distribution function
    //!param x value 0<=X<=1
    double arcsine_pdf(double x){
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
    double beta_pdf(double x, double a, double b){
        if (x > 1.0 || x < 0.0 || a <= 0.0 || b <= 0.0){
            std::cerr << "Beta PDF: Incorrect parameters\n";
            exit(2);
        }
        return betaPDF(x, a, b);
    }
    
    
    //!Logit Normal probability distribution function
    //!param x  Value
    //!param mu Mean
    //!param sigma  Std. deviation
    double logit_normal_pdf(double x, double mu, double sigma){
        if (x >= 1.0 || x <= 0.0 || sigma <= 0.0){
            std::cerr << "Logit Normal PDF: Incorrect parameters\n";
            exit(2);
        }
        return (1/(sigma*sqrtf(2*PI)))*exp(-1*pow((logit(x)-mu),2)/(2*pow(sigma,2)))*(1/(x*(1-x)));
    }
    
    
    //!Continuous Uniform probability distribution function
	//!param x Value or Position to calculate
    //!param a Minimum position
    //!param b Maximum position
    double continuous_uniform_pdf(double x, double a, double b ){
        if (b<=a){
            std::cerr << "Continuous uniform PDF: Incorrect parameters\n";
            exit(2);
        }
        else if (x < a || x > b){
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
    double kumaraswamy_pdf(double x, double a, double b){
        if (a <= 0 || b <= 0 || x < 0 || x > 1.0 ){
            std::cerr << "Kumaraswamy PDF:  Incorrect parameters\n";
            exit(2);
        }
        return a * b * powf(x, a-1) * powf((1 - powf(x, a)),b-1);
    }
    
    //!Raised cosine probability distribution function
    //!param mu
    //!param s   S>0
    //!param x  u-s<=x<=u+s
    double raised_cosine_pdf(double x, double mu, double s){
        if (s <= 0 || x < mu-s || x > mu+s){
            std::cerr << "Raised cosine PDF:  Incorrect parameters\n";
            exit(2);
        }
        return 1/(2*s)*(1+cosf(PI*(x-mu)/s));
    }
    
    //!Triangular probability distribution function
    //!param a  Real number
    //!param b  a<b
    //!param c  a<=c<=b
    //!param x  a<=x<=b
    double triangular_pdf(double x, double a, double b, double c){
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
    double truncated_normal_pdf(double x, double mu, double sd, double a, double b){
        if (sd<0 || x < a || x > b){
            std::cerr << "Truncated normal PDF:  Incorrect paramenters\n";
            exit(2);
        }
        
        double alpha = (a-mu)/sd;
        double beta  = (b-mu)/sd;
        double xi = (x-mu)/sd;
        
        double alpha_cdf = (0.5*(1+erf((alpha-mu)/sqrt(2*pow(sd,2)))));
        double beta_cdf = (0.5*(1+erf((beta-mu)/sqrt(2*pow(sd,2)))));
        double z = beta_cdf - alpha_cdf;
        
        return (0.5*(1+erf((xi-mu)/sqrt(2*pow(sd,2)))))/(sd * z);
    }
    
    
    //!U-quadratic probability distribution function
    //!param a
    //!param b
    //!param x
    double u_quadratic_pdf(double x, double a, double b){
        if (b <= a || x < a || x > b){
            std::cerr << "U-quadratic PDF: Incorrect parameters\n";
            exit(2);
        }
        
        double alpha = 12 / powf(b-a, 3);
        double beta = (b + a)/2;
        return alpha * powf(x-beta, 2);;
    }
    
    
    //!Wigner semicircle probability distribution function
    //!param r radius R>0
    //!param x value
    double wigner_semicircle_pdf(double x, double r){
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
    double beta_prime_pdf(double x, double a, double b){
        if (a<=0 || b<=0 || x<=0){
            std::cerr << "Beta prime PDF: Incorrect parameters\n";
            exit(2);
        }
        return (powf(x, a-1)*powf(1+x, -1*a-b))/ beta(a, b);
    }
    
    
    //!Chi probability distribution function
    //!param x Value x>=0
    //!param k degrees of freedom k>0
    double chi_pdf(double x, double k){
        if (k<=0){
            std::cerr << "Chi PDF: Incorrect parameters\n";
            exit(2);
        }
        return (powf(2, 1-(k/2))*powf(x, k-1)*exp(-1*(powf(x, 2)/2)))/gamma(k/2);
    }
    
    
    //!Chi-squared probability distribution function
    //!param k K>0
    //!param x x>0
    double chi_squared_pdf(double x, double k){
        if ( k < 0 || x < 0){
            std::cerr << "Chi-squared PDF: Incorrect parameters\n";
            exit(2);
        }
        
        return (powf(x, (k/2)-1)*exp(-1*x/2))/(powf(2, k/2)*gamma(k/2));
    }
    
    
    
    //!Inverse-Chi-squared probability distribution function
    //!param x Value x>0
    //!param v v>0
    double inverse_chi_squared_pdf(double x, double v){
        if ( x <= 0 || v <= 0 ){
            std::cerr << "Inverse Chi squared PDF: Incorrect parameters\n";
            exit(2);
        }
        
        return (powf(2, -1*v/2)*powf(x, -1*(v/2)-1)*exp(-1/(2*x)))/gamma(v/2);
    }
    
    //!Scaled Inverse Chi-squared probability distribution function
    //!param x Value x>0
    //!param v  v>0
    //!param sigma_sqrd  sigma_sqrd>0
    double scaled_inverse_chi_squared_pdf(double x, double v, double sigma_sqrd){
        if (x <= 0 || v<=0 || sigma_sqrd <= 0 ){
            std::cerr << "Scaled Inverse Chi-squared PDF: Incorrect paramenters\n";
            exit(2);
        }
        double numerator = powf(sigma_sqrd*(v/2),v/2);
        numerator *= exp((-1*v*sigma_sqrd)/(2*x));
        
        double denominator = gamma(v/2);
        denominator *= powf(x, 1+(v/2));
        
        return numerator/denominator;
    }
    
    
    //!Dagum probability distribution function
    //!param p Shape parameter p>0
    //!param a Shape parameter a>0
    //!param b Shape parameter b>0
    //!param x Value x>0
    double dagum_pdf(double x, double p, double a, double b){
        if ( p <= 0 || a <= 0 || b <= 0 || x <= 0){
            std::cerr << "Dagum PDF: Incorrect parameters\n";
            exit(2);
        }
        double pdf = (a*p)/x;
        pdf *= powf(x/b, a*p)/powf(powf(x/b, a)+1 , p+1);
        return pdf;
    }
    
    
    //!Exponential probability distribution function
    //!param lambda  Rate or Inverse scale lambda>0
    //!param x  Value x>=0
    double exponential_pdf(double x, double lambda){
        if ( lambda <= 0 || x < 0 ){
            std::cerr << "Exponential PDF: Incorrect parameters\n";
            exit(2);
        }
        
        return lambda * exp (-1 * lambda * x);
    }
    
    
    //!F-Distribution
    //!param x Value x>=0
    //!param d1 Degree of freedom d1>0
    //!param d2 Degree of freedom d2>0
    double f_pdf(double x, double d1, double d2){
        if (x < 0 || d1 <= 0 || d2 <= 0){
            std::cerr << "F Distribution PDF: Incorrect parameters\n";
            exit(2);
        }
        double numerator = sqrtf((powf(d1*x, d1)*powf(d2, d2))/powf(d1*x+d2, d1+d2));
        double denominator = x * beta(d1/2, d2/2);
        return numerator / denominator;
    }
    
    
    //!Fisher's z-distribution
    //!param x Value
    //!param d1 Degree of freedom d1>0
    //!param d2 Degree of freedom d2>0
    double fishers_z_pdf(double x, double d1, double d2){
        if (d1 <= 0 || d2 <= 0){
            std::cerr << "Fisher's z-Distribution PDF: Incorrect parameters\n";
            exit(2);
        }
        double numerator = 2*powf(d1, d1/2)*powf(d2, d2/2)*exp(d1*x);
        double denominator = beta(d1/2, d2/2) * powf(d1*exp(2*x)+d2, (d1+d2)/2);
        return numerator/denominator;
    }
    
    
    //!Folded Normal probability distribution function
    //!param x Value x>=0
    //!param mu Mean(location)
    //!param sigma_sqrd Scale
    double folded_normal_pdf(double x, double mu, double sigma_sqrd){
        if (x < 0 || sigma_sqrd <= 0 ){
            std::cerr << "Folded Normal PDF: Incorrect parameters\n";
            exit(2);
        }
        double l_pdf = 1/(sqrtf(sigma_sqrd*2 * PI));
        double l_temp = exp(-1* powf(-1*x-mu,2)/(2*sigma_sqrd));
        double r_temp = exp(-1* powf(x-mu,2)/(2*sigma_sqrd));
        
        return (l_pdf*l_temp) + (l_pdf *r_temp);
    }
    
    
    //!Frechet Probability distribution function
    //!param x Value X>m
    //!param alpha Shape parameter a>0
    //!param s  Scale parameter s>0 (default s=1)
    //!param m  Location minimum (default m=0)
    double frechet_pdf(double x, double alpha, double s, double m){
        if (alpha <= 0 || s <= 0 || x < m ){
            std::cerr << "Frechet PDF: Incorrect parameters\n";
            exit(2);
        }
        
        return (alpha/s)*powf((x-m)/s, -1-alpha) * exp(-1*powf((x-m)/s, -1*alpha));
    }
    
	//!Gamma probability distribution
	//!param x Value x>0
	//!param alpha Shape parameter a>0
	//!param beta Rate parameter b>0
	//!http://en.wikipedia.org/wiki/Gamma_distribution
	double gamma_pdf(double x, double alpha, double beta){
		if (x <= 0 ||  alpha <= 0 || beta <= 0){
			std::cerr << "Gamma PDF: Incorrect parameters\n";
            exit(2);
		}
		
		return (pow(beta,alpha)/gamma(alpha))*pow(x,alpha-1)*exp(-beta*x);
	}
	
	//!Inverse Gamma probability distribution
	//!param x Value x>0
	//!param alpha Shape parameter x>0
	//!param beta Scale parameter b>0
	double inv_gamma_pdf(double x, double alpha, double beta){
		if (x <= 0 ||  alpha <= 0 || beta <= 0){
			std::cerr << "Inverse Gamma PDF: Incorrect parameters\n";
            exit(2);
		}
		return (pow(beta,alpha)/gamma(alpha))*pow(x,-1*alpha-1)*exp(-beta/x);
	}
	
	//!Half Normal probability distribution
	//!param x Value x>0
	//!param sigma Standard Deviation sigma>0
	double half_normal_pdf(double x, double sigma){
		if (x <= 0 || sigma <= 0){
			std::cerr << "Half Normal PDF: Incorrect parameters\n";
            exit(2);
		}
		return sqrt(2)/(sigma*sqrt(PI)) * exp(-1* pow(x,2)/(2*pow(sigma,2)));
	}
	
	//!Inverse Gaussian probability distribution
	//!param x Value x>0
	//!param mu Average u>0
	//!param lambda Shape parameter l>0
	double inv_gaussian_pdf(double x, double mu, double lambda){
		if (x <= 0 || mu <= 0 || lambda <= 0){
			std::cerr << "Inverse Gaussian PDF: Incorrect parameters\n";
            exit(2);
		}
		return sqrt(lambda/(2*PI*pow(x, 3))) * exp((lambda*pow(x-mu,2))/(2*pow(mu,2)*x));
	}
	
	//!Levy probability distribution function
	//!param x Value  x>= mu and x<INFINITY
	//!param mu location parameter
	//!param scale scale parameter C>0
	double levy_pdf(double x, double mu, double scale){
		if (x < mu || x == INFINITY || scale <= 0){
			std::cerr << "Levy PDF: Incorrect parameters\n";
            exit(2);
		}
		return sqrt(scale/(2*PI)) * (exp((-1*scale)/(2*(x-mu)))/pow(x-mu,1.5));
	}
	
	//!Log Cauchy probability distribution function
	//!param x Value X>0 X<INFINITY;
	//!param mu Location
	//!param sigma scale parameter
	double log_cauchy_pdf(double x,double mu, double sigma){
		if (x <= 0 || sigma <= 0 || x== INFINITY){
			std::cerr << "Log Cauchy PDF: Incorrect parameters\n";
            exit(2);
		}
		return (1/(x*PI))*(sigma/(pow(log(x-mu),2)+pow(sigma,2)));
	}
	
	//!Log Laplace probability distribution function
	//!param x Value
	//!param mu parameter
	//!param b  parameter
	double log_laplace_pdf(double x, double mu, double b){
		if (x < mu){
			return (1/(b*x))*exp(-1*((mu-log(x))/b));
		}
		
		return (1/(b*x))*exp(-1*((log(x)-mu)/b));
	}
	
	//!Log logistic probability distribution function
	//!param x Value
	//!param alpha scale parameter
	//!param beta  shape parameter
	double log_logistic_pdf(double x, double a, double b){
		if (x<0 || a < 0 || b < 0 || x == INFINITY){
			std::cerr << "Log logistic PDF: Incorrect parameters\n";
			exit(2);
		}
		return ((b/a)*pow(x/a,b-1))/pow(1+pow(x/a,b),2);
	}
	
	
	//!Log Normal probability distribution function
	//!param x Value x>0 and x<INFINITY
	//!param mu Log scaled location parameter
	//!param sigma_sqrd Log scaled scaling parameter
	double log_normal_pdf(double x, double mu, double sigma_sqrd){
		if (x < 0 || x == INFINITY || sigma_sqrd <= 0){
			std::cerr << "Log Normal PDF: Incorrect parameters\n";
			exit(2);
		}
		return (1/(x*sqrt(2*PI*sigma_sqrd)))*exp(-1* (pow(log(x)-mu,2))/(2*sigma_sqrd));
	}
	
	
	//!Pareto probability distribution function
	//!param x Value
	//!param alpha shape parameter
	//!param x_m scale parameter
	double pareto_pdf(double x, double alpha, double x_m){
		if ( alpha <= 0 || x_m <= 0){
			std::cerr << "Pareto PDF: Incorrect parameters\n";
			exit(2);
		}
		
		if (x<x_m){
			return 1;
		}
		return pow(x_m/x, alpha);
	}
	
	//!Nakagami probability distribution function
	//!param x	Value
	//!param mu	shape parameter
	//!param w	spread parameter
	double nakagami_pdf(double x, double mu, double w){
		if (x<=0 || mu < 0.5 || w <=0){
			std::cerr << "Nakagami PDF: Incorrect parameters\n";
			exit(2);
		}
		
		return (pow(2*mu,2)/(gamma(mu)*pow(w, mu)))*pow(x,2*mu-1)*exp(-1*(mu/w)*pow(x,2));
	}
	
	
	//!Rayleigh probability distribution function
	//!param x Value x>=0 and x<INFINITY
	//!param sigma Mode
	double rayleigh_pdf(double x, double sigma){
		if (x<0 || x == INFINITY || sigma <= 0){
			std::cerr << "Rayleigh PDF: Incorrect parameters\n";
			exit(2);
		}
		return (x/pow(sigma, 2))*exp(-1*(pow(x,2))/(2*pow(sigma,2)));
	}
	
	//!Type 2 Gumbell Probability distribution function
	//!param x Value
	//!param a	parameter
	//!param b	shape parameter
	double gumbel_type_two_pdf(double x, double a, double b){
		if (x<=0 || x== INFINITY){
			std::cerr << "Type 2 Gumbel PDF: Incorrect parameters\n";
			exit(2);
		}
		return a*b*pow(x,-1*a-1) * exp(-1*b*pow(x,-1*a));
	}
	
	//!Weibull Probability distribution function
	//!param x	value
	//!param lambda	scale parameter
	//!param k	shape parameter
	double weibull_distribution(double x, double lambda, double k){
		if ( x==INFINITY || lambda <= 0 || k <= 0){
			std::cerr << "Weibull PDF: Incorrect parameters\n";
			exit(2);
		}
		
		if (x <0){
			return 0.0;
		}
		
		return (k/lambda)*pow(x/lambda,k-1)*exp(-1*(pow(x/lambda, k)));
	}
	
	
	
	
	/*-------------- Continuous with infinite Interval ---------------- */

	
	//!Cauchy probability distribution function
	//!param x	Value
	//!param x_o	location parameter (place of peak)
	//!param gamma	scale parameter gamma>0
	double cauchy_pdf(double x, double x_o, double gamma){
		if (gamma <= 0){
			std::cerr << "Cauchy PDF: Incorrect parameters\n";
			exit(2);
		}
		
		return (1/PI)* (gamma)/(pow(x-x_o, 2)+pow(gamma,2));
	}
	

	
	//!Gumbel Probability distribution function
	//!param x Value
	//!param mu	Location parameter
	//!param beta	Scale parameter
	double gumbel_pdf(double x, double mu, double beta){
		if (beta <= 0){
			std::cerr << "Gumbel PDF: Incorrect parameters\n";
			exit(2);
		}
		double z = (x-mu)/beta;
		return (1/beta)*exp(-1*(z)-exp(-1*z));
	}
	
//	//!Fisher's Z Probability distribution function
//	//!param x	Value
//	//!param d1	Degrees of freedom
//	//!param d2 Degrees of freedom
//	double fishers_z_pdf(double x, double d1, double d2){
//		if (d1 <= 0 || d2 <=0){
//			std::cerr << "Fisher's Z PDF: Incorrect parameters\n";
//			exit(2);
//		}
//		
//		return ((2*pow(d1,d1/2)*pow(d2,d2/2))/ beta(d1/2, d2/2)) * (exp(d1*x)/pow(d1*exp(2*x)+d2, (d1+d2)/2));
//	}
	
	//!Generalized Normal probability distribution function
	//!param x	Value
	//!param mu	Location parameter
	//!param alpha	scale parameter alpha>0
	//!param beta	shape parameter beta>0
	double generalized_normal_pdf(double x, double mu, double alpha, double beta){
		if ( alpha < 0 || beta < 0){
			std::cerr << "Generalized normal PDF: Incorrect parameters\n";
			exit(2);
		}
		return beta/(2*alpha*gamma(1/beta))*exp(-1*pow(abs(x-mu)/alpha, beta));
	}
	
	//!Hyperbolic secant probability distribution function
	//!param	x	Value
	double hyperbolic_secant_pdf(double x){
		return 0.5*pow(cosh((PI/2)*x),-1);
	}
	
	
	//!Laplace Probability distribution function
	//!param	x	Value
	//!param	mu	Location Parameter
	//!param	b	Shape parameter
	double laplace_pdf(double x, double mu, double b){
		return 1/(2*b) * exp( -1 * abs(x-mu)/b);
	}
	
	//!Logistic Probability distribution function
	//!param	x	Value
	//!param	mu	Location parameter
	//!param	s	Shape parameter
	double logistic_pdf(double x, double mu, double s){
		return (1/(4*s))*pow(cosh((x-mu)/s),-2);
	}
	
	//!Standard Normal probability distribution function
	//!param	x	Value
	double standard_normal_pdf(double x){
		return normal_pdf(x, 0, 1);
	}
	
	//!Normal probability distribution function
	//!param	x	Value
	//!param	mu	Location parameter(mean)
	//!param	sigma Scaling parameter (standard deviation)
	double normal_pdf(double x, double mu ,double sigma){
		if (sigma <= 0){
			std::cerr << "Normal PDF: Incorrect parameters\n";
			exit(2);
		}
		
		return (1/(sigma*sqrt(2*PI)))* exp(-1*pow(x-mu, 2)/(2*pow(sigma, 2)));
	}
	
	//!Student's t-probability distribution function
	//!param	x	Value
	//!param	v	Degrees of freedom
	double students_t_pdf(double x, double v){
		if (v<=0){
			std::cerr << "Student's t-PDF: Incorrect parameters\n";
			exit(2);
		}
		double i = gamma((v+1)/2)/(sqrt(v*PI)*gamma(v/2));
		return i * pow(1+(pow(x,2)/v),-1*(v+1)/2);
	}
	
	
	//!Gumbel Type-2 Probability distribution function
	//!param	x	Value
	//!param	a
	//!param	b	Shape
	double gumbel_type_one_pdf(double x, double a, double b){
		return a*b*exp(-1*(b*exp(-1*a*x))+a*x);
	}
	
	
	//!Generalized extreme value probability distribution function
	//!param x Value;
	//!param mu	location parameter
	//!param sigma	scale parameter
	//!param xi		shape parameter
	double generalized_extreme_value_pdf(double x, double mu, double sigma, double xi){
		if (sigma <= 0){
			std::cerr << "Generalized extreme value PDF: Incorrect parameters\n";
			exit(2);
		}
		
		if (xi > 0 && (x < mu - (sigma/xi) || x == INFINITY)){
			std::cerr << "Generalized extreme value PDF: Incorrect parameters\n";
			exit(2);
		}
		else if (xi == 0 && (x == -INFINITY || x == INFINITY)){
			std::cerr << "Generalized extreme value PDF: Incorrect parameters\n";
			exit(2);
		}
		else if (xi > 0 && (x == -INFINITY || x > mu - (sigma/xi))){
			std::cerr << "Generalized extreme value PDF: Incorrect parameters\n";
			exit(2);
		}
		
		double t = (xi == 0 ) ?  exp(-1*((x-mu)/sigma)) :  pow(1+((x-mu)/sigma)*xi, -1/xi);
		
		return (1/sigma)*pow(t,xi+1)*exp(-1*t);
	}
	
	
	//!Generalized Pareto probability distribution function
	//!param	x	Value
	//!param	mu	Location parameter
	//!param	sigma	Scale parameteer
	//!param	xi	Shape parameter
	double generalized_pareto_pdf(double x, double mu, double sigma, double xi){
		if (xi >= 0 && x < mu){
			std::cerr << "Generalized pareto PDF: Incorrect parameters\n";
			exit(2);
		}
		
		if (xi < 0 && ( x < mu || x > (mu-sigma/xi))){
			std::cerr << "Generalized pareto PDF: Incorrect parameters\n";
			exit(2);
		}
		
		if (xi ==0){
			return 1 - exp(-1*((x-mu)/sigma));
		}
		
		return 1 - pow(1+((xi*(x-mu))/sigma), -1/xi);
	}
	
	
	
	//!Dirichlet Distribution
	//!param	x	Vector of values
	//!param	alpha	concentration parameters where a_i>0
	double dirichlet_pdf(std::vector<double>& x, std::vector<double>& alpha){
		double gsum_alpha = gamma(sumVector(alpha));
		
		//Calc product of gamma transformed alpha values
		double b_n(gamma(alpha[0]));
		for (size_t i=1; i < alpha.size(); ++i){
			b_n*=gamma(alpha[i]);
		}
		
		double b = 1 / (b_n/gsum_alpha);  //Calculate B(alpha)
		
		//Calculate final value by finding product of xi^(ai-1)
		for (size_t i=0; i<alpha.size(); ++i){
			b*=pow(x[i],alpha[i]-1);
		}
		return b;
	}
	
	
	//!Multivariate Ewen's probability distribution function
	//!param	x	# of individual allele represented in sample
	//!param	theta	Parameter: details of evolutionary model
	double multivariate_ewens_pdf(std::vector<double>& x, double theta){
		
		double n = factorial(x.size());
		double d(theta);
		for (size_t i=0; i < x.size() -1 ; ++i){
			d *= theta + i;
		}
		
		double prob=1;
		for(size_t i=0; i < x.size() ; ++i){
			prob*= pow(theta,x[i]) / (pow(i,x[i]) * factorial(x[i]));
		}
		return (n/d) * prob;
	}
	
	
}
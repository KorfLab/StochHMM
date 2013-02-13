//statistics.cpp
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

#include "statistics.h"
namespace StochHMM{
    
    template <class T>
    T min(std::vector<T> &set){
        T min =set[0];
        for(long i=set.size()-1;i>0;i--){
            if (set[i]<min){ min=set[i];}
        }
        return min;
    }
    
    template <class T>
    T max(std::vector<T> &set){
        T max=set[0];
        for(long i=set.size()-1;i>0;i--){
            if (set[i]>max){ max=set[i];}
        }
        return min;
    }
    
    template <class T>
    T construct_histogram (std::vector<T> &set,int N_bins){
        T mini=min(set);
        T maxi=max(set);
        T delta=(maxi-mini+1)/N_bins;
        std::vector<T> bin (N_bins,0);
        for(int i=set.size()-1;i>=0;i--){
            T value=floor((set[i]-mini)/delta);
            bin(value)=bin(value)+1/(set.size()*delta);
        }
        return bin;
    }
    
    template <class T>
    T smooth_histogram (std::vector<T> histo, int intervals, int window_size, int iterations){
        std::vector<T> s_histo=histo;
        for (int i=1;i<=iterations;i++){
            for (int b=0;b<=intervals-window_size;b++){
                int c=b+floor((window_size-1)/2);
                s_histo[c]=0;
                for (int j=b;j<=b+window_size-1;i++){
                    s_histo[c]=s_histo[c]+histo[j]/window_size;
                }
            }
            for (int b=0;b<=((window_size-1)/2)-1;b++){
                s_histo=s_histo[floor((window_size-1)/2)];
            }
            for (int b=intervals-window_size+1+floor((window_size-1)/2); b<=intervals-1;i++){
                s_histo[b]=s_histo[intervals-window_size+floor((window_size-1) / 2)];
            }
            histo=s_histo;
        }
        T sum=0;
        for (int b=0;b<=intervals-1;b++){
            sum+=sum+s_histo[b];
        }
        for (int b=0;b<=intervals-1;b++){
            s_histo[b]/=sum;
        }
        return s_histo;
    }
    
    float entropy (std::vector<float> &set){
        float entropy=0;
        for(size_t i=0;i<set.size();i++){
            entropy+=set[i]*(log(set[i])/log(2));
        }
        return entropy*-1;
    }
    
    float rel_entropy (std::vector<float> &set, std::vector<float> &set2){
        float rel_entropy=0;
        if (set.size()!=set2.size()){
            std::cerr << "Distributions aren't the same size";
        }
        
        for(size_t i=0;i<set.size();i++){
            rel_entropy+=0.5*(set[i]* (log (set[i]/set2[i]) /log(2) )+set2[i]*(log(set2[i]/set[i])/log(2)));
        }
        return rel_entropy;
    }
    
    
    ///////////////////////////////////////////  Integration & Summation //////////////////////////////////////////////
    //Need to adapt to take up to three paramenter and pass them to the function.
    //If one parameter is given then only pass one to function and so on.
    
    //Fix: variable handed to function to include only std::vector<double>
    double _integrate (double (*funct)(double, std::vector<double>),double upper, double lower,std::vector<double> param){
        double mid=(lower+upper)/2.0;
        double h3=fabs(upper-lower)/6.0;
        return h3*(funct(lower,param)+4*funct(mid,param)+funct(upper,param));
    }
    
    double integrate (double (*funct)(double,std::vector<double>),double lower, double upper ,std::vector<double> param, double max_error=0.001, double sum=0){
        double mid=(upper+lower)/2.0;
        double left=_integrate(funct,lower, mid, param);
        double right=_integrate(funct,mid,upper, param);
        if (fabs(left+right-sum)<=15*max_error){
            return left+right +(left+right-sum)/15;
        }
        return integrate(funct,lower,mid,param,max_error/2,left) + integrate(funct,mid,upper,param,max_error/2,right);
    }
    
    
    
    //  Adaptive Simpson's numerical integration method
    double simpson (double (*funct)(double,double,double),double alpha, double beta, double lower, double upper){
        double mid=(lower+upper)/2.0;
        double h3=fabs(upper-lower)/6.0;
        return h3*(funct(lower,alpha,beta)+4*funct(mid,alpha,beta)+funct(upper,alpha,beta));
    }
    
    double adapt_simpson (double (*funct)(double,double,double),double alpha, double beta, double lower, double upper, double max_error, double sum){
        double mid=(upper+lower)/2.0;
        double left=simpson(funct,alpha,beta,lower, mid);
        
        double right=simpson(funct,alpha,beta,mid,upper);
        if (fabs(left+right-sum)<=15*max_error){
            return left+right +(left+right-sum)/15;
        }
        return adapt_simpson(funct,alpha,beta,lower,mid,max_error/2,left) + adapt_simpson(funct,alpha,beta,mid,upper,max_error/2,right);
    }
    
    
    /////////////////////////////////////// Distributions //////////////////////////////////////////////
    
    
    /* ---------------- Discrete Distributions ---------------- */
    /* ------------------- FINITE INTERVALS ------------------- */
    
    //!Discrete Uniform CDF
    //!param a Minimum position
    //!param b Maximum position
    //!param position Position to calculate
    float discrete_uniform_cdf(int a, int b, int position){
        if (position<a){
            return 0;
        }
        else if (position>=b){
            return 1;
        }
        else{
            return (float)(position-a+1)/(b-a+1);
        }
    }
    
    
    
    //!Binomial coefficient
    //!number of combinations possible given k unordered outcomes from n possibilities
    //!param n number of items
    //!param k number in set
    float bin_coef (int n, int k){
        float c=0;
        for(int i=k+1;i<=n;i++) {c+=log(i);}
        for(int j=1;j<=(n-k);j++) {c-=log(j);}
        return exp(c);
    }
    
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
    
    
    //!Binomial Cumulative Distribution Function
    //!param k Number of successes
    //!param n Number of trials
    //!param p Probability of success
    //!If we can get an incomplete beta function. We can update this function and
    //calculate without using the summation.
    float binomial_cdf(int k,int n, float p){
        std::vector<float> param (2,0.0);
        param[0]=n;
        param[1]=p;
        return summation(binomial_pdf,0,k,param);
    }
    
    
    
    //!Hypergeometric Cumulative Distribution Function
    //!param n  Number of draws from Population
    //!param N  Size of population
    //!param m  Number of successes in Population
    //!param k  Number of successes 
    double hypergeometric_cdf(int n, int N, int m, int k){
        if (n==N){
            return 1.0;
        }
        else{
            double prob=0;
            for(int i=1; i<=n;i++){
                prob+=(bin_coef(m, k)*bin_coef(N-m,i-k))/bin_coef(N, i);
            }
            return prob;
        }
    }
    
    
    
    //Wikipedia Windschitl approximation to gamma function
    double gamma_func (double x){
        if (x>0 && x<=0.5){
            return sin(PI*x)*gamma_func(1-x); }
        else {
            return sqrt((2*PI)/x)*pow(((x/exp(1))*sqrt(x*sinh(1/x)+1/810*pow(x,6))),x);
        }
    }
    
    double gamma_pdf (double x, double alpha, double beta) {
        if (x>0){
            return (pow(beta,alpha)/gamma_func(alpha))*pow(x,alpha-1)*exp(-beta*x);
        }
        else {
            return 0;
        }
    }
    
    double gamma_pdf (double x, std::vector<double> param){   //double x, double alpha, double beta) {
        return gamma_pdf(x,param[0],param[1]);
    }
    
    double gamma_cdf (double x, double alpha, double beta){
        std::vector<double> parameters (2,0.0);
        parameters[0]=alpha;
        parameters[1]=beta;
        return integrate(gamma_pdf,0.0,x,parameters);
    }
    
    double chi2_pdf(double x, double df){
        return gamma_pdf(x,df/2.0,0.5);
    }
    
    double chi2_cdf (double x, double df){
        return gamma_cdf(x,df/2.0,0.5);
    }
    
    double beta_pdf(double x, double alpha, double beta){
        return (gamma_func(alpha+beta)/(gamma_func(alpha)*gamma_func(beta)))*pow(x,alpha-1)*pow(1-x,beta-1);
    }
    
    double beta_pdf(double x, std::vector<double> param){
        return beta_pdf(x,param[0],param[1]);
    }
    
    double beta_cdf(double x, double alpha, double beta){
        std::vector<double> parameters (2,0.0);
        parameters[0]=alpha;
        parameters[1]=beta;
        return integrate(beta_pdf, 0.0, x,parameters);
    }
    
    double expon_pdf(double x, double beta){
        return gamma_pdf(x,1.0,beta);
    }
    
    double expon_cdf(double x, double beta){
        return gamma_cdf(x,1.0,beta);
    }
    
    double normal_pdf(double x, double mean, double variance){
        return (1/(sqrt(2*PI*variance)))*exp(-1*pow(x-mean,2)/(2*variance));
    }
    
    
    
   
    
    
    
    float summation (float (*funct)(int,std::vector<float>), int lower, int upper, std::vector<float> param){
        float sum=0;
        for(int i=lower;i<=upper;i++){
            sum+=funct(i,param);
        }
        return sum;
    }
    
    
    
   
    
    
    double multinomial_pdf(std::vector<int> r, int n, std::vector<double> p){
        double prob=1.0;
        double denom=1.0;
        if (r.size()==p.size()){return 0.0;}
        else {
            
            for(int i=0;i<r.size();i++){
                prob*=p[i];
                denom*=gamma_func(r[i]);
            }
        }
        
        return gamma_func((n+1)/denom)*prob;   ///check to make sure this is right;
    }
    
    //////////////////////// Need to fix everything below this point ///////////////////////////////////
    double _low_igamma (double x, std::vector<double> param){
        double value=param[0];
        return (pow(x,value-1)*exp(-x));
    }
    
    double low_igamma (double x, double alpha){
        std::vector<double> parameters (2,0.0);
        parameters[0]=alpha;
        return (integrate(_low_igamma,0.0,x,parameters));
    }
    
    double upper_igamma (double x , double alpha){
        return (1/gamma_func(alpha))-low_igamma(x,alpha);
    }
    
    double erf (double x) {
        return 1-(upper_igamma(0.5,pow(x,2))/sqrt(2));
    }
    
    double std_normal_cdf (double x){
        x/=sqrt(2);
        double erf_value=erf(x);
        return 0.5*(1.0+erf_value);
    }
}

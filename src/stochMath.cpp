//
//  stochMath.cpp

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

#include "stochMath.h"
namespace StochHMM{

    
    void generateUnknownIndices(std::vector<int>& results, int alphabetSize, int order ,int value){
        results.assign(alphabetSize,0);
        for (int i=0;i<alphabetSize;i++){
            results[i]= value + i * POWER[order-1][alphabetSize-1];
        }
        return;
    }

    
    //Linear interpolation of y value given two pair<x,y> and x value
    double interpolate(std::pair<double,double>& a, std::pair<double,double>& b, double& cx){
        //std::cout << a.first << "\t" << a.second <<std::endl;
        //std::cout << b.first << "\t" << b.second <<std::endl;
        
        return a.second+(b.second-a.second)*((cx-a.first)/(b.first-a.first));
    }
    
    //Linear extrapolation of y value given two pair<x,y> and x value
    double extrapolate(std::pair<double,double>& a, std::pair<double,double>& b, double& cx){
        return a.second+((cx-a.first)/(b.first-a.first))*(b.second-a.second);
    }
    
    //Shannon's entropy
    float entropy (std::vector<float> &set){
        float entropy=0;
        for(size_t i=0;i<set.size();i++){
            entropy+=set[i]*(log(set[i])/log(2));
        }
        return entropy*-1;
    }
	
	//Shannon's entropy
	double entropy (std::vector<double> &set){
        double entropy=0;
        for(size_t i=0;i<set.size();i++){
            entropy+=set[i]*(log(set[i])/log(2));
        }
        return entropy*-1;
    }
    
	//Shannon's Relative Entropy (Kullback-Liebler Convergence)
	//Normalized for A->B and B->A
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
	
	//Shannon's Relative Entropy (Kullback-Liebler Convergence)
	//Normalized for A->B and B->A
	double rel_entropy (std::vector<double> &set, std::vector<double> &set2){
        double rel_entropy=0;
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
    
    //Fix: variable handed to function to include only vector<double>
    double _integrate (double (*funct)(double, std::vector<double>&),double upper, double lower,std::vector<double>& param){
        double mid=(lower+upper)/2.0;
        double h3=fabs(upper-lower)/6.0;
        return h3*(funct(lower,param)+4*funct(mid,param)+funct(upper,param));
    }
    
    double integrate (double (*funct)(double,std::vector<double>&),double lower, double upper ,std::vector<double>& param, double max_error=0.000001, double sum=0){
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
    
    double summation (double (*funct)(int,std::vector<double>), int lower, int upper, std::vector<double> param){
        double sum=0;
        for(int i=lower;i<=upper;i++){
            sum+=funct(i,param);
        }
        return sum;
    }
    
    /////////////////////////////////////// FUNCTIONS //////////////////////////////////////////////
    
    double igamma_upper (double s, double x){
        return tgamma(s)-igamma_lower(s,x);
    }
    
    
    //incomplete lower gamma (a,x)
    ///http://wwwC/.rskey.org/CMS/index.php/the-library/11
    double igamma_lower (double a, double x){
        double constant = pow(x,a) * exp(-1 * x);
        double sum =0;
        for (int n=0;n<60;n++){
            double num = pow(x,(double)n);
            double denom=a;
            for(int m=1;m<=n;m++){
                denom*=a+m;
            }
            num/=denom;
            sum+=num;
        }
        return constant * sum;
    }
    
    
    double rgammap(double s, double x){
        return igamma_lower(s,x)/tgamma(s);
    }
    
    
    double rgammaq(double s, double x){
        return igamma_upper(s,x)/tgamma(s);
    }
    
    //Beta(a,b) = Beta function
    double beta(double a, double b){
        double value=(tgamma(a)*tgamma(b))/tgamma(a+b);
        return value;
    }
    
    //Incomplete beta function B(x,a,b)
    double ibeta(double x, double a, double b){
        return betaPDF(x,a,b) * exp(log(tgamma(a)) + log(tgamma(b)) - log(tgamma(a+b)));
    }
    
    
    //!Beta probability distribution function
    //!param x  Value 0<x<1
    //!param a  Shape parameter a>0
    //!param b  Shape parameter b>0
    float betaPDF(float x, float a, float b){
        float constant = 1/beta(a,b);
        constant*=pow(x,a-1) * pow(1-x,b-1);
        return constant;
    }
    
    
    //
    ////Incomplete beta function B(x,a,b)
    //double ibeta(double x, double a, double b){
    //	vector<double> parameters(2,0.0);
    //	parameters[0]=a;
    //	parameters[1]=b;
    //	return (integrate(_ibeta,0.0,x,parameters));
    //}
    //
    //double _ibeta(double x, vector<double>& parameters){
    //	double a=parameters[0];
    //	double b=parameters[1];
    //	return (pow(x,a-1)*pow(1-x,b-1));
    //}
    
    
    //Regularized incomplete beta  Ix(a,b)
    double ribeta(double x, double a, double b){
        return ibeta(x,a,b)/beta(a,b);
    }
    
    //Returns stirlings approximation of floor(X!)
    double factorial(double x){
        double f=sqrt((2*x+(1/3))*M_PI)*pow(x,x)*exp(-1*x);
        return floor(f);
    }
    
    //Calculate the binomial coefficient
    double bin_coef (double n, double r){
        double c=0;
        for(int i=r+1;i<=n;i++) {c+=log(i);}
        for(int j=1;j<=(n-r);j++) {c-=log(j);}
        return exp(c);
    }
    
    int bin_coef (int n, int r){
        float c=0;
        for(int i=r+1;i<=n;i++) {c+=log(i);}
        for(int j=1;j<=(n-r);j++) {c-=log(j);}
        return exp(c);
    }
    

}

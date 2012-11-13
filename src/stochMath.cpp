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

    
    double addLog(double& first, double& second){
        //cout << first <<endl;
        //cout << second <<endl;
        if (first==-INFINITY){
            return second;
        }
        else if (second==-INFINITY){
            return first;
        }
        else if (first>second){
            return first+log(1+exp(second-first));
        }
        else{
            return second+log(1+exp(first-second));
        }
    }
    
    //----------------------------------------------------------------------------//
    // Description: addVectorCombinatorial                                                              
    // Adds the values of vectors combinatorial
    // Example [ 1 4 ] + [ 2 7 ] = [ 2 6 8 11 ] 
    // 
    //----------------------------------------------------------------------------//
    void addVectorCombinatorial(std::vector<int>& result, std::vector<int>& lhs, std::vector<int>& rhs){
        
        if (result.size()>0){
            result.clear();
        }
        
        //If either vector is empty then return copy of the other
        if (lhs.size()==0){
            result.assign(rhs.begin(),rhs.end());
            return;
        }
        else if (rhs.size()==0){
            result.assign(lhs.begin(),lhs.end());
            return;
        }
        
        for(size_t i=0;i<lhs.size();i++){
            for(size_t j=0;j<rhs.size();j++){
                result.push_back(lhs[i]+rhs[j]);
            }
        }
        return;
    }
    
    
    void addVectorCombinatorial(std::vector<double>& result, std::vector<double>& lhs, std::vector<double>& rhs){
        if (result.size()>0){
            result.clear();  //Make sure result is empty
        }
        
        //If either vector is empty then return copy of the other
        if (lhs.size()==0){
            result.assign(rhs.begin(),rhs.end());
            return;
        }
        else if (rhs.size()==0){
            result.assign(lhs.begin(),lhs.end());
            return;
        }
        
        for(size_t i=0;i<lhs.size();i++){
            for(size_t j=0;j<rhs.size();j++){
                result.push_back(lhs[i]+rhs[j]);
            }
        }
        return;
    }
    
    
    void multiplyVectorCombinatorial(std::vector<double>& result, std::vector<double>&lhs, std::vector<double>&rhs){
        if (result.size()>0){
            result.clear();
        }
        
        //If either vector is empty then return copy of the other
        if (lhs.size()==0){
            result.assign(rhs.begin(),rhs.end());
            return;
        }
        else if (rhs.size()==0){
            result.assign(lhs.begin(),lhs.end());
            return;
        }
        
        for(size_t i=0;i<lhs.size();i++){
            for(size_t j=0;j<rhs.size();j++){
                result.push_back(lhs[i]*rhs[j]);
            }
        }
        return;
    }
    
    void multiplyVectorCombinatorial(std::vector<int>& result, std::vector<int>&lhs, std::vector<int>&rhs){
        if (result.size()>0){
            result.clear();
        }
        
        //If either vector is empty then return copy of the other
        if (lhs.size()==0){
            result.assign(rhs.begin(),rhs.end());
            return;
        }
        else if (rhs.size()==0){
            result.assign(lhs.begin(),lhs.end());
            return;
        }
        
        
        for(size_t i=0;i<lhs.size();i++){
            for(size_t j=0;j<rhs.size();j++){
                result.push_back(lhs[i]*rhs[j]);
            }
        }
        return;
    }
    
    
    void addValueToVector(std::vector<int>& vec, int value){
        for(size_t i=0;i<vec.size();i++){
            vec[i]+=value;
        }
        return;
    }
    
    void addValueToVector(std::vector<double>& vec, double value){
        for(size_t i=0;i<vec.size();i++){
            vec[i]+=value;
        }
        return;
    }
    
    
    void multiplyValueToVector(std::vector<double>& vec, double value){
        for(size_t i=0;i<vec.size();i++){
            vec[i]*=value;
        }
        return;
    }
    
    void multiplyValueToVector(std::vector<int>& vec, int value){
        for(size_t i=0;i<vec.size();i++){
            vec[i]*=value;
        }
        return;
    }
    
    
    void divideValueToVector(std::vector<int>& vec, int value){
        for(size_t i=0;i<vec.size();i++){
            vec[i]/=value;
        }
        return;
    }
    
    
    void divideValueToVector(std::vector<double>& vec, double value){
        for(size_t i=0;i<vec.size();i++){
            vec[i]/=value;
        }
        return;
    }
    
    
    
    
    void generateUnknownIndices(std::vector<int>& results, int alphabetSize, int order ,int value){
        results.assign(alphabetSize,0);
        for (int i=0;i<alphabetSize;i++){
            results[i]= value + i * POWER[order-1][alphabetSize-1];
        }
        return;
    }
    
    
    
    double sumVector(std::vector<double>& data){
        double sum=0;
        for(size_t i=0;i<data.size();i++){
            sum+=data[i];
        }
        return sum;
    }
    
    
    
    double minVector(std::vector<double>& data){
        return *min_element(data.begin(), data.end()); 
    }
    
    
    double maxVector(std::vector<double>& data){
        return *max_element(data.begin(), data.end());
    }
    
    
    double avgVector(std::vector<double>& data){
        return sumVector(data) / double(data.size());
    }
    
    void logVector(std::vector<double>& data){
        for(size_t i=0;i<data.size();i++){
            data[i]=log(data[i]);
        }
        return;
    }
    
    void expVector(std::vector<double>& data){
        for(size_t i=0;i<data.size();i++){
            data[i]=exp(data[i]);
        }
        return;
    }
    
    void probVector(std::vector<double>& data){
        double sum=sumVector(data);
        for(size_t iter=0;iter<data.size();iter++){
            if (sum==0){
                data[iter]=0;
            }
            else{
                data[iter]/=sum;
            }
        }
        return;
    }
    
	
	void logProbVector(std::vector<double>& data){
        double sum=sumVector(data);
        for(size_t iter=0;iter<data.size();iter++){
            if (sum==0){
                data[iter]=-INFINITY;
            }
            else{
                data[iter] = log(data[iter]/sum);
            }
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
		for(size_t i=set.size()-1;i>=0;i--){
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

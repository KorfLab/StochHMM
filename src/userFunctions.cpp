//
//  UserFunctions.cpp
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

#include "userFunctions.h"
namespace StochHMM{
    
	
	StateFuncs::StateFuncs(){
		_loadUnivariatePdf();
	}
    
    //!Assign a transition function to the StateFuncs class
    //!\param str Name of function
    //!\param ptrFunc  pt2StateFunc to use for StateFunc
    void StateFuncs::assignTransitionFunction(std::string& str, transitionFunc ptrFunc){
     
         if (transitionFunctions.count(str)==0){
             transitionFunctions[str]=ptrFunc;
         }
         else{
             std::cerr << "Function Name: " << str << " already exists.   You need to choose a new function name that doesn't already exist. For reference here are a list of names already assigned as external functions\nAssigned Names:\n";
         
             std::map<std::string,transitionFunc>::iterator it;
             
             for(it=transitionFunctions.begin();it!=transitionFunctions.end();it++){
                 std::cerr << "\t" << it->first <<std::endl;
             }
         }
     };
	
	
	void StateFuncs::assignTransitionFunction(const char* str, transitionFunc ptrFunc){
		
		std::string st(str);
		assignTransitionFunction(st, ptrFunc);
	};
    
    
    //!Assign a emission function to the StateFuncs class
    //!\param str Name of function
    //!\param ptrFunc  pt2StateFunc to use for StateFunc
    void StateFuncs::assignEmissionFunction(std::string& str, emissionFunc ptrFunc){
        
        if (emissionFunctions.count(str)==0){
            emissionFunctions[str]=ptrFunc;
        }
        else{
            std::cerr << "Function Name: " << str << " already exists.   You need to choose a new function name that doesn't already exist. For reference here are a list of names already assigned as external functions\nAssigned Names:\n";
            
            std::map<std::string,emissionFunc>::iterator it;
            
            for(it=emissionFunctions.begin();it!=emissionFunctions.end();it++){
                std::cerr << "\t" << it->first <<std::endl;
            }
        }
    }
	
	void StateFuncs::assignEmissionFunction(const char* str, emissionFunc ptrFunc){
		std::string st(str);
		assignEmissionFunction(st, ptrFunc);
	}
	
	
	
	//!Assign a Univariate Probability Distribution Function to the StateFuncs class
	//!\param str Name of function
	//!\param ptrFunc Pointer to pdfFunc to use for continuous emission
	void StateFuncs::assignPDFFunction(std::string& str, pdfFunc ptrFunc){
		if (pdfFunctions.count(str)==0){
			pdfFunctions[str]=ptrFunc;
		}
		else{
            std::cerr << "Function Name: " << str << " already exists.   You need to choose a new function name that doesn't already exist. For reference here are a list of names already assigned as external functions\nAssigned Names:\n";
            
            std::map<std::string,pdfFunc>::iterator it;
            
            for(it=pdfFunctions.begin();it!=pdfFunctions.end();it++){
                std::cerr << "\t" << it->first <<std::endl;
            }
        }
	}
	
	void StateFuncs::assignPDFFunction(const char* str, pdfFunc ptrFunc){
		std::string st(str);
		assignPDFFunction(st, ptrFunc);
	}
	
	//!Assign a Multivariate Probability Distribution Function to the StateFuncs class
	//!\param str Name of function
	//!\param ptrFunc Pointer to pdfFunc to use for continuous emission
	void StateFuncs::assignMultivariatePdfFunction(std::string& str, multiPdfFunc ptrFunc){
		if (multiPdfFunctions.count(str)==0){
			multiPdfFunctions[str]=ptrFunc;
		}
		else{
            std::cerr << "Function Name: " << str << " already exists.   You need to choose a new function name that doesn't already exist. For reference here are a list of names already assigned as external functions\nAssigned Names:\n";
            
            std::map<std::string,multiPdfFunc>::iterator it;
            
            for(it=multiPdfFunctions.begin();it!=multiPdfFunctions.end();it++){
                std::cerr << "\t" << it->first <<std::endl;
            }
        }
	}
	
	void StateFuncs::assignMultivariatePdfFunction(const char* str, multiPdfFunc ptrFunc){
		std::string st(str);
		assignMultivariatePdfFunction(st, ptrFunc);
	}
    
    
    
    //!Get pointer to function with given name
    //!\param name Name of the function 
    //!\return pt2StateFunc*
    transitionFunc* StateFuncs::getTransitionFunction(std::string& name){
        if (transitionFunctions.count(name)){
            return &transitionFunctions[name];
        }
        else{
            std::cerr << "Function named: " << name << " was not initialized. " <<std::endl;
            
            return NULL;
        }
    }
    
    
    //!Get pointer to function with given name
    //!\param name Name of the function 
    //!\return pt2StateFunc*
    emissionFunc* StateFuncs::getEmissionFunction(std::string& name){
        if (emissionFunctions.count(name)){
            return &emissionFunctions[name];
        }
        else{
            std::cerr << "Function named: " << name << " was not initialized. " <<std::endl;
            
            return NULL;
        }
    }
    
	
	//!Get pointer to Univariate probability distribution function with given name
    //!\param name Name of the function
    //!\return pdfFunc*
    pdfFunc* StateFuncs::getPDFFunction(std::string& name){
        if (pdfFunctions.count(name)){
            return &pdfFunctions[name];
        }
        else{
            std::cerr << "Function named: " << name << " was not initialized. " <<std::endl;
            
            return NULL;
        }
    }
	
	
	//!Get pointer to multivariate probability distribution function with given name
    //!\param name Name of the function
    //!\return pdfFunc*
    multiPdfFunc* StateFuncs::getMultivariatePdfFunction(std::string& name){
        if (multiPdfFunctions.count(name)){
            return &multiPdfFunctions[name];
        }
        else{
            std::cerr << "Function named: " << name << " was not initialized. " <<std::endl;
            
            return NULL;
        }
    }
	
	
	
	void StateFuncs::_loadUnivariatePdf(){
		//Discrete with finite support
		assignPDFFunction("BINOMIAL", static_cast<double (*)(const double , const std::vector<double>*)> (binomial_pdf));
		assignPDFFunction("BETA_BINOMIAL",static_cast<double (*)(const double, const std::vector<double>*)> (beta_binomial_pdf));
		assignPDFFunction("DEGENERATE",static_cast<double (*)(const double, const std::vector<double>*)> (degenerate_pdf));
		assignPDFFunction("DISCRETE_UNIFORM",static_cast<double (*)(const double, const std::vector<double>*)> (discrete_uniform_pdf));
		assignPDFFunction("HYPERGEOMETRIC", static_cast<double (*)(const double, const std::vector<double>*)> (hypergeometric_pdf));
		
		//Discrete with Infinite Support
		assignPDFFunction("BETA_NEGATIVE_BINOMIAL", static_cast<double (*)(const double, const std::vector<double>*)> (beta_negative_binomial_pdf));
		assignPDFFunction("MAXWELL_BOLTZMAN", static_cast<double (*)(const double, const std::vector<double>*)> (maxwell_boltzman_pdf));
		assignPDFFunction("GEOMETRIC", static_cast<double (*)(const double, const std::vector<double>*)> (geometric_pdf));
		assignPDFFunction("LOGARITHMIC", static_cast<double (*)(const double, const std::vector<double>*)> (logarithmic_pdf));
		assignPDFFunction("NEGATIVE_BINOMIAL", static_cast<double (*)(const double, const std::vector<double>*)> (negative_binomial_pdf));
		assignPDFFunction("POISSON", static_cast<double (*)(const double, const std::vector<double>*)> (poisson_pdf));
		assignPDFFunction("YULE_SIMON", static_cast<double (*)(const double, const std::vector<double>*)> (yule_simon_pdf));
		assignPDFFunction("ZIPF", static_cast<double (*)(const double, const std::vector<double>*)> (zipf_pdf));
		assignPDFFunction("ZIPF-MANDELBROT", static_cast<double (*)(const double, const std::vector<double>*)> (zipf_mandelbrot_pdf));
		
		//Continuous on bounded interval
		assignPDFFunction("ARCSINE", static_cast<double (*)(const double, const std::vector<double>*)> (arcsine_pdf));
		assignPDFFunction("BETA", static_cast<double (*)(const double, const std::vector<double>*)> (beta_pdf));
		assignPDFFunction("LOGIT_NORMAL", static_cast<double (*)(const double, const std::vector<double>*)> (logit_normal_pdf));
		assignPDFFunction("CONTINUOUS_UNIFORM", static_cast<double (*)(const double, const std::vector<double>*)> (continuous_uniform_pdf));
		assignPDFFunction("KUMARASWAMY", static_cast<double (*)(const double, const std::vector<double>*)> (kumaraswamy_pdf));
		assignPDFFunction("RAISED_COSINE", static_cast<double (*)(const double, const std::vector<double>*)> (raised_cosine_pdf));
		assignPDFFunction("TRIANGULAR", static_cast<double (*)(const double, const std::vector<double>*)> (triangular_pdf));
		assignPDFFunction("TRUNCATED_NORMAL", static_cast<double (*)(const double, const std::vector<double>*)> (truncated_normal_pdf));
		assignPDFFunction("U_QUADRATIC", static_cast<double (*)(const double, const std::vector<double>*)> (u_quadratic_pdf));
		assignPDFFunction("WIGNER_SEMICIRCLE", static_cast<double (*)(const double, const std::vector<double>*)> (wigner_semicircle_pdf));
		
		//Continous on semi-bounded interval
		assignPDFFunction("BETA_PRIME",  static_cast<double (*)(const double, const std::vector<double>*)> (beta_prime_pdf));
		assignPDFFunction("CHI",  static_cast<double (*)(const double, const std::vector<double>*)> (chi_pdf));
		assignPDFFunction("CHI_SQUARED",  static_cast<double (*)(const double, const std::vector<double>*)> (chi_squared_pdf));
		assignPDFFunction("INVERSE_CHI_SQUARED",  static_cast<double (*)(const double, const std::vector<double>*)> (inverse_chi_squared_pdf));
		assignPDFFunction("SCALED_INVERSE_CHI_SQUARED",  static_cast<double (*)(const double, const std::vector<double>*)> (scaled_inverse_chi_squared_pdf));
		assignPDFFunction("DAGUM",  static_cast<double (*)(const double, const std::vector<double>*)> (dagum_pdf));
		assignPDFFunction("EXPONENTIAL",  static_cast<double (*)(const double, const std::vector<double>*)> (exponential_pdf));
		assignPDFFunction("F_DIST",  static_cast<double (*)(const double, const std::vector<double>*)> (f_pdf));
		assignPDFFunction("FISHERS_Z",  static_cast<double (*)(const double, const std::vector<double>*)> (fishers_z_pdf));
		assignPDFFunction("FOLDED_NORMAL",  static_cast<double (*)(const double, const std::vector<double>*)> (folded_normal_pdf));
		assignPDFFunction("FRECHET",  static_cast<double (*)(const double, const std::vector<double>*)> (frechet_pdf));
		assignPDFFunction("GAMMA",  static_cast<double (*)(const double, const std::vector<double>*)> (gamma_pdf));
		assignPDFFunction("INVERSE_GAMMA", static_cast<double (*)(const double, const std::vector<double>*)> (inv_gamma_pdf));
		assignPDFFunction("HALF_NORMAL",  static_cast<double (*)(const double, const std::vector<double>*)> (half_normal_pdf));
		assignPDFFunction("INVERSE_GAUSSIAN",  static_cast<double (*)(const double, const std::vector<double>*)> (inv_gaussian_pdf));
		assignPDFFunction("LEVY",  static_cast<double (*)(const double, const std::vector<double>*)> (levy_pdf));
		assignPDFFunction("LOG_CAUCHY",  static_cast<double (*)(const double, const std::vector<double>*)> (log_cauchy_pdf));
		assignPDFFunction("LOG_LAPLACE",  static_cast<double (*)(const double, const std::vector<double>*)> (log_laplace_pdf));
		assignPDFFunction("LOG_LOGISTIC", static_cast<double (*)(const double, const std::vector<double>*)> (log_logistic_pdf));
		assignPDFFunction("LOG_NORMAL", static_cast<double (*)(const double, const std::vector<double>*)> (log_normal_pdf));
		assignPDFFunction("PARETO", static_cast<double (*)(const double, const std::vector<double>*)> (pareto_pdf));
		assignPDFFunction("NAKAGAMI",  static_cast<double (*)(const double, const std::vector<double>*)> (nakagami_pdf));
		assignPDFFunction("RAYLEIGH",  static_cast<double (*)(const double, const std::vector<double>*)> (rayleigh_pdf));
		assignPDFFunction("GUMBEL_TYPE_TWO",  static_cast<double (*)(const double, const std::vector<double>*)> (gumbel_type_two_pdf));
		assignPDFFunction("WEIBULL",  static_cast<double (*)(const double, const std::vector<double>*)> (weibull_distribution));
		
		
		//Continuous on unbounded interval
		assignPDFFunction("CAUCHY",  static_cast<double (*)(const double, const std::vector<double>*)> (cauchy_pdf));
		assignPDFFunction("GUMBEL",  static_cast<double (*)(const double, const std::vector<double>*)> (gumbel_pdf));
		assignPDFFunction("GENERALIZED_NORMAL",  static_cast<double (*)(const double, const std::vector<double>*)> (generalized_normal_pdf));
		assignPDFFunction("HYPERBOLIC_SECANT",  static_cast<double (*)(const double, const std::vector<double>*)> (hyperbolic_secant_pdf));
		assignPDFFunction("LAPLACE", static_cast<double (*)(const double, const std::vector<double>*)> (laplace_pdf));
		assignPDFFunction("LOGISTIC",  static_cast<double (*)(const double, const std::vector<double>*)> (logistic_pdf));
		assignPDFFunction("STANDARD_NORMAL",  static_cast<double (*)(const double, const std::vector<double>*)> (standard_normal_pdf));
		assignPDFFunction("NORMAL",  static_cast<double (*)(const double, const std::vector<double>*)> (normal_pdf));
		assignPDFFunction("STUDENT_T",  static_cast<double (*)(const double, const std::vector<double>*)> (students_t_pdf));
		assignPDFFunction("GUMBEL_TYPE_ONE",  static_cast<double (*)(const double, const std::vector<double>*)> (gumbel_type_one_pdf));
		assignPDFFunction("GENERALIZED_EXTREME_VALUE",  static_cast<double (*)(const double, const std::vector<double>*)> (generalized_extreme_value_pdf));
		assignPDFFunction("GENERALIZED_PARETO",  static_cast<double (*)(const double, const std::vector<double>*)> (generalized_pareto_pdf));
						  
	}
	
    
    
    //!Assign a function to the TrackFuncs class
    //!\param str Name of function
    //!\param ptrFunc  pt2TrackFunc to use for TrackFunc
    void TrackFuncs::assignFunction(std::string& str, pt2TrackFunc ptrFunc){
        
        if (functions.count(str)==0){
            functions[str]=ptrFunc;
        }
        else{
            std::cerr << "Function Name: " << str << " already exists.   You need to choose a new function name that doesn't already exist. For reference here are a list of names already assigned as external functions\nAssigned Names:\n";
            
            std::map<std::string,pt2TrackFunc>::iterator it;
            
            for(it=functions.begin();it!=functions.end();it++){
                std::cerr << "\t" << it->first <<std::endl;
            }
        }
    };

    
       
    //!Get pointer to function with given name
    //!\param name Name of the function 
    //!\return pt2TrackFunc*
    pt2TrackFunc* TrackFuncs::getFunction(std::string& name){
        if (functions.count(name)){
            return &functions[name];
        }
        else{
            std::cerr << "Function named: " << name << " was not initialized. " <<std::endl;
            
            return NULL;
        }
    }

}
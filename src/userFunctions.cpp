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
    
    
    //!Assign a emission function to the StateFuncs class
    //!\param str Name of function
    //!\param ptrFunc  pt2StateFunc to use for StateFunc
    void StateFuncs::assignEmmissionFunction(std::string& str, emissionFunc ptrFunc){
        
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
	
	
	
	//!Assign a Probability Distribution Function to the StateFuncs class
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
    
	
	//!Get pointer to probability distribution function with given name
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
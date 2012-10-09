//
//  stateTemplate.h
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

#ifndef StochHMM_stateTemplate_cpp
#define StochHMM_stateTemplate_cpp


#include <vector>
#include <string>
#include <iostream>
#include <map>
#include <algorithm>
#include <set>
#include "text.h"
#include <stdlib.h>
namespace StochHMM{
    
    //!class stencil: contains individial model template and keeps track of coordinates for replacement
    class stencil{
    public:
        stencil();
        stencil(std::string&); //! Construct stencil from string
        
        //MUTATORS
        bool parse(std::string&);
        
        //ACCESSOR
        std::string getTemplate(std::string&, std::string&, std::map<std::string,std::string>&); //!Returns the template with all user and identifiers filled out
        std::string CheckParameter(std::string&);
        
        inline void print(){std::cout << stringify() <<std::endl;};//!Prints template string
        inline std::string stringify(){return states;}; //!Returns template string
        
    private:
        std::string states; //!Template string
        std::set<std::string> UserParameterKeys; //!All UserDefinedVariables in Model Template
    };
    
    
    //!class templates: contains a map of all templates
    /*!class templates: contains a map of all templates, each should have a different identifier.  */
    class templates{
    public:
        templates();
        templates(std::string&); //!Create templates from file
        ~templates();
        
        //MUTATOR
        bool parse(std::string&); //!hand it a string with multiple stencles and have it parse out each one to be added to the templates
        
        void addTemplate(std::string&);        
        
        //ACCESSOR
        inline void print(){std::cout << stringify() <<std::endl;};
        std::string& stringify();
        
        
        //!returns string representation of states with all parameters filled out
        std::string getTemplate(std::string&, std::string&, std::map<std::string,std::string>&); 
        //std::string returnAllFilledTemplates(std::string&,std::map<std::string,std::string>& UserVals);
        //!returns all filled out templates for ID give
        
    private:
        std::map<std::string,stencil*> temps; //!Map of all stencils
    };
}
#endif

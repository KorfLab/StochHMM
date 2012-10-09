//
//  weight.h
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

#ifndef weight_cpp
#define weight_cpp
#include <iostream>
#include <vector>
#include <algorithm>
#include <map>
#include <stdlib.h>
#include "text.h"
#include "stochMath.h"
namespace StochHMM{


bool compXval(const std::pair<double,double>&, const std::pair<double,double>&);

//!weight class contains scaling factors for external scores

//!weight class contains either an 
//!1. Absolute scaling value: meaning the returned value multiplied by constant
//!2. Score Mapping: x and y values describing the transformation from x=returned score
//! to y the scaled score.  If the returned value isn't equal to a specific y the score 
//! is interpolated or extrapolated.
//! If the return score x<MIN then y=scaled MIN
//! If the return score x>MAX then y=scaled MAX
//! Scores provided to weight should be log(Prob)
//! It will return a log(Prob) scaled score


//Scaling weight given some score we map it to some probability or other weighed score
class weight{
public:
    weight();
    ~weight();
    
    //ACCESSORS
    
    /*!
     Checks to see if weight is defined as a constant
     \returns bool TRUE if weight is defined as a constant, else FALSE
     */
    inline bool isAbsolute(){return absolute;};
    
    inline double getAbsolute(){ return absoluteValue;};
    
    inline std::string& getName(){return name;};
    
    /*!
     Returns F(x) given some score x
     \param score The value returned by some function
     \return double F(x) which interpolated or extrapolated as necessary
     */
    double getWeightedScore(double); //!Given a score (x) it will map the score to F(x)

    /*!
     Returns the number of mapped values = equal to size of xVals handed to setWeights(...)
     \returns size_t value of vector size;
     */
    inline size_t size(){return vals->size();};
    
    /*!
     Prints out the weight type to text. Includes the x and F(x) vectors and MIN and MAX
     */
    inline void print(){std::cout << stringify() << std::endl;};
    
    
    /*!
     Converts weights to string
     */
    std::string stringify();
    
    //MUTATORS
    
    /*!
     Sets the scaling factor to a constant value
     \param val double that is returned x->Constant
    */
    inline void setAbsolute(double val){absoluteValue=val; absolute=true; return;};
    
    
    /*!
     Setup the score mapping from x->F(x) 2 vectors of equal length
     \param xVals vector of doubles that contains the scores for x
     \param yVals vector of doubles that contains the corresponding scores F(x)
     \param xMIN pair of doubles defines MIN x->F(x)
     \param xMAX pair of doubles defines MAX x->F(x)
     */
    void setWeights(std::vector<double>& ,std::vector<double>& ,std::pair<double,double>* , std::pair<double,double>* );
    
    
    /*!
     Setup the score mapping from x->F(x) 2 vectors of equal length
     \param xVals vector of doubles that contains the scores for x
     \param yVals vector of doubles that contains the corresponding scores F(x)
     */
    void setWeights(std::vector<double>&,std::vector<double>&);
    
    
    /*!
     Sets the MAX value for x to MAX F(MAX(x));
     \param x The MAX x value to allow
     \param y The value to return for x>=MAX
     */
    void setMaxWeight(double&,double&);
    
    
    /*!
     Sets the MIN value for x to F(MIN(x))
     \param x The MIN x value to allow
     \param y The value to return for x<=MIN
     */
    void setMinWeight(double&,double&);
    
    inline void setName(std::string& nm){name=nm;};
    
    
    //! Parse weight from string
    bool parse(const std::string&);
    
private:
    bool absolute; //! Stores whether Weight is a constant value
    double absoluteValue; //! Stores the constant value
    std::string name;
    
    std::pair<double,double>* maxValue; //!Stores the MAX value of x and the corresponding F(x)
    std::pair<double,double>* minValue; //!Stores the MIN value of x and the corresponding F(x)
    
    std::vector<std::pair<double,double> >* vals;  //!Values for mapping x->F(x)
};
    
    
    
    //!Weights class keeps track of multiple weight types
    //!Each weight can be indexed by a name 
    class weights{
    public:
        //MUTATORS
        
        void addWeight(weight*);
        
        //ACCESSORS
        weight* operator[](std::string&);
        inline size_t count(std::string& name){return wts.count(name);};
        
        void print();
        std::string stringify();
        
    private:
        std::map<std::string, weight*> wts;
    };
    


}
#endif

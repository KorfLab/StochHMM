//
//  externDefinitions.h
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

#ifndef StochHMM_externDefinitions_cpp
#define StochHMM_externDefinitions_cpp

#include <vector>
#include <string>
#include <iostream>
#include <math.h>
#include <map>
#include <set>
#include "text.h"
#include <stdlib.h>
#include "stateInfo.h"

namespace StochHMM{
    
    
    
    class ExDef;
    class weightDef;
    
    //! \class ExDefSequence
    //! Contains external definitions for the sequence
    //! An external definition is how a particular position of the sequence is either
    //! weighted toward a state or in some cases is absolutely produced by a given state
    //! Each position can contain an external definition either absolute or 
    //! states weighted by some value;
    class ExDefSequence{
    public:
        //Constructor
		ExDefSequence(){};
		ExDefSequence(size_t);
        
        //Copy Constructor
        ExDefSequence(const ExDefSequence&);
        
        //Destructor
        ~ExDefSequence();
        
        //Copy Operator
        ExDefSequence& operator=(const ExDefSequence&);
        
        friend class sequences;
        
        bool parse(std::ifstream&,stateInfo&);
        
        //bool getExDef(std::ifstream&,model*);
        
        //ACCESSOR
        bool defined(size_t); // Is ExDef defined at position
        
        bool isAbsolute(size_t); // Is ExDef Absolute at position
        size_t getAbsState(size_t); // Absolute State at position
        
        bool isWeighted(size_t); // Is position weighted
        bool isWeighted(size_t,size_t); // Is position and state weighted
        
        double getWeight(size_t,size_t); // Get Weighted Value for position and state
        
        void print();  // print ExDefSequence to stdout
        std::string stringify(); // Get string representation of ExDefSequence
        
    private:
        std::vector<ExDef*> defs;
        //bool enabled;
        bool _parseAbsDef(stringList& ln,stateInfo&);
        bool _parseWeightDef(stringList& ln,stateInfo&);
    };


    //! ExDef defines absolute external definition for model.   
    //! The model must pass through this state
    class ExDef{
    public:
        
        //Constructor
        ExDef();
        
        friend class ExDefSequence;
        
        inline size_t getState(){return weightedState;};
        inline bool isAbsolute(){return absolute;};
        
        inline void setState(size_t stIter){absolute=true; weightedState=stIter;};
        virtual inline void assignWeight(size_t stateIter,double val){weightedState=stateIter;absolute=true;};
        virtual inline double getWeight(size_t stateIter){if ( stateIter == weightedState){ return 0;} else { return -INFINITY;}};
        virtual std::string stringify();
        
    protected:
        bool absolute; //is it absolute
        //state* st;
        size_t weightedState;  //index iterator to state
    };


    //! weightDef defines weighted external definitions.
    //! Model weights the states at the position with a value
    class weightDef:public ExDef{
    public:
        weightDef(size_t);
        
        friend class ExDefSequence;
        
        inline double getWeight(size_t stateIter){return weights[stateIter];};
        void assignWeight(size_t,double);
        std::string stringify();
    private:
        std::vector<double> weights;
    };

}
#endif

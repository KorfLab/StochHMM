//
//  basicTrellis.h


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

#ifndef BASICTRELLIS_H
#define BASICTRELLIS_H

#include <vector>
#include <iostream>
#include <set>
#include <list>
#include "hmm.h"
#include "sequences.h"
#include "traceback_path.h"
#include <stdint.h>

namespace StochHMM{


    //!\class basicTrellis
    //! Base class for simpleTrellis, stochasticTrellis, and nth Trellis
    //! Contains the base functions that are generalized for all trellis and virtual functions for those function that are specific to simple, stochastic, or nth
    class basicTrellis{
    public:
        basicTrellis(model*,sequences*);
        ~basicTrellis();
        
        void resetCounters();
        
        bool moveToNext();
        bool moveToPrevious();
        
        void initForward();
        void initViterbi();
        void initNthViterbi(size_t);

        void initForwardViterbi();
        void initBackward();
        
        
        double getEndingTransition(size_t);
        
        double getTransition();
        
        int traceback_length();
        
        double exFuncTraceback(transitionFuncParam*);
        
        //Virtual functions defined in derived classes
        //Dependent upon Cell type
        virtual void print()=0;
        virtual std::string stringify()=0;
                
        virtual void calcViterbi()=0;
        virtual void calcNthViterbi(size_t)=0;
        virtual void calcForwardViterbi()=0;
        virtual void calcForward()=0;
        virtual void calcBackward()=0;
        
        virtual void calcEndViterbi()=0;
        virtual void calcNthEndViterbi(size_t)=0;
        virtual void calcEndForward()=0;
        virtual void calcEndForwardViterbi()=0;
        virtual void calcPosterior()=0;
        
        virtual void setViterbi(double&)=0;
        virtual void setNthViterbi(double&)=0;
        virtual void setForward(double&)=0;
        virtual void setBackward(double&)=0;
        virtual void setPtr(int)=0;
        
        virtual void setViterbi(size_t, size_t, double&)=0;
        virtual void setForward(size_t, size_t, double&)=0;
        virtual void setBackward(size_t, size_t, double&)=0;
        virtual void setPtr(size_t, size_t, int)=0;
        
        virtual int getPtr(void)=0;
        virtual double& getViterbi()=0;
        virtual double& getForward()=0;
        virtual double& getBackward()=0;
        
        virtual int getPtr(size_t,size_t)=0;
        virtual double& getViterbi(size_t, size_t)=0;
        virtual double& getForward(size_t, size_t)=0;
        virtual double& getBackward(size_t, size_t)=0;
        //virtual void finalize()=0;
        
        virtual void traceback(traceback_path&)=0;
        virtual void traceback(traceback_path&,size_t)=0;
        
        void tracePosterior(traceback_path&);
        void traceStochPosterior(traceback_path&);
        
        virtual void traceStochViterbi(traceback_path&)=0;
        virtual void traceStochForward(traceback_path&)=0;
    protected:
        
        model* hmm; //model
        sequences* seq; //sequence
        
        double probabilityOfSequence;  //P(x) of the sequence

        double currentProbability;
        
        //Trellis positions
        size_t sequencePosition;
        
        bool externalDefinitionSet;
        bool externalWeighted;
        
        //int currentState;   //State we are going to 
        //int previousState;  //What state we are coming from
        
        
        
        //std::set<int>* previous;
        //std::set<int>* next;
        //std::set<int>::iterator it;
        
        
        //Experimental
        state* currentStatePtr;
        int currentState;
        
        state* previousStatePtr;
        int previousState;
        
        std::set<state*>* currentStates;
        std::set<state*>::iterator currentStatesIterator;
        std::set<state*>* _nextStates;
        
        std::set<state*>* previousStates;
        std::set<state*>* DefinedStates;
        std::set<state*>::iterator iterStates;

        //Posterior Table
        std::vector<std::vector<double> >* posterior;
        std::vector<int>* posterior_pointer;

    private:
    };

}
#endif  /*BASICTRELLIS_H*/

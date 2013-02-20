//
//  externalFuncs.h
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

#ifndef StochHMM_externalFuncs_cpp
#define StochHMM_externalFuncs_cpp

#include <iostream>
#include <map>
#include <string>
#include <vector>
#include "weight.h"
#include "userFunctions.h"
#include "stochTypes.h"
#include "track.h"
#include "text.h"
#include "sequences.h"
#include <stdlib.h>

namespace StochHMM{
    
    /*! \file externalFuncs.h
     \brief Define additional function information for call to user function when calculating the emission or transition value
     At a state the function may be called and passed information.   The function may need traceback data or other sequence data.   Such options are defined by the externalFuncs class
     
     User-provided weight can be passed to externalFuncs to weight the score when passing to 
     emission or transition
     
     Functions must be defined before model import, then will be passed to the externaFuncs class, when created.
     
     Defined in the model definition files
     
     \sa userFunctions.h and userFunctions.cpp
     */
    
    
    
    class transitionFuncParam{
    public:
        //MUTATORS
        transitionFuncParam();
        
        //transitionFuncParam(stringList&, tracks&, weights*, StateFuncs*);
        bool parse(stringList&, tracks&, weights*, StateFuncs*);

        
        
        //External Function
        void setTransFunc(std::string&, transitionFunc*, weight*);
        void setTransTB(tracebackIdentifier, std::string&, combineIdentifier, std::string&);
        
        //ACCESSORS
        
        //!Check to see if a traceback is defined as necessary
        //! \return true if a traceback is necessary
        //! \return false if a traceback is not required
        inline bool isTracebackDefined(){return transFuncTraceback;};
        
        //!Returns what how extensive of a a traceback is required
        //! \return tracbackIdentifier  Type of traceback required 
        //! \sa enum tracebackIdentifier
        inline tracebackIdentifier getTracebackType(){return transFuncTracebackIdentifier;};
        
        //!Return where the traceback is to go back to before stopping
        //! \return std::string State/Label/GFF tag to traceback to
        inline std::string& getTracebackName(){return transFuncTracebackString;};
        
        //!Return what is suppose to be combined (State/Label/GFF)
        //! \return std::string
        inline std::string& getCombineName(){return transFuncCombineString;};
        
        //!Get the Combine type to apply to the traceback
        //!How the traceback is to be combined (editing of traceback)
        //! \return combineIdentifier
        //! \sa enum combineIdentifier
        inline combineIdentifier getCombineType(){return transFuncCombineIdentifier;};
        
        //!Get track name to pass to the function (as defined in the model)
        //! \return std::string Track Name as defined in the model to be passed to the function
        inline std::string& getTrackName(){return trackName;};
        
        inline track* getTrack(){return transFuncTrack;};
        
        double evaluate(const std::string*,size_t, const std::string*, size_t);
        
        void print(); //! prints the string representation of the externFuncs to stdout
        std::string stringify();  //! Return string representation of the externFuncs definition in the model 
        
    private:
        transitionFunc* transFunc;
        std::string trackName;  //! What track to pass to the function
        track* transFuncTrack;     //! Pointer to track function uses
        
        std::string transFuncName; //! < Name of the external function as used in the model
        
        weight* transFuncScaling;  //! < weighting information for values 
        
        //Traceback Information
        bool transFuncTraceback;  //! defines whether a traceback in necessary for the function
        
        tracebackIdentifier transFuncTracebackIdentifier;  //!< What to traceback until
        std::string transFuncTracebackString; //! Contains name of traceback 
        
        combineIdentifier transFuncCombineIdentifier; //! < Contains information on how to combine the traceback information
        std::string transFuncCombineString;  //! Contains what name of what to combine
        
    };
    
    
    
    
    class emissionFuncParam{
    public:
        //MUTATORS
        
        emissionFuncParam(std::string&, StateFuncs*,track*);
        emissionFuncParam();
        
        bool parse(stringList&, tracks&, weights*, StateFuncs*);

        
        
        //External Function
        void setEmissionFunc(std::string&, emissionFunc*, weight*);
        
        //ACCESSORS
        
        //!Get track name to pass to the function (as defined in the model)
        //! \return std::string Track Name as defined in the model to be passed to the function
        inline std::string& getTrackName(){return trackName;};
        
        inline std::string& getName(){return emissionFuncName;}
        
        inline track* getTrack(){return emissionFuncTrack;};
        
        double evaluate(const std::string*,size_t);
        double evaluate(sequences&, size_t);
		double evaluate(sequence&, size_t);
        
        void print(); //! prints the string representation of the externFuncs to stdout
        std::string stringify();  //! Return string representation of the externFuncs definition in the model 
        
    private:
        emissionFunc* emissionFunction;
        std::string trackName;  //! What track to pass to the function
        size_t trackNumber;
        track* emissionFuncTrack;     //! Pointer to track function uses
        
        std::string emissionFuncName; //! < Name of the external function as used in the model
        
        weight* emissionFuncScaling;  //! < weighting information for values 
        
    };

    

}
#endif

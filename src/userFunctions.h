//
//  UserFunctions.h
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

#ifndef StochHMM_userFunctions_cpp
#define StochHMM_userFunctions_cpp

#include <iostream>
#include <map>
#include <string>
#include <vector>
namespace StochHMM{

    
    //!\typedef Pointer to function that takes int, string ptr, string ptr
    //typedef double  (*pt2StateFunc) (const std::string*, const std::string*, size_t);
    typedef double  (*transitionFunc) (const std::string*, const size_t, const std::string*, const size_t);
    
    typedef double  (*emissionFunc) (const std::string*, const size_t);
    

    //!\typedef Pointer to Function that takes a string and returns a vector<float>
    typedef std::vector<double>* (*pt2TrackFunc) (const std::string*);
    
    //!\typedef Pointer to function that takes a string and returns a double
    typedef double (*pt2Attrib) (const std::string*);

    //The map stores the pointers to the different functions by const char* word.
    //!Allows the user to create and integrate their own functions into the model
    //!By specifying the name in the model, and assigning the ptr in the externalFuncs class.
    //!At any State or transition their function will be called and return float
    //!The float will then be applied in the trellis added* to transition or emission.   (**Adding in log space)

    //!\class StateFuncs
    //!Stores pointers to user functions used by the State's Emissions and Transitions
    class StateFuncs{
    public:
        
        void assignTransitionFunction(std::string&, transitionFunc);
        void assignEmmissionFunction(std::string&, emissionFunc);
        
        transitionFunc* getTransitionFunction(std::string&);
        emissionFunc* getEmissionFunction(std::string&);
        
    private:
        std::map<std::string, transitionFunc> transitionFunctions;
        std::map<std::string, emissionFunc> emissionFunctions;
    };


    //!Stores pointer to user functions used to create Real Number Tracks
    class TrackFuncs{
    public:
        void assignFunction(std::string&, pt2TrackFunc);
        pt2TrackFunc* getFunction(std::string&);
    private:
        std::map<std::string, pt2TrackFunc> functions;
    };

}

#endif

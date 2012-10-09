//  trellis.h
// 

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


#ifndef TRELLIS_H
#define TRELLIS_H
#include <vector>
#include <exception>
#include <stdexcept>
#include "simpleTrellis.h"
#include "stochasticTrellis.h"
#include "traceback_path.h"
#include "stochTypes.h"
#include <stdlib.h>
namespace StochHMM{

    //! \class trellis
    //! Description:  Trellis is a Handle class containing the trellis functions for trellis' can be accessed using -> operator
    //! Made to generalize the interface between trells for decoding algorithms Standard, Stochastic, and Nth viterbi
    
    class trellis{
    public:
        //trellis(int);
        
        trellis(model*,sequences*,trellisType type);  //State_size, sequence length
        
        void viterbi();     //Performs Viterbi algorithm
        void nthViterbi(size_t);
        void forward();     //Performs forward algorithm
        void forwardViterbi(); //Performs forward and viterbi algorithms
        void backward();    //Performs backward algorithm
        void posterior();   //Performs forward, backward algorithms
        void decodeAll();   //Performs forward, backward, and viterbi
        
        //void simpleNth();
        
        traceback_path* traceback();
        traceback_path* traceback(size_t);
        traceback_path* posteriorTraceback();
        
        traceback_path* stochasticPosterior();
        
        traceback_path* stochasticViterbiTraceback();
        traceback_path* stochasticForwardTraceback();
        
        multiTraceback* stochasticTraceback(int,decodingType);
        multiTraceback* stochasticViterbiTraceback(int);
        multiTraceback* stochasticForwardTraceback(int);
        
        void print();
        std::string* stringify();
        
    private:
        model* hmm; //Should I store the HMM reference here???
        sequences* seq;
        
        //Type of trellis
        bool simple;
        bool stoch;
        
        //Decoding that has been done on trellis
        bool forwardCompleted;
        bool viterbiCompleted;
        bool nthViterbiCompleted;
        bool backwardCompleted;
        bool posteriorCompleted;
        
        bool _initSimple();
        bool _initStochastic();
        
        traceback_path path;
        traceback_path stochPath;
        multiTraceback paths;
        
        basicTrellis* operator->() {return trell;};
        //basicTrellis& operator*  () const {if (trell) return *trell;
        
        basicTrellis* trell;  
    };
    
}

#endif /*NEWTRELLIS_H*/

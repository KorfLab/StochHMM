//
//  index.h

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

#ifndef StochHMM_index_cpp
#define StochHMM_index_cpp

#include <vector>
#include <algorithm>
#include <stdlib.h>
#include "stochMath.h"
namespace StochHMM{

    /* \class Index 
     \brief Class used for calculating index of word in the emission and transitions
     This class is used to seemlessly work with ambiguious character and their index in the emission or lexical transition tables.  If there is no ambiguous characters the index can be calculated as a int.  However, if there is ambiguous character then it has multiple indices.  
     */
    class Index{
    public:
        Index();
        Index(const Index&);
        
        ~Index();
        
        /*! \fn bool isAmbiguous()
         \brief Returns whether there is an ambiguous character
         */
        inline bool isAmbiguous(){return ambiguous;};
        
        
        void setAmbiguous(const std::vector<size_t>&); //! Set the Index as ambiguous
       
        
        Index& operator= (const Index&); //! Assign Index
        Index operator+ (const Index&); //! Add Index
        void operator += (const Index&); //! Add and assign Index
        void operator += (const size_t); //!Add int to Index
        void operator *= (const size_t); //!Multiply int to Index
        Index operator* (const size_t); //!Multiply int to Index
        size_t size();  //!Get Number of indices
        
        size_t operator[](const size_t);
        
    private:
        size_t index;
        std::vector<size_t>* amb;
        bool ambiguous;
    };


}
#endif

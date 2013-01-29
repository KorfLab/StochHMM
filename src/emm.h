//emm.h
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

#ifndef EMM_H
#define EMM_H
    
    
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <math.h>
#include "track.h"
#include "index.h"
#include "externalFuncs.h"
#include "weight.h"
#include "sequences.h"
#include "lexicalTable.h"
#include <stdlib.h>
namespace StochHMM{


	/*! Emissions for model
	 Contains the emission definition. Each emissions contains the probability, the log(p(x), and counts
	 Counts are used for calculating lower order emissions from higher order.  This is only applicable
	 at the beginning of the sequence.  
	 Each emission depends on some track or function, an emission can have multiple tracks.
	 Or in other words output a single character from each track it is associated with.
	 Tracks can be either alphabetic, real numbers values.
	 Emissions can also call an external function that is user defined.
	 If ambiguity is defined in the alphabet, the emission score can be defined as such
	 If ambiguity is not defined the returned value will be -INFINITY
	 */

	class emm{
	public:
		emm(); //!Constructs an empty emission
		emm(std::string&); //!Constructs emission from a string;
		
		~emm();
		
		friend class state;
		friend class model;
		
		//MUTATORS
		bool parse(std::string&, tracks&, weights*, StateFuncs* );
		
		//!Set the emission to a Real Number
		inline void setRealNumber(){real_number=true;};
		
		//!Set the emission to be the complement 1-P of given value
		inline void setComplement(){complement=true;};
		
		void setLexicalFunction(emissionFunc*);
		
		//ACCESSORS
		
		bool isReal();
		
		//!Check to see if emission will return the complement (1-P) value of emission
		inline bool isComplement(){return complement;};
		
		double get_emission(sequences& , size_t );
		
		//! Get the external Functions defined for the emission
		//! \return externalFuncs*
		inline emissionFuncParam* getExtFunction(){return tagFunc;};
		
		//! Print the string representation of the emission to stdout
		inline void print(){std::cout << stringify()<<std::endl;};
		
		std::string stringify();
		
		inline lexicalTable* getTables(){return &scores;};
		inline bool isSimple(){
			if (!function && tagFunc==NULL){return true;}
			return false;
		}
		
		inline bool isComplex(){
			if (function || tagFunc){return true;}
			return false;
		}
		
	private:
		
		//size_t track_size;
		bool real_number;
		bool continuous;
		bool complement;
		
		track* realTrack;
		
		//Lexical Scoring Tables
		lexicalTable scores;
		
		//Lexical Function Only
		bool function;
		emissionFuncParam* lexFunc;
		
		//Continuous Univariate Distribution
		pdfFunc* pdf;
		std::string pdfName;
		std::vector<double>* dist_parameters;
		
		
		//TODO: Implement Continuous Multivariate Distributions
		//Continuous Multivariate Distribution
		//multiPdfFunc* multiPdf;
		//std::string multiPdfName;
		//std::vector<double>* dist_parameters;
		//std::vector<tracks>* tracks
		
		//TODO:  Implement the external Function capabilities
		emissionFuncParam* tagFunc;
		
		//Private Methods
		bool _processTags(std::string&, tracks&, weights*, StateFuncs*);
			
	};

}
#endif /*EMM_H*/



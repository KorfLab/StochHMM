//
//  sequence.h
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

#ifndef SEQUENCE_H
#define SEQUENCE_H

#include <vector>
#include <stdlib.h>
#include <string>
#include <iostream>
#include <math.h>
#include <fstream>
#include <algorithm>
#include "text.h"
#include "track.h"
#include "stateInfo.h"
#include "externDefinitions.h"
#include "index.h"

//!  \file 

namespace StochHMM{
    //! \class sequence
    //! Contains individual sequence information and functions to deal with importing and digitizing the sequence
	//! Sequence can be either real numbers (double values)  or sequence(characters or words) discrete values
	//! class sequence supports 255 discrete values.
    class sequence{
    public:
        
        //Constructors
        
        sequence();
        sequence(bool);  //True if Real number track, False if alpha
						 //sequence(trackType);
        sequence(std::vector<double>*,track*);
        sequence(std::string&, track*);
        sequence(char* , track*);
		
        ~sequence();
        
        //Copy Constructors
        sequence(const sequence&);
        sequence& operator= (const sequence&);
        
        friend class sequences;
        friend class sequenceStream;
        
        //ACCESSOR
        
        //!Get reference to undigitized sequence
		//!If sequence hasn't been undigitized then it will undigitize it and
		//!store the result.   (Only undigitizes the sequence once, then passes
		//!reference to undigitized sequence)
        inline std::string* getUndigitized(){
            if (!undigitized.empty() || seq->empty()){
                return &undigitized;
            }
            else {
                undigitized = undigitize();
                return &undigitized;
            }
        }
        
        //!Get the size of the sequence
        inline size_t getLength(){return length;};//Returns length of sequence
        
        //!Get the attribute value for the sequence
        //!Selection of model may use this value to determine which model to use
        //! \sa setAttrib
        inline double getAttrib(){return attrib;}; //Returns the Attribute value for the sequence
        
        //!Get pointer to ExDefSequence for the sequence
        //! \return ExDefSequence* 
        inline ExDefSequence* getExDef(){return external;};
        
        //!Check to see if exDef is defined for the sequence
        //! \return true if ExDefSequence is defined for sequence
        //! \return false if no External definition exists for sequence
        inline bool exDefDefined(){if (external){return true;} return false;};
        
        double realValue(size_t);  // Returns Sequence Value at position
        uint8_t  seqValue (size_t);  // Returns Digitized Value at position
        //char   charValue(size_t);  // Returns Alpha Character Value at position
        
        //!Get the size of the sequence
        //! \return size_t size of the sequence
        inline size_t size(){if (realSeq){return real->size();} else {return seq->size();}};  // Returns size of sequence
        
        //! Get the pointer to the track that is defined for the sequence;
        //! \return pointer to track
        inline track* getTrack(){return seqtrk;};
		
		inline void setTrack(track* tr){
			seqtrk = tr;
			return;
		}
        
        
        //! Print the string represntation of the sequence to stdout
        //! Prints the digitized version
        inline void print(){std::cout << stringify() << std::endl;}; //Print sequence to stdout
        std::string stringify(); // Get sequence as string
		
		
		//! Undigitize the sequence
		//! If the sequence has not been digitized then it will return directly
		//! If the sequence has been digitized then it will undigitize it and return it
		//! \return character or word sequence from fasta
        std::string undigitize();
		
        //MUTATOR
        //!Set the sequence attribute value
        //!\param attr Value of attributes for sequence;
        inline void setAttrib(double attr){attrib=attr;}; //!Set the attribute value
        
        //!Set the header of the sequence
        //!\param head Header of the sequence
        inline void setHeader(std::string& head){header=head;};
        
        void setSeq(std::string&,track*);
        void setRealSeq(std::vector<double>*,track*);
		
		inline bool getFasta(std::ifstream& file){return getFasta(file,NULL,NULL);}
        inline bool getFasta(std::ifstream& file, track* trk){ return getFasta(file,trk,NULL);}
		bool getFasta(std::ifstream&, track*, stateInfo*);
        
		
		bool getMaskedFasta(std::ifstream&, track*);
        bool getFastq(std::ifstream&, track*);
        
        inline bool getReal (std::ifstream& file){return getReal(file,NULL,NULL);}
		inline bool getReal (std::ifstream& file, track* trk){ return getReal(file,trk,NULL);}
		bool getReal (std::ifstream&, track*, stateInfo*);
		
        int  getMaxMask(){return max_mask;}
        int  getMask(size_t);
    
        std::string getSymbol(size_t) const;
        
        void get_index(size_t position, int order, std::pair<Index, Index>& word_index);
        
		
		//! Returns the header of the sequence as a std::string
        inline std::string getHeader() { return header; }
        
        bool reverseComplement();
        bool complement();
        bool reverse();
        
		//!Converts sequence digital based on track alphabet
        bool digitize();
		
		//! Shuffles the sequence using std::random_shuffle
		void shuffle();
		
		inline std::vector<uint8_t>* getDigitalSeq(){return seq;}
        
        inline uint8_t operator[](size_t index){return (*seq)[index];}
		
		
		//!Empty Sequence
		void clear();
		
		
        //void getNext (std::ifstream&, track*);
        
        
        //bool _checkSequence(); //!Check the sequence adheres to the track alphabet

    private:
        bool realSeq; //If Real number sequence
        std::string header; // Header from the sequence
        
        double attrib; //Attribute value (Could be %GC or whatever user defines)
        size_t length; //Lenght of the Sequence
        
        track* seqtrk; //Ptr to track describing alphabet and type
        
        ExDefSequence* external; //External definitions
                                 //Stores defined states for given sequence 

        // FIXME:: DIGITIZED SEQUENCES STORED AS SHORT.  NEED TO STANDARDIZE BOTH TRACK AND SEQUENCE CLASS (Track stores as (int) but sequence stores as short.
        std::vector<uint8_t>* seq; // Digitized Sequence
        std::vector<double>* real; // Real Number Sequence
        std::vector<int>* mask; //Stores State masking information for training
        int max_mask;  //Maximum mask number
        
        
        std::string undigitized;  //Undigitized sequence
        
        bool _digitize();  //Digitize the sequence
    };
	
	
	
	//!Randomly generate a sequence based on Probabilities of each character
	//! \param freq  Reference to std::vector<double> that contains frequencies of alphabet corresponding to alphabet in track
	//! \param length  Length of sequence to generate
	//! \param tr Pointer to StochHMM::track where alphabet and ambiguous characters are defined
    sequence random_sequence(std::vector<double>& freq, size_t length, track* tr);
//	sequence random_sequence(emm*);    
//  sequence translate();
    
}
#endif /*SEQUENCE_H*/
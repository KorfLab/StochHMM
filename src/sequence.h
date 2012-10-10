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
#include <string>
#include <iostream>
#include <math.h>
#include <fstream>
#include <algorithm>
#include "text.h"
#include "track.h"
#include "externDefinitions.h"
#include "index.h"
#include <stdlib.h>

//!  \file 

namespace StochHMM{
    //! \class sequence
    //! Contains individual sequence information and functions to deal with importing and digitizing the sequence
    class sequence{
    public:
        
        //Constructors
        
        sequence();
        sequence(bool);  //True if Real number track, False if alpha
        sequence(trackType);
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
        std::string* getUndigitized(){
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
        short  seqValue (size_t);  // Returns Digitized Value at position
        //char   charValue(size_t);  // Returns Alpha Character Value at position
        
        //!Get the size of the sequence
        //! \return size_t size of the sequence
        inline size_t size(){if (realSeq){return real->size();} else {return seq->size();}};  // Returns size of sequence
        
        //! Get the pointer to the track that is defined for the sequence;
        //! \return pointer to track
        inline track* getTrack(){return seqtrk;};
        
        
        //! Print the string represntation of the sequence to stdout
        //! Prints the digitized version
        inline void print(){std::cout << stringify() << std::endl;}; //Print sequence to stdout
        std::string stringify(); // Get sequence as string
        std::string undigitize(); //Undigitize the sequences based on alphabet
        
        //MUTATOR
        //!Set the sequence attribute value
        //!\param attr Value of attributes for sequence;
        inline void setAttrib(double attr){attrib=attr;}; //!Set the attribute value
        
        //!Set the header of the sequence
        //!\param head Header of the sequence
        inline void setHeader(std::string& head){header=head;};
        
        void setSeq(std::string&,track*);
        void setRealSeq(std::vector<double>*,track*);

        bool getFasta(std::ifstream&, track*);
        bool getMaskedFasta(std::ifstream&, track*);
        bool getFastq(std::ifstream&, track*);
        
        
        bool getReal (std::ifstream&, track*);
        int  getMaxMask(){return max_mask;}
        int getMask(size_t);
    
        std::string getSymbol(size_t) const;
        
        void get_index(size_t position, int order, std::pair<Index, Index>& word_index);
        
        inline std::string getHeader() { return header; }
        
        bool reverseComplement();
        bool complement();
        bool reverse();
        
        bool digitize();
        
        
        //void getNext (std::ifstream&, track*);
        
        // Digitizes the character sequence using the track alphabet
        
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
        std::vector<short>* seq; // Digitized Sequence
        std::vector<double>* real; // Real Number Sequence
        std::vector<int>* mask; //Stores State masking information for training
        int max_mask;
        
        
        std::string undigitized;  //Undigitized sequence
        
        bool _digitize();  //Digitize the sequence
        

    };
    
    
    
}
#endif /*SEQUENCE_H*/
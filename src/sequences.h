//
//  sequences.h
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
#ifndef StochHMM_sequences_h
#define StochHMM_sequences_h

#include <vector>
#include <string>
#include <iostream>
#include <math.h>
#include <fstream>
#include <stdlib.h>
#include "text.h"
#include "track.h"
#include "externDefinitions.h"
#include "sequence.h"

namespace StochHMM {
    //TODO: Need to fix so it can handle unrelated sequences
    //Adding sequences error prone if index isn't set;
    

    //! \class sequences container hold the track sequence(s)
    //! Each sequence have to be the same length
    //! Created to pass multiple sequence tracks of different datasets to the HMM
    //! If you need multiple unrelated sequences, use std::vector<sequence> instead. 
    class sequences{
    public:
        //Constructors
        sequences();
        sequences(size_t sz);
        sequences(tracks* tr);
        
        //Copy Constructor
        sequences(const sequences&);
        
        
        sequences& operator= (const sequences&);
        
        ~sequences();
        
        
        //ACCESSOR
        
        //Sequence and Sequence Information
        double realValue(int , size_t); // Get value of Real track(i) in jth position
        short  seqValue( int , size_t); // Get digitized value for sequence track(i) in jth position

        sequence* getSeq(size_t);// Return sequence for Track (i)
        
        //!Get the attribute value for a particular sequence
        //! \param iter Sequence to get the attribute from
        //! \return double value of the attribute set for the sequence
        inline double getAttrib(size_t iter){
            if (iter<num_of_sequences && seq[iter]!=NULL){
                return seq[iter]->getAttrib();  //Returns attrib value for Track (i)
            }
            std::cerr << "No sequence defined at iterator " << iter << std::endl;
            exit(1);
        }
        
        //!Get the header for the first sequence
        //! \return std::string& The header for the first sequence in sequences
        //!
        inline std::string& getHeader(){
            if (num_of_sequences>0 && seq[0]!=NULL){
                return seq[0]->header;
            }
            std::cerr << "No sequence defined." << std::endl;
            exit(1);
        } 
        
        //TODO: fix if iter is not defined
        //!Get the header for the ith sequence
        //! \param iter size_t iterator for ith sequence
        //! \return std::string& The header for the ith sequence
        inline std::string& getHeader(size_t iter){
            if (iter<num_of_sequences){
                if (seq[iter]!=NULL){
                    return seq[iter]->header;
                }
            }
            std::cerr << "No sequence defined at iterator " << iter << std::endl;
            exit(1);
        }
                
        //TODO: need to fix so returns a reference to the sequence.
        //!Get the undigitized ith sequence from sequences
        //! \param iter  size_t iterator for ith sequence
        //! \return std::string of undigitized sequence at ith position
        inline std::string* getUndigitized(size_t iter){
            if (iter>seq.size()){
                std::cerr << "getUndigitized(size_t iter) called where iter is out of range\n";
                return NULL;
            }
            else{
                return seq[iter]->getUndigitized();
            }
        };
        
        // Sizes
        
        //!Get the number of sequence type in sequences
        //! \return size_t value of size
        inline size_t size(){return num_of_sequences;} //Get number of sequences
        
        //!Get the length of the ith sequence
        //! \param iter size_t iterator
        //! \return size_t value of length of sequence at iter
        inline size_t getLength(size_t iter){
            if (iter<num_of_sequences && seq[iter]!=NULL){
                return seq[iter]->getLength(); //Get lenght of sequence in position (i)
            }
            return 0;
        } 
        
        //!Get the length of sequences in general
        //!All of the sequence(s) should be the same length 
        //!\return size_t value of length of all sequences
        inline size_t getLength(){
            if (related_sequences){
                return length; // Get length of related sequence
            }
            return 0;
        }; 
        
        
        bool exDefDefined(size_t); // Is External definition set for state(i)
        bool exDefDefined(size_t,size_t);// Is External definitiion set for state (i) and position (j);
        double getWeight(size_t,int); //! Get Weight value for state at position
        bool exDefDefined();//! Is External definition defined 
        
        //!Print the string representation of digitized sequencs to the stdout
        inline void print(){std::cout<< stringify() << std::endl;}; //Print sequences to stdout
        
        std::string stringify();  //! Get string of sequences
        std::string undigitize(); //! Get sequence based on alphabet
        
        //MUTATOR
        
        //inline void addSeq(sequence* sq){seq.push_back(sq);}; //Add sequence to sequences
        
        void addSeq(sequence* sq);
        
        //! Add sequence in the track position
        //! \param 
        void addSeq(sequence*,size_t);
        
        //! Add sequence for track
        void addSeq(sequence*,track*);
        
        //! Set Length of sequences
        void setLength(size_t len);
        
        //! Set external definition
        //! \param ex  External definition to assign to the sequences
        inline void setExDef(ExDefSequence* ex){external=ex;};
        
        inline bool isSameSize(){
            return same_length;
        }
		
		sequence& operator[](size_t index){return *seq[index];}
        
    private:
        //EXTERNAL DEFINITIONS
        ExDefSequence* external;
        
        std::vector<sequence*> seq; 
        
        size_t length;  //Length of the sequence(s)
        size_t num_of_sequences;
        
        bool related_sequences;
        bool same_length;
    };

}

#endif

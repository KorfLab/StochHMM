//
//  sequences.cpp
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

#include "sequences.h"


namespace StochHMM {
        
    
    
    sequences::sequences(){
        length=-1;
        external=NULL;
        related_sequences=false;
        num_of_sequences=0;
        same_length=true;
    }
    
    
    //!Create an empty sequences
    //! \param sz Number of sequences
    sequences::sequences(size_t sz):seq(sz,NULL){
        length=-1;
        external=NULL;
        related_sequences=true;
        num_of_sequences=sz;
        same_length=true;
    }

    
    //!Create a sequences data t
    //! \param sz Number of sequences
    sequences::sequences(tracks* tr):seq(tr->size(),NULL){
        length=-1;
        external=NULL;
        related_sequences=true;
        num_of_sequences=tr->size();
        same_length=true;
    }
    
    
    sequences::sequences(const sequences& rhs){
        external = (rhs.external==NULL) ? NULL : new(std::nothrow) ExDefSequence(*rhs.external);
      
        if (rhs.external != NULL && external==NULL){
            std::cerr << "OUT OF MEMORY\nFile" << __FILE__ << "Line:\t"<< __LINE__ << std::endl;
            exit(1);
        }
        
        length = rhs.length;
        same_length=rhs.same_length;
        num_of_sequences=rhs.num_of_sequences;
        related_sequences=rhs.related_sequences;
        
        for(size_t i=0;i<seq.size();i++){
            sequence* temp=NULL;
            if (rhs.seq[i]!=NULL){
                temp = new(std::nothrow) sequence(*rhs.seq[i]);
                
                if (temp == NULL){
                    std::cerr << "OUT OF MEMORY\nFile" << __FILE__ << "Line:\t"<< __LINE__ << std::endl;
                    exit(1);
                }
            }            
            seq.push_back(temp);
        }
        
    }
    
    
    //!Destroy sequences
    sequences::~sequences(){
        for(size_t i=0;i<seq.size();i++){
            delete seq[i];
            seq[i]=NULL;
        }
        delete external;
        external = NULL;
    }

    
    //! Assignment Operator
    sequences& sequences::operator=(const sequences & rhs){
        external = (rhs.external==NULL) ? NULL : new(std::nothrow) ExDefSequence(*rhs.external);
        
        if (rhs.external != NULL && external==NULL){
            std::cerr << "OUT OF MEMORY\nFile" << __FILE__ << "Line:\t"<< __LINE__ << std::endl;
            exit(1);
        }
        
        length = rhs.length;
        same_length=rhs.same_length;
        num_of_sequences=rhs.num_of_sequences;
        related_sequences=rhs.related_sequences;
        
        for(size_t i=0;i<seq.size();i++){
            sequence* temp=NULL;
            if (rhs.seq[i]!=NULL){
                temp = new(std::nothrow) sequence(*rhs.seq[i]);
                
                if (temp == NULL){
                    std::cerr << "OUT OF MEMORY\nFile" << __FILE__ << "Line:\t"<< __LINE__ << std::endl;
                    exit(1);
                }
            }            
            seq.push_back(temp);
        }
        
        return *this;
    }
    
    
    //!Get the value from a real Number sequence for track trck at position 
    //! \param trck  Sequence track to use
    //! \param position Position in sequence to get value from
    //! \return double value of the real sequence at the position
    double sequences::realValue(int trck,size_t position){
        if (seq[trck]->realSeq){
            return seq[trck]->realValue(position);
        }
        return -INFINITY;
        
    }
    
    //TODO: fix if the sequence isn't a alpha sequence
    //!Get the digitized value from the sequence at trck at position
    //! \param trck Sequence track to use
    //! \param position Position in sequence to get the value from
    //! \return short digitized value of the sequence based on track type 
    short sequences::seqValue(int trck, size_t position){
        return seq[trck]->seqValue(position);
    }
    
    //!Get pointer to ith sequence from sequences
    //! \param iter Iterator to use for extracting sequence;
    //! \return sequence* pointer to sequence
    //! \return NULL if no sequence exists at iter
    sequence* sequences::getSeq(size_t iter){
        if(iter<seq.size()){
            return seq[iter];
        }
        
        return NULL;
    }
    
    //!Get std:string representation of all the digitized sequence(s) in sequences
    //! \return std::string representation of all string (digitized)
    std::string sequences::stringify(){
        std::string tmp;
        for(int i=0; i<size(); i++){
            
            
            if (seq[i]==NULL){
                tmp+= ">TRACK: " + int_to_string(i) + ":\t" ;
                tmp+= "<<EMPTY>>\n" ;
            }
            else{
                track* trk = seq[i]->getTrack();
                tmp+= ">" + trk->getName();
                tmp+= seq[i]->stringify();
            }
        }
        return tmp;
    }
    
    //!Get std::string representation of all the undigitized sequence(s) in sequences
    //! \return std::string representation of all string (digitized)
    std::string sequences::undigitize(){
        std::string output;
        for(int i=0;i<size();i++){
            if (seq[i]==NULL){
                output+= ">TRACK: " + int_to_string(i) + ":\t" ;
                output+= "<<EMPTY>>\n" ;
            }
            else{
                track* trk = seq[i]->getTrack();
                output+= ">" + trk->getName();
                output+= seq[i]->undigitize();
            }
        }
        return output;
    }
    
    //!Check to see if ther is an external definition defined for a certain position within the sequences
    //! \param pos Position within the sequence
    //! \return True if there is an external definition defined for the position
    bool sequences::exDefDefined(size_t pos){
        if (external && external->defs[pos]){
            return true;
        }
        else{
            return false;
        }
    }
    
    //!Check to see if there is an external definition defined for a certain state at a certain position
    //! \param pos Position within the sequence
    //! \return True if there is an external definition fo the the state and position
    bool sequences::exDefDefined(size_t pos, size_t stateIter){
        if (external && external->defs[pos]){
            if (external->defined(stateIter)){
                return true;
            }
        }
        return false;
    }
    
    
    //!Check to see if there are any external definitions defined
    //! \return true if there are external definitions defined
    bool sequences::exDefDefined(){
        return external;
    }
    
    //!Get the weight for the state at a certain position in the sequence
    //! \param position Position with the sequence
    //! \param stateIter integer iterator of the state
    double sequences::getWeight(size_t position, int stateIter){
        return external->getWeight(position, stateIter);
    }
    
    //!Add a sequence to the sequences
    //!Sequence is added to a certain position based on the track used by the sequence
    
    //! \param sq Pointer to sequence
    void sequences::addSeq(sequence* sq){
        if (related_sequences){
            track* trk = sq->getTrack();
            seq[trk->getIndex()]=sq;
            setLength(sq->getLength());
            return;
        }
        else{
            seq.push_back(sq);
            setLength(sq->getLength());
            num_of_sequences++;
        }
        
    }
    
    //!Add a sequence to the sequences at a certain position
    //!Sequence is added to a certain position based on the iterator
    //! \param sq Pointer to sequence
    //! \param iter Position in sequences to add sequence
    void sequences::addSeq(sequence* sq, size_t iter){
        if (seq.size()>iter){  //If position already exists just add it
            seq[iter]=sq;
        }
        else{ //If the position doesn't exist extend the vector, then add
            seq.resize(iter+1,NULL);
            seq[iter]=sq;
        }
        
        setLength(sq->getLength());
        
        if (!related_sequences){
            num_of_sequences=iter;
        }
        
        return;
    }
    
    //!Add a sequence to the sequences given a certain track
    //!Sequence is added to a certain position based on the track
    //! \param sq Pointer to sequence
    //! \param tr Track to use when adding sequence
    void sequences::addSeq(sequence* sq,track* tr){
        if (tr!=NULL){
            if (related_sequences){
                size_t index = tr->getIndex();
                addSeq(sq,index);
            }
            else{
                addSeq(sq);
            }
        }
        else{
            if (related_sequences){
                std::cerr << "Track is undefined for related sequences.  Unable to add related sequences if track is not defined\n";
                exit(1);
            }
            else{
                addSeq(sq);
            }
        }
        return;
    }
    
    //! Set the length of the sequence(s) in sequence
    //! Because all the sequences should be the same size
    //! If there size differs when adding a sequence
    //! \exception sDifferentSizeSequences thrown if the sizes differ
    void sequences::setLength(size_t len){
        if (length == std::numeric_limits<size_t>::max()){
            length=len;
        }
        
        if (len==length){
            return;
        }
        else{
            same_length=false;
            if (related_sequences){
                std::cerr << "Sequences have different lengths. Sequences should all have the same length because they are suppose to be related (from different datasets).  For multiple unrelated sequence types, use a different structure.\n";
                exit(20);
            }
        }
        return;
    }
}

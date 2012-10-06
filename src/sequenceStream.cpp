//
//  sequenceStream.cpp
//  StochHMM
//
//  Created by Ken Yu on 7/17/12.
//  Copyright 2012 University of California, Davis. All rights reserved.
//

#include "sequenceStream.h"

namespace StochHMM {
    sequenceStream::sequenceStream(): sequence(){
        bufferSize=0;
        retainSize=0;
        readingFile = false;
    }
    
    sequenceStream::sequenceStream(bool realTrack): sequence(realTrack){
        bufferSize=0;
        retainSize=0;
        readingFile = false;
    }
    
    sequenceStream::sequenceStream(std::vector<double>*vec, track* tr ): sequence(vec, tr){
        bufferSize=0;
        retainSize=0;
        readingFile = false;
    }
    
    sequenceStream::sequenceStream(char* seq, track* tr ): sequence(seq, tr){
        bufferSize=0;
        retainSize=0;
        readingFile = false;
    }
    
    sequenceStream::sequenceStream(std::string& sq, track* tr ): sequence(sq, tr){
        bufferSize=0;
        retainSize=0;
        readingFile = false;
    }
    
    sequenceStream::sequenceStream(size_t buff, size_t ret): sequence(){
        bufferSize=buff;
        retainSize=ret;
        readingFile = false;
    }
    
    sequenceStream::sequenceStream(size_t buff, size_t ret, bool realTrack): sequence(realTrack){
        bufferSize=buff;
        retainSize=ret;
        readingFile = false;
    }

    sequenceStream::sequenceStream(size_t buff, size_t ret, std::vector<double>*vec, track* tr ): sequence(vec, tr){
        bufferSize=buff;
        retainSize=ret;
        readingFile = false;
    }
    
    sequenceStream::sequenceStream(size_t buff, size_t ret, char* seq, track* tr ): sequence(seq, tr){
        bufferSize=buff;
        retainSize=ret;
        readingFile = false;
    }
    
    sequenceStream::sequenceStream(size_t buff, size_t ret, std::string& sq, track* tr ): sequence(sq, tr){
        bufferSize=buff;
        retainSize=ret;
        readingFile = false;
    }

//    sequenceStream::~sequenceStream(){
//        
//    }
    
//    sequenceStream::sequenceStream(std::ifstream& file, size_t buff, size_t retain) {
//        bufferSize=buff;
//        retainSize=retain;
//    }
    
    bool sequenceStream::getFasta(std::ifstream &file, track* trk) {
        
        resetSeq();
        seqtrk=trk;
        
        if (readingFile == false) {
            //Find next header mark 
            while(file.peek() != '>'){
                std::string temp;
                getline(file,temp,'\n');
                
                if (!file.good()){
                    return false;
                }
            }
            getline(file,header,'\n');
            //std::cout << header << std::endl;
            readingFile = true;
        }
        

        bool success;
        
        size_t fillBuffer=0;
        
        //Sequence always begins with whatever was retained
        fillBuffer+=retain.size();
        undigitized+=retain;
        retain.clear();
        
        //The following clears up the previousSeq when the "getline" size 
        //is 2 or more times the buffer size
        fillBuffer+=previousSeq.size();
        if (fillBuffer >= bufferSize) {
            size_t overhanging, charsToKeep;
            overhanging = fillBuffer - bufferSize;
            if (overhanging == 0) {
                undigitized+=previousSeq;
                previousSeq.clear();
            }
            else{
                charsToKeep = previousSeq.size() - overhanging;
                undigitized+=previousSeq.substr(0,charsToKeep);
                std::string temp = previousSeq;
                previousSeq.clear();
                previousSeq+=temp.substr(charsToKeep);
            }
            
            retain+=undigitized.substr(undigitized.size()-retainSize);
            
            //
            //std::cout << "1: " << undigitized << std::endl;
            //
            success = _digitize();
            
            length=seq->size();
            
            return success;
        }
        
        //The remaining line sequence (if any) after the last buffer filled up goes to the next buffer
        undigitized+=previousSeq;
        previousSeq.clear();
        
        //For the case where the buffer filled right before the next header or eof;
        //the "nl_peek" portion inside the loop would be skipped, and when overhanging!=0,
        //there is still some sequence to deal with in the final line
        char nl_peek = file.peek();
        if (nl_peek =='>' || nl_peek == EOF) {
            if (undigitized.size() >= 1) {
                //
               // std::cout << "2: " << undigitized << std::endl;
                //
                if (!_digitize()) {
                    std::cerr << "sequence was not digitized" << std::endl;
                }
            }
            
            success = false;
            readingFile = false;
            length=seq->size();
            return success;
        }
        
        std::string line;
        
        while(getline(file,line,'\n')) {
            fillBuffer+=line.size();
            
            //When the buffer fills up, the entire line may fill it up 
            //or some of the line is left overhanging
            if (fillBuffer >= bufferSize) {
                size_t overhanging, charsToKeep;
                overhanging = fillBuffer - bufferSize;
                if (overhanging == 0) {
                    undigitized+=line;
                }
                else{
                    //Overhanging is stored in previousSeq and accounted for in the next calling
                    //of this function
                    charsToKeep = line.size() - overhanging;
                    undigitized+=line.substr(0,charsToKeep);
                    previousSeq+=line.substr(charsToKeep);

                }
                
                retain+=undigitized.substr(undigitized.size()-retainSize);
                
                //
               // std::cout << "3: " << undigitized << std::endl;
                //
                
                success = _digitize();
                
                break;
            }
            
            undigitized+=line;
            
            //For the case where the buffer has not filled up but the file has reached the next header or eof
            char nl_peek = file.peek();
            if (nl_peek =='>' || nl_peek == EOF) {
                //
               // std::cout << "4: " << undigitized << std::endl;
                //
                if (!_digitize()) {
                    std::cerr << "sequence was not digitized" << std::endl;
                }
                success = false;
                readingFile = false;
                break;
            }

        }
        
        length = seq->size();
        return success;
    }
    
    void sequenceStream::resetSeq() {
        if (realSeq) {
            real=new(std::nothrow) std::vector<double>;
            
            if (real==NULL){
                std::cerr << "OUT OF MEMORY\nFile" << __FILE__ << "Line:\t"<< __LINE__ << std::endl;
                exit(1);
            }
            
            seq=NULL;
        }
        else{
            real=NULL;
            seq=new(std::nothrow) std::vector<short>;
            
            if (seq==NULL){
                std::cerr << "OUT OF MEMORY\nFile" << __FILE__ << "Line:\t"<< __LINE__ << std::endl;
                exit(1);
            }
        }
    }

}
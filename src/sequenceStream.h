//
//  sequenceStream.h
//  StochHMM
//
//  Created by Ken Yu on 7/17/12.
//  Copyright 2012 University of California, Davis. All rights reserved.
//
#ifndef sequenceStream_H
#define sequenceStream_H

#include "sequence.h"
#include "track.h"
#include <iostream>
#include <fstream>
#include <string>

namespace StochHMM {
    
    class sequenceStream: public sequence {
    public:
        
        sequenceStream();
        sequenceStream(bool);  //True if Real number track, False if alpha
        sequenceStream(std::vector<double>*,track*);
        sequenceStream(char* , track*);
        sequenceStream(std::string&, track*);

        
        sequenceStream(size_t, size_t);
        sequenceStream(size_t, size_t, bool);  //True if Real number track, False if alpha
        sequenceStream(size_t, size_t, std::vector<double>*,track*);
        sequenceStream(size_t, size_t, std::string&, track*);
        sequenceStream(size_t, size_t, char* , track*);
        
        //~sequenceStream();
        
        bool getFasta(std::ifstream&, track*);
        
        inline void setBuffer (size_t buff) { bufferSize = buff; }
        inline void setRetain (size_t ret) { retainSize = ret; }
        
    private:
        size_t bufferSize;
        size_t retainSize;
        
        //!What is left in "getline" after the buffer has been filled
        std::string previousSeq;
        
        //!The retain sequence
        std::string retain;
        
        //!Keeps track of whether the sequence under the same header is being read
        bool readingFile;
        
        //!Reset seq or realseq each time getfasta is called
        void resetSeq ();
    };
}
#endif /*sequenceStream_H*/
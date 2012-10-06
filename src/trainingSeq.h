//
//  trainingSeq.h
//  StochHMM
//
//  Created by Paul Lott on 3/5/12.
//  Copyright (c) 2012 University of California, Davis. All rights reserved.
//

#ifndef StochHMM_trainingSeq_cpp
#define StochHMM_trainingSeq_cpp

#include <string>
#include <iostream>
#include <stdint.h>
#include <queue>
#include <vector>
#include <map>
#include <fstream>
#include <algorithm>
#include <stdio.h>
#include <set>
#include "text.h"

using namespace StochHMM;

#define MAX_BUFFER 1000000;
#define MIN_BUFFER 100;

class stateInfo;
class trainingEmission;

class trainingSeq{
public:
    sequence(string&);
    bool next();
    void reverseComplement();
    int operator[](int pos){return seq[pos];};
    bool importMask();
    int getMask(int i){return mask[i];};
    inline size_t size(){return seq.size();};
    inline string& getDNA(){return DNA;};
    inline string& getHeader(){return header;};
    inline vector<int8_t>& getSequence(){return seq;};
    inline void setMasked(){masked=true;};
    inline bool isMasked(){return masked;};
private:
    std::string header;
    queue<int> seq;
};


class trainingSeqs{
public:
    inline int fileCount(){return seqFiles.size();};
    inline bool eof(size_t iter){return seqFiles[iter]->eof();}
    inline bool good(size_t iter){return seqFiles[iter]->good();}
    inline bool is_open(size_t iter){return seqFiles[iter]->is_open();}
    
    
    inline int queueSize(size_t iter){return seqs.size();}
    bool openFiles(std::string&);
    bool openFiles(char**);
    
    bool importMask(std::string&);
    bool importMask(char**);
    
    void determineAlphabet(size_t iter);
    void setAlphabet(size_t, stringList&);
    
private:
    vector<std::ifstream*> seqFiles;
    std::ifstream* maskFile;
    
    bool masked;
    queue<int> stateMask;
    
    std::vector<std::map<std::string, int> > alphaIndex;
    std::vector<std::vector<std::string> > indexAlpha;
    
    std::vector<stateInfo> stateEmm;
    std::map<int,map<int,int> > transitions //from state <toState,count>
    
    //State definitions of emission/order/joint emission
    
    std::vector<trainingSeq> seqs;
};

class stateInfo{
public:
    
private:
    std::vector<trainingEmission> *emissions;
};

class trainingEmission{
public:
    
private:
    std::set<std::pair<int,int> > seqTracks;  //Track and order information
    //if more than one track then it is a joint
}



#endif

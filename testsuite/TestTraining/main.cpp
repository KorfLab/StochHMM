//
//  main.cpp
//  TestTraining
//
//  Created by Ken Yu on 6/6/12.
//  Copyright 2012 __MyCompanyName__. All rights reserved.
//

#include <string>
#include <iostream>
#include <fstream>
#include "externDefinitions.h"
#include "sequence.h"
#include "sequences.h"
#include "sequenceStream.h"
#include "lexicalTable.h"
#include "Counter.h"
#include "track.h"

using std::string;

bool real_num_track = false;

int period = 1;
StochHMM::Counter count(period);
//StochHMM::sequence seq(real_num_track); 

int main (int argc, const char * argv[])
{
    StochHMM::track* trk = new StochHMM::track();
    
    trk->addAlphabetChar("A");
    trk->addAlphabetChar("C");
    trk->addAlphabetChar("G");
    trk->addAlphabetChar("T");
    
    trk->addComplement("A","T");
    trk->addComplement("C","G");
    trk->addComplement("G","C");
    trk->addComplement("T","A");
    
    count.setTrack(trk);
    size_t ord = 1;
    count.setOrder(ord);
    
    //string line;
    std::ifstream infasta ("S288C_reference_sequence_R64-1-1_20110203.fsa", std::ios::in);
    //std::ifstream infasta ("NC_000913.fa", std::ios::in);
    //std::ifstream infasta ("test_mask.fa", std::ios::in);
    //std::ifstream infasta ("test_stream.fa", std::ios::in);
    
    if (!infasta.is_open()) {
        std::cout << "Error opening file" << std::endl;
    }
    
    //getline(infasta, line);
    //std::cout << line << std::endl;
    
    //size_t num_seqs = 6628;
    //std::vector<StochHMM::sequence*> seqs;
    
    int upstream = 1, downstream = 0;
    
    //while (!infasta.eof() && infasta.good()) {
        //StochHMM::sequence seq(real_num_track);
        //StochHMM::sequence *seq = new StochHMM::sequence(real_num_track);
        //seq.getFasta(infasta,trk);
        //seq->getFasta(infasta,trk);
        //seq->getMaskedFasta(infasta,trk);
        
        //std::cout << seq.getHeader() << std::endl;
        //std::cout << seq.getLength() << std::endl;
        //seq.complement();
        //seq->complement();
        //seq->print();
        
        //seqs.push_back(seq);
        
        //count.countPWM(seq, upstream, downstream);
    //}
    

    
    size_t bufferSize = 10000000;
    size_t retainSize = 1;
    StochHMM::sequenceStream stream(bufferSize, retainSize, real_num_track);
    
    while (!infasta.eof() && infasta.good()) {
        while (stream.getFasta(infasta,trk)) {
            //stream.complement();
            //stream.reverse();
            //stream.reverseComplement();
            count.countGeneral(stream, upstream, downstream);
        }
        //stream.complement();
        //stream.reverse();
        //stream.reverseComplement();
        count.countGeneral(stream, upstream, downstream);
    }
    infasta.close();
    
    //count.countMask(seqs, upstream, downstream);
    
    //exit(1);

    //count.countGeneral(seq, 100, 100);

    count.printTable();
    
    return 0;
}


//
//  main.cpp
//  Test_newAlphabet
//
//  Created by Paul Lott on 11/13/12.
//  Copyright (c) 2012 Korf Lab, Genome Center, UC Davis, Davis, CA. All rights reserved.
//


#include <iostream>
#include <string>
#include <bitset>
#include "mem_trellis.h"
#include "mem_viterbi.h"
#include "hmm.h"
#include "state.h"
#include "sequence.h"
#include "seqTracks.h"
#include "options.h"

using namespace StochHMM;

seqTracks jobs;

int main(int argc, const char * argv[])
{
    model hmm;
    //Check that model argument is defined and import the model
    std::string model_file = "Lettuce_Final.hmm";
    std::string seq_file = "Test.fa";
    
	hmm.import(model_file);
	
	sequence temp("ACGACGTACGTNNNK",hmm.getTrack(0));
	
	for(size_t i=0;i<hmm.state_size();i++){
		std::bitset<STATE_MAX>* temp = hmm[i]->getTo();
		for(size_t j=0;j<1024;j++){
			if (temp->test(j)){
				std::cout << "State:\t" << i << "\t to " << j << std::endl;
				std::cout << "State:\t" << hmm[i]->getName() << "\t to " << hmm[j]->getName() << std::endl;
			}
			
		}
	}
		
	//temp.print();
    
    //jobs.loadSeqs(hmm, seq_file, FASTA);
    
    //seqJob *job=jobs.getJob();
	
	//sequences* temp = job->getSeqs();
	//temp->print();
    
    //mem_trellis test_trellis;
    
    //viterbi_three(&test_trellis, job->getModel(), job->getSeqs());
	//viterbi_two(&test_trellis, job->getModel(), job->getSeqs());
	//viterbi_four(&test_trellis, job->getModel(), job->getSeqs());
    
    return 0;
}


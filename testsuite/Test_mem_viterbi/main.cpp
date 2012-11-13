//
//  main.cpp
//  Test_mem_viterbi
//
//  Created by Paul Lott on 11/2/12.
//  Copyright (c) 2012 Korf Lab, Genome Center, UC Davis, Davis, CA. All rights reserved.
//

#include <iostream>
#include <string>
#include "mem_trellis.h"
#include "mem_viterbi.h"
#include "hmm.h"
#include "sequence.h"
#include "seqTracks.h"
#include "trellis.h"
#include "options.h"

using namespace StochHMM;

seqTracks jobs;

int main(int argc, const char * argv[])
{
    model hmm;
    //Check that model argument is defined and import the model
    std::string model_file = "Lettuce_Final.hmm";
    std::string seq_file = "Test_no_N.fa";
    
	hmm.import(model_file);
    
    jobs.loadSeqs(hmm, seq_file, FASTA);
    
    seqJob *job=jobs.getJob();
    
    mem_trellis test_trellis;
    
    viterbi_three(&test_trellis, job->getModel(), job->getSeqs());
	//viterbi_two(&test_trellis, job->getModel(), job->getSeqs());
	//viterbi_four(&test_trellis, job->getModel(), job->getSeqs());
    
    return 0;
}


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
#include <time.h>
//#include "mem_trellis.h"
//#include "mem_viterbi.h"
#include "new_trellis.h"
#include "hmm.h"
#include "state.h"
#include "sequence.h"
#include "seqTracks.h"
#include "options.h"
#include "traceback_path.h"

using namespace StochHMM;

seqTracks jobs;

double fair(const double val, const std::vector<double>* param){
	return 0.167;
}

double loaded(const double val, const std::vector<double>* param){
	if (val == 6.0){
		return 0.5;
	}
	
	return 0.1;
}


int main(int argc, const char * argv[])
{
    srand(time(NULL));
	model hmm;
	StateFuncs functions;
	functions.assignPDFFunction("FAIR", *fair);
	functions.assignPDFFunction("LOADED", *loaded);
	
	
	
    //Check that model argument is defined and import the model
//    std::string model_file = "Lettuce_Final.hmm";
//    std::string seq_file = "Test.fa";
//
//	std::string model_file = "Dice.hmm";
//	std::string seq_file = "Dice.fa";

//	std::string model_file = "model_V.txt";
//	std::string seq_file = "TestTCR1.fa";

//	std::string model_file = "3_16Eddy.hmm";
//	std::string seq_file = "3_17Eddy.fa";

//	std::string model_file = "Dice_explicit.hmm";
//	std::string seq_file = "Dice_short.fa";
	
	std::string model_file = "Dice_Continuous.hmm";
	std::string seq_file   = "Dice_real.fa";
	
//	std::string model_file = "Dice_amb_wo_definition.hmm";
//	std::string seq_file   = "Dice_amb_short.fa";
	
	hmm.import(model_file,&functions);

//	sequence temp("ACGACGTACGTNNNK",hmm.getTrack(0));
//	
//	track* tr = hmm.getTrack(0);
//	std::cout << "Max Size:" << tr->getAlphaMax() << std::endl;
//	
//	uint8_t word[3];
//	tr->convertIndexToDigital(63, 3, word);
//	
//	std::cout << (int) word[2] << (int) word[1] << (int) word[0] << std::endl;
//	
//	for(size_t i=0;i<hmm.state_size();i++){
//		std::bitset<STATE_MAX>* temp = hmm[i]->getTo();
//		for(size_t j=0;j<1024;j++){
//			if (temp->test(j)){
//				std::cout << "State:\t" << i << "\t to " << j << std::endl;
//				std::cout << "State:\t" << hmm[i]->getName() << "\t to " << hmm[j]->getName() << std::endl;
//			}
//		}
//	}
		
	//temp.print();
    
    jobs.loadSeqs(hmm, seq_file, FASTA);
    
    seqJob *job=jobs.getJob();
	
	//sequences* temp = job->getSeqs();
	//temp->print();
    
    trellis test_trellis(job->getModel(),job->getSeqs());
	
	hmm.print();
    
    //viterbi_three(&test_trellis, job->getModel(), job->getSeqs());
	//viterbi_two(&test_trellis, job->getModel(), job->getSeqs());
	//viterbi_four(&test_trellis, job->getModel(), job->getSeqs());
    
	
	//test_trellis.stochastic_viterbi();
	
	//test_trellis.viterbi();
	
	test_trellis.simple_viterbi();
	
	//test_trellis.simple_posterior();
	//test_trellis.simple_forward();
	//test_trellis.simple_backward();
	
	traceback_path test(job->getModel());
	test_trellis.traceback(test);
	
	//test_trellis.stochastic_traceback(test);
	
	test.print_label();
	
	//test_trellis.forward();
	//test_trellis.backward();
	
    return 0;
}


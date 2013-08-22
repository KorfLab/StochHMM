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
//#include "trellis.h"
//#include "hmm.h"
//#include "state.h"
//#include "sequence.h"
//#include "seqTracks.h"
//#include "options.h"
//#include "traceback_path.h"
#include "StochHMMlib.h"

using namespace StochHMM;

seqTracks jobs;

double fair(const double val, const std::vector<double>* param){
	return log(0.167);
}

double loaded(const double val, const std::vector<double>* param){
	if (val == 6.0){
		return log(0.5);
	}
	
	return log(0.1);
}


int main(int argc, const char * argv[])
{	
	srand(time(NULL));
	model hmm;
	StateFuncs functions;
	functions.assignPDFFunction("FAIR", *fair);
	functions.assignPDFFunction("LOADED", *loaded);
	
	
	
    //Check that model argument is defined and import the model	
	std::string model_file = "Dice_Continuous.hmm";
	std::string seq_file   = "Dice_real.fa";

	//Import model (Pass State functions to model import)
	hmm.import(model_file,&functions);
	//hmm.print();

    //Import sequences
    jobs.loadSeqs(hmm, seq_file, FASTA);
    
    seqJob *job=jobs.getJob();
    
	//Create Trellis
    trellis test_trellis(job->getModel(),job->getSeqs());
    
	//Perform Viterbi
	test_trellis.simple_viterbi();
	
	//Get traceback and print the labels
	traceback_path test(job->getModel());
	test_trellis.traceback(test);
	test.print_label();
	
    return 0;
}


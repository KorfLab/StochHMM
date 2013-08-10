//
//  main.cpp
//  M82_PEN_HET
//
//  Created by Paul Lott on 8/10/13.
//  Copyright (c) 2013 Korf Lab, Genome Center, UC Davis, Davis, CA. All rights reserved.
//

#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <string>
#include <iomanip>
#include <time.h>
#include <fstream>
#include "StochHMMlib.h"

void print_posterior(StochHMM::trellis&);

double m82_m82(const double pos, const std::vector<double>* values){
	size_t position = static_cast<size_t>(pos);
	std::cout << "m82_m82 transition call:" << position << "\t" << (*values)[position] << std::endl;
	return 0;
}

double m82_pen(const double pos, const std::vector<double>* values){
	size_t position = static_cast<size_t>(pos);
	std::cout << "m82_pen transition call:" << position << "\t" << (*values)[position] << std::endl;
	return 0;
}

double m82_het(const double pos, const std::vector<double>* values){
	size_t position = static_cast<size_t>(pos);
	std::cout << "m82_het transition call:" << position << "\t" << (*values)[position] << std::endl;
	return 0;
}

double pen_m82(const double pos, const std::vector<double>* values){
	size_t position = static_cast<size_t>(pos);
	std::cout << "pen_m82 transition call:" << position << "\t" << (*values)[position] << std::endl;
	return 0;
}

double pen_pen(const double pos, const std::vector<double>* values){
	size_t position = static_cast<size_t>(pos);
	std::cout << "pen_pen transition call:" << position << "\t" << (*values)[position] << std::endl;
	return 0;
}

double pen_het(const double pos, const std::vector<double>* values){
	size_t position = static_cast<size_t>(pos);
	std::cout << "pen_het transition call:" << position << "\t" << (*values)[position] << std::endl;
	return 0;
}

double het_m82(const double pos, const std::vector<double>* values){
	size_t position = static_cast<size_t>(pos);
	std::cout << "het_m82 transition call:" << position << "\t" << (*values)[position] << std::endl;
	return 0;
}

double het_pen(const double pos, const std::vector<double>* values){
	size_t position = static_cast<size_t>(pos);
	std::cout << "het_pen transition call:" << position << "\t" << (*values)[position] << std::endl;
	return 0;
}

double het_het(const double pos, const std::vector<double>* values){
	size_t position = static_cast<size_t>(pos);
	std::cout << "het_het transition call:" << position << "\t" << (*values)[position] << std::endl;
	return 0;
}


int main(int argc, const char * argv[])
{
	StochHMM::StateFuncs myFunctions;
	myFunctions.assignPDFFunction("M82_M82", *m82_m82);
	myFunctions.assignPDFFunction("M82_PEN", *m82_pen);
	myFunctions.assignPDFFunction("M82_HET", *m82_het);
	myFunctions.assignPDFFunction("PEN_M82", *pen_m82);
	myFunctions.assignPDFFunction("PEN_PEN", *pen_pen);
	myFunctions.assignPDFFunction("PEN_HET", *pen_het);
	myFunctions.assignPDFFunction("HET_M82", *het_m82);
	myFunctions.assignPDFFunction("HET_PEN", *het_pen);
	myFunctions.assignPDFFunction("HET_HET", *het_het);
	
	std::string model_file = "M82_PEN_HET.hmm";
	StochHMM::model hmm;
	hmm.import(model_file, &myFunctions);
	hmm.print();
	
	StochHMM::seqTracks jobs;
	
	std::string sequence_file = "test.fa";
	jobs.loadSeqs(hmm, sequence_file, StochHMM::FASTA);
	StochHMM::seqJob *job=jobs.getJob();
	job->getSeqs()->print();
	StochHMM::trellis trell(&hmm,job->getSeqs());
	trell.posterior();
	print_posterior(trell);

	
    return 0;
}

//Print the posterior probabilities for each state at each position
//Each state is in separate column
//Each row is on different row
void print_posterior(StochHMM::trellis& trell){
	StochHMM::model* hmm = trell.getModel();
	StochHMM::double_2D* table = trell.getPosteriorTable();
	size_t state_size = hmm->state_size();
	char cstr[200];
	
	std::string output;
	output+="Posterior Probabilities Table\n";
	output+="Model:\t" + hmm->getName() + "\n";
	output+="Sequence:\t" + trell.getSeq()->getHeader() + "\n";
	sprintf(cstr, "Probability of Sequence from Forward: Natural Log'd\t%f\n",trell.getForwardProbability());
	output+= cstr;
	sprintf(cstr, "Probability of Sequence from Backward:Natural Log'd\t%f\n",trell.getBackwardProbability());
	output+= cstr;
	output+= "Position";
	for(size_t i=0;i< state_size; ++i){
		output+= "\t" + hmm->getStateName(i);
	}
	output+="\n";
	
	std::cout <<  output;
	
	
	for(size_t position = 0; position < table->size(); ++position){
		sprintf(cstr, "%ld", position+1);
		output= cstr;
		for (size_t st = 0 ; st < state_size ; st++){
			float val  = exp((*table)[position][st]);
			if (val<= 0.001){
				output+="\t0";
			}
			else if (val == 1.0){
				output+="\t1";
			}
			else{
				sprintf(cstr,"\t%.3f", exp((*table)[position][st]));
				output+= cstr;
			}
			
		}
		output+="\n";
		std::cout << output;
	}
	
	std::cout << std::endl;
	
	return;
	
}



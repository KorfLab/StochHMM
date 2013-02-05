//
//  baum_welch.cpp
//  StochHMM
//
//  Created by Paul Lott on 2/4/13.
//  Copyright (c) 2013 Korf Lab, Genome Center, UC Davis, Davis, CA. All rights reserved.
//

#include "new_trellis.h"

namespace StochHMM {
	
	//TODO:  Need to implement functions to allow Baum-Welch to update the model.
	
	
	void trellis::naive_baum_welch(model* h, sequences* sqs){
		
	}
	
	
	void trellis::naive_baum_welch(){
		if (naive_forward_score == NULL){
			naive_forward();
		}
		
		if (naive_backward_score == NULL){
			naive_backward();
		}
		
		naive_baum_welch_score = new (std::nothrow) double_3D(seq_size, std::vector<std::vector<double> >(state_size, std::vector<double>(state_size, -INFINITY)));
		
		if (naive_baum_welch_score == NULL){
			std::cerr << "Can't allocate memory. OUT OF MEMORY\t" << __FUNCTION__ << std::endl;
			exit(2);
		}
		
		
		double sum(-INFINITY);
		
		for(size_t position = 0; position < seq_size-1; position++){ // Time(t)
			sum = (-INFINITY);
			for (size_t previous = 0; previous < state_size ; previous++){ // state(i)
				for (size_t current = 0; current < state_size; current++){ // state(j)
					(*naive_baum_welch_score)[position][previous][current] = (*naive_forward_score)[position][previous] + getTransition(hmm->getState(previous), current, position) + (*hmm)[current]->get_emission_prob(*seqs, position+1) + (*naive_backward_score)[position+1][current];
					sum = addLog((*naive_baum_welch_score)[position][previous][current], sum);
				}
			}
			for (size_t previous = 0; previous < state_size ; previous++){ // state(i)
				for (size_t current = 0; current < state_size; current++){ // state(j)
					(*naive_baum_welch_score)[position][previous][current] -= sum;
				}
			}
		}
		return;
	}
	
	
	
	void trellis::update_transitions(){
		
		std::cout << "Transitions to Start:\n";
		double updated(-INFINITY);
		for(size_t st = 0 ; st < state_size ; st++){
			updated = ((*naive_backward_score)[0][st] + (*naive_forward_score)[0][st])-ending_forward_prob;
			std::cout << hmm->getStateName(st) << "\t" << exp(updated) << std::endl;
		}
		
		float_2D numerator(state_size, std::vector<float>(state_size,-INFINITY));
		float_2D denominator(state_size, std::vector<float>(state_size,-INFINITY));
		
		for (size_t position = 0; position < seq_size-1; position++){
			for (size_t i = 0; i < state_size; i++){
				for (size_t j = 0; j < state_size; j++){
					numerator[i][j]= addLog((double)numerator[i][j], (*naive_baum_welch_score)[position][i][j]);
					denominator[i][j] = addLog((double)denominator[i][j], ((*naive_forward_score)[position][i] + (*naive_backward_score)[position][i])-ending_forward_prob);
				}
			}
		}
		
		for (size_t i = 0; i < state_size; i++){
			std::cout << hmm->getStateName(i) << "\t";
			for (size_t j = 0; j < state_size; j++){
				std::cout << exp((numerator[i][j]-denominator[i][j])) << "\t";
			}
			std::cout << std::endl;
		}
	}
	
	void trellis::update_emissions(){
		
	}

	
	
}
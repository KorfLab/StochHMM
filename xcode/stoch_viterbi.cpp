//
//  stoch_viterbi.cpp
//  StochHMM
//
//  Created by Paul Lott on 2/4/13.
//  Copyright (c) 2013 Korf Lab, Genome Center, UC Davis, Davis, CA. All rights reserved.
//

#include "new_trellis.h"

namespace StochHMM{
	void trellis::stochastic_viterbi(){
		scoring_previous = new (std::nothrow) std::vector<double> (state_size,-INFINITY);
		scoring_current  = new (std::nothrow) std::vector<double> (state_size,-INFINITY);
		stochastic_table = new (std::nothrow) stochTable(seq_size);
		
		std::bitset<STATE_MAX> next_states;
		std::bitset<STATE_MAX> current_states;
		
		double  viterbi_temp(-INFINITY);
		double  emission(-INFINITY);
		bool	exDef_position(false);
		
		state* init = hmm->getInitial();
		
		std::bitset<STATE_MAX>* initial_to = hmm->getInitialTo();
		std::bitset<STATE_MAX>* from_trans(NULL);
		
		//Calculate Viterbi from transitions from INIT (initial) state
		for(size_t i = 0; i < state_size; ++i){
			if ((*initial_to)[i]){  //if the bitset is set (meaning there is a transition to this state), calculate the viterbi
				
				viterbi_temp = (*hmm)[i]->get_emission_prob(*seqs,0) + getTransition(init, i, 0);;
				
				if (viterbi_temp > -INFINITY){
					if ((*scoring_current)[i] < viterbi_temp){
						(*scoring_current)[i] = viterbi_temp;
					}
					next_states |= (*(*hmm)[i]->getTo());
				}
			}
		}
		
		//		for(size_t i=0; i < state_size; ++i){
		//			std::cout << "Position: 0" << std::endl;
		//			std::cout << exp((*viterbi_current)[i]) << std::endl;
		//		}
		//
		
		for(size_t position = 1; position < seq_size ; ++position ){
			
			//Swap current and previous viterbi scores
			scoring_previous->assign(state_size,-INFINITY);
			swap_ptr = scoring_previous;
			scoring_previous = scoring_current;
			scoring_current = swap_ptr;
			
			//Swap current_states and next states sets
			
			current_states.reset();
			current_states |= next_states;
			next_states.reset();
			
			if (exDef_defined){
				exDef_position = seqs->exDefDefined(position);
			}
			
			//			std::cout << "\nPosition:\t" << position << "\t";
			//			std::cout << "Letter:\t" << seqs->seqValue(0, position) << std::endl;
			
			for (size_t i = 0; i < state_size; ++i){ //i is current state that emits value
				if (!current_states[i]){
					continue;
				}
				
				//current_state = (*hmm)[i];
				//emission = current_state->get_emission(*seqs,position);
				emission = (*hmm)[i]->get_emission_prob(*seqs, position);
				
				//				std::cout << "State Emission:\t" << i << "\t" << exp(emission) << std::endl;
				
				if (exDef_defined && exDef_position){
					emission += seqs->getWeight(position, i);
				}
				
				from_trans = (*hmm)[i]->getFrom();
				
				for (size_t j = 0; j < state_size ; ++j) {  //j is previous state
					if (!(*from_trans)[j]){
						continue;
					}
					
					if ((*scoring_previous)[j] != -INFINITY){
						viterbi_temp = getTransition((*hmm)[j], i , position) + emission + (*scoring_previous)[j];
						
						
						//						std::cout << "Temp Viterbi:\tTransFrom: "<< j << "\tto\t" << i << "\t" << viterbi_temp / log(2) << std::endl;
						//Save partial value to stochastic table
						stochastic_table->push(position-1,i,j,viterbi_temp);
						
						if (viterbi_temp > (*scoring_current)[i]){
							
							(*scoring_current)[i] = viterbi_temp;
						}
						
						next_states |= (*(*hmm)[i]->getTo());
					}
				}
			}
		}
		
		//TODO:  Calculate ending and set the final viterbi and traceback pointer
		//Swap current and previous viterbi scores
		scoring_previous->assign(state_size,-INFINITY);
		swap_ptr = scoring_previous;
		scoring_previous = scoring_current;
		scoring_current = swap_ptr;
		
		for(size_t i = 0; i < state_size ;++i){
			if ((*scoring_previous)[i] > -INFINITY){
				viterbi_temp = (*scoring_previous)[i] + (*hmm)[i]->getEndTrans();
				
				stochastic_table->push(seq_size-1,SIZE_T_MAX, i,viterbi_temp);
				
				if (viterbi_temp > ending_viterbi_score){
					ending_viterbi_score = viterbi_temp;
					ending_viterbi_tb = i;
				}
			}
		}
		
		stochastic_table->finalize();
		stochastic_table->print();
		
		delete scoring_previous;
		delete scoring_current;
		scoring_current = NULL;
		scoring_previous = NULL;
	}
}

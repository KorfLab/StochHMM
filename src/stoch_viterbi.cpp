//
//  stoch_viterbi.cpp
//  StochHMM
//
//  Created by Paul Lott on 2/4/13.
//  Copyright (c) 2013 Korf Lab, Genome Center, UC Davis, Davis, CA. All rights reserved.
//

#include "trellis.h"

namespace StochHMM{
	
	void trellis::stochastic_viterbi(model* h, sequences *sqs){
		hmm = h;
		seqs = sqs;
		seq_size		= seqs->getLength();
		state_size		= hmm->state_size();
		exDef_defined	= seqs->exDefDefined();
		
		
		stochastic_viterbi();
		
	}
	
	void trellis::stochastic_viterbi(){
		if (hmm->isBasic()){
			simple_stochastic_viterbi();
		}
	}
	
	
	void trellis::simple_stochastic_viterbi(model* h, sequences* sqs){
		hmm = h;
		seqs = sqs;
		seq_size		= seqs->getLength();
		state_size		= hmm->state_size();
		exDef_defined	= seqs->exDefDefined();
		
		simple_stochastic_viterbi();
	}
	
	void trellis::simple_stochastic_viterbi(){
		scoring_previous = new (std::nothrow) std::vector<double> (state_size,-INFINITY);
		scoring_current  = new (std::nothrow) std::vector<double> (state_size,-INFINITY);
		traceback_table = new int_2D(seq_size,std::vector<int16_t> (state_size,-1));
		stochastic_table = new (std::nothrow) stochTable(seq_size);
		
		std::bitset<STATE_MAX> next_states;
		std::bitset<STATE_MAX> current_states;
		
		double  viterbi_temp(-INFINITY);
		double  emission(-INFINITY);
		bool	exDef_position(false);
		ending_viterbi_score = -INFINITY;
		
		state* init = hmm->getInitial();
		
		std::bitset<STATE_MAX>* initial_to = hmm->getInitialTo();
		std::bitset<STATE_MAX>* from_trans(NULL);
		
		//Calculate Viterbi from transitions from INIT (initial) state
		for(size_t st = 0; st < state_size; ++st){
			if ((*initial_to)[st]){  //if the bitset is set (meaning there is a transition to this state), calculate the viterbi
				
				viterbi_temp = (*hmm)[st]->get_emission_prob(*seqs,0) + getTransition(init, st, 0);;
				
				if (viterbi_temp > -INFINITY){
					if ((*scoring_current)[st] < viterbi_temp){
						(*scoring_current)[st] = viterbi_temp;
					}
					next_states |= (*(*hmm)[st]->getTo());
				}
			}
		}
		
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
			
			for (size_t st_current = 0; st_current < state_size; ++st_current){ //i is current state that emits value
				
				//Check to see if current state is valid
				if (!current_states[st_current]){
					continue;
				}
				
				emission = (*hmm)[st_current]->get_emission_prob(*seqs, position);
				
				
				if (exDef_defined && exDef_position){
					emission += seqs->getWeight(position, st_current);
				}
				
				if (emission == -INFINITY){
					continue;
				}
				
				//Get list of states that are valid previous states
				from_trans = (*hmm)[st_current]->getFrom();
				
				for (size_t st_previous = 0; st_previous < state_size ; ++st_previous) {  //j is previous state
					if (!(*from_trans)[st_previous]){
						continue;
					}
					
					if ((*scoring_previous)[st_previous] != -INFINITY){
						viterbi_temp = getTransition((*hmm)[st_previous], st_current , position) + emission + (*scoring_previous)[st_previous];
						
						if (viterbi_temp == -INFINITY){
							continue;
						}
						
						//Save partial value to stochastic table
						stochastic_table->push(position-1,st_current,st_previous,viterbi_temp);
						
						if (viterbi_temp > (*scoring_current)[st_current]){
							(*scoring_current)[st_current] = viterbi_temp;
							(*traceback_table)[position][st_current] = st_previous;
						}
						
						next_states |= (*(*hmm)[st_current]->getTo());
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
		
		for(size_t st_previous = 0; st_previous < state_size ;++st_previous){
			if ((*scoring_previous)[st_previous] > -INFINITY){
				viterbi_temp = (*scoring_previous)[st_previous] + (*hmm)[st_previous]->getEndTrans();
				if (viterbi_temp ==  -INFINITY){
					continue;
				}
				stochastic_table->push(seq_size-1,SIZE_T_MAX, st_previous,viterbi_temp);
				
				if (viterbi_temp > ending_viterbi_score){
					ending_viterbi_score = viterbi_temp;
					ending_viterbi_tb = st_previous;
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
	
	void trellis::naive_stochastic_viterbi(model* h, sequences* sqs){
		hmm = h;
		seqs = sqs;
		seq_size		= seqs->getLength();
		state_size		= hmm->state_size();
		exDef_defined	= seqs->exDefDefined();
		
		naive_stochastic_viterbi();
	}
	
	void trellis::naive_stochastic_viterbi(){
		traceback_table = new(std::nothrow) int_2D(seq_size, std::vector<int16_t>(state_size,-1));
		dbl_viterbi_score = new (std::nothrow) double_2D(seq_size, std::vector<double>(state_size,-INFINITY));
		stochastic_table = new (std::nothrow) stochTable(seq_size);
		
		double emission(-INFINITY);
		double viterbi_temp(-INFINITY);
		double trans(-INFINITY);
		double previous(-INFINITY);
		bool	exDef_position(false);
		ending_viterbi_tb = -1;
		ending_viterbi_score = -INFINITY;
		
		state* init = hmm->getInitial();
		
		//Calculate from Initial states
		for(size_t st = 0; st < state_size; ++st){
			viterbi_temp = (*hmm)[st]->get_emission_prob(*seqs,0) +  getTransition(init, st, 0);
			(*dbl_viterbi_score)[0][st]=viterbi_temp;
			(*traceback_table)[0][st]=-1;
		}
		
		//Calculate Forward for all states
		for (size_t position = 1 ; position < seq_size ; ++position){
			
			if (exDef_defined){
				exDef_position = seqs->exDefDefined(position);
			}
			
			for (size_t st_current = 0; st_current < state_size; ++st_current){
				//Calc emissions
				emission = (*hmm)[st_current]->get_emission_prob(*seqs, position);
				
				if (exDef_defined && exDef_position){
					emission += seqs->getWeight(position, st_current);
				}
				
				if (emission == -INFINITY){
					continue;
				}
				
				for (size_t st_previous = 0; st_previous < state_size; ++st_previous){
					previous = (*dbl_viterbi_score)[position-1][st_previous];
					if (previous == -INFINITY){
						continue;
					}
					
					trans = getTransition(hmm->getState(st_previous), st_current, position);
					
					if (trans !=-INFINITY){
						viterbi_temp = emission + trans + previous;
						
						//Save partial value to stochastic table
						stochastic_table->push(position-1,st_current,st_previous,viterbi_temp);
						
						if (viterbi_temp > (*dbl_viterbi_score)[position][st_current]){
							(*dbl_viterbi_score)[position][st_current] = viterbi_temp;
							(*traceback_table)[position][st_current] = previous;
						}
					}
				}
			}
		}
		
		
		//Calculate Ending Transition
		ending_viterbi_tb = -1;
		for (size_t st_previous = 0; st_previous < state_size; ++st_previous){
			
			if ((*hmm)[st_previous]->getEndTrans() != -INFINITY){
				
				if ((*dbl_viterbi_score)[seq_size-1][st_previous] != -INFINITY){
					viterbi_temp = (*dbl_viterbi_score)[seq_size-1][st_previous] + (*hmm)[st_previous]->getEndTrans();
					stochastic_table->push(seq_size-1,SIZE_T_MAX, st_previous,viterbi_temp);
					
					if (viterbi_temp > ending_viterbi_score){
						ending_viterbi_score = viterbi_temp;
						ending_viterbi_tb = st_previous;
					}
				}
			}
			
		}
		
		stochastic_table->finalize();
		stochastic_table->print();
		
		return;
	}

	
	
}

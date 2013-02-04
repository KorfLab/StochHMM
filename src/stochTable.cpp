//
//  stoch_table.cpp
//  StochHMM
//
//  Created by Paul Lott on 2/4/13.
//  Copyright (c) 2013 Korf Lab, Genome Center, UC Davis, Davis, CA. All rights reserved.
//

#include "stochTable.h"

namespace StochHMM {

	stochTable::stochTable(size_t seq_size){
		
		state_val = new(std::nothrow) std::vector<stoch_value>;
		state_val->reserve(seq_size);
		
		position = new(std::nothrow) std::vector<size_t>(seq_size,0);
		if (state_val == NULL){
			std::cerr << "Cannot allocate stochTable- OUT OF MEMORY" << std::endl;
			delete state_val;
			delete position;
			state_val = NULL;
			position = NULL;
			exit(2);
		}
		last_position = 0;
		return;
	}

	stochTable::~stochTable(){
		delete state_val;
		delete position;
		state_val = NULL;
		position = NULL;
	}

	void stochTable::push(size_t pos, size_t st, size_t st_to, float value){
		if (pos != last_position){
			(*position)[last_position] = state_val->size()-1;
			last_position = pos;
		}
		
		state_val->push_back( stoch_value(st,st_to,value));
	}

	void stochTable::finalize(){
		(*position)[last_position] = state_val->size()-1;
		
		last_position = 0;
		
		double sum(-INFINITY);
		int current_state(-1);
		int state_start(-1);
		
		//For each position in the sequence we want to normalize the viterbi values
		//for the state from previous states
		for(size_t i=0;i<position->size();++i){
			
			//Need to sum the values.  If we see another state then we'll need to
			//normalize the previous states that contributed to the sum.
			//Then set the sum for the current state.
			for (size_t j = (i==0) ? 0 : (*position)[i-1]+1; j <= (*position)[i] ; ++j ){
				
				//If this is a different state then we'll need to apply the sum
				//if this is the first state then we'll just set the current state
				//and set the sum value.  Also need to keep track of which state
				// is the first state so when we normalize we only normalize the
				// the previous state
				if ((*state_val)[j].state_id != current_state){
					
					//Normalize
					for(size_t k = state_start; k < j ;++k){
						(*state_val)[k].prob = exp((*state_val)[k].prob - sum);
					}
					
					//Set state and sum value and starting state
					current_state = (*state_val)[j].state_id;
					sum = (*state_val)[j].prob;
					state_start = j;
				}
				else{
					sum = addLog(sum,(double) (*state_val)[j].prob);
				}
			}
			
			//Normalize the last states seen (b/c) we'll exit out of for loop before
			//we have done these
			for(size_t k = state_start; k <= (*position)[i]; ++k){
				(*state_val)[k].prob = exp((*state_val)[k].prob - sum);
			}
			
			//Set values for next loop
			state_start = (*position)[i]+1;
			current_state = -1;
		}
		
		//Assign previous cell relative to (position) value.  Used when traceback
		//That way function won't have to search for correct value;
		for (size_t i =1; i < position->size(); ++i){
			for (size_t j = (*position)[i-1]+1; j <= (*position)[i] ; ++j){
				(*state_val)[j].prev_cell = get_state_position(i-1, (*state_val)[j].state_prev);
			}
		}
		
		return;
	}

	uint16_t stochTable::get_state_position(size_t pos, uint16_t state){
		size_t start_val = (pos == 0) ? 0 : (*position)[pos-1]+1;
		for (size_t i = start_val ; i <= (*position)[pos] ;  ++i){
			if ((*state_val)[i].state_id == state){
				return i-start_val;
			}
		}
		
		return UINT16_MAX;
	}

	void stochTable::print(){
		last_position = 0;
		
		for(size_t i=0;i<position->size();++i){
			for (size_t j = (i==0) ? 0 : (*position)[i-1]+1 ; j <= (*position)[i] ; ++j ){
				std::cout << (*state_val)[j].state_id <<":"<< (*state_val)[j].state_prev << " : " << (*state_val)[j].prob << "\t";
			}
			std::cout << std::endl;
		}
		return;
	}
	
	void stochTable::traceback(traceback_path& path){
		
		double random((double)rand()/((double)(RAND_MAX)+(double)(1)));
		double cumulative_prob(0.0);
		std::cout << random << std::endl;
		
		size_t offset(0);
		uint16_t state_prev;
		
		//Get traceback from END state
		for(size_t i = (*position)[position->size()-2]+1 ; i <= position->back() ; ++i){
			cumulative_prob += (*state_val)[i].prob;
			if (random < cumulative_prob){
				state_prev = (*state_val)[i].state_prev;
				path.push_back(state_prev);
				offset = (*state_val)[i].prev_cell;
				std::cout << "Chose:\t" << state_prev << std::endl;
				break;
			}
		}
		
		for(size_t i = position->size()-2; i != SIZE_T_MAX ; --i){
			random = (double)rand()/((double)(RAND_MAX)+(double)(1));
			std::cout << random << std::endl;
			cumulative_prob = 0;
			
			for (size_t cells = (i == 0) ? 0 + offset : (*position)[i-1]+offset+1; cells < (*position)[i] ; ++cells){
				if ((*state_val)[cells].state_id != state_prev ){
					std::cout << "Houston, We have an error!" << std::endl;
				}
                
				cumulative_prob+=(*state_val)[cells].prob;
				
                if (random <= cumulative_prob){
                    state_prev = (*state_val)[cells].state_prev;
					path.push_back(state_prev);
					offset = (*state_val)[cells].prev_cell;
					std::cout << "Chose:\t" << state_prev << std::endl;
					break;
                }
            }
		}
		return;
	}


}


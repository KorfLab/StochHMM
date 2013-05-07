//
//  stoch_table.cpp
//  StochHMM
//
//  Created by Paul Lott on 2/4/13.
//  Copyright (c) 2013 Korf Lab, Genome Center, UC Davis, Davis, CA. All rights reserved.
//

#include "stochTable.h"

namespace StochHMM {
	
	//!  stochTable constructor
	//! \param [in] seq_size  Size of sequence
	stochTable::stochTable(size_t seq_size):state_val(NULL){
		
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
	
	//! Pushes the information for traceback pointer onto the table
	//! \param pos Position of current traceback pointer in sequence
	//! \param st  Current State
	//! \param st_to Previous State
	//! \param value Unnormalized value of traceback pointer
	//! \sa stochTable::finalize
	void stochTable::push(size_t pos, size_t st, size_t st_to, float value){
		if (pos != last_position){
			(*position)[last_position] = state_val->size()-1;
			last_position = pos;
		}
		
		state_val->push_back( stoch_value(st,st_to,value));
	}
	
	//! Finalized stochTable
	//! Normalizes the traceback pointer values and calculates the previous cell
	//! iterator within the previous position segment.  This speeds up the
	//! the referencing necessary for tracebacks across the traceback table.
	void stochTable::finalize(){
		(*position)[last_position] = state_val->size()-1;
		
		last_position = 0;
		
		double sum(-INFINITY);
		uint16_t current_state(SIZE_MAX);
		size_t state_start(-1);
		
		//For each position in the sequence we want to normalize the viterbi/forward values
		//for the state from previous states
		for(size_t i=0;i<position->size()-1;++i){
			
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
			state_start = (*position)[i] + 1;
			current_state = SIZE_MAX;
		}
		
		
		//Process and normalize the ending cell
		//Calculate ending state sum
		sum = -INFINITY;
		for (size_t j=state_start; j < state_val->size();j++){
			sum=addLog(sum,(double)(*state_val)[j].prob);
		}
		
		//Normalize the ending cell probability
		for (size_t j=state_start; j < state_val->size();j++){
			(*state_val)[j].prob = exp((*state_val)[j].prob - sum);
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
	
	//! Returns the states iterator position within the array
	//! \param pos Position within the sequence
	//! \param state Iterator of state that you want
	//! \return Position within the table where state can be found
	size_t stochTable::get_state_position(size_t pos, uint16_t state){
		size_t start_val = (pos == 0) ? 0 : (*position)[pos-1]+1;
		for (size_t i = start_val ; i <= (*position)[pos] ;  ++i){
			if ((*state_val)[i].state_id == state){
				return i-start_val;
			}
		}
		
		return UINT16_MAX;
	}
	
	//! Prints the stochTable to stdout
	void stochTable::print(){
		std::cout << stringify();
		return;
	}
	
	//! Stringify the stochTable
	std::string stochTable::stringify(){
		std::stringstream str;
		
		last_position = 0;
		
		for(size_t i=0;i<position->size();++i){
			for (size_t j = (i==0) ? 0 : (*position)[i-1]+1 ; j <= (*position)[i] ; ++j ){
				str << (*state_val)[j].state_id <<":"<< (*state_val)[j].state_prev << " : " << (*state_val)[j].prob << "\t";
			}
			str << std::endl;
		}
		return str.str();
	}
	
	
	//TODO: Test stochastic tracebacks
	//! Traceback through the table using the traceback probabilities
	//! \param[out] path Reference to traceback_path
	void stochTable::traceback(traceback_path& path){
		
		double random((double)rand()/((double)(RAND_MAX)+(double)(1)));
		double cumulative_prob(0.0);
		//std::cout << random << std::endl;
		
		size_t offset(0);
		uint16_t state_prev(UINT16_MAX);
		
		//Get traceback from END state
		for(size_t i = (*position)[position->size()-2]+1 ; i <= position->back() ; ++i){
			cumulative_prob += (*state_val)[i].prob;
			if (random <= cumulative_prob){
				state_prev = (*state_val)[i].state_prev;
				path.push_back(state_prev);
				offset = (*state_val)[i].prev_cell;
				//std::cout << "Chose:\t" << state_prev << std::endl;
				break;
			}
		}
		
		//For the rest of the table traceback to the beginning of the sequence
		for(size_t i = position->size()-2; i != SIZE_MAX ; --i){
			random = (double)rand()/((double)(RAND_MAX)+(double)(1));
			//std::cout << random << std::endl;
			cumulative_prob = 0;
			
			for (size_t cells = (i == 0) ? 0 + offset : (*position)[i-1]+offset+1; cells <= (*position)[i] ; ++cells){
//				if ((*state_val)[cells].state_id != state_prev ){
//					std::cout << "Houston, We have an error!" << std::endl;
//				}
                
				cumulative_prob+=(*state_val)[cells].prob;
				
                if (random <= cumulative_prob){
                    state_prev = (*state_val)[cells].state_prev;
					path.push_back(state_prev);
					offset = (*state_val)[cells].prev_cell;
					//std::cout << "Chose:\t" << state_prev << std::endl;
					break;
                }
            }
		}
		return;
	}
	
	
	
	
//	void alt_stochTable::push(size_t pos, size_t state, size_t state_to, float val){
//		(*table)[state][pos].push_back(stoch_val(state_to,val));
//		return;
//	}
//	
//	void alt_stochTable::push_ending(size_t state_to, float val){
//		ending.push_back(stoch_val(state_to,val));
//		return;
//	}
//	
//
//	//! Traceback through the table using the traceback probabilities
//	//! \param[out] path Reference to traceback_path
//	void alt_stochTable::traceback(traceback_path& path){
//		
//		double random((double)rand()/((double)(RAND_MAX)+(double)(1)));
//		double cumulative_prob(0.0);
//		
//		uint16_t state_prev(UINT16_MAX);
//		
//		//Get traceback from END state
//		for(size_t i = 0; i< ending.size(); ++i){
//			cumulative_prob += (ending)[i].prob;
//			if (random <= cumulative_prob){
//				state_prev = ending[i].previous_state;
//				path.push_back(state_prev);
//				break;
//			}
//		}
//		
//		//For the rest of the table traceback to the beginning of the sequence
//		for(size_t position = table[state_prev].size()-1; position != SIZE_MAX ; --position){
//			random = (double)rand()/((double)(RAND_MAX)+(double)(1));
//			cumulative_prob = 0;
//			
//			for (size_t st = 0; st < (*table)[state_prev][position].size(); st++){
//                
//				cumulative_prob+=(*table)[state_prev][position][st].prob;
//				
//                if (random <= cumulative_prob){
//                    state_prev = (*table)[state_prev][position][st].previous_state;
//					path.push_back(state_prev);
//					break;
//                }
//            }
//		}
//		return;
//	}
//	
//	alt_stochTable::alt_stochTable(size_t st, size_t length){
//		states = st;
//		seq_length = length;
//		
//		table = new (std::nothrow) std::vector<sparseArray<std::vector<stoch_val> > >(states);
//		
//		for (size_t i=0;i<states;i++){
//			(*table)[i].resize(seq_length);
//		}
//	}
//	
//	alt_stochTable::~alt_stochTable(){
//		delete table;
//		table=NULL;
//	}
//	
//	//! Finalized stochTable
//	//! Normalizes the traceback pointer values and calculates the previous cell
//	//! iterator within the previous position segment.  This speeds up the
//	//! the referencing necessary for tracebacks across the traceback table.
//	void alt_stochTable::finalize(){		
//		
//		double sum(-INFINITY);		
//		
//		//For each position in the sequence we want to normalize the viterbi/forward values
//		//for the state from previous states
//		for(size_t position=0; position < seq_length;++position){
//			
//			for(size_t state=0; state < states; ++state){
//				if ((*table)[state].defined(position)){
//					sum = -INFINITY;
//					for(size_t value=0; value < (*table)[state][position].size();value++){
//						sum = addLog(sum,(double) (*table)[state][position][value].prob);
//					}
//					
//					for(size_t value=0; value < (*table)[state][position].size();value++){
//						(*table)[state][position][value].prob = exp((*table)[state][position][value].prob-sum);
//					}
//				}
//			}
//		}
//		
//		sum = -INFINITY;
//		for(size_t i=0;i<ending.size();i++){
//			sum = addLog(sum,(double) ending[i].prob);
//		}
//		
//		for(size_t i=0;i<ending.size();i++){
//			ending[i].prob=exp(ending[i].prob-sum);
//		}
//		
//		return;
//	}
//	
//	std::string alt_stochTable::stringify(){
//		
//		std::stringstream str;		
//		for(size_t position = 0; position < seq_length; position++){
//			str << position ;
//			for(size_t state = 0; state < states; state++){
//				if ((*table)[state].defined(position)){
//					for (size_t val = 0; val < (*table)[state][position].size();val++){
//						str << "\t" << state  << ":" << (*table)[state][position][val].previous_state;
//					}
//				}
//			}
//			str << "\n";
//		}
//		
//		return str.str();
//	}
//	
//	void alt_stochTable::print(){
//		std::cout << stringify() << std::endl;
//		return;
//	}
	
	
	alt_simple_stochTable::alt_simple_stochTable(size_t st, size_t length): states(st), seq_length(length){
		table = new (std::nothrow) std::vector<std::vector<std::vector<stoch_val> > > ;
		table->resize(seq_length,std::vector<std::vector<stoch_val> > (states, std::vector<stoch_val>()));
	}
	
	alt_simple_stochTable::~alt_simple_stochTable(){
		delete table;
	}
	
	void alt_simple_stochTable::push(size_t pos, size_t state, size_t state_to, float val){
		(*table)[pos][state].push_back(stoch_val(state_to,val));
		return;
	}
	
	void alt_simple_stochTable::push_ending(size_t state_to, float val){
		ending.push_back(stoch_val(state_to,val));
		return;
	}
	
	//! Finalized stochTable
	//! Normalizes the traceback pointer values and calculates the previous cell
	//! iterator within the previous position segment.  This speeds up the
	//! the referencing necessary for tracebacks across the traceback table.
	void alt_simple_stochTable::finalize(){
		
		double sum(-INFINITY);
		
		//For each position in the sequence we want to normalize the viterbi/forward values
		//for the state from previous states
		for(size_t position=0; position < seq_length;++position){
			
			for(size_t state=0; state < states; ++state){
				if ((*table)[position][state].size() > 0 ){
					sum = -INFINITY;
					for(size_t value=0; value < (*table)[position][state].size(); value++){
						sum = addLog(sum,(double) (*table)[position][state][value].prob);
					}
					
					for(size_t value=0; value < (*table)[position][state].size();value++){
						(*table)[position][state][value].prob = exp((*table)[position][state][value].prob-sum);
					}
				}
			}
		}
		
		sum = -INFINITY;
		for(size_t i=0;i<ending.size();i++){
			sum = addLog(sum,(double) ending[i].prob);
		}
		
		for(size_t i=0;i<ending.size();i++){
			ending[i].prob=exp(ending[i].prob-sum);
		}
		
		return;
	}
	
	//! Traceback through the table using the traceback probabilities
	//! \param[out] path Reference to traceback_path
	void alt_simple_stochTable::traceback(traceback_path& path){
		
		double random((double)rand()/((double)(RAND_MAX)+(double)(1)));
		double cumulative_prob(0.0);
		
		uint16_t state_prev(UINT16_MAX);
		
		//Get traceback from END state
		for(size_t i = 0; i< ending.size(); ++i){
			cumulative_prob += (ending)[i].prob;
			if (random <= cumulative_prob){
				state_prev = ending[i].previous_state;
				path.push_back(state_prev);
				break;
			}
		}
		
		//For the rest of the table traceback to the beginning of the sequence
		for(size_t position = table[state_prev].size()-1; position != SIZE_MAX ; --position){
			random = (double)rand()/((double)(RAND_MAX)+(double)(1));
			cumulative_prob = 0;
			
			for (size_t st = 0; st < (*table)[position][state_prev].size(); st++){
				
				cumulative_prob+=(*table)[position][state_prev][st].prob;
			
				if (random <= cumulative_prob){
					state_prev = (*table)[position][state_prev][st].previous_state;
					path.push_back(state_prev);
					break;
				}
            }
		}
		return;
	}
	
	std::string alt_simple_stochTable::stringify(){

		std::stringstream str;
		for(size_t position = 0; position < seq_length-1; position++){
			str << position ;
			for(size_t state = 0; state < states; state++){
				for (size_t val = 0; val < (*table)[position][state].size();val++){
					str << "\t" << state  << ":" << (*table)[position][state][val].previous_state << " : " << (*table)[position][state][val].prob ;
				}
			}
			str << "\n";
		}

		return str.str();
	}

	void alt_simple_stochTable::print(){
		std::cout << stringify() << std::endl;
		return;
	}


}


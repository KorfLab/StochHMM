//
//  viterbi.cpp
//  StochHMM
//
//  Created by Paul Lott on 2/4/13.
//  Copyright (c) 2013 Korf Lab, Genome Center, UC Davis, Davis, CA. All rights reserved.
//

#include "new_trellis.h"

namespace StochHMM {
	
	void trellis::viterbi(model* h, sequences* sqs){
		//Initialize the table
		hmm = h;
		seqs = sqs;
		seq_size		= seqs->getLength();
		state_size		= hmm->state_size();
		exDef_defined	= seqs->exDefDefined();
		
		//TODO: determine which model and chose the type of algorithm to use;
		viterbi();
	}
	
	void trellis::viterbi(){
		if (hmm->isBasic()){
			simple_viterbi();
		}
	}
	
	void trellis::simple_viterbi(model* h, sequences* sqs){
		//Initialize the table
		hmm = h;
		seqs = sqs;
		seq_size		= seqs->getLength();
		state_size		= hmm->state_size();
		exDef_defined	= seqs->exDefDefined();
		
		simple_viterbi();
	}
	
	void trellis::simple_viterbi(){
		
		if (!hmm->isBasic()){
			std::cerr << "Model isn't a simple/basic HMM.  Use complex algorithms\n";
			return;
		}
		
		//Initialize the traceback table
		if (traceback_table != NULL){
			delete traceback_table;
		}
		
		traceback_table = new int_2D(seq_size,std::vector<int16_t> (state_size,-1));
		scoring_previous = new (std::nothrow) std::vector<double> (state_size,-INFINITY);
		scoring_current  = new (std::nothrow) std::vector<double> (state_size,-INFINITY);
		
		if (scoring_previous == NULL || scoring_current == NULL || traceback_table == NULL){
			std::cerr << "Can't allocate Viterbi score and traceback table. OUT OF MEMORY" << std::endl;
			exit(2);
		}
		
		
		std::bitset<STATE_MAX> next_states;
		std::bitset<STATE_MAX> current_states;
		
		double  viterbi_temp(-INFINITY);
		double  emission(-INFINITY);
		bool	exDef_position(false);
		ending_viterbi_tb = -1;
		ending_viterbi_score = -INFINITY;
		
		state* init = hmm->getInitial();
		
		std::bitset<STATE_MAX>* initial_to = hmm->getInitialTo();
		std::bitset<STATE_MAX>* from_trans(NULL);
		
		//Calculate Viterbi from transitions from INIT (initial) state
		for(size_t st = 0; st < state_size; ++st){
			if ((*initial_to)[st]){  //if the bitset is set (meaning there is a transition to this state), calculate the viterbi
				
				viterbi_temp = (*hmm)[st]->get_emission_prob(*seqs,0) + getTransition(init, st, 0);
				
				if (viterbi_temp > -INFINITY){
					if ((*scoring_current)[st] < viterbi_temp){
						(*scoring_current)[st] = viterbi_temp;
					}
					next_states |= (*(*hmm)[st]->getTo());
				}
			}
		}
		
		//Each position in the sequence
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
			
			//Current states
			for (size_t st_current = 0; st_current < state_size; ++st_current){ //Current state that emits value
				
				//Check to see if current state is valid
				if (!current_states[st_current]){
					continue;
				}
				
				//Get emission of current state
				emission = (*hmm)[st_current]->get_emission_prob(*seqs, position);
				
				
				if (exDef_defined && exDef_position){
					emission += seqs->getWeight(position, st_current);
				}
				
				if (emission == -INFINITY){
					continue;
				}
				
				//Get list of states that are valid previous states
				from_trans = (*hmm)[st_current]->getFrom();
				
				for (size_t st_previous = 0; st_previous < state_size ; ++st_previous) {  //for previous states
					if (!(*from_trans)[st_previous]){
						continue;
					}
					
					//Check that previous state has transition to current state
					//and that the previous viterbi score is not -INFINITY
					if ((*scoring_previous)[st_previous] != -INFINITY){
						viterbi_temp = getTransition((*hmm)[st_previous], st_current , position) + emission + (*scoring_previous)[st_previous];
						
						
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
		
		
		//Calculate ending viterbi score and traceback from END state
		for(size_t st_previous = 0; st_previous < state_size ;++st_previous){
			if ((*scoring_previous)[st_previous] > -INFINITY){
				viterbi_temp = (*scoring_previous)[st_previous] + (*hmm)[st_previous]->getEndTrans();
				
				if (viterbi_temp > ending_viterbi_score){
					ending_viterbi_score = viterbi_temp;
					ending_viterbi_tb = st_previous;
				}
			}
		}
		
		delete scoring_previous;
		delete scoring_current;
		scoring_previous = NULL;
		scoring_current  = NULL;
	}
	
	
	void trellis::naive_viterbi(model* h, sequences* sqs){
		//Initialize the table
		hmm = h;
		seqs = sqs;
		seq_size		= seqs->getLength();
		state_size		= hmm->state_size();
		exDef_defined	= seqs->exDefDefined();
		
		//TODO: determine which model and chose the type of algorithm to use;
		naive_viterbi();

	}
	
	
	void trellis::naive_viterbi(){
		traceback_table = new(std::nothrow) int_2D(seq_size, std::vector<int16_t>(state_size,-1));
		dbl_viterbi_score = new (std::nothrow) double_2D(seq_size, std::vector<double>(state_size,-INFINITY));
		
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
						
						if (viterbi_temp > (*dbl_viterbi_score)[position][st_current]){
							(*dbl_viterbi_score)[position][st_current] = viterbi_temp;
							(*traceback_table)[position][st_current] = st_previous;
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
					if (viterbi_temp > ending_viterbi_score){
						ending_viterbi_score = viterbi_temp;
						ending_viterbi_tb = st_previous;
					}
				}
			}
		}
		return;
	}
	
	
	
		
	
	
//	void trellis::viterbi(){
//		
//		//Initialize the traceback table
//		if (traceback_table != NULL){
//			delete traceback_table;
//		}
//		
//		traceback_table = new int_2D(seq_size,std::vector<int16_t> (state_size,-1));
//		scoring_previous = new (std::nothrow) std::vector<double> (state_size,-INFINITY);
//		scoring_current  = new (std::nothrow) std::vector<double> (state_size,-INFINITY);
//		
//		if (scoring_previous == NULL || scoring_current == NULL || traceback_table == NULL){
//			std::cerr << "Can't allocate Viterbi score and traceback table. OUT OF MEMORY" << std::endl;
//			exit(2);
//		}
//		
//		
//		std::bitset<STATE_MAX> next_states;
//		std::bitset<STATE_MAX> current_states;
//		
//		double  viterbi_temp(-INFINITY);
//		double  emission(-INFINITY);
//		bool	exDef_position(false);
//		
//		
//		//If model is not a basic model, then we need to initialize the explicit duration vector
//		//bool extend_duration keeps track of whether the transition to same state was selected.
//		if (!hmm->isBasic()){
//			explicit_duration_current = new(std::nothrow) std::vector<size_t>(state_size,0);
//			explicit_duration_previous= new(std::nothrow) std::vector<size_t>(state_size,0);
//		}
//		bool extend_duration(false);
//		std::vector<bool>* duration = hmm->get_explicit();
//		
//		state* init = hmm->getInitial();
//		
//		std::bitset<STATE_MAX>* initial_to = hmm->getInitialTo();
//		std::bitset<STATE_MAX>* from_trans(NULL);
//		
//		//Calculate Viterbi from transitions from INIT (initial) state
//		for(size_t i = 0; i < state_size; ++i){
//			if ((*initial_to)[i]){  //if the bitset is set (meaning there is a transition to this state), calculate the viterbi
//				
//				viterbi_temp = (*hmm)[i]->get_emission_prob(*seqs,0) + getTransition(init, i, 0);
//				
//				if (viterbi_temp > -INFINITY){
//					if ((*scoring_current)[i] < viterbi_temp){
//						(*scoring_current)[i] = viterbi_temp;
//					}
//					next_states |= (*(*hmm)[i]->getTo());
//				}
//			}
//		}
//		
//		//		for(size_t i=0; i < state_size; ++i){
//		//			std::cout << "Position: 0" << std::endl;
//		//			std::cout << exp((*viterbi_current)[i]) << std::endl;
//		//		}
//		//
//		
//		for(size_t position = 1; position < seq_size ; ++position ){
//			
//			//Swap current and previous viterbi scores
//			scoring_previous->assign(state_size,-INFINITY);
//			swap_ptr = scoring_previous;
//			scoring_previous = scoring_current;
//			scoring_current = swap_ptr;
//			
//			//Swap current_states and next states sets
//			
//			current_states.reset();
//			current_states |= next_states;
//			next_states.reset();
//			
//			if (exDef_defined){
//				exDef_position = seqs->exDefDefined(position);
//			}
//			
//			//TODO: Check use of external definitions below.
//			
//			std::cout << "\nPosition:\t" << position << "\n";
//			//			std::cout << "Letter:\t" << seqs->seqValue(0, position) << std::endl;
//			
//			for (size_t i = 0; i < state_size; ++i){ //i is current state that emits value
//				if (!current_states[i]){
//					continue;
//				}
//				
//				//current_state = (*hmm)[i];
//				//emission = current_state->get_emission(*seqs,position);
//				emission = (*hmm)[i]->get_emission_prob(*seqs, position);
//				
//				
//				//				std::cout << "State Emission:\t" << i << "\t" << exp(emission) << std::endl;
//				
//				if (exDef_defined && exDef_position){
//					emission += seqs->getWeight(position, i);
//				}
//				
//				from_trans = (*hmm)[i]->getFrom();
//				
//				for (size_t j = 0; j < state_size ; ++j) {  //j is previous state
//					if (!(*from_trans)[j]){
//						continue;
//					}
//					
//					if ((*scoring_previous)[j] != -INFINITY){
//						viterbi_temp = getTransition((*hmm)[j], i , position) + emission + (*scoring_previous)[j];
//						
//						std::cout << exp(getTransition((*hmm)[j],i,position)) << std::endl;
//						
//						//						std::cout << "Temp Viterbi:\tTransFrom: "<< j << "\tto\t" << i << "\t" << viterbi_temp / log(2) << std::endl;
//						
//						
//						if (viterbi_temp > (*scoring_current)[i]){
//							//If transition is from same to same then if it is
//							//explicit duration we'll need to change
//							extend_duration = (i==j) ? true : false;
//							
//							(*scoring_current)[i] = viterbi_temp;
//							(*traceback_table)[position][i] = j;
//						}
//						
//						next_states |= (*(*hmm)[i]->getTo());
//					}
//				}
//				
//				//If explicit durations vector defined, and transition from-to same
//				//then we'll increment the value. Otherwise, set to zero;
//				if (explicit_duration_current){
//					if (extend_duration && (*duration)[i]){
//						(*explicit_duration_current)[i]=(*explicit_duration_previous)[i]+1;
//						extend_duration=false;
//					}
//					else{
//						(*explicit_duration_current)[i]=0;
//					}
//				}
//			}
//			
//			if(explicit_duration_current){
//				swap_ptr_duration = explicit_duration_previous;
//				explicit_duration_previous = explicit_duration_current;
//				explicit_duration_current= swap_ptr_duration;
//				explicit_duration_current->assign(state_size,0);
//			}
//			
//		}
//		
//		//TODO:  Calculate ending and set the final viterbi and traceback pointer
//		//Swap current and previous viterbi scores
//		scoring_previous->assign(state_size,-INFINITY);
//		swap_ptr = scoring_previous;
//		scoring_previous = scoring_current;
//		scoring_current = swap_ptr;
//		
//		for(size_t i = 0; i < state_size ;++i){
//			if ((*scoring_previous)[i] > -INFINITY){
//				viterbi_temp = (*scoring_previous)[i] + (*hmm)[i]->getEndTrans();
//				
//				if (viterbi_temp > ending_viterbi_score){
//					ending_viterbi_score = viterbi_temp;
//					ending_viterbi_tb = i;
//				}
//			}
//		}
//		
//		delete scoring_previous;
//		delete scoring_current;
//		scoring_previous = NULL;
//		scoring_current = NULL;
//	}
	
	
	//TODO:  Need to test and finalize these algorithms perform
	/*	Need to simplify the basic calls so that it checks model and chooses the
	 algorithm to perform
	 
	 Need to establish duration algorithms for forward/backward that first to viterbi
	 to calculate the transition probability.
	 
	 */
	
	
	//! Sparse Complex Viterbi
	//! Stores the transition duration probabilities in a hashmap (memory efficient, slower)
	//! The duratio probabilities can then be used in forward and backward algorithms.
	void trellis::sparse_complex_viterbi(){
		
		//Initialize the traceback table
        if (traceback_table != NULL){
            delete traceback_table;
        }
        
        traceback_table = new int_2D(seq_size,std::vector<int16_t> (state_size,1));
		scoring_previous = new (std::nothrow) std::vector<double> (state_size,-INFINITY);
        scoring_current  = new (std::nothrow) std::vector<double> (state_size,-INFINITY);
		
		if (scoring_previous == NULL || scoring_current == NULL || traceback_table == NULL){
			std::cerr << "Can't allocate Viterbi score and traceback table. OUT OF MEMORY" << std::endl;
			exit(2);
		}
		
		
        std::bitset<STATE_MAX> next_states;
        std::bitset<STATE_MAX> current_states;
		
        double  viterbi_temp(-INFINITY);
        double  emission(-INFINITY);
        bool	exDef_position(false);
		double  transition_prob(-INFINITY);
		ending_viterbi_score = -INFINITY;
		ending_viterbi_tb = -1;
		
		
		//If model is not a basic model, then we need to initialize the explicit duration vector
		//bool extend_duration keeps track of whether the transition to same state was selected.
		if (!hmm->isBasic()){
			explicit_duration_current = new(std::nothrow) std::vector<size_t>(state_size,0);
			explicit_duration_previous= new(std::nothrow) std::vector<size_t>(state_size,0);
		}
		
		bool extend_duration(false);
		
		// Get list of States with explicit duration
		std::vector<bool>* duration = hmm->get_explicit();
		
		//Initialize Duration storage table
		complex_transitions = new std::vector<std::vector< std::map<uint16_t,double>* >* > (state_size,NULL);
		for(size_t i=0; i < state_size; i++){
			if((*duration)[i]){
				(*complex_transitions)[i] = new std::vector<std::map<uint16_t,double>* > (seq_size, NULL);
			}
		}
        
        state* init = hmm->getInitial();
		
		std::bitset<STATE_MAX>* initial_to = hmm->getInitialTo();
		std::bitset<STATE_MAX>* from_trans(NULL);
		
        //Calculate Viterbi from transitions from INIT (initial) state
        for(size_t i = 0; i < state_size; ++i){
            if ((*initial_to)[i]){  //if the bitset is set (meaning there is a transition to this state), calculate the viterbi
				
				//Transitions here are guarenteed to be standard from the initial state
				viterbi_temp = (*hmm)[i]->get_emission_prob(*seqs,0) + getTransition(init, i, 0);
                
				if (viterbi_temp > -INFINITY){
                    if ((*scoring_current)[i] < viterbi_temp){
                        (*scoring_current)[i] = viterbi_temp;
                    }
					next_states |= (*(*hmm)[i]->getTo());
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
			
			//TODO: Check use of external definitions below.
			
			std::cout << "\nPosition:\t" << position << "\n";
			//			std::cout << "Letter:\t" << seqs->seqValue(0, position) << std::endl;
			
            for (size_t st_current = 0; st_current < state_size; ++st_current){ //i is current state that emits value
                if (!current_states[st_current]){
                    continue;
                }
				
                emission = (*hmm)[st_current]->get_emission_prob(*seqs, position);
				
				
				//Check External definitions
				if (exDef_defined && exDef_position){
                    emission += seqs->getWeight(position, st_current);
                }
                
				//Get list of state that transition to current state
				from_trans = (*hmm)[st_current]->getFrom();
				
                for (size_t st_previous = 0; st_previous < state_size ; ++st_previous) {  //j is previous state
                    if (!(*from_trans)[st_previous]){ //if transition not possible next
                        continue;
                    }
					
                    if ((*scoring_previous)[st_previous] != -INFINITY){
						
						transition_prob = getTransition((*hmm)[st_previous], st_current , position);
						
						if ((*duration)[st_previous]){
							if ((*(*complex_transitions)[st_current])[position] == NULL){
								(*(*complex_transitions)[st_current])[position] = new std::map<uint16_t,double>;
							}
							
							(*(*(*complex_transitions)[st_current])[position])[st_previous] = transition_prob;
						}
						
                        viterbi_temp = transition_prob + emission + (*scoring_previous)[st_previous];
                        
						
						if (viterbi_temp > (*scoring_current)[st_current]){
							//If transition is from same to same then if it is
							//explicit duration we'll need to change
							extend_duration = (st_current==st_previous) ? true : false;
							
                            (*scoring_current)[st_current] = viterbi_temp;
                            (*traceback_table)[position][st_current] = st_previous;
                        }
						
						next_states |= (*(*hmm)[st_current]->getTo());
                    }
                }
				
				//If explicit durations vector defined, and transition from-to same
				//then we'll increment the value. Otherwise, set to zero;
				if (explicit_duration_current){
					if (extend_duration && (*duration)[st_current]){
						(*explicit_duration_current)[st_current]=(*explicit_duration_previous)[st_current]+1;
						extend_duration=false;
					}
					else{
						(*explicit_duration_current)[st_current]=0;
					}
				}
            }
			
			if(explicit_duration_current){
				swap_ptr_duration = explicit_duration_previous;
				explicit_duration_previous = explicit_duration_current;
				explicit_duration_current= swap_ptr_duration;
				explicit_duration_current->assign(state_size,0);
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
                
                if (viterbi_temp > ending_viterbi_score){
                    ending_viterbi_score = viterbi_temp;
                    ending_viterbi_tb = st_previous;
                }
            }
        }
        
        delete scoring_previous;
        delete scoring_current;
		scoring_previous = NULL;
		scoring_current = NULL;
	}
	
	
	void trellis::fast_complex_viterbi(){
		
		//Initialize the traceback table
        if (traceback_table != NULL){
            delete traceback_table;
        }
        
        traceback_table = new int_2D(seq_size,std::vector<int16_t> (state_size,-1));
		scoring_previous = new (std::nothrow) std::vector<double> (state_size,-INFINITY);
        scoring_current  = new (std::nothrow) std::vector<double> (state_size,-INFINITY);
		
		if (scoring_previous == NULL || scoring_current == NULL || traceback_table == NULL){
			std::cerr << "Can't allocate Viterbi score and traceback table. OUT OF MEMORY" << std::endl;
			exit(2);
		}
		
		
        std::bitset<STATE_MAX> next_states;
        std::bitset<STATE_MAX> current_states;
		
        double  viterbi_temp(-INFINITY);
        double  emission(-INFINITY);
        bool	exDef_position(false);
		
		
		//If model is not a basic model, then we need to initialize the explicit duration vector
		//bool extend_duration keeps track of whether the transition to same state was selected.
		if (!hmm->isBasic()){
			explicit_duration_current = new(std::nothrow) std::vector<size_t>(state_size,0);
			explicit_duration_previous= new(std::nothrow) std::vector<size_t>(state_size,0);
		}
		bool extend_duration(false);
		std::vector<bool>* duration = hmm->get_explicit();
        
        state* init = hmm->getInitial();
		
		std::bitset<STATE_MAX>* initial_to = hmm->getInitialTo();
		std::bitset<STATE_MAX>* from_trans(NULL);
		
        //Calculate Viterbi from transitions from INIT (initial) state
        for(size_t i = 0; i < state_size; ++i){
            if ((*initial_to)[i]){  //if the bitset is set (meaning there is a transition to this state), calculate the viterbi
				
				viterbi_temp = (*hmm)[i]->get_emission_prob(*seqs,0) + getTransition(init, i, 0);
                
				if (viterbi_temp > -INFINITY){
                    if ((*scoring_current)[i] < viterbi_temp){
                        (*scoring_current)[i] = viterbi_temp;
                    }
					next_states |= (*(*hmm)[i]->getTo());
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
			
			//TODO: Check use of external definitions below.
			
			std::cout << "\nPosition:\t" << position << "\n";
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
                        
						std::cout << exp(getTransition((*hmm)[j],i,position)) << std::endl;
						
						//						std::cout << "Temp Viterbi:\tTransFrom: "<< j << "\tto\t" << i << "\t" << viterbi_temp / log(2) << std::endl;
                        
						
						if (viterbi_temp > (*scoring_current)[i]){
							//If transition is from same to same then if it is
							//explicit duration we'll need to change
							extend_duration = (i==j) ? true : false;
							
                            (*scoring_current)[i] = viterbi_temp;
                            (*traceback_table)[position][i] = j;
                        }
						
						next_states |= (*(*hmm)[i]->getTo());
                    }
                }
				
				//If explicit durations vector defined, and transition from-to same
				//then we'll increment the value. Otherwise, set to zero;
				if (explicit_duration_current){
					if (extend_duration && (*duration)[i]){
						(*explicit_duration_current)[i]=(*explicit_duration_previous)[i]+1;
						extend_duration=false;
					}
					else{
						(*explicit_duration_current)[i]=0;
					}
				}
            }
			
			if(explicit_duration_current){
				swap_ptr_duration = explicit_duration_previous;
				explicit_duration_previous = explicit_duration_current;
				explicit_duration_current= swap_ptr_duration;
				explicit_duration_current->assign(state_size,0);
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
                
                if (viterbi_temp > ending_viterbi_score){
                    ending_viterbi_score = viterbi_temp;
                    ending_viterbi_tb = i;
                }
            }
        }
        
        delete scoring_previous;
        delete scoring_current;
		scoring_previous = NULL;
		scoring_current = NULL;
	}
	

}


//
//  forward.cpp
//  StochHMM
//
//  Created by Paul Lott on 2/4/13.
//  Copyright (c) 2013 Korf Lab, Genome Center, UC Davis, Davis, CA. All rights reserved.
//

#include "new_trellis.h"

namespace StochHMM {
	
	
	void trellis::forward(model* h, sequences* sqs){
		//Initialize the table
		hmm = h;
		seqs = sqs;
		seq_size		= seqs->getLength();
		state_size		= hmm->state_size();
		exDef_defined	= seqs->exDefDefined();
		
        forward();
    }
	
	
	
	void trellis::forward(){
		if (hmm->isBasic()){
			simple_forward();
		}
	}

	
	void trellis::simple_forward(model* h, sequences* sqs){
		hmm = h;
		seqs = sqs;
		seq_size		= seqs->getLength();
		state_size		= hmm->state_size();
		exDef_defined	= seqs->exDefDefined();
		
        simple_forward();
	}
	
	void trellis::simple_forward(){
		
		if (!hmm->isBasic()){
			std::cerr << "Model isn't a simple/basic HMM.  Use complex algorithms\n";
			return;
		}
		
		forward_score	= new (std::nothrow) float_2D(seq_size, std::vector<float>(state_size,-INFINITY));
		scoring_current = new (std::nothrow) std::vector<double> (state_size,-INFINITY);
		scoring_previous= new (std::nothrow) std::vector<double> (state_size,-INFINITY);
		
		if (scoring_current == NULL || scoring_previous == NULL || forward_score == NULL){
			std::cerr << "Can't allocate forward score table. OUT OF MEMORY" << std::endl;
			exit(2);
		}
		
        std::bitset<STATE_MAX> next_states;
        std::bitset<STATE_MAX> current_states;
		
        double  forward_temp(-INFINITY);
        double  emission(-INFINITY);
        bool	exDef_position(false);
        
        state* init = hmm->getInitial();
		
		std::bitset<STATE_MAX>* initial_to = hmm->getInitialTo();
		std::bitset<STATE_MAX>* from_trans(NULL);
		
		
		//		std::cout << "Position: 0" << std::endl;
        //Calculate Viterbi from transitions from INIT (initial) state
        for(size_t st = 0; st < state_size; ++st){
            if ((*initial_to)[st]){  //if the bitset is set (meaning there is a transition to this state), calculate the viterbi
				
				forward_temp = (*hmm)[st]->get_emission_prob(*seqs,0) + getTransition(init, st, 0);
                
				if (forward_temp > -INFINITY){
                    
					(*forward_score)[0][st] = forward_temp;
					(*scoring_current)[st] = forward_temp;
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
			
			//			std::cout << "\nPosition: " << position << std::endl;
            
            for (size_t st_current = 0; st_current < state_size; ++st_current){ //i is current state that emits value
                if (!current_states[st_current]){
                    continue;
                }
                
                emission = (*hmm)[st_current]->get_emission_prob(*seqs, position);
				
				if (exDef_defined && exDef_position){
                    emission += seqs->getWeight(position, st_current);
                }
                
				from_trans = (*hmm)[st_current]->getFrom();
				
                for (size_t previous = 0; previous < state_size ; ++previous) {  //j is previous state
                    if (!(*from_trans)[previous]){
                        continue;
                    }
					
					if ((*scoring_previous)[previous] != -INFINITY){
                        forward_temp = (*scoring_previous)[previous] + emission + getTransition((*hmm)[previous], st_current , position);
						
						if ((*scoring_current)[st_current] == -INFINITY){
							(*scoring_current)[st_current] = forward_temp;
							(*forward_score)[position][st_current] = forward_temp;
						}
						else{
							(*scoring_current)[st_current] = addLog(forward_temp, (*scoring_current)[st_current]);
							(*forward_score)[position][st_current] = (*scoring_current)[st_current];						}
						
						next_states |= (*(*hmm)[st_current]->getTo());
                    }
                }
				//				std::cout << "State: " << current <<"\t" << exp((*forward_score)[position][current]) << std::endl;
            }
		}
		
        //Swap current and previous scores
        scoring_previous->assign(state_size,-INFINITY);
        swap_ptr = scoring_previous;
        scoring_previous = scoring_current;
        scoring_current = swap_ptr;
		
        ending_forward_prob = -INFINITY;
        for(size_t st_previous = 0; st_previous < state_size ;++st_previous){
            if ((*scoring_previous)[st_previous] != -INFINITY){
                forward_temp = (*scoring_previous)[st_previous] + (*hmm)[st_previous]->getEndTrans();
                
                if (forward_temp > -INFINITY){
					if (ending_forward_prob == -INFINITY){
						ending_forward_prob = forward_temp;
					}
					else{
						ending_forward_prob = addLog(ending_forward_prob,forward_temp);
					}
                }
            }
        }
		
		delete scoring_previous;
		delete scoring_current;
		scoring_previous = NULL;
		scoring_current = NULL;
		
	}
	
	
	void trellis::naive_forward(){
		dbl_forward_score = new (std::nothrow) double_2D(seq_size, std::vector<double>(state_size,-INFINITY));
		
		if (dbl_forward_score == NULL){
			std::cerr << "Can't allocate Forward table and traceback table. OUT OF MEMORY\t" << __FUNCTION__ << std::endl;
			exit(2);
		}
		
		double emission(-INFINITY);
		double forward_temp(-INFINITY);
		double trans(-INFINITY);
		double previous(-INFINITY);
		bool	exDef_position(false);

		
		state* init = hmm->getInitial();
		
		//Calculate from Initial states
		for(size_t st = 0; st < state_size; ++st){
			forward_temp = (*hmm)[st]->get_emission_prob(*seqs,0) +  getTransition(init, st, 0);
			(*dbl_forward_score)[0][st]=forward_temp;
		}
		
		//Calculate Forward for all states
		for (size_t position = 1 ; position < seq_size ; ++position){
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
					previous = (*dbl_forward_score)[position-1][st_previous];
					if (previous == -INFINITY){
						continue;
					}
					
					trans = getTransition(hmm->getState(st_previous), st_current, position);
					
					
					if (trans !=-INFINITY){
						forward_temp = previous + emission + trans;
						
						if ((*dbl_forward_score)[position][st_current]==-INFINITY){
							(*dbl_forward_score)[position][st_current] = forward_temp;
						}
						else{
							(*dbl_forward_score)[position][st_current] = addLog(forward_temp, (*dbl_forward_score)[position][st_current]);
						}
					}
				}
			}
		}
		
		
		//Calculate Ending Transition
		ending_forward_prob = -INFINITY;
		for (size_t st_previous = 0; st_previous < state_size; ++st_previous){
			if ((*hmm)[st_previous]->getEndTrans() != -INFINITY){
				if ((*dbl_forward_score)[seq_size-1][st_previous] != -INFINITY){
					ending_forward_prob = addLog((*dbl_forward_score)[seq_size-1][st_previous] + (*hmm)[st_previous]->getEndTrans(), ending_forward_prob);
				}
			}
		}
		return;
	}
	
	
		
//	//Fix: Need to fix so it calcuates value in double and stores in float
//	void trellis::forward(){
//		
//		//Initialize forward score table
//		forward_score = new (std::nothrow) float_2D(seq_size, std::vector<float>(state_size,-INFINITY));
//		
//		if (forward_score == NULL){
//			std::cerr << "Can't allocate forward score table. OUT OF MEMORY" << std::endl;
//			exit(2);
//		}
//		
//        std::bitset<STATE_MAX> next_states;
//        std::bitset<STATE_MAX> current_states;
//		
//        double  forward_temp(-INFINITY);
//        double  emission(-INFINITY);
//        bool	exDef_position(false);
//        
//        state* init = hmm->getInitial();
//		
//		std::bitset<STATE_MAX>* initial_to = hmm->getInitialTo();
//		std::bitset<STATE_MAX>* from_trans(NULL);
//		
//		
//		//		std::cout << "Position: 0" << std::endl;
//        //Calculate Viterbi from transitions from INIT (initial) state
//        for(size_t i = 0; i < state_size; ++i){
//            if ((*initial_to)[i]){  //if the bitset is set (meaning there is a transition to this state), calculate the viterbi
//				
//				forward_temp = (*hmm)[i]->get_emission_prob(*seqs,0) + getTransition(init, i, 0);
//                
//				if (forward_temp > -INFINITY){
//                    
//					(*forward_score)[0][i] = forward_temp;
//					//					std::cout << "State: " << i << "\t" << exp(forward_temp) << std::endl;
//					next_states |= (*(*hmm)[i]->getTo());
//                }
//            }
//        }
//        
//        
//        for(size_t position = 1; position < seq_size ; ++position ){
//            
//            //Swap current_states and next states sets
//			
//			current_states.reset();
//            current_states |= next_states;
//            next_states.reset();
//            
//            if (exDef_defined){
//                exDef_position = seqs->exDefDefined(position);
//            }
//			
//			//			std::cout << "\nPosition: " << position << std::endl;
//            
//            for (size_t i = 0; i < state_size; ++i){ //i is current state that emits value
//                if (!current_states[i]){
//                    continue;
//                }
//                
//                emission = (*hmm)[i]->get_emission_prob(*seqs, position);
//				
//				if (exDef_defined && exDef_position){
//                    emission += seqs->getWeight(position, i);
//                }
//                
//				from_trans = (*hmm)[i]->getFrom();
//				
//                for (size_t j = 0; j < state_size ; ++j) {  //j is previous state
//                    if (!(*from_trans)[j]){
//                        continue;
//                    }
//					
//                    if ((*forward_score)[position-1][j] != INFINITY){
//                        forward_temp = getTransition((*hmm)[j], i , position) + emission + (*forward_score)[position-1][j];
//						
//						if ((*forward_score)[position][i] == -INFINITY){
//							(*forward_score)[position][i] = forward_temp;
//						}
//						else{
//							(*forward_score)[position][i] = addLog((double)forward_temp, (double)(*forward_score)[position][i]);
//						}
//						
//						next_states |= (*(*hmm)[i]->getTo());
//                    }
//                }
//				//				std::cout << "State: " << i <<"\t" << exp((*forward_score)[position][i]) << std::endl;
//            }
//            
//        }
//        
//        for(size_t i = 0; i < state_size ;++i){
//            if ((*forward_score)[seq_size-1][i] > -INFINITY){
//                forward_temp = (*forward_score)[seq_size-1][i] + (*hmm)[i]->getEndTrans();
//                
//                if (forward_temp > -INFINITY){
//					if (ending_forward_prob == -INFINITY){
//						ending_forward_prob = forward_temp;
//					}
//					else{
//						ending_forward_prob = addLog(ending_forward_prob,forward_temp);
//					}
//                }
//            }
//        }
//		
//		//		std::cout << exp(ending_posterior) << std::endl;
//	}
	
	

}


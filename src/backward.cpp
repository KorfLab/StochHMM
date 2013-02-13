//
//  backward.cpp
//  StochHMM
//
//  Created by Paul Lott on 2/4/13.
//  Copyright (c) 2013 Korf Lab, Genome Center, UC Davis, Davis, CA. All rights reserved.
//

#include "trellis.h"

namespace StochHMM{
	
	void trellis::backward(model* h, sequences* sqs){
		//Initialize the table
		hmm = h;
		seqs = sqs;
		seq_size		= seqs->getLength();
		state_size		= hmm->state_size();
		exDef_defined	= seqs->exDefDefined();
		
        backward();
	}
	
	
	void trellis::backward(){
		if (hmm->isBasic()){
			simple_backward();
		}
	}
	
	
	//! Naive Backward
	//! Implements the backward algorithm (simple coded, little or no optimizations)
	//! Stores the scores as doubles in table. Scoring table accessible from trellis -> getNaiveBackward();
	void trellis::naive_backward(){
		dbl_backward_score = new double_2D(seq_size, std::vector<double>(state_size,-INFINITY));
		
		double emission(-INFINITY);
		double backward_temp(-INFINITY);
		double trans(-INFINITY);
		double previous(-INFINITY);
		
		//Initialize Backward score with transition from Ending Cell
		for(size_t st_current = 0 ; st_current < state_size; st_current++){
			trans = (*hmm)[st_current]->getEndTrans();
			if (trans != -INFINITY){
				(*dbl_backward_score)[seq_size-1][st_current] = trans;
			}
		}
		
		
		typedef std::numeric_limits< double > dbl;
		std::cout << std::setprecision(dbl::digits10);
		
		for (size_t position = seq_size-2; position != SIZE_T_MAX ; position--){
			
			for(size_t st_current = 0; st_current < state_size ; st_current++){
				
				
				for (size_t st_previous = 0; st_previous < state_size; st_previous++){
					
					
					previous = (*dbl_backward_score)[position+1][st_previous];
					
					if (previous == -INFINITY){
						continue;
					}
					
					trans = getTransition(hmm->getState(st_current), st_previous, position);
					emission = (*hmm)[st_previous]->get_emission_prob(*seqs, position+1);
					
					
					if (trans == -INFINITY || emission == -INFINITY){
						continue;
					}
					
					backward_temp = previous + emission + trans;
					
					if ((*dbl_backward_score)[position][st_current] == -INFINITY){
						(*dbl_backward_score)[position][st_current] = backward_temp;
					}
					else{
						(*dbl_backward_score)[position][st_current] = addLog( backward_temp, (*dbl_backward_score)[position][st_current]);
					}
				}
			}
		}
		
		
		//Calculate Final Probability
		ending_backward_prob = -INFINITY;
		state* init = hmm->getInitial();
		for (size_t st_current=0; st_current < state_size; st_current++ ){
			if ((*dbl_backward_score)[0][st_current] == -INFINITY){
				continue;
			}
			
			previous = (*dbl_backward_score)[0][st_current];
			
			trans = getTransition(init, st_current, 0);
			
			emission = (*hmm)[st_current]->get_emission_prob(*seqs, 0);
			
			if (trans == -INFINITY || emission == -INFINITY){
				continue;
			}
			backward_temp = previous + emission + trans;
			
			if (ending_backward_prob == -INFINITY){
				ending_backward_prob = backward_temp;
			}
			else{
				ending_backward_prob = addLog( backward_temp, ending_backward_prob);
			}
		}
	}
	
	
	void trellis::simple_backward(model* h, sequences* sqs){
		hmm = h;
		seqs = sqs;
		seq_size		= seqs->getLength();
		state_size		= hmm->state_size();
		exDef_defined	= seqs->exDefDefined();
		
        simple_backward();
	}
	
	
	//Performs the backward algorithm using the model
	void trellis::simple_backward(){
		
		if (!hmm->isBasic()){
			std::cerr << "Model isn't a simple/basic HMM.  Use complex algorithms\n";
			return;
		}
		
		//Allocate backward score table
		backward_score = new (std::nothrow) float_2D(seq_size, std::vector<float>(state_size,-INFINITY));
		
		//Allocate scoring vectors
		scoring_previous = new (std::nothrow) std::vector<double> (state_size,-INFINITY);
        scoring_current  = new (std::nothrow) std::vector<double> (state_size,-INFINITY);
		
		if (scoring_previous == NULL || scoring_current == NULL || backward_score == NULL){
			std::cerr << "Can't allocate Backward score table. OUT OF MEMORY" << std::endl;
			exit(2);
		}
		
		std::bitset<STATE_MAX> next_states;
		std::bitset<STATE_MAX> current_states;
		
		double  backward_temp(-INFINITY);
		double  emission(-INFINITY);
		bool	exDef_position(false);
		
		std::bitset<STATE_MAX>* ending_from = hmm->getEndingFrom();
		std::bitset<STATE_MAX>* from_trans(NULL);
		
		
		//Calculate initial Backward from ending state
		for(size_t st_current = 0; st_current < state_size; ++st_current){
			if ((*ending_from)[st_current]){  //if the bitset is set (meaning there is a transition to this state)
				
				backward_temp = (*hmm)[st_current]->getEndTrans();
				
				if (backward_temp > -INFINITY){
					(*backward_score)[seq_size-1][st_current] = backward_temp;
					(*scoring_current)[st_current] = backward_temp;
					next_states[st_current] = 1;
				}
			}
		}
		
		
		for(size_t position = seq_size-2; position != SIZE_T_MAX ; --position ){
			
			//Swap current_states and next states sets
			current_states.reset();
			current_states |= next_states;
			next_states.reset();
			
			//Swap current and previous viterbi scores
            scoring_previous->assign(state_size,-INFINITY);
            swap_ptr = scoring_previous;
			scoring_previous = scoring_current;
			scoring_current = swap_ptr;
			
			
			if (exDef_defined){
				exDef_position = seqs->exDefDefined(position);
			}
						
			for (size_t st_previous = 0; st_previous < state_size; ++st_previous){ //i is previous state that emits value
				if (!current_states[st_previous]){
					continue;
				}
				
				emission = (*hmm)[st_previous]->get_emission_prob(*seqs, position+1);
				
				if (exDef_defined && exDef_position){
					emission += seqs->getWeight(position+1, st_previous);
				}
				
				if (emission == -INFINITY){
					continue;
				}
				
				from_trans = (*hmm)[st_previous]->getFrom();
				
				for (size_t st_current = 0; st_current < state_size ; ++st_current) {  //j is current state
					if (!(*from_trans)[st_current]){
						continue;
					}
					
					
					if ((*scoring_previous)[st_previous] != -INFINITY){
						
						backward_temp = (*scoring_previous)[st_previous] + emission +  getTransition((*hmm)[st_current], st_previous , position);
						
						if ((*scoring_current)[st_current] == -INFINITY){
							(*scoring_current)[st_current] = backward_temp;
							(*backward_score)[position][st_current] = backward_temp;
						}
						else{
							(*scoring_current)[st_current] = addLog(backward_temp, (*scoring_current)[st_current]);
							(*backward_score)[position][st_current] = (*scoring_current)[st_current];
						}
						
						next_states[st_current] = 1;
					}
				}
			}
		}
		
		ending_backward_prob = -INFINITY;
		state* init = hmm->getInitial();
		for(size_t i = 0; i < state_size ;++i){
			if ((*scoring_current)[i] != -INFINITY){
				backward_temp = (*scoring_current)[i] + (*hmm)[i]->get_emission_prob(*seqs,0) + getTransition(init, i, 0);
				if (backward_temp > -INFINITY){
					if (ending_backward_prob == -INFINITY){
						ending_backward_prob = backward_temp;
					}
					else{
						ending_backward_prob = addLog(ending_backward_prob,backward_temp);
					}
				}
			}
		}
		
		
		delete scoring_previous;
		delete scoring_current;
		scoring_previous = NULL;
		scoring_current = NULL;
		
	}
	
	
	
	
	
//	//TODO: Fix calculation in double not float (store in float)
//	void trellis::backward(){
//		//Initialize forward score table
//		backward_score = new (std::nothrow) float_2D(seq_size, std::vector<float>(state_size,-INFINITY));
//		if (backward_score == NULL){
//			std::cerr << "Can't allocate Backward score table. OUT OF MEMORY" << std::endl;
//			exit(2);
//		}
//		
//        std::bitset<STATE_MAX> next_states;
//        std::bitset<STATE_MAX> current_states;
//		
//        double  backward_temp(-INFINITY);
//        double  emission(-INFINITY);
//        bool	exDef_position(false);
//		
//		std::bitset<STATE_MAX>* ending_from = hmm->getEndingFrom();
//		std::bitset<STATE_MAX>* from_trans(NULL);
//		
//		
//		//		std::cout << "Position: 3" << std::endl;
//        //Calculate initial Backward from ending state
//        for(size_t i = 0; i < state_size; ++i){
//            if ((*ending_from)[i]){  //if the bitset is set (meaning there is a transition to this state), calculate the viterbi
//				
//				backward_temp = (*hmm)[i]->getEndTrans();
//                
//				if (backward_temp > -INFINITY){
//                    
//					(*backward_score)[seq_size-1][i] = backward_temp;
//					//					std::cout << "State: " << i << "\t" << exp(backward_temp) << std::endl;
//					next_states |= (*(*hmm)[i]->getFrom());
//                }
//            }
//        }
//        
//        
//        for(size_t position = seq_size-1; position > 0 ; --position ){
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
//                    if ((*backward_score)[position-1][j] != INFINITY){
//						
//						//						double temp_trans = getTransition((*hmm)[j], i , position-1);
//						//						double temp_score = (*backward_score)[position][i];
//						//
//						//						std::cout << "\nTransition from " << j << " to " << i << "\t" << exp(temp_trans) << std::endl;
//						//						std::cout << "Previous Score: " << exp(temp_score) << std::endl;
//						//						std::cout << "Emission: " << exp(emission) << std::endl;
//						//						backward_temp = temp_trans + emission + temp_score;
//						//						std::cout << "Temp Score: " << exp(backward_temp) << std::endl;
//						
//                        backward_temp = getTransition((*hmm)[j], i , position-1) + emission + (*backward_score)[position][i];
//						
//						if ((*backward_score)[position-1][j] == -INFINITY){
//							(*backward_score)[position-1][j] = backward_temp;
//						}
//						else{
//							(*backward_score)[position-1][j] = addLog((double)backward_temp, (double)(*backward_score)[position-1][j]);
//						}
//						
//						next_states |= (*(*hmm)[i]->getFrom());
//                    }
//                }
//				//				std::cout << "State: " << i <<"\t" << exp((*backward_score)[position][i]) << std::endl;
//            }
//            
//        }
//		
//		//		std::cout << exp((*backward_score)[0][0]) << std::endl;
//		//		std::cout << exp((*backward_score)[0][1]) << std::endl;
//		
//		double backward_posterior = -INFINITY;
//        state* init = hmm->getInitial();
//        for(size_t i = 0; i < state_size ;++i){
//            if ((*backward_score)[0][i] > -INFINITY){
//                
//				backward_temp = (*backward_score)[0][i] + (*hmm)[i]->get_emission_prob(*seqs,0) + getTransition(init, i, 0);
//                
//                if (backward_temp > -INFINITY){
//					if (backward_posterior == -INFINITY){
//						backward_posterior = backward_temp;
//					}
//					else{
//						backward_posterior = addLog(backward_posterior,backward_temp);
//					}
//                }
//            }
//        }
//		
//		//		std::cout << exp(backward_posterior) << std::endl;
//	}
//	

}
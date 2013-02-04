//
//  forward_viterbi.cpp
//  StochHMM
//
//  Created by Paul Lott on 2/4/13.
//  Copyright (c) 2013 Korf Lab, Genome Center, UC Davis, Davis, CA. All rights reserved.
//

#include "new_trellis.h"

namespace StochHMM {
	
	void trellis::forward_viterbi(){
		
		//Initialize the traceback table
        if (traceback_table != NULL){
            delete traceback_table;
        }
        
        traceback_table = new int_2D(seq_size,std::vector<int16_t> (state_size,1));
		
		//Storing Viterbi Scores
		scoring_previous = new (std::nothrow) std::vector<double> (state_size,-INFINITY);
        scoring_current  = new (std::nothrow) std::vector<double> (state_size,-INFINITY);
		
		//Storing Forward Scores
		alt_scoring_current  = new (std::nothrow) std::vector<double> (state_size,-INFINITY);
		alt_scoring_previous = new (std::nothrow) std::vector<double> (state_size,-INFINITY);
		
		
		if (scoring_previous == NULL || scoring_current == NULL || alt_scoring_current == NULL ||
			alt_scoring_previous == NULL || traceback_table == NULL){
			std::cerr << "Can't allocate Viterbi score and traceback table. OUT OF MEMORY" << std::endl;
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
        for(size_t i = 0; i < state_size; ++i){
            if ((*initial_to)[i]){  //if the bitset is set (meaning there is a transition to this state), calculate the viterbi
				
				forward_temp = (*hmm)[i]->get_emission_prob(*seqs,0) + getTransition(init, i, 0);
                
				if (forward_temp > -INFINITY){
                    
					(*forward_score)[0][i] = forward_temp;
					//					std::cout << "State: " << i << "\t" << exp(forward_temp) << std::endl;
					next_states |= (*(*hmm)[i]->getTo());
                }
            }
        }
		
		for(size_t position = 1; position < seq_size ; ++position ){
            
            //Swap current_states and next states sets
			
			current_states.reset();
            current_states |= next_states;
            next_states.reset();
            
            if (exDef_defined){
                exDef_position = seqs->exDefDefined(position);
            }
			
			//			std::cout << "\nPosition: " << position << std::endl;
            
            for (size_t i = 0; i < state_size; ++i){ //i is current state that emits value
                if (!current_states[i]){
                    continue;
                }
                
                emission = (*hmm)[i]->get_emission_prob(*seqs, position);
				
				if (exDef_defined && exDef_position){
                    emission += seqs->getWeight(position, i);
                }
                
				from_trans = (*hmm)[i]->getFrom();
				
                for (size_t j = 0; j < state_size ; ++j) {  //j is previous state
                    if (!(*from_trans)[j]){
                        continue;
                    }
					
                    if ((*forward_score)[position-1][j] != INFINITY){
                        forward_temp = getTransition((*hmm)[j], i , position) + emission + (*forward_score)[position-1][j];
						
						if ((*forward_score)[position][i] == -INFINITY){
							(*forward_score)[position][i] = forward_temp;
						}
						else{
							(*forward_score)[position][i] = addLog((double)forward_temp, (double)(*forward_score)[position][i]);
						}
						
						next_states |= (*(*hmm)[i]->getTo());
                    }
                }
				//				std::cout << "State: " << i <<"\t" << exp((*forward_score)[position][i]) << std::endl;
            }
            
        }
        
        for(size_t i = 0; i < state_size ;++i){
            if ((*forward_score)[seq_size-1][i] > -INFINITY){
                forward_temp = (*forward_score)[seq_size-1][i] + (*hmm)[i]->getEndTrans();
                
                if (forward_temp > -INFINITY){
					if (ending_posterior == -INFINITY){
						ending_posterior = forward_temp;
					}
					else{
						ending_posterior = addLog(ending_posterior,forward_temp);
					}
                }
            }
        }
		//		std::cout << exp(ending_posterior) << std::endl;
	}
	
	
	
	void trellis::simple_forward_viterbi(model *h, sequences *sqs){
		
	}
	
	void trellis::simple_forward_viterbi(){
		
		if (!hmm->isBasic()){
			std::cerr << "Model isn't a simple/basic HMM.  Use complex algorithms\n";
			return;
		}
		
		//Initialize the traceback table
        if (traceback_table != NULL || forward_score != NULL){
            delete traceback_table;
			delete forward_score;
			traceback_table = NULL;
			forward_score = NULL;
        }
        
        traceback_table = new int_2D(seq_size,std::vector<int16_t> (state_size,-1));
		forward_score	= new (std::nothrow) float_2D(seq_size, std::vector<float>(state_size,-INFINITY));
		scoring_current = new (std::nothrow) std::vector<double> (state_size,-INFINITY);
		scoring_previous= new (std::nothrow) std::vector<double> (state_size,-INFINITY);
		alt_scoring_current = new (std::nothrow) std::vector<double> (state_size,-INFINITY);
		alt_scoring_previous= new (std::nothrow) std::vector<double> (state_size,-INFINITY);
		
		if (scoring_current == NULL || scoring_previous == NULL ||
			alt_scoring_current == NULL || alt_scoring_previous == NULL ||
			traceback_table == NULL || forward_score == NULL){
			std::cerr << "Can't allocate forward score table. OUT OF MEMORY" << std::endl;
			exit(2);
		}
		
        std::bitset<STATE_MAX> next_states;
        std::bitset<STATE_MAX> current_states;
		
        double  forward_temp(-INFINITY);
		double  viterbi_temp(-INFINITY);
        double  emission(-INFINITY);
        bool	exDef_position(false);
        
        state* init = hmm->getInitial();
		
		std::bitset<STATE_MAX>* initial_to = hmm->getInitialTo();
		std::bitset<STATE_MAX>* from_trans(NULL);
		
		
        //Calculate Viterbi and Forward from transitions from INIT (initial) state
        for(size_t i = 0; i < state_size; ++i){
            if ((*initial_to)[i]){  //if the bitset is set (meaning there is a transition to this state), calculate the viterbi
				
				forward_temp = (*hmm)[i]->get_emission_prob(*seqs,0) + getTransition(init, i, 0);
				
                
				if (forward_temp > -INFINITY){
                    
					(*forward_score)[0][i] = forward_temp;
					(*alt_scoring_current)[i] = forward_temp;
					(*scoring_current)[i] = forward_temp;
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
			
			//Swap current and previous viterbi scores
            alt_scoring_previous->assign(state_size,-INFINITY);
            swap_ptr = alt_scoring_previous;
			alt_scoring_previous = alt_scoring_current;
			alt_scoring_current = swap_ptr;
			
			//Swap current_states and next states sets
			current_states.reset();
            current_states |= next_states;
            next_states.reset();
            
            if (exDef_defined){
                exDef_position = seqs->exDefDefined(position);
            }
			
            //Current State
            for (size_t current = 0; current < state_size; ++current){ //current state that emits value
				
                if (!current_states[current]){
                    continue;
                }
                
                emission = (*hmm)[current]->get_emission_prob(*seqs, position);
				
				if (exDef_defined && exDef_position){
                    emission += seqs->getWeight(position, current);
                }
                
				from_trans = (*hmm)[current]->getFrom();
				
                for (size_t previous = 0; previous < state_size ; ++previous) {  //j is previous state
					
					//Check that previous state has transition to current state
					//and that the previous viterbi score is not -INFINITY
					if ((*from_trans)[previous] && (*scoring_previous)[previous] != -INFINITY){
                        
						viterbi_temp = getTransition((*hmm)[previous], current , position) + emission;
						
						forward_temp = viterbi_temp + (*alt_scoring_previous)[previous];
						viterbi_temp += (*scoring_previous)[previous];
						
						if (viterbi_temp > (*scoring_current)[current]){
							(*scoring_current)[current] = viterbi_temp;
							(*traceback_table)[position][current] = previous;
						}
						
						if ((*alt_scoring_current)[current] == -INFINITY){
							(*alt_scoring_current)[current] = forward_temp;
							(*forward_score)[position][current] = forward_temp;
						}
						else{
							(*alt_scoring_current)[current] = addLog(forward_temp, (*alt_scoring_current)[current]);
							(*forward_score)[position][current] = (*alt_scoring_current)[current];						}
						
						next_states |= (*(*hmm)[current]->getTo());
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
		
		//Swap current and previous scores
        alt_scoring_previous->assign(state_size,-INFINITY);
        swap_ptr = alt_scoring_previous;
        alt_scoring_previous = alt_scoring_current;
        alt_scoring_current = swap_ptr;
		
		
        ending_posterior = -INFINITY;
        for(size_t i = 0; i < state_size ;++i){
			//Ending Forward Transition
            if ((*alt_scoring_previous)[i] != -INFINITY){
                forward_temp = (*alt_scoring_previous)[i] + (*hmm)[i]->getEndTrans();
                
                if (forward_temp > -INFINITY){
					if (ending_posterior == -INFINITY){
						ending_posterior = forward_temp;
					}
					else{
						ending_posterior = addLog(ending_posterior,forward_temp);
					}
                }
            }
			
			
			//Ending Viterbi Transition
			if ((*scoring_previous)[i] > -INFINITY){
				viterbi_temp = (*scoring_previous)[i] + (*hmm)[i]->getEndTrans();
				if (viterbi_temp > ending_viterbi_score){
                    ending_viterbi_score = viterbi_temp;
                    ending_viterbi_tb = i;
                }
			}
			
        }
		
		//		std::cout << exp(ending_posterior) << std::endl;
		delete scoring_previous;
		delete scoring_current;
		delete alt_scoring_current;
		delete alt_scoring_previous;
		scoring_previous = NULL;
		scoring_current = NULL;
		alt_scoring_current = NULL;
		alt_scoring_previous = NULL;
	}

}

//
//  posterior.cpp
//  StochHMM
//
//  Created by Paul Lott on 2/4/13.
//  Copyright (c) 2013 Korf Lab, Genome Center, UC Davis, Davis, CA. All rights reserved.
//

#include "new_trellis.h"

namespace StochHMM {
	
	void trellis::posterior(){
		if (hmm->isBasic()){
			simple_posterior();
		}
	}
	
	void trellis::posterior(model *h, sequences* sqs){
		hmm = h;
		seqs = sqs;
		seq_size		= seqs->getLength();
		state_size		= hmm->state_size();
		exDef_defined	= seqs->exDefDefined();
		
		if (posterior_score!=NULL){
			delete posterior_score;
			posterior_score = NULL;
		}
		ending_backward_prob = -INFINITY;
		ending_forward_prob  = -INFINITY;
		
		posterior();
	}
	
	void trellis::simple_posterior(model* h, sequences* sqs){
		hmm = h;
		seqs = sqs;
		seq_size		= seqs->getLength();
		state_size		= hmm->state_size();
		exDef_defined	= seqs->exDefDefined();
		
		if (posterior_score!=NULL){
			delete posterior_score;
			posterior_score = NULL;
		}
		ending_backward_prob = -INFINITY;
		ending_forward_prob  = -INFINITY;
		
		if (hmm->isBasic()){
			simple_posterior();
		}		
	}

	
	void trellis::simple_posterior(){
		
		if (!hmm->isBasic()){
			std::cerr << "Model isn't a simple/basic HMM.  Use complex algorithms\n";
			return;
		}
		
		posterior_score = new (std::nothrow) float_2D(seq_size, std::vector<float>(state_size,-INFINITY));
		scoring_current = new (std::nothrow) std::vector<double> (state_size,-INFINITY);
		scoring_previous= new (std::nothrow) std::vector<double> (state_size,-INFINITY);
		
		if (scoring_current == NULL || scoring_previous == NULL || posterior_score == NULL){
			std::cerr << "Can't allocate Posterior score table. OUT OF MEMORY" << std::endl;
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
		
		
        //Calculate Forward from transitions from INIT (initial) state
        for(size_t i = 0; i < state_size; ++i){
            if ((*initial_to)[i]){  //if the bitset is set (meaning there is a transition to this state), calculate the viterbi
				forward_temp = (*hmm)[i]->get_emission_prob(*seqs,0) +  getTransition(init, i, 0);
				
				if (forward_temp > -INFINITY){
					(*scoring_current)[i] = forward_temp;
					next_states |= (*(*hmm)[i]->getTo());
                }
            }
        }
        
        
        for(size_t position = 1; position < seq_size ; ++position ){
			(*posterior_score)[position-1].assign(scoring_current->begin(), scoring_current->end());
            
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
			            
            for (size_t current = 0; current < state_size; ++current){ //i is current state that emits value
                if (!current_states[current]){
                    continue;
                }
                
                emission = (*hmm)[current]->get_emission_prob(*seqs, position);
				
				if (exDef_defined && exDef_position){
                    emission += seqs->getWeight(position, current);
                }
				
				if (emission == -INFINITY){
					continue;
				}
                
				from_trans = (*hmm)[current]->getFrom();
				
                for (size_t previous = 0; previous < state_size ; ++previous) {  //j is previous state
                    if (!(*from_trans)[previous]){
                        continue;
                    }
					
					
					
					if ((*scoring_previous)[previous] != -INFINITY){
                        forward_temp = (*scoring_previous)[previous] + emission + getTransition((*hmm)[previous], current , position);
						
						if ((*scoring_current)[current] == -INFINITY){
							(*scoring_current)[current] = forward_temp;
						}
						else{
							(*scoring_current)[current] = addLog(forward_temp, (*scoring_current)[current]);
						}
						
						next_states |= (*(*hmm)[current]->getTo());
                    }
                }
            }
		}
		
		(*posterior_score)[seq_size-1].assign(scoring_current->begin(), scoring_current->end());
		
        //Swap current and previous scores
        scoring_previous->assign(state_size,-INFINITY);
        swap_ptr = scoring_previous;
        scoring_previous = scoring_current;
        scoring_current = swap_ptr;
		
		
        ending_forward_prob = -INFINITY;
        for(size_t i = 0; i < state_size ;++i){
            if ((*scoring_previous)[i] != -INFINITY){
                forward_temp = (*scoring_previous)[i] + (*hmm)[i]->getEndTrans();
				
                if (forward_temp > -INFINITY){
					if (ending_forward_prob == -INFINITY){
						ending_forward_prob = forward_temp;
					}
					else{
						ending_forward_prob = addLog(forward_temp, ending_forward_prob);
					}
                }
            }
        }

		// Perform Backward Algorithm
		double  backward_temp(-INFINITY);		
		scoring_previous->assign(state_size, -INFINITY);
		scoring_current->assign(state_size, -INFINITY);
		std::vector<float> posterior_sum(seq_size,-INFINITY);

		
		std::bitset<STATE_MAX>* ending_from = hmm->getEndingFrom();

		//Calculate initial Backward from ending state
		for(size_t st_current = 0; st_current < state_size; ++st_current){
			if ((*ending_from)[st_current]){  //if the bitset is set (meaning there is a transition to this state)

				backward_temp = (*hmm)[st_current]->getEndTrans();

				if (backward_temp > -INFINITY){
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

			for (size_t i=0;i<state_size;++i){
				if ((*posterior_score)[position+1][i] != -INFINITY && (*scoring_current)[i]!= -INFINITY){
					(*posterior_score)[position+1][i] = ((double)(*posterior_score)[position+1][i] + (double)(*scoring_current)[i]) - ending_forward_prob;
					if ((*posterior_score)[position+1][i] > -7.6009){  //Above significant value;
						posterior_sum[position+1] = addLog(posterior_sum[position+1], (*posterior_score)[position+1][i]);
					}
				}
			}

			//Swap current and previous viterbi scores
            scoring_previous->assign(state_size,-INFINITY);
            swap_ptr = scoring_previous;
			scoring_previous = scoring_current;
			scoring_current = swap_ptr;



			if (exDef_defined){
				exDef_position = seqs->exDefDefined(position);
			}


			for (size_t st_previous	= 0; st_previous < state_size; ++st_previous){ //i is current state that emits value
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

				for (size_t st_current = 0; st_current < state_size ; ++st_current) {  //j is previous state
					if (!(*from_trans)[st_current]){
						continue;
					}

					if ((*scoring_previous)[st_previous] != -INFINITY){
						
						backward_temp =(*scoring_previous)[st_previous] + emission + getTransition((*hmm)[st_current], st_previous , position);

						if ((*scoring_current)[st_current] == -INFINITY){
							(*scoring_current)[st_current] = backward_temp;
						}
						else{
							(*scoring_current)[st_current] = addLog(backward_temp, (*scoring_current)[st_current]);
						}

						next_states[st_current] = 1;
					}
				}
			}
		}

		for (size_t i=0;i<state_size;++i){
			(*posterior_score)[0][i] = ((double)(*posterior_score)[0][i] + (double)(*scoring_current)[i]) - ending_forward_prob;
			posterior_sum[0] = addLog(posterior_sum[0], (*posterior_score)[0][i]);
		}
		


		ending_backward_prob = -INFINITY;
		init = hmm->getInitial();
		for(size_t i = 0; i < state_size ;++i){

			if ((*scoring_current)[i] != -INFINITY){
				backward_temp = (*scoring_current)[i] + (*hmm)[i]->get_emission_prob(*seqs,0) + getTransition(init, i, 0);
				
				if (backward_temp > -INFINITY){
					if (ending_backward_prob == -INFINITY){
						ending_backward_prob = backward_temp;
					}
					else{
						ending_backward_prob = addLog(backward_temp, ending_backward_prob);
					}
				}
			}
		}
		
		if (abs(ending_backward_prob - ending_forward_prob) > 0.0000001){
			std::cerr << "Ending sequence probabilities calculated by Forward and Backward algorithm are different.  They should be the same.\t" << __FUNCTION__ << std::endl;
		}
		
		
		
		for (size_t i=0;i<seq_size;i++){
			for(size_t j=0;j<state_size;j++){
				if ((*posterior_score)[i][j] == -INFINITY){
					continue;
				}
//				std::cout << "Sum:\t" << exp(posterior_sum[i]) << std::endl;
//				std::cout << "Score:\t" << exp((*posterior_score)[i][j]) << std::endl;
				
				if ((*posterior_score)[i][j] > -7.6009){  //Above significant value;
					(*posterior_score)[i][j] -= posterior_sum[i];
				}
				else{
					(*posterior_score)[i][j] = -INFINITY;
				}
			}
		}
		
		
		delete scoring_previous;
		delete scoring_current;
		scoring_previous = NULL;
		scoring_current = NULL;
	}
	
	
	void trellis::traceback_posterior(traceback_path& path){
		if (posterior_score == NULL){
			std::cerr << __FUNCTION__ << " called before trellis::posterior was completed\n";
			exit(2);
		}
		
		double max(-INFINITY);
		int max_ptr(SIZE_T_MAX);
		for(size_t position=seq_size-1; position > 0 ;position--){
			max = -INFINITY;
			max_ptr = SIZE_T_MAX ;
			
            for (size_t st=0; st < state_size; st++){
                if ((*posterior_score)[position][st] > max){
					max = (*posterior_score)[position][st];
					max_ptr = st;
				}
                
			}
            path.push_back(max_ptr);
        }
		return;
	}
	
	void trellis::traceback_stoch_posterior(traceback_path& path){
		for (size_t position =seq_size-1; position != SIZE_MAX; --position){
            double random=((double)rand()/((double)(RAND_MAX)+(double)(1)));
            double cumulative_prob(0.0);
            
			for (size_t st = 0; st < state_size ; ++st){
                cumulative_prob+=exp((*posterior_score)[position][st]);
                if (random < cumulative_prob){
                    path.push_back(st);
                }
            }
        }
		return;
	}
	
	void trellis::traceback_stoch_posterior(multiTraceback& paths, size_t reps){
		for (size_t i = 0; i < reps; i++){
			traceback_path path(hmm);
			traceback_stoch_posterior(path);
			paths.assign(path);
		}
		return;
	}
	
	//	void trellis::simple_posterior(){
	//
	//		if (!hmm->isBasic()){
	//			std::cerr << "Model isn't a simple/basic HMM.  Use complex algorithms\n";
	//			return;
	//		}
	//
	//		if (forward_score == NULL){
	//			simple_forward();
	//		}
	//
	//		if (backward_score == NULL){
	//			simple_backward();
	//		}
	//
	//		if (abs(ending_backward_prob - ending_forward_prob) > 0.0000001){
	//			std::cerr << "Ending Forward and Backward Probabilities are not equal.\n Check the model.\t" << __FUNCTION__ << std::endl;
	//			exit(2);
	//		}
	//
	//		posterior_score = new(std::nothrow) float_2D(seq_size, std::vector<float>(state_size,-INFINITY));
	//
	//		for(size_t position = 0;position < seq_size; ++position){
	//			for(size_t state = 0; state < state_size; ++state){
	//				(*posterior_score)[position][state] = ((*forward_score)[position][state] + (*backward_score)[position][state]) - ending_forward_prob;
	//			}
	//		}
	//		
	//		return;
	//	}
	
	
}
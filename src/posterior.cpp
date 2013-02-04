//
//  posterior.cpp
//  StochHMM
//
//  Created by Paul Lott on 2/4/13.
//  Copyright (c) 2013 Korf Lab, Genome Center, UC Davis, Davis, CA. All rights reserved.
//

#include "new_trellis.h"

namespace StochHMM {
	
	void trellis::simple_posterior(model* h, sequences* sqs){
		hmm = h;
		seqs = sqs;
		seq_size		= seqs->getLength();
		state_size		= hmm->state_size();
		exDef_defined	= seqs->exDefDefined();
		
        simple_posterior();
	}
	
	void trellis::simple_posterior(){
		
		if (!hmm->isBasic()){
			std::cerr << "Model isn't a simple/basic HMM.  Use complex algorithms\n";
			return;
		}
		
		if (forward_score == NULL){
			simple_forward();
		}
		
		if (backward_score == NULL){
			simple_backward();
		}
		
		posterior_score = new(std::nothrow) float_2D(seq_size, std::vector<float>(state_size,-INFINITY));
		
		for(size_t position = 0;position < seq_size; ++position){
			for(size_t state = 0; state < state_size; ++state){
				(*posterior_score)[position][state] = ((*forward_score)[position][state] + (*backward_score)[position][state]) - ending_posterior;
			}
		}
		
		//		for(size_t position = 0;position < seq_size; ++position){
		//			for(size_t state = 0; state < state_size; ++state){
		//				std::cout << exp((*posterior_score)[position][state]);
		//			}
		//			std::cout << std::endl;
		//		}
		
		//		std::cout << "here" << std::endl;
		
		return;
	}
	
	
	void trellis::simple_posterior_second(){
		
		typedef std::numeric_limits< double > dbl;
		
		std::cout.precision(dbl::digits10);
		
		if (!hmm->isBasic()){
			std::cerr << "Model isn't a simple/basic HMM.  Use complex algorithms\n";
			return;
		}
		
		posterior2Score = new (std::nothrow) std::vector<std::vector<double> >(seq_size, std::vector<double>(state_size,-INFINITY));
		scoring_current = new (std::nothrow) std::vector<double> (state_size,-INFINITY);
		scoring_previous= new (std::nothrow) std::vector<double> (state_size,-INFINITY);
		
		if (scoring_current == NULL || scoring_previous == NULL || posterior2Score == NULL){
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
		
		
		//		std::cout << "Position: 0" << std::endl;
        //Calculate Forward from transitions from INIT (initial) state
        for(size_t i = 0; i < state_size; ++i){
            if ((*initial_to)[i]){  //if the bitset is set (meaning there is a transition to this state), calculate the viterbi
				std::cout << "0\t-1\t" << i << "\t" << (*hmm)[i]->get_emission_prob(*seqs,0) << "\t" << getTransition(init, i, 0) << std::endl;
				forward_temp = (*hmm)[i]->get_emission_prob(*seqs,0) +  getTransition(init, i, 0);
				
				if (forward_temp > -INFINITY){
					(*scoring_current)[i] = forward_temp;
					next_states |= (*(*hmm)[i]->getTo());
                }
            }
        }
        
        
        for(size_t position = 1; position < seq_size ; ++position ){
			(*posterior2Score)[position-1].assign(scoring_current->begin(), scoring_current->end());
            
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
						
						std::cout <<position << "\t" << previous << "\t" << current << "\t" <<(*hmm)[current]->get_emission_prob(*seqs, position) << "\t" <<getTransition((*hmm)[previous], current , position) <<"\t" << (*scoring_previous)[previous] << std::endl;
						
                        forward_temp = getTransition((*hmm)[previous], current , position) + emission + (*scoring_previous)[previous];
						
						if ((*scoring_current)[current] == -INFINITY){
							(*scoring_current)[current] = forward_temp;
						}
						else{
							(*scoring_current)[current] = addLog(forward_temp, (*scoring_current)[current]);
						}
						
						next_states |= (*(*hmm)[current]->getTo());
                    }
                }
				//				std::cout << "State: " << current <<"\t" << exp((*forward_score)[position][current]) << std::endl;
            }
		}
		
		(*posterior2Score)[seq_size-1].assign(scoring_current->begin(), scoring_current->end());
		
        //Swap current and previous scores
        scoring_previous->assign(state_size,-INFINITY);
        swap_ptr = scoring_previous;
        scoring_previous = scoring_current;
        scoring_current = swap_ptr;
		
		
        ending_posterior = -INFINITY;
		
        for(size_t i = 0; i < state_size ;++i){
            if ((*scoring_previous)[i] != -INFINITY){
                forward_temp = (*scoring_previous)[i] + (*hmm)[i]->getEndTrans();
				
				std::cout << "Ending\t" << i << "\t" << (*hmm)[i]->getEndTrans() << "\t" << (*scoring_previous)[i] << std::endl;
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
		
		for(size_t i=0;i<seq_size;i++){
			std::cout << i+1 << "\t";
			for (size_t j=0; j<state_size; j++){
				std::cout << (*posterior2Score)[i][j] << "\t";
				
			}
			std::cout << std::endl;
		}
		
		std::cout << "FORWARD  PROB:\t" << ending_posterior << std::endl;
		
		//		std::cout << std::endl << std::endl;
		//
		//
		//
		//		//		std::cout << exp(ending_posterior) << std::endl;
		//		scoring_previous->assign(state_size, -INFINITY);
		//		scoring_current->assign(state_size, -INFINITY);
		//		Lscoring_previous->assign(state_size, -INFINITY);
		//		Lscoring_current->assign(state_size, -INFINITY);
		//
		////		for (size_t i=0; i < seq_size; i++){
		////			std::cout << i+1 << "\t";
		////			for(size_t j=0; j < state_size; j++){
		////				std::cout << exp((*posterior2Score)[i][j]) << "\t";
		////			}
		////			std::cout << std::endl;
		////		}
		//
		//
		//		//Backward
		//
		//		double  backward_temp(-INFINITY);
		//		long double Lbackward_temp(-INFINITY);
		//
		//		std::bitset<STATE_MAX>* ending_from = hmm->getEndingFrom();
		//
		//		//Calculate initial Backward from ending state
		//		for(size_t i = 0; i < state_size; ++i){
		//			if ((*ending_from)[i]){  //if the bitset is set (meaning there is a transition to this state), calculate the viterbi
		//
		//				backward_temp = (*hmm)[i]->getEndTrans();
		//				Lbackward_temp = (*hmm)[i]->getEndTrans();
		//
		//				if (backward_temp > -INFINITY){
		//					(*scoring_current)[i] = backward_temp;
		//					(*Lscoring_current)[i] = Lbackward_temp;
		//					next_states[i] = 1;
		//				}
		//			}
		//		}
		//
		//
		//
		//		for(size_t position = seq_size-1; position > 0 ; --position ){
		//
		//			//Swap current_states and next states sets
		//			current_states.reset();
		//			current_states |= next_states;
		//			next_states.reset();
		//
		//			for (size_t i=0;i<state_size;++i){
		////				std::cout << "Backward:\t" << (*scoring_current)[state_size] << std::endl;
		////				std::cout << "Forward:\t"  << (*posterior2Score)[position][i] << std::endl;
		//
		//				(*posterior2Score)[position][i] += (*scoring_current)[i];
		//				Lposterior2Score[position][i] += (*Lscoring_current)[i];
		//
		//				(*posterior2Score)[position][i] -= ending_posterior;
		//				Lposterior2Score[position][i] -= Lending_posterior;
		////				\std::cout << "Posterior:\t" << exp((*posterior2Score)[position][i]) << std::endl;
		//			}
		////			(*backward2Score)[position] = (*scoring_current);
		//
		//			//Swap current and previous viterbi scores
		//            scoring_previous->assign(state_size,-INFINITY);
		//            swap_ptr = scoring_previous;
		//			scoring_previous = scoring_current;
		//			scoring_current = swap_ptr;
		//
		//			Lscoring_previous->assign(state_size,-INFINITY);
		//            Lptr = Lscoring_previous;
		//			Lscoring_previous = Lscoring_current;
		//			Lscoring_current = Lptr;
		//
		//
		//			if (exDef_defined){
		//				exDef_position = seqs->exDefDefined(position);
		//			}
		//
		//			//			std::cout << "\nPosition: " << position << std::endl;
		//
		//			for (size_t i = 0; i < state_size; ++i){ //i is current state that emits value
		//				if (!current_states[i]){
		//					continue;
		//				}
		//
		//				emission = (*hmm)[i]->get_emission_prob(*seqs, position);
		//				Lemission = (*hmm)[i]->get_emission_prob(*seqs, position);
		//
		//				if (exDef_defined && exDef_position){
		//					emission += seqs->getWeight(position, i);
		//				}
		//
		//				if (emission == -INFINITY){
		//					continue;
		//				}
		//
		//				from_trans = (*hmm)[i]->getFrom();
		//
		//				for (size_t j = 0; j < state_size ; ++j) {  //j is previous state
		//					if (!(*from_trans)[j]){
		//						continue;
		//					}
		//
		//					//if ((*backward_score)[position-1][j] != -INFINITY){
		//					if ((*scoring_previous)[i] != -INFINITY){
		//
		//						std::cout << position << "\t" << j <<"\t"<< i << "\t" << emission << "\t" << getTransition((*hmm)[j], i , position-1) <<"\t" << (*scoring_previous)[i] << std::endl;
		//
		//						//backward_temp = getTransition((*hmm)[j], i , position-1) + emission + (*backward_score)[position][i];
		//						backward_temp = getTransition((*hmm)[j], i , position-1) + emission + (*scoring_previous)[i];
		//						Lbackward_temp = getTransition((*hmm)[j], i , position-1) + Lemission + (*Lscoring_previous)[i];
		//
		//						if ((*scoring_current)[j] == -INFINITY){
		//							(*scoring_current)[j] = backward_temp;
		//							(*Lscoring_current)[j] = Lbackward_temp;
		//						}
		//						else{
		//							(*scoring_current)[j] = addLog(backward_temp, (*scoring_current)[j]);
		//							(*Lscoring_current)[j] = addLog(Lbackward_temp, (*Lscoring_current)[j]);
		//						}
		//
		//						next_states |= (*(*hmm)[i]->getFrom());
		//					}
		//				}
		//				//				std::cout << "State: " << i <<"\t" << exp((*backward_score)[position][i]) << std::endl;
		//			}
		//
		//		}
		//
		//		for (size_t i=0;i<state_size;++i){
		//			(*posterior2Score)[0][i] += (*scoring_current)[i];
		//			(*posterior2Score)[0][i] -= ending_posterior;
		//
		//			Lposterior2Score[0][i] += (*Lscoring_current)[i];
		//			Lposterior2Score[0][i] -= Lending_posterior;
		//
		//		}
		//
		////		(*backward2Score)[0] = (*scoring_current);
		//
		////
		////		for (size_t i=0; i < seq_size; i++){
		////			std::cout << i+1 << "\t";
		////			for(size_t j=0; j < state_size; j++){
		////				std::cout << exp((*backward2Score)[i][j]) << "\t";
		////			}
		////			std::cout << std::endl;
		////		}
		//
		//		double backward_posterior = -INFINITY;
		//		long double Lbackward_posterior = -INFINITY;
		//		init = hmm->getInitial();
		//		for(size_t i = 0; i < state_size ;++i){
		//			if (init->getTrans(i) == NULL){
		//				continue;
		//			}
		//
		//			if ((*scoring_current)[i] != -INFINITY){
		//				std::cout << "-1\t" << "-1\t" << i <<"\t"<< (*scoring_current)[i] <<"\t" << (*hmm)[i]->get_emission_prob(*seqs,0) <<  "\t" << getTransition(init, i, 0) << std::endl;
		//
		//				//if ((*backward_score)[0][i] > -INFINITY){
		//
		//				//backward_temp = (*backward_score)[0][i] + (*hmm)[i]->get_emission_prob(*seqs,0) + getTransition(init, i, 0);
		//				backward_temp = (*scoring_current)[i] + (*hmm)[i]->get_emission_prob(*seqs,0) + getTransition(init, i, 0);
		//				Lbackward_temp = (*Lscoring_current)[i] + (*hmm)[i]->get_emission_prob(*seqs,0) + getTransition(init, i, 0);
		//				if (backward_temp > -INFINITY){
		//					if (backward_posterior == -INFINITY){
		//						backward_posterior = backward_temp;
		//						Lbackward_posterior = Lbackward_temp;
		//					}
		//					else{
		//						backward_posterior = addLog(backward_posterior,backward_temp);
		//						Lbackward_posterior = addLog(Lbackward_posterior,Lbackward_temp);
		//					}
		//				}
		//			}
		//		}
		//
		//		std::cout << "FORWARD  PROB:\t" << ending_posterior << std::endl;
		//		std::cout << "BACKWARD PROB:\t" << backward_posterior << std::endl;
		//		std::cout << "FORWARD  PROB:\t" << exp(ending_posterior) << std::endl;
		//		std::cout << "BACKWARD PROB:\t" << exp(backward_posterior) << std::endl;
		//
		//		typedef std::numeric_limits< long double > ldbl;
		//		
		//		std::cout.precision(ldbl::digits10);
		//		std::cout << "L_FORWARD  PROB:\t" << Lending_posterior << std::endl;
		//		std::cout << "L_BACKWARD PROB:\t" << Lbackward_posterior << std::endl;
		//
		//		
		delete scoring_previous;
		delete scoring_current;
		scoring_previous = NULL;
		scoring_current = NULL;
	}
}
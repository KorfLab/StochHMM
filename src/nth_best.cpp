//
//  nth_best.cpp
//  StochHMM
//
//  Created by Paul Lott on 2/6/13.
//  Copyright (c) 2013 Korf Lab, Genome Center, UC Davis, Davis, CA. All rights reserved.
//

#include "trellis.h"

namespace StochHMM{
	

	
	
	//!Sort the viterbi scores in the nth trellis cells
    void sort_scores(std::vector<nthScore>& nth_scores){
        sort(nth_scores.begin(), nth_scores.end(), _vec_sort );
        return;
    }
	
    //Sort vector of pairs using the first value in the pair
    bool _vec_sort(const nthScore& i, const nthScore& j){
        return (i.score > j.score);
    }

	
	//TODO:  Currently adapt naive algorithm to optimized algorithm
	
	//As written it's slower than the naive and may use as much memory....
	void trellis::simple_nth_viterbi(size_t n){
		if (!hmm->isBasic()){
			std::cerr << "Model isn't a simple/basic HMM.  Use complex algorithms\n";
			return;
		}
		
		//Initialize the traceback table
		if (nth_traceback_table != NULL){delete nth_traceback_table;}
		if (ending_nth_viterbi	!= NULL){delete ending_nth_viterbi;}
		if (nth_scoring_previous!= NULL){delete nth_scoring_previous;}
		if (nth_scoring_current	!= NULL){delete nth_scoring_current;}
		
		nth_scoring_previous = new (std::nothrow) std::vector<std::vector<nthScore> > (state_size);
		nth_scoring_current  = new (std::nothrow) std::vector<std::vector<nthScore> > (state_size);
		nth_traceback_table  = new (std::nothrow) std::vector<nthTrace>(seq_size);
		ending_nth_viterbi = new(std::nothrow) std::vector<nthScore>;
		
		if (nth_scoring_previous == NULL || nth_scoring_current == NULL || nth_traceback_table == NULL){
			std::cerr << "Can't allocate Viterbi score and traceback table. OUT OF MEMORY" << std::endl;
			exit(2);
		}
		
		
		std::bitset<STATE_MAX> next_states;
		std::bitset<STATE_MAX> current_states;
		
		double  viterbi_temp(-INFINITY);
		double  emission(-INFINITY);
		double	trans(-INFINITY);
		bool	exDef_position(false);
		
		
		state* init = hmm->getInitial();
		
		std::bitset<STATE_MAX>* initial_to = hmm->getInitialTo();
		std::bitset<STATE_MAX>* from_trans(NULL);
		
		//Calculate Viterbi from transitions from INIT (initial) state
		for(size_t st = 0; st < state_size; ++st){
			if ((*initial_to)[st]){  //if the bitset is set (meaning there is a transition to this state), calculate the viterbi
				
				viterbi_temp = (*hmm)[st]->get_emission_prob(*seqs,0) + getTransition(init, st, 0);
				
				if (viterbi_temp > -INFINITY){
					(*nth_scoring_current)[st].push_back(nthScore(-1,-1,viterbi_temp));
					
					next_states |= (*(*hmm)[st]->getTo());
				}
			}
		}
		
		//Each position in the sequence
		for(size_t position = 1; position < seq_size ; ++position ){
			
			//Swap current and previous viterbi scores
			nth_scoring_previous->assign(state_size,std::vector<nthScore>());
			nth_swap_ptr = nth_scoring_previous;
			nth_scoring_previous = nth_scoring_current;
			nth_scoring_current = nth_swap_ptr;
			
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
					if (!(*nth_scoring_previous)[st_previous].empty()){
						
						trans = getTransition((*hmm)[st_previous], st_current, position);
						if (trans== -INFINITY){
							continue;
						}
						
						for(size_t vit_previous=0; vit_previous < (*nth_scoring_previous)[st_previous].size(); vit_previous++){
							viterbi_temp = emission + trans + (*nth_scoring_previous)[st_previous][vit_previous].score;
							(*nth_scoring_current)[st_current].push_back(nthScore(st_previous,vit_previous,viterbi_temp));
						}
						
						next_states |= (*(*hmm)[st_current]->getTo());
					}
				}
				
				sort_scores((*nth_scoring_current)[st_current]);
				if ((*nth_scoring_current)[st_current].size() > nth_size){
					(*nth_scoring_current)[st_current].resize(nth_size);
				}
				
				for (size_t i=0; i < (*nth_scoring_current)[st_current].size(); i++){
					(*nth_traceback_table)[position].assign(st_current, i , (*nth_scoring_current)[st_current][i].st_tb, (*nth_scoring_current)[st_current][i].score_tb);
				}
			}
		}
		
		//TODO:  Calculate ending and set the final viterbi and traceback pointer
		//Swap current and previous viterbi scores
		nth_scoring_previous->assign(state_size,std::vector<nthScore>());
		nth_swap_ptr = nth_scoring_previous;
		nth_scoring_previous = nth_scoring_current;
		nth_scoring_current = nth_swap_ptr;
		
		
		//Calculate ending viterbi score and traceback from END state
		for(size_t st_previous = 0; st_previous < state_size ;++st_previous){
			if (!(*nth_scoring_previous)[st_previous].empty()){
				
				trans = (*hmm)[st_previous]->getEndTrans();
				
				if (trans == -INFINITY){
					continue;
				}
				
				for(size_t vit_previous=0; vit_previous < (*nth_scoring_previous)[st_previous].size(); vit_previous++){
					viterbi_temp = trans + (*nth_scoring_previous)[st_previous][vit_previous].score;
					ending_nth_viterbi->push_back(nthScore(st_previous,vit_previous,viterbi_temp));
				}
			}
		}
		
		sort_scores(*ending_nth_viterbi);
		if (ending_nth_viterbi->size() > nth_size){
			ending_nth_viterbi->resize(nth_size);
		}
		
		delete nth_scoring_previous;
		delete nth_scoring_current;
		nth_scoring_previous = NULL;
		nth_scoring_current  = NULL;
		
		std::cout << "Ending:\t" << (*ending_nth_viterbi)[0].score << std::endl;
		std::cout << "Ending:\t" << (*ending_nth_viterbi)[1].score << std::endl;
		std::cout << "Ending:\t" << (*ending_nth_viterbi)[2].score << std::endl;
	}
	
	void trellis::naive_nth_viterbi(size_t n){
		nth_size = n;
		
		if (naive_nth_scores != NULL){delete naive_nth_scores; naive_nth_scores=NULL;}
		if (ending_nth_viterbi != NULL){delete ending_nth_viterbi; ending_nth_viterbi = NULL;}
		
		naive_nth_scores = new (std::nothrow) std::vector<std::vector<std::vector<nthScore >* > >(seq_size, std::vector<std::vector<nthScore>* >(state_size,NULL));
		
		std::vector<nthScore>* temp_scores;
		
		double emission(-INFINITY);
		double viterbi_temp(-INFINITY);
		double trans(-INFINITY);
		bool	exDef_position(false);
		
		state* init = hmm->getInitial();
		
		//Calculate from Initial states
		for(size_t st = 0; st < state_size; ++st){
			viterbi_temp = (*hmm)[st]->get_emission_prob(*seqs,0) +  getTransition(init, st, 0);
			
			if (viterbi_temp == -INFINITY){
				continue;
			}
			
			temp_scores = new (std::nothrow) std::vector<nthScore>;
			temp_scores->push_back(nthScore(-1,-1,viterbi_temp));
			//nth_table[0][st]=temp_scores;
			(*naive_nth_scores)[0][st]=temp_scores;
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
				
				temp_scores = new (std::nothrow) std::vector<nthScore>;
				for (size_t st_previous = 0; st_previous < state_size; ++st_previous){
					trans = getTransition(hmm->getState(st_previous), st_current, position);
					if (trans== -INFINITY){
						continue;
					}
					
					if ((*naive_nth_scores)[position-1][st_previous] == NULL){
						continue;
					}
					
					for(size_t vit_previous=0; vit_previous < (*naive_nth_scores)[position-1][st_previous]->size(); vit_previous++){
						viterbi_temp = emission + trans + (*(*naive_nth_scores)[position-1][st_previous])[vit_previous].score;
						temp_scores->push_back(nthScore(st_previous,vit_previous,viterbi_temp));
					}
				}
				
				sort_scores(*temp_scores);
				if (temp_scores->size() > nth_size){
					temp_scores->resize(nth_size);
				}
				(*naive_nth_scores)[position][st_current]=temp_scores;				
			}
		}
		
		//Calculate Ending Transition
		temp_scores = new(std::nothrow) std::vector<nthScore>;
		for (size_t st_previous = 0; st_previous < state_size; ++st_previous){
			
			//Check validity of transition to end for the state
			if ((*hmm)[st_previous]->getEndTrans() != -INFINITY){
				
				//Check to see that previous viterbi scores are defined for the state
				if ((*naive_nth_scores)[seq_size-1][st_previous] != NULL){
					
					//Calculate viterbi score for each score in previous cell
					for (size_t prev_scores=0; prev_scores < (*naive_nth_scores)[seq_size-1][st_previous]->size(); prev_scores++){

						viterbi_temp = (*(*naive_nth_scores)[seq_size-1][st_previous])[prev_scores].score + (*hmm)[st_previous]->getEndTrans();
						temp_scores->push_back(nthScore(st_previous, prev_scores, viterbi_temp));
					}

				}
			}
		}

		sort_scores(*temp_scores);
		if (temp_scores->size() > nth_size){
			temp_scores->resize(nth_size);
		}
		
		ending_nth_viterbi = temp_scores;
		
		std::cout << "Ending:\t" << (*ending_nth_viterbi)[0].score << std::endl;
		std::cout << "Ending:\t" << (*ending_nth_viterbi)[1].score << std::endl;
		std::cout << "Ending:\t" << (*ending_nth_viterbi)[2].score << std::endl;
		
		
		return;
	}

}
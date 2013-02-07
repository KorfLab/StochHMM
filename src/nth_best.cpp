//
//  nth_best.cpp
//  StochHMM
//
//  Created by Paul Lott on 2/6/13.
//  Copyright (c) 2013 Korf Lab, Genome Center, UC Davis, Davis, CA. All rights reserved.
//

#include "new_trellis.h"

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

	
	void trellis::nth_viterbi(){
		return;
	}
	
	void trellis::naive_nth_viterbi(size_t n){
		
		naive_nth_scores = new (std::nothrow) std::vector<std::vector<std::vector<nthScore >* > >(seq_size, std::vector<std::vector<nthScore>* >(state_size,NULL));
		
		//std::vector<std::vector<std::vector<nthScore >* > > nth_table(seq_size, std::vector<std::vector<nthScore>* >(state_size,NULL));
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
					
					for(size_t vit_previous=0; vit_previous < (*naive_nth_scores)[position-1][st_previous]->size(); vit_previous++){
						viterbi_temp = emission + trans + (*(*naive_nth_scores)[position-1][st_previous])[vit_previous].score;
						temp_scores->push_back(nthScore(st_previous,vit_previous,viterbi_temp));
					}
					
//					for(size_t vit_previous=0; vit_previous < nth_table[position-1][st_previous]->size(); vit_previous++){
//						viterbi_temp = emission + trans + (*nth_table[position-1][st_previous])[vit_previous].score;
//						temp_scores->push_back(nthScore(st_previous,vit_previous,viterbi_temp));
//					}

					
				}
				
				sort_scores(*temp_scores);
				if (temp_scores->size() > n){
					temp_scores->resize(n);
				}
//				nth_table[position][st_current]=temp_scores;
				(*naive_nth_scores)[position][st_current]=temp_scores;				
			}
		}
		
		//Calculate Ending Transition
		temp_scores = new(std::nothrow) std::vector<nthScore>;
		for (size_t st_previous = 0; st_previous < state_size; ++st_previous){
			
			//Check validity of transition to end for the state
			if ((*hmm)[st_previous]->getEndTrans() != -INFINITY){
				
				//Check to see that previous viterbi scores are defined for the state
//				if (nth_table[seq_size-1][st_previous] != NULL){
				if ((*naive_nth_scores)[seq_size-1][st_previous] != NULL){

					
//					//Calculate viterbi score for each score in previous cell
//					for (size_t prev_scores=0; prev_scores < nth_table[seq_size-1][st_previous]->size(); prev_scores++){
//						
//						viterbi_temp = (*nth_table[seq_size-1][st_previous])[prev_scores].score + (*hmm)[st_previous]->getEndTrans();
//						temp_scores->push_back(nthScore(st_previous, prev_scores, viterbi_temp));
//					}
					//Calculate viterbi score for each score in previous cell
					for (size_t prev_scores=0; prev_scores < (*naive_nth_scores)[seq_size-1][st_previous]->size(); prev_scores++){

						viterbi_temp = (*(*naive_nth_scores)[seq_size-1][st_previous])[prev_scores].score + (*hmm)[st_previous]->getEndTrans();
						temp_scores->push_back(nthScore(st_previous, prev_scores, viterbi_temp));
					}

				}
			}
		}

		sort_scores(*temp_scores);
		if (temp_scores->size() > n){
			temp_scores->resize(n);
		}
		
		ending_nth_viterbi = temp_scores;
		
		std::cout << "Ending:\t" << (*ending_nth_viterbi)[0].score << std::endl;
		
		
		return;
	}

}
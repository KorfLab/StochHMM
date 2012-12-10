//
//  new_trellis.cpp
//  StochHMM
//
//  Created by Paul Lott on 11/13/12.
//  Copyright (c) 2012 Korf Lab, Genome Center, UC Davis, Davis, CA. All rights reserved.
//

#include "new_trellis.h"


namespace StochHMM {
	trellis::trellis(){
		hmm=NULL;
		seqs=NULL;
		
		type = SIMPLE;
		store_values=false;
		
		traceback_table=NULL;
		stoch_traceback=NULL;
		nth_traceback=NULL;
		viterbi_score=NULL;
		forward_score=NULL;
		backward_score=NULL;
		posterior_score=NULL;
		ending_stoch_tb=NULL;
		ending_nth_viterbi=NULL;
	}
	
	trellis::trellis(model* h, sequences* sqs){
		hmm=h;
		seqs=sqs;
		
		seq_size		= seqs->getLength();
		state_size		= hmm->state_size();
		exDef_defined	= seqs->exDefDefined();
		
		type = SIMPLE;
		store_values=false;
		
		traceback_table=NULL;
		stoch_traceback=NULL;
		nth_traceback=NULL;
		viterbi_score=NULL;
		forward_score=NULL;
		backward_score=NULL;
		posterior_score=NULL;
		ending_stoch_tb=NULL;
		ending_nth_viterbi=NULL;
	}
	
	
	
	trellis::~trellis(){
		delete traceback_table;
		delete stoch_traceback;
		delete nth_traceback;
		delete viterbi_score;
		delete forward_score;
		delete backward_score;
		delete posterior_score;
		delete ending_stoch_tb;
		delete ending_nth_viterbi;
	}
	
	void trellis::viterbi(){
		
		viterbi_previous = new std::vector<double> (state_size,-INFINITY);
        viterbi_current  = new std::vector<double> (state_size,-INFINITY);
		
		//Initialize the traceback table
        if (traceback_table != NULL){
            delete traceback_table;
        }
        
        traceback_table = new int_2D(seq_size,std::vector<uint16_t> (state_size));
		
		//std::vector<bool>* empty = new std::vector<bool>(state_size,false);
        //std::vector<bool>* next_states = new std::vector<bool>(state_size, false);
        //std::vector<bool>* current_states = new std::vector<bool>(state_size,false);
        //std::vector<bool>* temp_states;
		
        std::bitset<STATE_MAX> next_states;
        std::bitset<STATE_MAX> current_states;
        
		
        double  viterbi_temp(-INFINITY);
        double  emission(-INFINITY);
        bool	exDef_position(false);
        
        state* init = hmm->getInitial();
        //state* ending  = hmm->getEnding();
		
		std::bitset<STATE_MAX>* initial_to = hmm->getInitialTo();
		std::bitset<STATE_MAX>* from_trans(NULL);
		
        //Calculate Viterbi from transitions from INIT (initial) state
        for(size_t i = 0; i < state_size; ++i){
            if ((*initial_to)[i]){  //if the bitset is set (meaning there is a transition to this state), calculate the viterbi
				
				viterbi_temp = (*hmm)[i]->get_emission_prob(*seqs,0) + init->getTrans(i)->getTransition(0,NULL);
                
				if (viterbi_temp > -INFINITY){
                    if ((*viterbi_current)[i] < viterbi_temp){
                        (*viterbi_current)[i] = viterbi_temp;
                    }
					next_states[i] = 1;
                }
            }
        }
        
        
        for(size_t position = 1; position < seq_size ; ++position ){
            
            //Swap current and previous viterbi scores
            viterbi_previous->assign(state_size,-INFINITY);
            swap_ptr = viterbi_previous;
			viterbi_previous = viterbi_current;
			viterbi_current = swap_ptr;
            
            //Swap current_states and next states sets
			
			current_states.reset();
            current_states |= next_states;
            next_states.reset();
            
            if (exDef_defined){
                exDef_position = seqs->exDefDefined(position);
            }
            
            for (size_t i = 0; i < state_size; ++i){ //i is current state that emits value
                if (!current_states[i]){
                    continue;
                }
                
                //current_state = (*hmm)[i];
                //emission = current_state->get_emission(*seqs,position);
                emission = (*hmm)[i]->get_emission_prob(*seqs, position);
				
				if (exDef_defined && exDef_position){
                    emission += seqs->getWeight(position, i);
                }
                
				from_trans = (*hmm)[i]->getFrom();
				
                for (size_t j = 0; j < state_size ; ++j) {  //j is previous state
                    if (!(*from_trans)[j]){
                        continue;
                    }
					
                    if ((*viterbi_previous)[j] != INFINITY){
                        viterbi_temp = (*hmm)[j]->get_transition_prob(*seqs, i ,NULL) + emission + (*viterbi_previous)[j];
                        //std::cout << viterbi_temp << " ";
                        if (viterbi_temp > (*viterbi_current)[i]){
                            (*viterbi_current)[i] = viterbi_temp;
                            (*traceback_table)[position][i] = j;
                        }
						
						next_states |= (*(*hmm)[i]->getTo());
                    }
                }
            }
            
        }
        
        //TODO:  Calculate ending and set the final viterbi and traceback pointer
        //Swap current and previous viterbi scores
        viterbi_previous->assign(state_size,-INFINITY);
        swap_ptr = viterbi_previous;
        viterbi_previous = viterbi_current;
        viterbi_current = swap_ptr;
        
        for(size_t i = 0; i < state_size ;++i){
            if ((*viterbi_previous)[i] > -INFINITY){
                viterbi_temp = (*viterbi_previous)[i] + (*hmm)[i]->getEndTrans();
                
                if (viterbi_temp > -INFINITY){
                    ending_viterbi_score = viterbi_temp;
                    ending_viterbi_tb = i;
                    std::cout << viterbi_temp << std::endl;
                }
            }
        }
        
        
        delete viterbi_previous;
        delete viterbi_current;
	}
	
	void trellis::viterbi(model* h, sequences* sqs){
		//Initialize the table
		hmm = h;
		seqs = sqs;
		seq_size		= seqs->getLength();
		state_size		= hmm->state_size();
		exDef_defined	= seqs->exDefDefined();

        viterbi();
    }
	
	
}




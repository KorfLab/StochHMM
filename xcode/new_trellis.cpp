//
//  new_trellis.cpp
//  StochHMM
//
//  Created by Paul Lott on 11/13/12.
//  Copyright (c) 2012 Korf Lab, Genome Center, UC Davis, Davis, CA. All rights reserved.
//

#include "new_trellis.h"


namespace StochHMM {
	trellis::mem_trellis(){
		hmm=NULL;
		seqs=NULL;
		
		type(SIMPLE);
		store_values=false;
		
		traceback=NULL;
		stoch_traceback=NULL;
		nth_traceback=NULL;
		viterbi_score=NULL;
		forward_score=NULL;
		backward_score=NULL;
		posterior=NULL;
		ending_stoch_tb=NULL;
		ending_nth_viterbi_tb=NULL;
	}
	
	trellis::trellis(model* h, sequences* sqs){
		hmm=h;
		seqs=sq;
		
		type(SIMPLE);
		store_values=false;
		
		traceback=NULL;
		stoch_traceback=NULL;
		nth_traceback=NULL;
		viterbi_score=NULL;
		forward_score=NULL;
		backward_score=NULL;
		posterior=NULL;
		ending_stoch_tb=NULL;
		ending_nth_viterbi_tb=NULL;
	}
	
	
	
	trellis::~mem_trellis(){
		delete traceback;
		delete stoch_traceback;
		delete nth_traceback;
		delete viterbi_score;
		delete forward_score;
		delete backward_score;
		delete posterior;
		delete ending_stoch_tb;
		delete ending_nth_viterbi;
	}
	
	void trellis::init_table(){
		seq_size = seqs->seq_size;
		state_size = hmm->state_size;
		exDef_defined = seqs->exDefDefined();
		
		
		previous_calc_cells = new std::vector<double> (state_size,-INFINITY);
        current_calc_cells  = new std::vector<double> (state_size,-INFINITY);
	}
	
	void trellis::viterbi(model* h, sequences* sqs){
		//Initialize the table
		hmm = h;
		seqs = sqs;
		init_table();
        
		
        clock_t start = clock();
		
        double  viterbi_temp(-INFINITY);
        double  emission(-INFINITY);
        bool exDef_position(false);
        
        size_t current_state_iter;
        state* current_state;
        state* previous_state = hmm->getInitial();
        size_t previous_state_iter;
        
        
//        std::vector<std::vector<bool> >* to_trans = process_to_trans_new(hmm);
//        std::vector<std::vector<bool> >* from_trans = process_from_trans_new(hmm);
//        std::vector<bool> initial(state_size,false);
//        process_initial_new(hmm,initial);
//        
//        std::vector<bool>* empty = new std::vector<bool>(state_size,false);
//        std::vector<bool>* next_states = new std::vector<bool>(state_size, false);
//        std::vector<bool>* current_states = new std::vector<bool>(state_size,false);
//        std::vector<bool>* temp_states;
        
        state* init = hmm->getInitial();
        state* ending  = hmm->getEnding();
        
        //Initialize the traceback table
        if (trell->traceback!= NULL){
            delete trell->traceback;
        }
        
        trell->traceback = new two_int_table(seq_size,std::vector<uint16_t> (state_size));
        
        clock_t stop = clock();
        std::cout << "Initialization Time:\t" << (double)(stop-start)/CLOCKS_PER_SEC << std::endl;
        start = clock();
        
        //Calculate Viterbi from INIT state
        for(size_t i = 0; i < state_size; ++i){
            if (initial[i]){
                //current_state = (*hmm)[i];
                //viterbi_temp = current_state->get_emission(*seqs,0) + init->getTrans(i)->getTransition(0,NULL);
				
				viterbi_temp = new_emm[i]->get_emission(0) + init->getTrans(i)->getTransition(0,NULL);
                
                //Emission:  current_state->get_emission(*seq,0);
                //Transition: (previous_state->getTrans(current_state_iter))->getTransition(0,NULL);
                
                
                if (viterbi_temp > -INFINITY){
                    if ((*viterbi_current)[i] < viterbi_temp){
                        (*viterbi_current)[i] = viterbi_temp;
                    }
					for (size_t j=0;j<state_size;++j){
						(*next_states)[j] = (*next_states)[j] || (*to_trans)[i][j];
					}
                    
                }
            }
        }
        
        
        for(size_t position = 1; position < seq_size ; ++position ){
            
            //Swap current and previous viterbi scores
            viterbi_previous->assign(state_size,-INFINITY);
            viterbi_swap = viterbi_previous;
            viterbi_previous = viterbi_current;
            viterbi_current = viterbi_swap;
            
            //Swap current_states and next states sets
            temp_states = current_states;
            current_states = next_states;
			next_states = temp_states;
            next_states->assign(state_size,false);
            
            if (exDef_defined){
                exDef_position = seqs->exDefDefined(position);
            }
            
            for (size_t i = 0; i < state_size; ++i){
                if (!(*current_states)[i]){
                    continue;
                }
                
                //current_state = (*hmm)[i];
                //emission = current_state->get_emission(*seqs,position);
                emission = new_emm[i]->get_emission(position);
				if (exDef_defined && exDef_position){
                    emission += seqs->getWeight(position, i);
                }
                
                for (size_t j = 0; j < state_size ; ++j) {
                    if (!(*from_trans)[i][j]){
                        continue;
                    }
                    
                    previous_state = (*hmm)[j];
                    
                    if ((*viterbi_previous)[j] != INFINITY){
                        viterbi_temp = previous_state->get_trans(*seqs, i, 0) + emission + (*viterbi_previous)[j];
                        //std::cout << viterbi_temp << " ";
                        if (viterbi_temp > (*viterbi_current)[i]){
                            (*viterbi_current)[i] = viterbi_temp;
                            (*(*trell).traceback)[position][i] = j;
                        }
						for (size_t j=0;j<state_size;++j){
							(*next_states)[j] = (*next_states)[j] || (*to_trans)[i][j];
						}
                    }
                }
            }
            
        }
        
        //TODO:  Calculate ending and set the final viterbi and traceback pointer
        //Swap current and previous viterbi scores
        viterbi_previous->assign(state_size,-INFINITY);
        viterbi_swap = viterbi_previous;
        viterbi_previous = viterbi_current;
        viterbi_current = viterbi_swap;
        
        for(size_t i = 0; i < state_size ;++i){
            if ((*viterbi_previous)[i] > -INFINITY){
                previous_state = (*hmm)[i];
                viterbi_temp = (*viterbi_previous)[i] + previous_state->getEndTrans();
                
                if (viterbi_temp > -INFINITY){
                    trell->ending_viterbi_score = viterbi_temp;
                    trell->ending_viterbi_tb = i;
                    std::cout << viterbi_temp << std::endl;
                }
            }
        }
        
        stop = clock();
        std::cout << "Trellis Time:\t" << (double)(stop-start)/CLOCKS_PER_SEC << std::endl;
        start = clock();
        
        
        delete viterbi_previous;
        delete viterbi_current;
        
        stop = clock();
        std::cout << "Cleanup Time:\t" << (double)(stop-start)/CLOCKS_PER_SEC << std::endl;
    }
	
	
}




//
//  mem_viterbi.cpp
//  StochHMM
//
//  Created by Paul Lott on 11/2/12.
//  Copyright (c) 2012 Korf Lab, Genome Center, UC Davis, Davis, CA. All rights reserved.
//

#include "mem_viterbi.h"

namespace StochHMM {
    
    void viterbi(mem_trellis* trell, model* hmm, sequences* seqs){
        
        clock_t start = clock();
        
        size_t  seq_size    = seqs->getLength();
        size_t  state_size  = hmm->state_size();
        std::vector<double>* viterbi_previous = new std::vector<double> (state_size,-INFINITY);
        std::vector<double>* viterbi_current  = new std::vector<double> (state_size,-INFINITY);
        std::vector<double>* viterbi_swap;
        double  viterbi_temp;
        double  emission;

        bool exDef_defined = seqs->exDefDefined();
        bool exDef_position(false);
        
        size_t current_state_iter;
        state* current_state;
        state* previous_state = hmm->getInitial();
        size_t previous_state_iter;
        
        std::set<state*>* next_states    = new std::set<state*>;
        std::set<state*>* defined_states = new std::set<state*>;
        std::set<state*>* current_states = new std::set<state*>;
        std::set<state*>* previous_states= new std::set<state*>(previous_state->getToBegin(), previous_state->getToEnd());
        std::set<state*>* temp_states = NULL;

        std::set<state*>::iterator state_set_it;
        std::vector<state*>::iterator state_vec_it;
    
        //Initialize the traceback table
        if (trell->traceback!= NULL){
            delete trell->traceback;
        }
        
        trell->traceback = new two_int_table(seq_size,std::vector<uint16_t> (state_size));
        
        clock_t stop = clock();
        std::cout << "Initialization Time:\t" << (double)(stop-start)/CLOCKS_PER_SEC << std::endl;
        start = clock();
        
        //Calculate Viterbi from INIT state
        for(state_set_it = previous_states->begin(); state_set_it != previous_states->end(); ++state_set_it){
            current_state = (*state_set_it);
            current_state_iter = current_state->getIterator();
            viterbi_temp = current_state->get_emission(*seqs,0) + (previous_state->getTrans(current_state_iter))->getTransition(0,NULL);
            
            //Emission:  current_state->get_emission(*seq,0);
            //Transition: (previous_state->getTrans(current_state_iter))->getTransition(0,NULL);
            
            
            if (viterbi_temp > -INFINITY){
                if ((*viterbi_current)[current_state_iter] < viterbi_temp){
                    (*viterbi_current)[current_state_iter] = viterbi_temp;
                }
                
                defined_states->insert(current_state);
                next_states->insert(current_state->getToBegin(),current_state->getToEnd());
            }
            
        }
        
        
        for(size_t position = 1; position < seq_size ; ++position ){
            
            //Swap current and previous viterbi scores
            viterbi_previous->assign(state_size,-INFINITY);
            viterbi_swap = viterbi_previous;
            viterbi_previous = viterbi_current;
            viterbi_current = viterbi_swap;
            
//            //TODO: remove these b/c don't use them
//            //Swap previous and defined states sets
//            previous_states->clear();
//            temp_states = previous_states;
//            previous_states = defined_states;
//            defined_states = temp_states;
            
            //Swap current_states and next states sets
            current_states->clear();
            temp_states = current_states;
            current_states = next_states;
            next_states = temp_states;
            
            if (current_states->size() == 0){
                std::cerr << "No valid path through the sequence using the given model\n";
                return;
            }
            
            if (exDef_defined){
                exDef_position = seqs->exDefDefined(position);
            }
            
            for (state_set_it = current_states->begin(); state_set_it != current_states->end(); ++state_set_it){
                current_state = (*state_set_it);
                current_state_iter = current_state->getIterator();
                emission = current_state->get_emission(*seqs,position);
                if (exDef_defined && exDef_position){
                    emission += seqs->getWeight(position, current_state_iter);
                }
                
                for (state_vec_it=current_state->getFromBegin(); state_vec_it != current_state->getFromEnd(); ++state_vec_it) {
                    previous_state = (*state_vec_it);
                    previous_state_iter = previous_state->getIterator();
                    
                    if ((*viterbi_previous)[previous_state_iter] != INFINITY){
                        viterbi_temp = previous_state->get_trans(*seqs, current_state_iter, 0) + emission + (*viterbi_previous)[previous_state_iter];
                        //std::cout << viterbi_temp << " ";
                        if (viterbi_temp > (*viterbi_current)[current_state_iter]){
                            (*viterbi_current)[current_state_iter] = viterbi_temp;
                            (*(*trell).traceback)[position][current_state_iter] = previous_state_iter;
                            next_states->insert(current_state->getToBegin(), current_state->getToEnd());
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
            if ((*viterbi_previous)[i] >- INFINITY){
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
        delete next_states;
        delete defined_states;
        delete previous_states;
        
        stop = clock();
        std::cout << "Cleanup Time:\t" << (double)(stop-start)/CLOCKS_PER_SEC << std::endl;
    }
    
    
    std::vector<std::bitset<400> >* process_to_trans(model *hmm){
        size_t state_size = hmm->state_size();
        
        std::vector<std::bitset<400> >* trans = new std::vector<std::bitset<400> > (state_size);
        
        for(size_t i = 0; i < state_size ; ++i){
            state* current = (*hmm)[i];
            std::vector<state*>* to_vector = current->getTo();
            for(size_t j = 0; j < to_vector->size(); ++j){
                size_t iterator = (*to_vector)[j]->getIterator();
                (*trans)[i][j]=1;
            }
        }
        
        return trans;
    }
    
    std::vector<std::bitset<400> >* process_from_trans(model *hmm){
        size_t state_size = hmm->state_size();
        
        std::vector<std::bitset<400> >* trans = new std::vector<std::bitset<400> > (state_size);
        
        for(size_t i = 0; i < state_size ; ++i){
            state* current = (*hmm)[i];
            std::vector<state*>* to_vector = current->getFrom();
            for(size_t j = 0; j < to_vector->size(); ++j){
                size_t iterator = (*to_vector)[j]->getIterator();
                (*trans)[i][j]=1;
            }
        }
        
        return trans;
    }
    
    void process_initial(model* hmm, std::bitset<400>& ending){
        state* initial = hmm->getInitial();
        std::vector<state*>* to_vector = initial->getTo();
        for(size_t j = 0; j < to_vector->size(); ++j){
            size_t iterator = (*to_vector)[j]->getIterator();
            ending[j]=1;
        }
        return;
    }
    
    void viterbi_two(mem_trellis* trell, model* hmm, sequences* seqs){
        
       
        
        clock_t start = clock();
        
        size_t  seq_size    = seqs->getLength();
        size_t  state_size  = hmm->state_size();
        std::vector<double>* viterbi_previous = new std::vector<double> (state_size,-INFINITY);
        std::vector<double>* viterbi_current  = new std::vector<double> (state_size,-INFINITY);
        std::vector<double>* viterbi_swap;
        double  viterbi_temp;
        double  emission;
        
        bool exDef_defined = seqs->exDefDefined();
        bool exDef_position(false);
        
        size_t current_state_iter;
        state* current_state;
        state* previous_state = hmm->getInitial();
        size_t previous_state_iter;
        
        
        std::vector<std::bitset<400> >* to_trans = process_to_trans(hmm);
        std::vector<std::bitset<400> >* from_trans = process_from_trans(hmm);
        std::bitset<400> initial;
        process_initial(hmm,initial);
        
        std::bitset<400> empty;
        std::bitset<400> next_states;
        std::bitset<400> current_states;
        std::bitset<400> temp_states;
        
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
            if (initial.test(i)){
                current_state = (*hmm)[i];
                viterbi_temp = current_state->get_emission(*seqs,0) + init->getTrans(i)->getTransition(0,NULL);
                
                //Emission:  current_state->get_emission(*seq,0);
                //Transition: (previous_state->getTrans(current_state_iter))->getTransition(0,NULL);
                
                
                if (viterbi_temp > -INFINITY){
                    if ((*viterbi_current)[i] < viterbi_temp){
                        (*viterbi_current)[i] = viterbi_temp;
                    }
                    next_states |= (*to_trans)[i];
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
            
            current_states.reset();
            current_states |= next_states;
            next_states.reset();
            
            if (current_states.none()){
                std::cerr << "No valid path through the sequence using the given model\n";
                return;
            }
            
            if (exDef_defined){
                exDef_position = seqs->exDefDefined(position);
            }
            
            for (size_t i = 0; i < state_size; ++i){
                if (!current_states.test(i)){
                    continue;
                }
                
                current_state = (*hmm)[i];
                emission = current_state->get_emission(*seqs,position);
                if (exDef_defined && exDef_position){
                    emission += seqs->getWeight(position, i);
                }
                
                for (size_t j = 0; j < state_size ; ++j) {
                    if (!(*from_trans)[i].test(j)){
                        continue;
                    }
                    
                    previous_state = (*hmm)[j];
                    
                    if ((*viterbi_previous)[j] != INFINITY){
                        viterbi_temp = previous_state->get_trans(*seqs, i, 0) + emission + (*viterbi_previous)[j];
                        //std::cout << viterbi_temp << " ";
                        if (viterbi_temp > (*viterbi_current)[i]){
                            (*viterbi_current)[i] = viterbi_temp;
                            (*(*trell).traceback)[position][i] = j;
                            next_states |= (*to_trans)[i];
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
            if ((*viterbi_previous)[i] >- INFINITY){
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
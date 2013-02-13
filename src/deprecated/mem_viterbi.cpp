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
                    //std::cout << viterbi_temp << std::endl;
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
    
    
    std::vector<std::bitset<200> >* process_to_trans(model *hmm){
        size_t state_size = hmm->state_size();
        
        std::vector<std::bitset<200> >* trans = new std::vector<std::bitset<200> > (state_size);
        
        for(size_t i = 0; i < state_size ; ++i){
            state* current = (*hmm)[i];
            std::vector<state*>* to_vector = current->getTo();
            for(size_t j = 0; j < to_vector->size(); ++j){
                size_t iterator = (*to_vector)[j]->getIterator();
                (*trans)[i][iterator]=1;
            }
            //std::cout << "State " << i << "\t" << (*trans)[i].to_string() << std::endl;
        }
        
        return trans;
    }
    
    std::vector<std::bitset<200> >* process_from_trans(model *hmm){
        size_t state_size = hmm->state_size();
        
        std::vector<std::bitset<200> >* trans = new std::vector<std::bitset<200> > (state_size);
        
        for(size_t i = 0; i < state_size ; ++i){
            state* current = (*hmm)[i];
            std::vector<state*>* to_vector = current->getFrom();
            for(size_t j = 0; j < to_vector->size(); ++j){
                size_t iterator = (*to_vector)[j]->getIterator();
                (*trans)[i][iterator]=1;
            }
        }
        
        return trans;
    }
    
    void process_initial(model* hmm, std::bitset<200>& ending){
        state* initial = hmm->getInitial();
        std::vector<state*>* to_vector = initial->getTo();
        for(size_t j = 0; j < to_vector->size(); ++j){
            size_t iterator = (*to_vector)[j]->getIterator();
            ending[iterator]=1;
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
        
        
        std::vector<std::bitset<200> >* to_trans = process_to_trans(hmm);
        std::vector<std::bitset<200> >* from_trans = process_from_trans(hmm);
        std::bitset<200> initial;
        process_initial(hmm,initial);
        
        std::bitset<200> empty;
        std::bitset<200> next_states;
        std::bitset<200> current_states;
        std::bitset<200> temp_states;
        
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
                if (!current_states[i]){
                    continue;
                }
                
                current_state = (*hmm)[i];
                emission = current_state->get_emission(*seqs,position);
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
                        next_states |= (*to_trans)[i];
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
                    //std::cout << viterbi_temp << std::endl;
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
    
    
    void convert_seq(sequence* sq, std::vector<short>& new_seq){
        for(size_t i = 0; i< sq->size(); ++i){
            short value = sq->seqValue(i);
            new_seq[i] = (value > 0) ? value : 4;
        }
        return;
    }
    
    
    
    emissions::emissions(emm* emiss,std::vector<short>* sq):seq(sq){
        lexicalTable* scores = emiss->getTables();
        order = scores->getOrder(0);
        
        zero      = new zero_order(5,0);
		first     = (order>=1) ? new first_order(5,*zero)	: NULL;
        second    = (order>=2) ? new second_order(5,*first)	: NULL;
        third     = (order>=3) ? new third_order(5,*second)	: NULL;
        fourth    = (order>=4) ? new fourth_order(5,*third)	: NULL;
        fifth     = (order>=5) ? new fifth_order(5,*fourth)	: NULL;
        parse_emission(emiss);
    }
    
    void emissions::parse_emission(emm* emiss){
        lexicalTable* scores = emiss->getTables();
        std::vector<std::vector<double> >* counts = scores->getCountsTable();
        
		//Transfer values from old emission class
        for(size_t row=0;row < counts->size(); ++row){
            size_t* indices = indexToindices(row,order);
            for(size_t col = 0; col < (*counts)[row].size(); ++col){
                double value = (*counts)[row][col];
				//std::cout <<"Column:"<< col << "\tValue:" << value << std::endl;
                switch (order){
                    case 5:
                        (*fifth)[indices[0]][indices[1]][indices[2]][indices[3]][indices[4]][col]=value;
                        (*fourth)[indices[1]][indices[2]][indices[3]][indices[4]][col]+=value;
                        (*third)[indices[2]][indices[3]][indices[4]][col]+=value;
                        (*second)[indices[3]][indices[4]][col]+=value;
                        (*first)[indices[4]][col]+=value;
                        (*zero)[col]+=value;
                        break;
                    case 4:
                        (*fourth)[indices[0]][indices[1]][indices[2]][indices[3]][col]+=value;
                        (*third)[indices[1]][indices[2]][indices[3]][col]+=value;
                        (*second)[indices[2]][indices[3]][col]+=value;
                        (*first)[indices[3]][col]+=value;
                        (*zero)[col]+=value;
                        break;
                    case 3:
                        (*third)[indices[0]][indices[1]][indices[2]][col]+=value;
                        (*second)[indices[1]][indices[2]][col]+=value;
                        (*first)[indices[2]][col]+=value;
                        (*zero)[col]+=value;
                        break;
                    case 2:
                        (*second)[indices[0]][indices[1]][col]+=value;
                        (*first)[indices[1]][col]+=value;
                        (*zero)[col]+=value;
                        break;
                    case 1:
                        (*first)[indices[0]][col]+=value;
                        (*zero)[col]+=value;
                        break;
                    default:
                        (*zero)[col]+=value;
                        break;
                }
            }
            delete indices;
        }
		
		//Convert vectors to Probabilities and Calculate emission of N
		for(size_t row=0;row< counts->size(); ++row){
            size_t* indices = indexToindices(row,order);
            for(size_t col = 0; col < (*counts)[row].size(); ++col){
                switch (order){
                    case 5:
                        probVector((*fifth)[indices[0]][indices[1]][indices[2]][indices[3]][indices[4]]);
						probVector((*fourth)[indices[1]][indices[2]][indices[3]][indices[4]]);
						probVector((*third)[indices[2]][indices[3]][indices[4]]);
						probVector((*second)[indices[3]][indices[4]]);
						probVector((*first)[indices[4]]);
						probVector(*zero);
						
						(*fifth)[indices[0]][indices[1]][indices[2]][indices[3]][indices[4]][4]=sumVector((*fifth)[indices[0]][indices[1]][indices[2]][indices[3]][indices[4]])/4;
						(*fourth)[indices[1]][indices[2]][indices[3]][indices[4]][4]=sumVector((*fourth)[indices[1]][indices[2]][indices[3]][indices[4]])/4;
						(*third)[indices[2]][indices[3]][indices[4]][4]=sumVector((*third)[indices[2]][indices[3]][indices[4]])/4;
						(*second)[indices[3]][indices[4]][4]=sumVector((*second)[indices[3]][indices[4]])/4;
						(*first)[indices[4]][4]=sumVector((*first)[indices[4]])/4;
						(*zero)[4]=sumVector((*zero))/4;
                        break;
                    case 4:
                        probVector((*fourth)[indices[0]][indices[1]][indices[2]][indices[3]]);
                        probVector((*third)[indices[1]][indices[2]][indices[3]]);
                        probVector((*second)[indices[2]][indices[3]]);
                        probVector((*first)[indices[3]]);
                        probVector(*zero);
						
						(*fourth)[indices[0]][indices[1]][indices[2]][indices[3]][4]=sumVector((*fourth)[indices[0]][indices[1]][indices[2]][indices[3]])/4;
                        (*third)[indices[1]][indices[2]][indices[3]][4]=sumVector((*third)[indices[1]][indices[2]][indices[3]])/4;
                        (*second)[indices[2]][indices[3]][4]=sumVector((*second)[indices[2]][indices[3]])/4;
                        (*first)[indices[3]][4]=sumVector((*first)[indices[3]])/4;
                        (*zero)[4]=sumVector((*zero))/4;
                        break;
                    case 3:
                        probVector((*third)[indices[0]][indices[1]][indices[2]]);
                        probVector((*second)[indices[1]][indices[2]]);
                        probVector((*first)[indices[2]]);
                        probVector(*zero);
						
						(*third)[indices[0]][indices[1]][indices[2]][4]=sumVector((*third)[indices[0]][indices[1]][indices[2]])/4;
                        (*second)[indices[1]][indices[2]][4]=sumVector((*second)[indices[1]][indices[2]])/4;
                        (*first)[indices[2]][4]=sumVector((*first)[indices[2]])/4;
                        (*zero)[4]=sumVector(*zero)/4;
                        break;
                    case 2:
						probVector((*second)[indices[0]][indices[1]]);
                        probVector((*first)[indices[1]]);
                        probVector(*zero);
						
						(*second)[indices[0]][indices[1]][4]=sumVector((*second)[indices[0]][indices[1]])/4;
                        (*first)[indices[1]][4]=sumVector((*first)[indices[1]])/4;
                        (*zero)[4]=sumVector(*zero)/4;
                        break;
                    case 1:
                        probVector((*first)[indices[0]]);
                        probVector(*zero);
						
						(*first)[indices[0]][4]=sumVector((*first)[indices[0]])/4;
                        (*zero)[4]=sumVector(*zero)/4;
                        break;
                    default:
                        probVector(*zero);
						
						(*zero)[4]=sumVector(*zero)/4;
                        break;
                }
            }
            delete indices;
        }
		
		//Fill out Ambiguous characters in higher order
		calc_ambig();
		
        
        return;
    }
	
	
	void emissions::calc_ambig(){
		
		//Create ambiguous character N;
		std::vector<size_t> ambiguous;
		ambiguous.push_back(0);
		ambiguous.push_back(1);
		ambiguous.push_back(2);
		ambiguous.push_back(3);
		
		for(unsigned int i=1; i <= order ; --i){
			
			std::vector<std::vector<size_t> > words;
			permute(i,words);
			for(size_t j=0;j<words.size();j++){
				
				std::vector<std::vector<size_t> > new_words;
				expand_amb(0,ambiguous,words[j],new_words);
				
				for(size_t k = 0; k < new_words.size(); i++){
					std::vector<std::vector<size_t> > unambig_words;
					std::vector<size_t>& word = new_words[i];
					expand_amb(5,ambiguous,word,unambig_words);
					
					add_word(word,unambig_words);
					
				}
			}
		}
		return;
	}
	
	void emissions::add_word(std::vector<size_t>& ambig, std::vector<std::vector<size_t> >& words){
		
		size_t word_size = ambig.size();
		size_t num_of_words = words.size();
		for(size_t i = 0; i < num_of_words; i++){
			if (word_size == 1){
				addVector((*first)[ambig[0]],(*first)[words[i][0]]);
			}
			
			if (word_size == 2){
				addVector((*second)[ambig[0]][ambig[1]],(*second)[words[i][0]][words[i][1]]);
			}
			
			if (word_size == 3){
				addVector((*third)[ambig[0]][ambig[1]][ambig[2]],(*third)[words[i][0]][words[i][1]][words[i][2]]);
			}
			
			if (word_size == 4){
				addVector((*fourth)[ambig[0]][ambig[1]][ambig[2]][ambig[3]],(*fourth)[words[i][0]][words[i][1]][words[i][2]][words[i][3]]);
			}
			
			if (word_size == 5){
				addVector((*fifth)[ambig[0]][ambig[1]][ambig[2]][ambig[3]][ambig[4]],(*fifth)[words[i][0]][words[i][1]][words[i][2]][words[i][3]][words[i][4]]);
			}
		}
		
		if (word_size == 1){
			divideValueToVector((*first)[ambig[0]],num_of_words);
		}
		
		if (word_size == 2){
			divideValueToVector((*second)[ambig[0]][ambig[1]],num_of_words);
		}
		
		if (word_size == 3){
			divideValueToVector((*third)[ambig[0]][ambig[1]][ambig[2]],num_of_words);
		}
		
		if (word_size == 4){
			divideValueToVector((*fourth)[ambig[0]][ambig[1]][ambig[2]][ambig[3]],num_of_words);
		}
		
		if (word_size == 5){
			divideValueToVector((*fifth)[ambig[0]][ambig[1]][ambig[2]][ambig[3]][ambig[4]],num_of_words);
		}
		
		return;
	}
	
	void permute(size_t word_length, std::vector<std::vector<size_t> >& words){
		
		for(size_t i = 1; i < word_length ;++i){
			std::vector<size_t > word(word_length,0);
			for(size_t j = 0; j < i; ++j){
				word[j]=4;
			}
			sort(word.begin(),word.end());
			
			do{
				words.push_back(word);
			} while( next_permutation(word.begin(),word.end()));
			
		}
		return;
	}
	
	void expand_amb(size_t value, std::vector<size_t>& ambiguous, std::vector<size_t>& word, std::vector<std::vector<size_t> >& ret_val){
		_expand_amb(value, 0ULL, ambiguous, word, ret_val);
		return;
	}
	
	void _expand_amb(size_t value, size_t position, std::vector<size_t>& ambiguous,std::vector<size_t>& word, std::vector<std::vector<size_t> >& ret_val){
		
		if (position == word.size()){
			ret_val.push_back(word);
			return;
		}
		//Words that have a expanding character. We'll expand and pass the expanded
		//word to the evaluate the next position
		//Add all the return values to the return array
		else if (word[position] == value){
			for(size_t i=0; i< ambiguous.size(); ++i){
				std::vector<size_t> temp(word);
				temp[position]=ambiguous[i];
				std::vector<std::vector<size_t> > tempor;
				_expand_amb(value,position+1,ambiguous,temp,tempor);
				ret_val.insert(ret_val.end(), tempor.begin(), tempor.end());
			}
		}
		else{ //Doesn't have expanding character so evaluate next position
			std::vector<std::vector<size_t> > tempor;
			_expand_amb(value,position+1,ambiguous,word,tempor);
			ret_val.insert(ret_val.end(), tempor.begin(), tempor.end());
		}
		return;		
	}
	
	
	
    
    size_t* indexToindices(size_t index, size_t length){
        size_t position = length;
        size_t* output = new size_t[length];
        for(size_t i=0; i<length ;++i){
            size_t dreg = POWER[position-1][3];
            float value=(float)index/dreg;
            if (value<1){
                output[i]=0;
            }
            else if (value<2){
                output[i]=1;
                index-=dreg;}
            else if (value<3){
                output[i]=2;
                index-=2*dreg;}
            else {
                output[i]=3;
                index-=3*dreg;
            }
            position--;
        }
        
        return output;
    }
	
	double emissions::get_emission(size_t position){
		if (order==0 || position == 0){
			return (*zero)[(*seq)[position]];
		}
		else if (order == 1 || position == 1){
			return (*first)[(*seq)[position-1]][(*seq)[position]];
		}
		else if (order == 2 || position == 2){
			return (*second)[(*seq)[position-2]][(*seq)[position-1]][(*seq)[position]];
		}
		else if (order == 3 || position == 3){
			return (*third)[(*seq)[position-3]][(*seq)[position-2]][(*seq)[position-1]][(*seq)[position]];
		}
		else if (order == 4 || position ==4){
			return (*fourth)[(*seq)[position-4]][(*seq)[position-3]][(*seq)[position-2]][(*seq)[position-1]][(*seq)[position]];
		}
		else if (order==5){
			return (*fifth)[(*seq)[position-5]][(*seq)[position-4]][(*seq)[position-3]][(*seq)[position-2]][(*seq)[position-1]][(*seq)[position]];
		}
	}


    
    
    void viterbi_three(mem_trellis* trell, model* hmm, sequences* seqs){
        
        size_t  seq_size    = seqs->getLength();
        size_t  state_size  = hmm->state_size();

        //Get New Sequence
        sequence* seq = seqs->getSeq(0);
        std::vector<short> new_seq (seq_size);
		
        for(size_t i=0;i < seq_size; ++i){
            short letter = seq->seqValue(i);
            if (letter > 0){
                new_seq[i]=letter;
            }
            else{
                new_seq[i]=4;
            }
        }
        
        //Get New Emission tables;
        std::vector<emissions*> new_emm (state_size);
        for(size_t i = 0; i < state_size ; ++i){
			//std::cout << i << std::endl;
            state* temp_state = (*hmm)[i];
            emissions* temp = new emissions(temp_state->getEmission(0),&new_seq);
            new_emm[i]=temp;
        }
        
        
        
        clock_t start = clock();
        
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
        
        
        std::vector<std::bitset<200> >* to_trans = process_to_trans(hmm);
        std::vector<std::bitset<200> >* from_trans = process_from_trans(hmm);
        std::bitset<200> initial;
        process_initial(hmm,initial);
        
        std::bitset<200> empty;
        std::bitset<200> next_states;
        std::bitset<200> current_states;
        std::bitset<200> temp_states;
        
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
                if (!current_states[i]){
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
                        next_states |= (*to_trans)[i];
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
	
	
	
	
	std::vector<std::vector<bool> >* process_to_trans_new(model *hmm){
        size_t state_size = hmm->state_size();
        
		std::vector<bool> temp (state_size,false);
        std::vector<std::vector<bool> >* trans = new std::vector<std::vector<bool> > (state_size,temp);
        
        for(size_t i = 0; i < state_size ; ++i){
            state* current = (*hmm)[i];
            std::vector<state*>* to_vector = current->getTo();
            for(size_t j = 0; j < to_vector->size(); ++j){
                size_t iterator = (*to_vector)[j]->getIterator();
                (*trans)[i][iterator]=1;
            }
            //std::cout << "State " << i << "\t" << (*trans)[i].to_string() << std::endl;
        }
        
        return trans;
    }
    
    std::vector<std::vector<bool> >* process_from_trans_new(model *hmm){
        size_t state_size = hmm->state_size();
        
        std::vector<bool> temp (state_size,false);
        std::vector<std::vector<bool> >* trans = new std::vector<std::vector<bool> > (state_size,temp);
        
        for(size_t i = 0; i < state_size ; ++i){
            state* current = (*hmm)[i];
            std::vector<state*>* to_vector = current->getFrom();
            for(size_t j = 0; j < to_vector->size(); ++j){
                size_t iterator = (*to_vector)[j]->getIterator();
                (*trans)[i][iterator]=1;
            }
        }
        
        return trans;
    }
    
    void process_initial_new(model* hmm, std::vector<bool>& ending){
        state* initial = hmm->getInitial();
        std::vector<state*>* to_vector = initial->getTo();
        for(size_t j = 0; j < to_vector->size(); ++j){
            size_t iterator = (*to_vector)[j]->getIterator();
            ending[iterator]=1;
        }
        return;
    }
	
	
	
	
	
	void viterbi_four(mem_trellis* trell, model* hmm, sequences* seqs){
        
        size_t  seq_size    = seqs->getLength();
        size_t  state_size  = hmm->state_size();
		
        //Get New Sequence
        sequence* seq = seqs->getSeq(0);
        std::vector<short> new_seq (seq_size);
		
        for(size_t i=0;i < seq_size; ++i){
            short letter = seq->seqValue(i);
            if (letter > 0){
                new_seq[i]=letter;
            }
            else{
                new_seq[i]=4;
            }
        }
        
        //Get New Emission tables;
        std::vector<emissions*> new_emm (state_size);
        for(size_t i = 0; i < state_size ; ++i){
			//std::cout << i << std::endl;
            state* temp_state = (*hmm)[i];
            emissions* temp = new emissions(temp_state->getEmission(0),&new_seq);
            new_emm[i]=temp;
        }
        
        
        
        clock_t start = clock();
        
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
        
        
        std::vector<std::vector<bool> >* to_trans = process_to_trans_new(hmm);
        std::vector<std::vector<bool> >* from_trans = process_from_trans_new(hmm);
        std::vector<bool> initial(state_size,false);
        process_initial_new(hmm,initial);
        
        std::vector<bool>* empty = new std::vector<bool>(state_size,false);
        std::vector<bool>* next_states = new std::vector<bool>(state_size, false);
        std::vector<bool>* current_states = new std::vector<bool>(state_size,false);
        std::vector<bool>* temp_states;
        
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
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
		
		seq_size=0;
		state_size=0;
		
		type = SIMPLE;
		store_values=false;
		exDef_defined=false;
		
		traceback_table	= NULL;
		stochastic_table= NULL;
		nth_traceback	= NULL;
		viterbi_score	= NULL;
		forward_score	= NULL;
		backward_score	= NULL;
		posterior_score	= NULL;
		
		ending_viterbi_score = -INFINITY;
		ending_viterbi_tb = 0;
		ending_posterior = -INFINITY;
		
		scoring_current = NULL;
		scoring_previous= NULL;
		alt_scoring_current = NULL;
		alt_scoring_previous = NULL;
	}
	
	trellis::trellis(model* h, sequences* sqs){
		hmm=h;
		seqs=sqs;
		
		seq_size		= seqs->getLength();
		state_size		= hmm->state_size();
		
		type = SIMPLE;
		store_values=false;
		exDef_defined	= seqs->exDefDefined();
		
		traceback_table		= NULL;
		stochastic_table	= NULL;
		nth_traceback		= NULL;
		viterbi_score		= NULL;
		forward_score		= NULL;
		backward_score		= NULL;
		posterior_score		= NULL;
		
		ending_viterbi_score = -INFINITY;
		ending_viterbi_tb = 0;
		ending_posterior = -INFINITY;

		scoring_current = NULL;
		scoring_previous= NULL;
		alt_scoring_current = NULL;
		alt_scoring_previous = NULL;

	}
	
	
	
	trellis::~trellis(){
		delete traceback_table;
		delete stochastic_table;
		delete nth_traceback;
		delete viterbi_score;
		delete forward_score;
		delete backward_score;
		delete posterior_score;
		
		delete scoring_previous;
		delete scoring_current;
		delete alt_scoring_previous;
		delete alt_scoring_current;
	}
	
	void trellis::reset(){
		hmm=NULL;
		seqs=NULL;
		state_size=0;
		seq_size=0;
		type= SIMPLE;
		store_values = false;
		exDef_defined = false;
		
		delete traceback_table;
		delete stochastic_table;
		delete nth_traceback;
		
		delete viterbi_score;
		delete forward_score;
		delete backward_score;
		delete posterior_score;
		
		ending_viterbi_score = -INFINITY;
		ending_viterbi_tb = 0;
		ending_posterior = -INFINITY;
		
		delete scoring_current;
		delete scoring_previous;
		
		delete alt_scoring_current;
		delete alt_scoring_previous;
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
	
	
	void trellis::simple_viterbi(model* h, sequences* sqs){
		//Initialize the table
		hmm = h;
		seqs = sqs;
		seq_size		= seqs->getLength();
		state_size		= hmm->state_size();
		exDef_defined	= seqs->exDefDefined();
		
        simple_viterbi();
	}
	
	void trellis::simple_viterbi(){
		
		if (!hmm->isBasic()){
			std::cerr << "Model isn't a simple/basic HMM.  Use complex algorithms\n";
			return;
		}
		
		//Initialize the traceback table
        if (traceback_table != NULL){
            delete traceback_table;
        }
        
        traceback_table = new int_2D(seq_size,std::vector<uint16_t> (state_size));
		scoring_previous = new (std::nothrow) std::vector<double> (state_size,-INFINITY);
        scoring_current  = new (std::nothrow) std::vector<double> (state_size,-INFINITY);
		
		if (scoring_previous == NULL || scoring_current == NULL || traceback_table == NULL){
			std::cerr << "Can't allocate Viterbi score and traceback table. OUT OF MEMORY" << std::endl;
			exit(2);
		}
		
		
        std::bitset<STATE_MAX> next_states;
        std::bitset<STATE_MAX> current_states;
		
        double  viterbi_temp(-INFINITY);
        double  emission(-INFINITY);
        bool	exDef_position(false);
        
        state* init = hmm->getInitial();
		
		std::bitset<STATE_MAX>* initial_to = hmm->getInitialTo();
		std::bitset<STATE_MAX>* from_trans(NULL);
		
        //Calculate Viterbi from transitions from INIT (initial) state
        for(size_t i = 0; i < state_size; ++i){
            if ((*initial_to)[i]){  //if the bitset is set (meaning there is a transition to this state), calculate the viterbi
				
				viterbi_temp = (*hmm)[i]->get_emission_prob(*seqs,0) + getTransition(init, i, 0);
                
				if (viterbi_temp > -INFINITY){
                    if ((*scoring_current)[i] < viterbi_temp){
                        (*scoring_current)[i] = viterbi_temp;
                    }
					next_states[i] = 1;
                }
            }
        }
		
        //Each position in the sequence
        for(size_t position = 1; position < seq_size ; ++position ){
            
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
			
			//Current states
            for (size_t current = 0; current < state_size; ++current){ //Current state that emits value
				
				//Check to see if current state is valid
				if (!current_states[current]){
                    continue;
                }
				
				//Get emission of current state
                emission = (*hmm)[current]->get_emission_prob(*seqs, position);
				
								
				if (exDef_defined && exDef_position){
                    emission += seqs->getWeight(position, current);
                }
                
				//Get list of states that are valid previous states
				from_trans = (*hmm)[current]->getFrom();
				
				
                for (size_t prev = 0; prev < state_size ; ++prev) {  //for previous states
                    
					//Check that previous state has transition to current state
					//and that the previous viterbi score is not -INFINITY
                    if ((*from_trans)[prev] && (*scoring_previous)[prev] != -INFINITY){
                        viterbi_temp = getTransition((*hmm)[prev], current , position) + emission + (*scoring_previous)[prev];
                        
						
						if (viterbi_temp > (*scoring_current)[current]){
                            (*scoring_current)[current] = viterbi_temp;
                            (*traceback_table)[position][current] = prev;
                        }
						
						next_states |= (*(*hmm)[current]->getTo());
                    }
                }
            }
        }
        
        //TODO:  Calculate ending and set the final viterbi and traceback pointer
        //Swap current and previous viterbi scores
        scoring_previous->assign(state_size,-INFINITY);
        swap_ptr = scoring_previous;
        scoring_previous = scoring_current;
        scoring_current = swap_ptr;
        
        for(size_t i = 0; i < state_size ;++i){
            if ((*scoring_previous)[i] > -INFINITY){
                viterbi_temp = (*scoring_previous)[i] + (*hmm)[i]->getEndTrans();
                
                if (viterbi_temp > ending_viterbi_score){
                    ending_viterbi_score = viterbi_temp;
                    ending_viterbi_tb = i;
                }
            }
        }
        
        delete scoring_previous;
        delete scoring_current;
		scoring_previous = NULL;
		scoring_current  = NULL;
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
		for(size_t i = 0; i < state_size; ++i){
			if ((*ending_from)[i]){  //if the bitset is set (meaning there is a transition to this state), calculate the viterbi
				
				backward_temp = (*hmm)[i]->getEndTrans();
				
				if (backward_temp > -INFINITY){
					(*backward_score)[seq_size-1][i] = backward_temp;
					(*scoring_current)[i] = backward_temp;
					next_states[i] = 1;
				}
			}
		}
		
		
		for(size_t position = seq_size-1; position > 0 ; --position ){
			
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
					
					//if ((*backward_score)[position-1][j] != -INFINITY){
					if ((*scoring_previous)[j] != -INFINITY){
						
						//backward_temp = getTransition((*hmm)[j], i , position-1) + emission + (*backward_score)[position][i];
						backward_temp = getTransition((*hmm)[j], i , position-1) + emission + (*scoring_previous)[i];

						
						if ((*backward_score)[position-1][j] == -INFINITY){
							(*scoring_current)[j] = backward_temp;
							(*backward_score)[position-1][j] = backward_temp;
						}
						else{
							(*scoring_current)[j] = addLog(backward_temp, (*scoring_current)[j]);
							(*backward_score)[position-1][j] = (*scoring_current)[j];
							
							//(*backward_score)[position-1][j] = addLog((double)backward_temp, (double)(*backward_score)[position-1][j]);
						}
						
						next_states |= (*(*hmm)[i]->getFrom());
					}
				}
//				std::cout << "State: " << i <<"\t" << exp((*backward_score)[position][i]) << std::endl;
			}
			
		}
		
		ending_posterior = -INFINITY;
		double backward_posterior = -INFINITY;
		state* init = hmm->getInitial();
		for(size_t i = 0; i < state_size ;++i){
			if ((*scoring_current)[i] != -INFINITY){
				
				//if ((*backward_score)[0][i] > -INFINITY){
				
				//backward_temp = (*backward_score)[0][i] + (*hmm)[i]->get_emission_prob(*seqs,0) + getTransition(init, i, 0);
				backward_temp = (*scoring_current)[i] + (*hmm)[i]->get_emission_prob(*seqs,0) + getTransition(init, i, 0);
				if (backward_temp > -INFINITY){
					if (backward_posterior == -INFINITY){
						backward_posterior = backward_temp;
					}
					else{
						backward_posterior = addLog(backward_posterior,backward_temp);
					}
				}
			}
		}
		
		ending_posterior = backward_posterior;
		
		delete scoring_previous;
		delete scoring_current;
		scoring_previous = NULL;
		scoring_current = NULL;
		
//		std::cout << exp(backward_posterior) << std::endl;
	}
	
	
	void trellis::simple_forward(model* h, sequences* sqs){
		hmm = h;
		seqs = sqs;
		seq_size		= seqs->getLength();
		state_size		= hmm->state_size();
		exDef_defined	= seqs->exDefDefined();
		
        simple_forward();
	}
	
	void trellis::simple_forward(){
		
		if (!hmm->isBasic()){
			std::cerr << "Model isn't a simple/basic HMM.  Use complex algorithms\n";
			return;
		}

		forward_score	= new (std::nothrow) float_2D(seq_size, std::vector<float>(state_size,-INFINITY));
		scoring_current = new (std::nothrow) std::vector<double> (state_size,-INFINITY);
		scoring_previous= new (std::nothrow) std::vector<double> (state_size,-INFINITY);
		
		if (scoring_current == NULL || scoring_previous == NULL || forward_score == NULL){
			std::cerr << "Can't allocate forward score table. OUT OF MEMORY" << std::endl;
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
					(*scoring_current)[i] = forward_temp;
					next_states[i] = 1;
                }
            }
        }
        
        
        for(size_t position = 1; position < seq_size ; ++position ){
            
            
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
                
				from_trans = (*hmm)[current]->getFrom();
				
                for (size_t previous = 0; previous < state_size ; ++previous) {  //j is previous state
                    if (!(*from_trans)[previous]){
                        continue;
                    }
					
					if ((*scoring_previous)[previous] != -INFINITY){
                        forward_temp = getTransition((*hmm)[previous], current , position) + emission + (*scoring_previous)[previous];
						
						if ((*scoring_current)[current] == -INFINITY){
							(*scoring_current)[current] = forward_temp;
							(*forward_score)[position][current] = forward_temp;
						}
						else{
							(*scoring_current)[current] = addLog(forward_temp, (*scoring_current)[current]);
							(*forward_score)[position][current] = (*scoring_current)[current];						}
						
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

        ending_posterior = -INFINITY;
        for(size_t i = 0; i < state_size ;++i){
            if ((*scoring_previous)[i] != -INFINITY){
                forward_temp = (*scoring_previous)[i] + (*hmm)[i]->getEndTrans();
                
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
		delete scoring_previous;
		delete scoring_current;
		scoring_previous = NULL;
		scoring_current = NULL;

	}
	
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
        
        traceback_table = new int_2D(seq_size,std::vector<uint16_t> (state_size));
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
					next_states[i] = 1;
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
		
	
	void trellis::fast_complex_viterbi(){
		
		//Initialize the traceback table
        if (traceback_table != NULL){
            delete traceback_table;
        }
        
        traceback_table = new int_2D(seq_size,std::vector<uint16_t> (state_size));
		scoring_previous = new (std::nothrow) std::vector<double> (state_size,-INFINITY);
        scoring_current  = new (std::nothrow) std::vector<double> (state_size,-INFINITY);
		
		if (scoring_previous == NULL || scoring_current == NULL || traceback_table == NULL){
			std::cerr << "Can't allocate Viterbi score and traceback table. OUT OF MEMORY" << std::endl;
			exit(2);
		}
		
		
        std::bitset<STATE_MAX> next_states;
        std::bitset<STATE_MAX> current_states;
		
        double  viterbi_temp(-INFINITY);
        double  emission(-INFINITY);
        bool	exDef_position(false);
		
		
		//If model is not a basic model, then we need to initialize the explicit duration vector
		//bool extend_duration keeps track of whether the transition to same state was selected.
		if (!hmm->isBasic()){
			explicit_duration_current = new(std::nothrow) std::vector<size_t>(state_size,0);
			explicit_duration_previous= new(std::nothrow) std::vector<size_t>(state_size,0);
		}
		bool extend_duration(false);
		std::vector<bool>* duration = hmm->get_explicit();
        
        state* init = hmm->getInitial();
		
		std::bitset<STATE_MAX>* initial_to = hmm->getInitialTo();
		std::bitset<STATE_MAX>* from_trans(NULL);
		
        //Calculate Viterbi from transitions from INIT (initial) state
        for(size_t i = 0; i < state_size; ++i){
            if ((*initial_to)[i]){  //if the bitset is set (meaning there is a transition to this state), calculate the viterbi
				
				viterbi_temp = (*hmm)[i]->get_emission_prob(*seqs,0) + getTransition(init, i, 0);
                
				if (viterbi_temp > -INFINITY){
                    if ((*scoring_current)[i] < viterbi_temp){
                        (*scoring_current)[i] = viterbi_temp;
                    }
					next_states[i] = 1;
                }
            }
        }
		
        
        for(size_t position = 1; position < seq_size ; ++position ){
            
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
			
			//TODO: Check use of external definitions below.
			
			std::cout << "\nPosition:\t" << position << "\n";
			//			std::cout << "Letter:\t" << seqs->seqValue(0, position) << std::endl;
			
            for (size_t i = 0; i < state_size; ++i){ //i is current state that emits value
                if (!current_states[i]){
                    continue;
                }
				
                //current_state = (*hmm)[i];
                //emission = current_state->get_emission(*seqs,position);
                emission = (*hmm)[i]->get_emission_prob(*seqs, position);
				
				
				//				std::cout << "State Emission:\t" << i << "\t" << exp(emission) << std::endl;
				
				if (exDef_defined && exDef_position){
                    emission += seqs->getWeight(position, i);
                }
                
				from_trans = (*hmm)[i]->getFrom();
				
                for (size_t j = 0; j < state_size ; ++j) {  //j is previous state
                    if (!(*from_trans)[j]){
                        continue;
                    }
					
                    if ((*scoring_previous)[j] != -INFINITY){
                        viterbi_temp = getTransition((*hmm)[j], i , position) + emission + (*scoring_previous)[j];
                        
						std::cout << exp(getTransition((*hmm)[j],i,position)) << std::endl;
						
						//						std::cout << "Temp Viterbi:\tTransFrom: "<< j << "\tto\t" << i << "\t" << viterbi_temp / log(2) << std::endl;
                        
						
						if (viterbi_temp > (*scoring_current)[i]){
							//If transition is from same to same then if it is
							//explicit duration we'll need to change
							extend_duration = (i==j) ? true : false;
							
                            (*scoring_current)[i] = viterbi_temp;
                            (*traceback_table)[position][i] = j;
                        }
						
						next_states |= (*(*hmm)[i]->getTo());
                    }
                }
				
				//If explicit durations vector defined, and transition from-to same
				//then we'll increment the value. Otherwise, set to zero;
				if (explicit_duration_current){
					if (extend_duration && (*duration)[i]){
						(*explicit_duration_current)[i]=(*explicit_duration_previous)[i]+1;
						extend_duration=false;
					}
					else{
						(*explicit_duration_current)[i]=0;
					}
				}
            }
			
			if(explicit_duration_current){
				swap_ptr_duration = explicit_duration_previous;
				explicit_duration_previous = explicit_duration_current;
				explicit_duration_current= swap_ptr_duration;
				explicit_duration_current->assign(state_size,0);
			}
            
        }
        
        //TODO:  Calculate ending and set the final viterbi and traceback pointer
        //Swap current and previous viterbi scores
        scoring_previous->assign(state_size,-INFINITY);
        swap_ptr = scoring_previous;
        scoring_previous = scoring_current;
        scoring_current = swap_ptr;
        
        for(size_t i = 0; i < state_size ;++i){
            if ((*scoring_previous)[i] > -INFINITY){
                viterbi_temp = (*scoring_previous)[i] + (*hmm)[i]->getEndTrans();
                
                if (viterbi_temp > ending_viterbi_score){
                    ending_viterbi_score = viterbi_temp;
                    ending_viterbi_tb = i;
                }
            }
        }
        
        delete scoring_previous;
        delete scoring_current;
		scoring_previous = NULL;
		scoring_current = NULL;
	}
	
	void trellis::viterbi(){
		
		//Initialize the traceback table
        if (traceback_table != NULL){
            delete traceback_table;
        }
        
        traceback_table = new int_2D(seq_size,std::vector<uint16_t> (state_size));
		scoring_previous = new (std::nothrow) std::vector<double> (state_size,-INFINITY);
        scoring_current  = new (std::nothrow) std::vector<double> (state_size,-INFINITY);
		
		if (scoring_previous == NULL || scoring_current == NULL || traceback_table == NULL){
			std::cerr << "Can't allocate Viterbi score and traceback table. OUT OF MEMORY" << std::endl;
			exit(2);
		}

		
        std::bitset<STATE_MAX> next_states;
        std::bitset<STATE_MAX> current_states;
		
        double  viterbi_temp(-INFINITY);
        double  emission(-INFINITY);
        bool	exDef_position(false);
		
		
		//If model is not a basic model, then we need to initialize the explicit duration vector
		//bool extend_duration keeps track of whether the transition to same state was selected.
		if (!hmm->isBasic()){
			explicit_duration_current = new(std::nothrow) std::vector<size_t>(state_size,0);
			explicit_duration_previous= new(std::nothrow) std::vector<size_t>(state_size,0);
		}
		bool extend_duration(false);
		std::vector<bool>* duration = hmm->get_explicit();
        
        state* init = hmm->getInitial();
		
		std::bitset<STATE_MAX>* initial_to = hmm->getInitialTo();
		std::bitset<STATE_MAX>* from_trans(NULL);
		
        //Calculate Viterbi from transitions from INIT (initial) state
        for(size_t i = 0; i < state_size; ++i){
            if ((*initial_to)[i]){  //if the bitset is set (meaning there is a transition to this state), calculate the viterbi
				
				viterbi_temp = (*hmm)[i]->get_emission_prob(*seqs,0) + getTransition(init, i, 0);
                
				if (viterbi_temp > -INFINITY){
                    if ((*scoring_current)[i] < viterbi_temp){
                        (*scoring_current)[i] = viterbi_temp;
                    }
					next_states[i] = 1;
                }
            }
        }
		
//		for(size_t i=0; i < state_size; ++i){
//			std::cout << "Position: 0" << std::endl;
//			std::cout << exp((*viterbi_current)[i]) << std::endl;
//		}
//        
        
        for(size_t position = 1; position < seq_size ; ++position ){
            
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
			
			//TODO: Check use of external definitions below.
			
			std::cout << "\nPosition:\t" << position << "\n";
//			std::cout << "Letter:\t" << seqs->seqValue(0, position) << std::endl;

            for (size_t i = 0; i < state_size; ++i){ //i is current state that emits value
                if (!current_states[i]){
                    continue;
                }
                				
                //current_state = (*hmm)[i];
                //emission = current_state->get_emission(*seqs,position);
                emission = (*hmm)[i]->get_emission_prob(*seqs, position);
				
				
//				std::cout << "State Emission:\t" << i << "\t" << exp(emission) << std::endl;
				
				if (exDef_defined && exDef_position){
                    emission += seqs->getWeight(position, i);
                }
                
				from_trans = (*hmm)[i]->getFrom();
				
                for (size_t j = 0; j < state_size ; ++j) {  //j is previous state
                    if (!(*from_trans)[j]){
                        continue;
                    }
					
                    if ((*scoring_previous)[j] != -INFINITY){
                        viterbi_temp = getTransition((*hmm)[j], i , position) + emission + (*scoring_previous)[j];
                        
						std::cout << exp(getTransition((*hmm)[j],i,position)) << std::endl;
						
//						std::cout << "Temp Viterbi:\tTransFrom: "<< j << "\tto\t" << i << "\t" << viterbi_temp / log(2) << std::endl;
                        
						
						if (viterbi_temp > (*scoring_current)[i]){
							//If transition is from same to same then if it is
							//explicit duration we'll need to change
							extend_duration = (i==j) ? true : false;
					
                            (*scoring_current)[i] = viterbi_temp;
                            (*traceback_table)[position][i] = j;
                        }
						
						next_states |= (*(*hmm)[i]->getTo());
                    }
                }
				
				//If explicit durations vector defined, and transition from-to same
				//then we'll increment the value. Otherwise, set to zero;
				if (explicit_duration_current){
					if (extend_duration && (*duration)[i]){
						(*explicit_duration_current)[i]=(*explicit_duration_previous)[i]+1;
						extend_duration=false;
					}
					else{
						(*explicit_duration_current)[i]=0;
					}
				}
            }
			
			if(explicit_duration_current){
				swap_ptr_duration = explicit_duration_previous;
				explicit_duration_previous = explicit_duration_current;
				explicit_duration_current= swap_ptr_duration;
				explicit_duration_current->assign(state_size,0);
			}
            
        }
        
        //TODO:  Calculate ending and set the final viterbi and traceback pointer
        //Swap current and previous viterbi scores
        scoring_previous->assign(state_size,-INFINITY);
        swap_ptr = scoring_previous;
        scoring_previous = scoring_current;
        scoring_current = swap_ptr;
        
        for(size_t i = 0; i < state_size ;++i){
            if ((*scoring_previous)[i] > -INFINITY){
                viterbi_temp = (*scoring_previous)[i] + (*hmm)[i]->getEndTrans();
                
                if (viterbi_temp > ending_viterbi_score){
                    ending_viterbi_score = viterbi_temp;
                    ending_viterbi_tb = i;
                }
            }
        }
        
        delete scoring_previous;
        delete scoring_current;
		scoring_previous = NULL;
		scoring_current = NULL;
	}
	
	
	void trellis::forward_viterbi(){
		
		//Initialize the traceback table
        if (traceback_table != NULL){
            delete traceback_table;
        }
        
        traceback_table = new int_2D(seq_size,std::vector<uint16_t> (state_size));
		
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
					next_states[i] = 1;
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
	
	
	//Fix: Need to fix so it calcuates value in double and stores in float
	void trellis::forward(){
		
		//Initialize forward score table
		forward_score = new (std::nothrow) float_2D(seq_size, std::vector<float>(state_size,-INFINITY));
		
		if (forward_score == NULL){
			std::cerr << "Can't allocate forward score table. OUT OF MEMORY" << std::endl;
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
					next_states[i] = 1;
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
	
	
	
	void trellis::forward(model* h, sequences* sqs){
		//Initialize the table
		hmm = h;
		seqs = sqs;
		seq_size		= seqs->getLength();
		state_size		= hmm->state_size();
		exDef_defined	= seqs->exDefDefined();
		
        forward();
    }
	
	
	void trellis::backward(model* h, sequences* sqs){
		//Initialize the table
		hmm = h;
		seqs = sqs;
		seq_size		= seqs->getLength();
		state_size		= hmm->state_size();
		exDef_defined	= seqs->exDefDefined();
		
        forward();
	}
	
	
	//TODO: Fix calculation in double not float (store in float)
	void trellis::backward(){
		//Initialize forward score table
		backward_score = new (std::nothrow) float_2D(seq_size, std::vector<float>(state_size,-INFINITY));
		if (backward_score == NULL){
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
		
		
//		std::cout << "Position: 3" << std::endl;
        //Calculate initial Backward from ending state
        for(size_t i = 0; i < state_size; ++i){
            if ((*ending_from)[i]){  //if the bitset is set (meaning there is a transition to this state), calculate the viterbi
				
				backward_temp = (*hmm)[i]->getEndTrans();
                
				if (backward_temp > -INFINITY){
                    
					(*backward_score)[seq_size-1][i] = backward_temp;
//					std::cout << "State: " << i << "\t" << exp(backward_temp) << std::endl;
					next_states[i] = 1;
                }
            }
        }
        
        
        for(size_t position = seq_size-1; position > 0 ; --position ){
            
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
					
                    if ((*backward_score)[position-1][j] != INFINITY){
						
//						double temp_trans = getTransition((*hmm)[j], i , position-1);
//						double temp_score = (*backward_score)[position][i];
//						
//						std::cout << "\nTransition from " << j << " to " << i << "\t" << exp(temp_trans) << std::endl;
//						std::cout << "Previous Score: " << exp(temp_score) << std::endl;
//						std::cout << "Emission: " << exp(emission) << std::endl;
//						backward_temp = temp_trans + emission + temp_score;
//						std::cout << "Temp Score: " << exp(backward_temp) << std::endl;
						
                        backward_temp = getTransition((*hmm)[j], i , position-1) + emission + (*backward_score)[position][i];
						
						if ((*backward_score)[position-1][j] == -INFINITY){
							(*backward_score)[position-1][j] = backward_temp;
						}
						else{
							(*backward_score)[position-1][j] = addLog((double)backward_temp, (double)(*backward_score)[position-1][j]);
						}
						
						next_states |= (*(*hmm)[i]->getFrom());
                    }
                }
//				std::cout << "State: " << i <<"\t" << exp((*backward_score)[position][i]) << std::endl;
            }
            
        }
		
//		std::cout << exp((*backward_score)[0][0]) << std::endl;
//		std::cout << exp((*backward_score)[0][1]) << std::endl;
		
		double backward_posterior = -INFINITY;
        state* init = hmm->getInitial();
        for(size_t i = 0; i < state_size ;++i){
            if ((*backward_score)[0][i] > -INFINITY){
                
				backward_temp = (*backward_score)[0][i] + (*hmm)[i]->get_emission_prob(*seqs,0) + getTransition(init, i, 0);
                
                if (backward_temp > -INFINITY){
					if (backward_posterior == -INFINITY){
						backward_posterior = backward_temp;
					}
					else{
						backward_posterior = addLog(backward_posterior,backward_temp);
					}
                }
            }
        }
		
//		std::cout << exp(backward_posterior) << std::endl;
	}
	
	
	stochTable::stochTable(size_t seq_size){
		
		state_val = new(std::nothrow) std::vector<stoch_value>;
		state_val->reserve(seq_size);
		
		position = new(std::nothrow) std::vector<size_t>(seq_size,0);
		if (state_val == NULL){
			std::cerr << "Cannot allocate stochTable- OUT OF MEMORY" << std::endl;
			delete state_val;
			delete position;
			state_val = NULL;
			position = NULL;
			exit(2);
		}
		last_position = 0;
		return;
	}
	
	stochTable::~stochTable(){
		delete state_val;
		delete position;
		state_val = NULL;
		position = NULL;
	}
	
	void stochTable::push(size_t pos, size_t st, size_t st_to, float value){
		if (pos != last_position){
			(*position)[last_position] = state_val->size()-1;
			last_position = pos;
		}
		
		state_val->push_back( stoch_value(st,st_to,value));
	}
	
	void stochTable::finalize(){
		(*position)[last_position] = state_val->size()-1;
		
		last_position = 0;
		
		double sum(-INFINITY);
		int current_state(-1);
		int state_start(-1);
		
		//For each position in the sequence we want to normalize the viterbi values
		//for the state from previous states
		for(size_t i=0;i<position->size();++i){
			
			//Need to sum the values.  If we see another state then we'll need to
			//normalize the previous states that contributed to the sum.
			//Then set the sum for the current state.
			for (size_t j = (i==0) ? 0 : (*position)[i-1]+1; j <= (*position)[i] ; ++j ){
				
				//If this is a different state then we'll need to apply the sum
				//if this is the first state then we'll just set the current state
				//and set the sum value.  Also need to keep track of which state
				// is the first state so when we normalize we only normalize the
				// the previous state
				if ((*state_val)[j].state_id != current_state){
					
					//Normalize
					for(size_t k = state_start; k < j ;++k){
						(*state_val)[k].prob = exp((*state_val)[k].prob - sum);
					}
					
					//Set state and sum value and starting state
					current_state = (*state_val)[j].state_id;
					sum = (*state_val)[j].prob;
					state_start = j;
				}
				else{
					sum = addLog(sum,(double) (*state_val)[j].prob);
				}
			}
			
			//Normalize the last states seen (b/c) we'll exit out of for loop before
			//we have done these
			for(size_t k = state_start; k <= (*position)[i]; ++k){
				(*state_val)[k].prob = exp((*state_val)[k].prob - sum);
			}
			
			//Set values for next loop
			state_start = (*position)[i]+1;
			current_state = -1;	
		}
		
		//Assign previous cell relative to (position) value.  Used when traceback
		//That way function won't have to search for correct value;
		for (size_t i =1; i < position->size(); ++i){
			for (size_t j = (*position)[i-1]+1; j <= (*position)[i] ; ++j){
				(*state_val)[j].prev_cell = get_state_position(i-1, (*state_val)[j].state_prev);
			}
		}
		
		return;
	}
	
	uint16_t stochTable::get_state_position(size_t pos, uint16_t state){
		size_t start_val = (pos == 0) ? 0 : (*position)[pos-1]+1;
		for (size_t i = start_val ; i <= (*position)[pos] ;  ++i){
			if ((*state_val)[i].state_id == state){
				return i-start_val;
			}
		}
		
		return UINT16_MAX;
	}
	
	void stochTable::print(){
		last_position = 0;
		
		for(size_t i=0;i<position->size();++i){
			for (size_t j = (i==0) ? 0 : (*position)[i-1]+1 ; j <= (*position)[i] ; ++j ){
				std::cout << (*state_val)[j].state_id <<":"<< (*state_val)[j].state_prev << " : " << (*state_val)[j].prob << "\t";
			}
			std::cout << std::endl;
		}
		return;
	}
	
	
	
	void trellis::stochastic_viterbi(){
		scoring_previous = new (std::nothrow) std::vector<double> (state_size,-INFINITY);
        scoring_current  = new (std::nothrow) std::vector<double> (state_size,-INFINITY);
		stochastic_table = new (std::nothrow) stochTable(seq_size);
        		
        std::bitset<STATE_MAX> next_states;
        std::bitset<STATE_MAX> current_states;
		
        double  viterbi_temp(-INFINITY);
        double  emission(-INFINITY);
        bool	exDef_position(false);
        
        state* init = hmm->getInitial();
		
		std::bitset<STATE_MAX>* initial_to = hmm->getInitialTo();
		std::bitset<STATE_MAX>* from_trans(NULL);
		
        //Calculate Viterbi from transitions from INIT (initial) state
        for(size_t i = 0; i < state_size; ++i){
            if ((*initial_to)[i]){  //if the bitset is set (meaning there is a transition to this state), calculate the viterbi
				
				viterbi_temp = (*hmm)[i]->get_emission_prob(*seqs,0) + getTransition(init, i, 0);;
                
				if (viterbi_temp > -INFINITY){										
                    if ((*scoring_current)[i] < viterbi_temp){
                        (*scoring_current)[i] = viterbi_temp;
                    }
					next_states[i] = 1;
                }
            }
        }
		
//		for(size_t i=0; i < state_size; ++i){
//			std::cout << "Position: 0" << std::endl;
//			std::cout << exp((*viterbi_current)[i]) << std::endl;
//		}
//
        
        for(size_t position = 1; position < seq_size ; ++position ){
            
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
			
			//			std::cout << "\nPosition:\t" << position << "\t";
			//			std::cout << "Letter:\t" << seqs->seqValue(0, position) << std::endl;
			
            for (size_t i = 0; i < state_size; ++i){ //i is current state that emits value
                if (!current_states[i]){
                    continue;
                }
				
                //current_state = (*hmm)[i];
                //emission = current_state->get_emission(*seqs,position);
                emission = (*hmm)[i]->get_emission_prob(*seqs, position);
				
				//				std::cout << "State Emission:\t" << i << "\t" << exp(emission) << std::endl;
				
				if (exDef_defined && exDef_position){
                    emission += seqs->getWeight(position, i);
                }
                
				from_trans = (*hmm)[i]->getFrom();
				
                for (size_t j = 0; j < state_size ; ++j) {  //j is previous state
                    if (!(*from_trans)[j]){
                        continue;
                    }
					
                    if ((*scoring_previous)[j] != -INFINITY){
                        viterbi_temp = getTransition((*hmm)[j], i , position) + emission + (*scoring_previous)[j];
                        
						
//						std::cout << "Temp Viterbi:\tTransFrom: "<< j << "\tto\t" << i << "\t" << viterbi_temp / log(2) << std::endl;
//Save partial value to stochastic table
                        stochastic_table->push(position-1,i,j,viterbi_temp);
						
						if (viterbi_temp > (*scoring_current)[i]){
							
                            (*scoring_current)[i] = viterbi_temp;
                        }
						
						next_states |= (*(*hmm)[i]->getTo());
                    }
                }
            }
        }
        
        //TODO:  Calculate ending and set the final viterbi and traceback pointer
        //Swap current and previous viterbi scores
        scoring_previous->assign(state_size,-INFINITY);
        swap_ptr = scoring_previous;
        scoring_previous = scoring_current;
        scoring_current = swap_ptr;
        
        for(size_t i = 0; i < state_size ;++i){
            if ((*scoring_previous)[i] > -INFINITY){
                viterbi_temp = (*scoring_previous)[i] + (*hmm)[i]->getEndTrans();
				
				stochastic_table->push(seq_size-1,SIZE_T_MAX, i,viterbi_temp);
                
                if (viterbi_temp > ending_viterbi_score){
                    ending_viterbi_score = viterbi_temp;
                    ending_viterbi_tb = i;
                }
            }
        }
		
		stochastic_table->finalize();
		stochastic_table->print();
        
        delete scoring_previous;
        delete scoring_current;
		scoring_current = NULL;
		scoring_previous = NULL;
	}
	
	
	
	double trellis::getTransition(state* st, size_t trans_to_state, size_t sequencePosition){
		        
        transition* trans = st->getTrans(trans_to_state);
		transType trans_type= trans->getTransitionType();
        
        double transition_prob=-INFINITY;
        
		
        if (trans_type == STANDARD ){  //if the transition type is standard then just return the standard probability
            transition_prob= trans->getTransition(0,NULL);
        }
        else if (trans_type == DURATION){
			
//            //TODO: Check traceback_length function
//			if ((*explicit_duration_current)[st->getIterator()] != 0 ){
//				transition_prob = trans->getTransition((*explicit_duration_current)[st->getIterator()]+1,NULL);
//			}
//			else{
				size_t size = get_explicit_duration_length(trans,sequencePosition, st->getIterator(), trans_to_state);
				transition_prob=trans->getTransition(size,NULL);
//				(*explicit_duration_current)[st->getIterator()]=size;
//			}
			
        }
        else if (trans_type == LEXICAL){
            transition_prob=trans->getTransition(sequencePosition, seqs);
        }
        
        //TODO:  Fix the exFuncTraceback(...)
        //Is external function define for the transition
        if (trans->FunctionDefined()){
//            transition_prob+=exFuncTraceback(trans->getExtFunction());
        }
        
        return transition_prob;		
	}
	
	size_t trellis::get_explicit_duration_length(transition* trans, size_t sequencePosition, size_t state_iter, size_t to_state){
		
		if ((*explicit_duration_previous)[state_iter]!=0){
			return (*explicit_duration_previous)[state_iter]+1;
		}
		
		
		//If it hasn't been defined then traceback until the ending parameter is reached
		
		size_t length(0);
	
		
		size_t tbState(state_iter);  //First traceback pointer to use in traceback_table
		
		//tracebackIdentifier traceback_identifier = previousState->transi[transitionTo].traceback_identifier;
		tracebackIdentifier traceback_identifier = trans->getTracebackIdentifier();
		
		//string identifier = previousState->transi[transitionTo].traceback_string;
		std::string identifier = trans->getTracebackString();
		
		
		for(size_t trellPos=sequencePosition-1 ; trellPos != SIZE_MAX ;trellPos--){
			//for(;trellisPos>=0;trellisPos--){
			length++;
			//state=trellis.trell[trellisPos][state].ptr;  //Get previous state traceback
			
			tbState = (*traceback_table)[trellPos][tbState];
			state* st = hmm->getState(tbState);
			
			//Check to see if stop conditions of traceback are met, if so break;
			if(traceback_identifier == START_INIT && tbState == -1) {break;}
			else if (traceback_identifier == DIFF_STATE  && state_iter != st->getIterator())	{ break;}
			else if (traceback_identifier == STATE_NAME  && identifier.compare(st->getName())==0)	{ break;}
			else if (traceback_identifier == STATE_LABEL && identifier.compare(st->getLabel())==0)	{ break;}
			else if (traceback_identifier == STATE_GFF   && identifier.compare(st->getGFF())==0)	{ break;}
			
		}
		return length+1;

	}
	
	//!Perform traceback through trellis
    //!\return [out] path trackback_path
    void trellis::traceback(traceback_path& path){
				
		if (seq_size==0){
			return;
		}
		
		if (ending_viterbi_score == -INFINITY){
			return;
		}
        else{
            path.setScore(ending_viterbi_score);
            path.push_back(ending_viterbi_tb);
            
			uint16_t pointer = ending_viterbi_tb;
			
            for( size_t position = seq_size -1 ; position>0 ; position--){
				pointer = (*traceback_table)[position][pointer];
				path.push_back(pointer);
            }
        }
        return;
    }
	
	void trellis::stochastic_traceback(traceback_path& path){
		
		stochastic_table->traceback(path);
		return;
	}
	
	void stochTable::traceback(traceback_path& path){
		
		double random((double)rand()/((double)(RAND_MAX)+(double)(1)));
		double cumulative_prob(0.0);
		std::cout << random << std::endl;
		
		size_t offset(0);
		uint16_t state_prev;
		
		//Get traceback from END state
		for(size_t i = (*position)[position->size()-2]+1 ; i <= position->back() ; ++i){
			cumulative_prob += (*state_val)[i].prob;
			if (random < cumulative_prob){
				state_prev = (*state_val)[i].state_prev;
				path.push_back(state_prev);
				offset = (*state_val)[i].prev_cell;
				std::cout << "Chose:\t" << state_prev << std::endl;
				break;
			}
		}
		
		for(size_t i = position->size()-2; i != SIZE_T_MAX ; --i){
			random = (double)rand()/((double)(RAND_MAX)+(double)(1));
			std::cout << random << std::endl;
			cumulative_prob = 0;
			
			for (size_t cells = (i == 0) ? 0 + offset : (*position)[i-1]+offset+1; cells < (*position)[i] ; ++cells){
				if ((*state_val)[cells].state_id != state_prev ){
					std::cout << "Houston, We have an error!" << std::endl;
				}
                
				cumulative_prob+=(*state_val)[cells].prob;
				
                if (random <= cumulative_prob){
                    state_prev = (*state_val)[cells].state_prev;
					path.push_back(state_prev);
					offset = (*state_val)[cells].prev_cell;
					std::cout << "Chose:\t" << state_prev << std::endl;
					break;
                }
            }
		}
		return;
	}
	
}




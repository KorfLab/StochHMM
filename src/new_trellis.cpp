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
		naive_forward_score = NULL;
		naive_viterbi_score = NULL;
		naive_backward_score= NULL;
		
		ending_viterbi_score = -INFINITY;
		ending_viterbi_tb = -1;
//		ending_posterior = -INFINITY;
		ending_forward_prob = -INFINITY;
		ending_backward_prob= -INFINITY;
		
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
		naive_forward_score = NULL;
		naive_viterbi_score = NULL;
		naive_backward_score= NULL;
		
		ending_viterbi_score = -INFINITY;
		ending_viterbi_tb = -1;
//		ending_posterior = -INFINITY;
		ending_forward_prob = -INFINITY;
		ending_backward_prob= -INFINITY;

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
		delete naive_forward_score;
		delete naive_viterbi_score;
		delete naive_backward_score;
		
		delete scoring_previous;
		delete scoring_current;
		delete alt_scoring_previous;
		delete alt_scoring_current;
		
		traceback_table		= NULL;
		stochastic_table	= NULL;
		nth_traceback		= NULL;
		viterbi_score		= NULL;
		forward_score		= NULL;
		backward_score		= NULL;
		posterior_score		= NULL;
		naive_forward_score	= NULL;
		naive_viterbi_score	= NULL;
		naive_backward_score= NULL;
		
		scoring_previous	= NULL;
		scoring_current		= NULL;
		alt_scoring_previous= NULL;
		alt_scoring_current	= NULL;
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
		
		delete scoring_current;
		delete scoring_previous;
		
		delete alt_scoring_current;
		delete alt_scoring_previous;
		
		delete naive_forward_score;
		delete naive_backward_score;
		delete naive_viterbi_score;
		
		traceback_table		= NULL;
		stochastic_table	= NULL;
		nth_traceback		= NULL;
		
		viterbi_score		= NULL;
		forward_score		= NULL;
		backward_score		= NULL;
		posterior_score		= NULL;
		
		scoring_current		= NULL;
		scoring_previous	= NULL;
		
		alt_scoring_current	= NULL;
		alt_scoring_previous= NULL;
		
		naive_forward_score	= NULL;
		naive_backward_score= NULL;
		naive_viterbi_score	= NULL;
		
		ending_viterbi_score = -INFINITY;
		ending_viterbi_tb = -1;
//		ending_posterior = -INFINITY;
		ending_forward_prob = -INFINITY;
		ending_backward_prob= -INFINITY;
		
		
	}
	
	
	double trellis::getTransition(state* st, size_t trans_to_state, size_t sequencePosition){
		double transition_prob(-INFINITY);
        transition* trans = st->getTrans(trans_to_state);
		if (trans==NULL){
			return transition_prob;
		}
				
		transType trans_type= trans->getTransitionType();
        
        
        
		
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

	

	
}




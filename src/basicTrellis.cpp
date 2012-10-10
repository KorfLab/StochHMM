//
//  basicTrellis.cpp

//Copyright (c) 2007-2012 Paul C Lott 
//University of California, Davis
//Genome and Biomedical Sciences Facility
//UC Davis Genome Center
//Ian Korf Lab
//Website: www.korflab.ucdavis.edu
//Email: lottpaul@gmail.com
//
//Permission is hereby granted, free of charge, to any person obtaining a copy of
//this software and associated documentation files (the "Software"), to deal in
//the Software without restriction, including without limitation the rights to
//use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of
//the Software, and to permit persons to whom the Software is furnished to do so,
//subject to the following conditions:
//
//The above copyright notice and this permission notice shall be included in all
//copies or substantial portions of the Software.
//
//THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
//IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS
//FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR
//COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER
//IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN
//CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

#include "basicTrellis.h"
namespace StochHMM{

    //!Create basicTrellis class
    //!\param modl  Pointer to model
    //!\param seqs  Pointer to sequences
    basicTrellis::basicTrellis(model* modl, sequences* seqs){
        
        //TODO: Check that modl or seqs are valid and not NULL
        hmm=modl;
        seq=seqs;
        
        posterior=NULL;
        

        probabilityOfSequence=-INFINITY;
        
        sequencePosition=0; //Start of sequence
        currentStatePtr=NULL;  //Undefined state
        previousStatePtr=NULL; //Undefined state
        
        currentState = -1;
        previousState = -1;
        
        //currentPrelim = 0;
        
        currentProbability=-INFINITY;
        DefinedStates=new(std::nothrow) std::set<state*>;
        previousStates=NULL;
        
        currentStates=NULL;
        _nextStates= new(std::nothrow) std::set<state*>;
        
        externalDefinitionSet=false;
        if (DefinedStates==NULL || _nextStates==NULL){
            std::cerr << "OUT OF MEMORY\nFile" << __FILE__ << "Line:\t"<< __LINE__ << std::endl;
            exit(1);
        }
    }
    
    
    //! Reset the counter information in basicTrellis                                                              
    //! The counter information keeps track of what state we're in and what position. We also keep track in next and what are the states to be evaluated next. Previous has the states that were valid in the previous position of the sequence
    void basicTrellis::resetCounters(){
        currentStatePtr=NULL;
        previousStatePtr=NULL;

        sequencePosition=0;
        currentStates=NULL;
        if (DefinedStates!=NULL){
            delete DefinedStates;
        }
        if (previousStates!=NULL){
            delete previousStates;
        }
        
        DefinedStates=new(std::nothrow) std::set<state*>;
        
        if (DefinedStates==NULL){
            std::cerr << "OUT OF MEMORY\nFile" << __FILE__ << "Line:\t"<< __LINE__ << std::endl;
            exit(1);
        }
        
        previousStates=NULL;
    }
    
    
    
    //! Move to next - moves to the next position in the trellis table.  If there are no more states to be evaluated in the current sequence position then we move to the next position and also change previous to next and previous to new set.
    //! If there is another state to be evaluated then we move to that state
    //! Check external definitions here and change Previous set accordingly for Absolute for Weighted we'll apply to the scores.
    bool basicTrellis::moveToNext(){
        
        // previous should only be null right after initialized
        if (previousStates!=NULL){
            delete previousStates;
            delete currentStates;
        }
        
        if (seq->exDefDefined()){
            externalDefinitionSet=seq->exDefDefined(sequencePosition);
        }
        
        previousStates=DefinedStates;
        DefinedStates=new(std::nothrow) std::set<state*>;
        
        
        currentStates = _nextStates;
        _nextStates = new(std::nothrow) std::set<state*>;
        
        if (DefinedStates==NULL || _nextStates==NULL){
            std::cerr << "OUT OF MEMORY\nFile" << __FILE__ << "Line:\t"<< __LINE__ << std::endl;
            exit(1);
        }
        
        
        sequencePosition++;
        
        
        //If the size of the set is zero of if we have gotten to the end of the sequence return false
        if (sequencePosition==seq->getLength()){
            return false;
        }    
        else if (currentStates->size()==0){
            
            std::cout << "ended before end.  Not valid grammar" << std::endl;
            return false;
        }
        
        return true;
    }
    
    
    
    //! Similar to moveToNext() but moving backward for backward algorithm                                                               
    // TODO: Test implementation
    bool basicTrellis::moveToPrevious(){
        
        // previous should only be null right after initialized
        if (previousStates!=NULL){
            delete previousStates;
            delete currentStates;
            
        }
            
        if (seq->exDefDefined()){
            externalDefinitionSet = seq->exDefDefined(sequencePosition-1);
        }
            
        previousStates=DefinedStates;
        DefinedStates=new(std::nothrow) std::set<state*>;
        
        currentStates = _nextStates;
        _nextStates = new(std::nothrow) std::set<state*>;
        
        if (DefinedStates==NULL|| _nextStates==NULL){
            std::cerr << "OUT OF MEMORY\nFile" << __FILE__ << "Line:\t"<< __LINE__ << std::endl;
            exit(1);
        }
        
        sequencePosition--;
        
        //If the size of the set is zero of if we have gotten to the end of the sequence return false
        if (sequencePosition==0){
            return false;
        }    
        else if (currentStates->size()==0){
            std::cout << "ended before end.  Not valid grammar" <<std::endl;
            return false;
        }
        
        return true;
    }
    
    
    
    //! Performs the initial calculations for the forward algorithm                                                              
    void basicTrellis::initForward(){
        sequencePosition=0;
        std::vector<state*>* init=hmm->getInitialTo();
        state* initial=hmm->getInitial();
        
        
        for(size_t x=0; x<init->size(); x++){
            currentStatePtr=(*init)[x];
            currentState = currentStatePtr->getIterator();
            
            double prelim = ((*initial).getTrans(currentState))->getTransition(0,NULL) + currentStatePtr->get_emission(*seq,sequencePosition);
            
            if (prelim>-INFINITY){
                DefinedStates->insert(currentStatePtr);
                _nextStates->insert(currentStatePtr->getToBegin(),currentStatePtr->getToEnd());
                setForward(prelim);
            }
            
            
        }
        
        
        return;
    }
    
    //!Performs the initial calculations for the viterbi algorithm                                                            
    void basicTrellis::initViterbi(){
        sequencePosition=0;
        previousStatePtr=NULL;
        std::vector<state*>* init = hmm->getInitialTo();
        state* initial=hmm->getInitial();
        
        for (size_t x=0; x<init->size();x++){
            currentStatePtr=(*init)[x];
            currentState = currentStatePtr->getIterator();
            
            double viterbiValue = ((*initial).getTrans(currentState))->getTransition(0,NULL) + currentStatePtr->get_emission(*seq,sequencePosition);
            
            if (viterbiValue>-INFINITY){
                DefinedStates->insert(currentStatePtr);
                _nextStates->insert(currentStatePtr->getToBegin(),currentStatePtr->getToEnd());
                setViterbi(viterbiValue);
            }
            
#ifdef VERBOSE
            std::cout << "Position:\t" << sequencePosition <<std::endl;
            std::cout << "State:\t" << currentState << std::endl;
            std::cout << "Previous State:\tSTART" << std::endl;
            std::cout << "Transition:\t" << ((*initial).getTrans(currentState))->getTransition(NULL,NULL) <<std::endl;
            std::cout << "Emission:\t" << hmm->getState(currentState)->get_emission(*seq,sequencePosition) << std::endl<<std::endl;
#endif
        }
        
        return;
    }
    
    
    //!Performs the initial calculations for the viterbi algorithm
    void basicTrellis::initNthViterbi(size_t n){
                
        sequencePosition=0;
        previousStatePtr=NULL;
        std::vector<state*>* init = hmm->getInitialTo();
        state* initial=hmm->getInitial();
        
        for (size_t x=0; x<init->size();x++){
            currentStatePtr=(*init)[x];
            currentState = currentStatePtr->getIterator();
            
            double viterbiValue = ((*initial).getTrans(currentState))->getTransition(0,NULL) + currentStatePtr->get_emission(*seq,sequencePosition);
            
            if (viterbiValue>-INFINITY){
                DefinedStates->insert(currentStatePtr);
                _nextStates->insert(currentStatePtr->getToBegin(),currentStatePtr->getToEnd());
                setNthViterbi(viterbiValue);
            }
            
#ifdef VERBOSE
            std::cout << "Position:\t" << sequencePosition <<std::endl;
            std::cout << "State:\t" << currentState << std::endl;
            std::cout << "Previous State:\tSTART" << std::endl;
            std::cout << "Transition:\t" << ((*initial).getTrans(currentState))->getTransition(0,NULL) <<std::endl;
            std::cout << "Emission:\t" << hmm->getState(currentState)->get_emission(*seq,sequencePosition) << std::endl<<std::endl;
#endif
        }
        
        return;
    }
    
    
//    //!Performs the initial calculations for the viterbi algorithm
//    void basicTrellis::initNthViterbi(size_t n){
//        
//        
//        sequencePosition=0;
//        previousStatePtr=NULL;
//        std::vector<state*>* init = hmm->getInitialTo();
//        state* initial=hmm->getInitial();
//        
//        for (unsigned int x=0; x<init->size();x++){
//            currentStatePtr=(*init)[x];
//            currentState = currentStatePtr->getIterator();
//            
//            double viterbiValue = ((*initial).getTrans(currentState))->getTransition(NULL,NULL) + currentStatePtr->get_emission(*seq,sequencePosition);
//            
//            if (viterbiValue>-INFINITY){
//                DefinedStates->insert(currentStatePtr);
//                _nextStates->insert(currentStatePtr->getToBegin(),currentStatePtr->getToEnd());
//                setViterbi(viterbiValue);
//            }
//            
//#ifdef VERBOSE
//            std::cout << "Position:\t" << sequencePosition <<std::endl;
//            std::cout << "State:\t" << currentState << std::endl;
//            std::cout << "Previous State:\tSTART" << std::endl;
//            std::cout << "Transition:\t" << ((*initial).getTrans(currentState))->getTransition(NULL,NULL) <<std::endl;
//            std::cout << "Emission:\t" << hmm->getState(currentState)->get_emission(*seq,sequencePosition) << std::endl<<std::endl;
//#endif
//        }
//        
//        return;
//    }
    
    
    
    //!Perform calculations from initial state for viterbi and forward                                                           
    void basicTrellis::initForwardViterbi(){
        sequencePosition=0;
        previousStatePtr=NULL;
        std::vector<state*>* init = hmm->getInitialTo();
        state* initial=hmm->getInitial();
        
        for (unsigned int x=0; x<init->size();x++){
            currentStatePtr=(*init)[x];
            currentState = currentStatePtr->getIterator();

            double viterbiValue = ((*initial).getTrans(currentState))->getTransition(0,NULL) + currentStatePtr->get_emission(*seq,sequencePosition);
            
            if (viterbiValue>-INFINITY){
                DefinedStates->insert(currentStatePtr);
                _nextStates->insert(currentStatePtr->getToBegin(),currentStatePtr->getToEnd());
                setViterbi(viterbiValue);
                setForward(viterbiValue);
            }
        }
        
        return;
    }
    
    //!Perform initial calculation for Backward algorithm from ending state                                                            
    void basicTrellis::initBackward(){
        
        sequencePosition=seq->getLength();
        
        std::vector<state*>* ending=hmm->getEndingFrom();
        for(unsigned int x=0;x<ending->size();x++){
            currentStatePtr=(*ending)[x];
            currentState = currentStatePtr->getIterator();
            
            double backScore = getEndingTransition(currentState);
            
            setBackward(sequencePosition-1,currentState,backScore);
            
            DefinedStates->insert(currentStatePtr);
            _nextStates->insert(currentStatePtr->getFromBegin(),currentStatePtr->getFromEnd());
            
        }
        
        return;
    }
    
    
    
    //!Get ending transition probability                                                          
    double basicTrellis::getEndingTransition(size_t transitionFrom){
        static state* currState = hmm->getState(transitionFrom);
        return (currState->getEnding())->getTransition(0,NULL);
        
        //static state* ending = hmm->getEnding();
        //return (ending->getTrans(transitionFrom))->getTransition(NULL,NULL);
        //return model->states[transitionFrom]->endi.log_trans;
    }



    
    //!Get transition probability                                                              
    //! \return double log(p(x)) of transition from state (st) to state (transitionTo)
    double basicTrellis::getTransition(){
        
        
        transition* trans=previousStatePtr->getTrans(currentState);
        
        transType trans_type= trans->getTransitionType();
        
        double transition_prob=-INFINITY;
        
           
        if (trans_type == STANDARD ){  //if the transition type is standard then just return the standard probability
            transition_prob= trans->getTransition(0,NULL);
        }
        else if (trans_type == DURATION){
            
            //TODO: Check traceback_length function
            int size=traceback_length();
            
            transition_prob=trans->getTransition(size,NULL);
            
        }
        else if (trans_type == LEXICAL){
            transition_prob=trans->getTransition(sequencePosition, seq);
        }
        
        
        //Is external function define for the transition
        if (trans->FunctionDefined()){
            transition_prob+=exFuncTraceback(trans->getExtFunction());
        }
        
        return transition_prob;
    }




    //!Traceback from given position and return the distance to the ending requirement
    //! The probability can then be extracted from the transition distribution
    int basicTrellis::traceback_length(){
        
        //sequencePosition, previousStatePtr, currentState, this

        
        int length=0; 
        
        int tbState=previousState;  //Starting state to use
        int starting_state=tbState;
        
        transition* trans = previousStatePtr->getTrans(currentState);
        
        //tracebackIdentifier traceback_identifier = previousState->transi[transitionTo].traceback_identifier;  
        tracebackIdentifier traceback_identifier = previousStatePtr->getTrans(currentState)->getTracebackIdentifier();
        
        //string identifier = previousState->transi[transitionTo].traceback_string;
        std::string identifier = trans->getTracebackString();
        
        
        for(size_t trellPos=sequencePosition-1 ; trellPos != SIZE_MAX ;trellPos--){
        //for(;trellisPos>=0;trellisPos--){
            length++;
            //state=trellis.trell[trellisPos][state].ptr;  //Get previous state traceback
            
            tbState=this->getPtr(trellPos,tbState);
            state* st = hmm->getState(tbState);
            
            //Check to see if stop conditions of traceback are met, if so break;
            if(traceback_identifier == START_INIT && tbState == -1) {break;}
            else if (traceback_identifier == DIFF_STATE  && starting_state != st->getIterator()) {  break;}
            else if (traceback_identifier == STATE_NAME  && identifier.compare(st->getName())==0){ break;}
            else if (traceback_identifier == STATE_LABEL && identifier.compare(st->getLabel())==0) {  break;}
            else if (traceback_identifier == STATE_GFF   && identifier.compare(st->getGFF())==0) {  break;}
            
        }        
        return length;                
    }



    //----------------------------------------------------------------------------//
    // Description: hmmerTraceback  performs a traceback and collects the sequence                                                               
    // at the given states and returns a vector sequence.   This can then be translated
    // and have HMMER search performed on the sequence
    // 
    //----------------------------------------------------------------------------//
    double basicTrellis::exFuncTraceback(transitionFuncParam* func){
        
        std::vector<int> tracebackPath;
        std::vector<std::string> tracebackString;

        
        int tbState=previousState;  //Starting state to use
        int starting_state=tbState;
        
        //How far to traceback
        tracebackIdentifier traceback_identifier = func->getTracebackType();
        const std::string& tracebackIdentifierName = func->getTracebackName();
        
        
        //What to combine
        combineIdentifier combineIdent = func->getCombineType();
        const std::string& combineIdentName = func->getCombineName();
        
        
        //Deterimine which track to use
        track* alphaTrack = func->getTrack();
        size_t trackIndex = alphaTrack->getIndex();
        if (!alphaTrack->isAlpha()){
            
            //TODO:  Throw Error
            
            return 0;
        }
        
        const sequence* sq = seq->getSeq(trackIndex);
        state* st = hmm->getState(tbState);
        
                
        for(size_t trellisPos = sequencePosition-1; trellisPos != SIZE_MAX ; trellisPos--){
            
            tracebackPath.push_back(tbState);
            
            if ((combineIdent == FULL) || 
                (combineIdent == STATENAME && combineIdentName.compare(st->getName())==0)||
                (combineIdent == STATELABEL && combineIdentName.compare(st->getLabel())==0)||
                (combineIdent == STATEGFF && combineIdentName.compare(st->getGFF())==0))
            {
                tracebackString.push_back(sq->getSymbol(trellisPos));
            }
        
            tbState=this->getPtr(tbState,trellisPos);
            state* st = hmm->getState(tbState);
            
            //Check to see if stop conditions of traceback are met, if so break;
            if(traceback_identifier == START_INIT && tbState == -1) {break;}
            else if (traceback_identifier == DIFF_STATE  && starting_state != st->getIterator()) {  break;}
            else if (traceback_identifier == STATE_NAME  && tracebackIdentifierName.compare(st->getName())==0){ break;}
            else if (traceback_identifier == STATE_LABEL && tracebackIdentifierName.compare(st->getLabel())==0) {  break;}
            else if (traceback_identifier == STATE_GFF   && tracebackIdentifierName.compare(st->getGFF())==0) {  break;}
        }        
        
        size_t length=tracebackPath.size();
        
        
        std::string CombinedString;
        size_t maxSymbolSize = alphaTrack->getAlphaMax();
        if (maxSymbolSize ==1){
            for(std::vector<std::string>::reverse_iterator rit = tracebackString.rbegin(); rit!=tracebackString.rend();++rit){
                CombinedString+=(*rit);
            }
        }
        
        std::string temp_string = seq->getUndigitized(trackIndex);
        double transitionValue = func->evaluate(&temp_string, sequencePosition, &CombinedString, length);
        
        return transitionValue;
    }
    
    
    //!Tracing posterior not currently supported on simpleTrellis
    void basicTrellis::tracePosterior(traceback_path& path){
        
        size_t seq_length = seq->getLength();
        if (seq_length==0){
            return;
        }
        
        for (size_t position = seq_length-1;position!=SIZE_MAX;--position){
            path.push_back((*posterior_pointer)[position]);
        }
        return;
    }
    
    void basicTrellis::traceStochPosterior(traceback_path& path){
        
        size_t states = hmm->state_size();
        size_t seq_length = seq->getLength();
        if (seq_length==0){
            return;
        }
        
        for (size_t position = seq_length-1;position!=SIZE_MAX;--position){
            
            double random=((double)rand()/((double)(RAND_MAX)+(double)(1)));
            double cumulative_prob(0.0);
            for (int state = 0; state<states ; ++state){
                cumulative_prob+=(*posterior)[position][state];
                if (random<cumulative_prob){
                    path.push_back(state);
                    break;
                }
            }
        }
        return;
    }

    
}




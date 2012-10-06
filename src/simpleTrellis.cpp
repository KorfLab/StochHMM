//
//  simpleTrellis.cpp

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

#include "simpleTrellis.h"
namespace StochHMM{

    //!Create simpleTrellis
    //!\param modl  Model file to use for decoding
    //!\param seqs  Sequences to use for decoding
    simpleTrellis::simpleTrellis(model* modl,sequences* seqs): basicTrellis(modl,seqs){
        size_t stateSize=modl->state_size();
        size_t seqSize=seq->getLength();
        
        simpleCell initial;
        //try {
            trell.assign(seqSize,std::vector<simpleCell>(stateSize,initial));
//        } catch (std::bad_alloc) {
//            error(sOutOfMemory);
//        }
        
        return;
    }
    
    //TODO: Create a function to save and import trellis file
    
    std::string simpleTrellis::stringify(){
        size_t state_size = hmm->state_size();
        size_t seq_size = seq->getLength();
        
        std::string output("");
        
        //Header
        output+= ",,,,,,,,,,Stochastic Viterbi Scores, Stochastic Forward Scores\n";
        output+= "Position, State Name, State Label, Emission, Transition, Forward, Backward, Posterior, Viterbi, Viterbi Traceback Pointer";
        
        std::string state_names("");
        for(size_t i=0;i<state_size;i++){
            state_names+="," + hmm->getStateName(i);
        }
        
        output+= state_names + state_names + "\n";
        
        for(size_t position = 0 ; position<seq_size;++position){
            for (size_t state_iter = 0 ; state_iter<state_size;++state_iter){
                simpleCell& cell = trell[position][state_iter];
                
                if (cell.forw!=-INFINITY || cell.viti!=-INFINITY || cell.back!=-INFINITY){
                    output += position + ",";
                    output += hmm->getStateName(state_iter) + ",";
                    output += hmm->getStateLabel(state_iter) + ",";
                    output += double_to_string(cell.emm) + ",";
                    
                    //Need to do state transitions
                    
                    output += double_to_string(cell.forw) + ",";
                    output += double_to_string(cell.back) + ",";
                    output += double_to_string(cell.viti) + ",";
                    if (cell.statePtr != NULL){
                         output += cell.statePtr->getName();
                    }
                   
                    
                }
            }
        }
        
        return output;
    }
    
    //!Print the trellis to stdout
    void simpleTrellis::print(){
        
        //Print Table
        for(int i=0;i<seq->getLength();i++){
            for(int j=0;j<hmm->state_size();j++){
                std::cout << "Position: " << i << "\t" << "State: " << j << std::endl;
                trell[i][j].print();
                std::cout << std::endl;
            }
        }
        return;
    }
    
    //!Calculate the viterbi and forward scores for a given cell                                                            
    void simpleTrellis::calcForwardViterbi(){
        
        for(currentStatesIterator = currentStates->begin(); currentStatesIterator!=currentStates->end();currentStatesIterator++){
            currentStatePtr = (*currentStatesIterator);
            currentState = currentStatePtr->getIterator();
            double emm;
            
            //If emission is not yet calculated for the cell
            if (!trell[sequencePosition][currentState].emmCalculated){
                
                //Calculate emission
                emm=currentStatePtr->get_emission(*seq, sequencePosition);
                
                //Check for external definition
                if (externalDefinitionSet){
                    emm+=seq->getWeight(sequencePosition, currentState);
                }
                
                //Store emission value in cell
                trell[sequencePosition][currentState].emmCalculated=true;
                trell[sequencePosition][currentState].emm=emm;
            }
            else{  //emm previously calculated
                
                //Get previously calculated emm value;
                emm = trell[sequencePosition][currentState].emm;
            }
            
            std::vector<state*>::iterator it;
            for(it=currentStatePtr->getFromBegin();it!=currentStatePtr->getFromEnd();it++){
                previousStatePtr=(*it);
                
                if (previousStates->count(previousStatePtr)){
                    previousState = previousStatePtr->getIterator();
                    
                    
                    double trans;
                    if (trell[sequencePosition][currentState].trans.count(previousStatePtr)){
                        trans = trell[sequencePosition][currentState].trans[previousStatePtr];
                    }
                    else{
                        trans = getTransition();
                        trell[sequencePosition][currentState].trans[previousStatePtr]=trans;
                        //trans = getTransition(sequencePosition, previousStatePtr, currentState, this);
                    }
                    
#ifdef VERBOSE
                    std::cout << "Position:\t" << sequencePosition <<std::endl;
                    std::cout << "State:\t" << currentState <<std::endl;
                    std::cout << "Previous State:\t"<<previousState << std::endl;
                    std::cout << "Transition:\t" << trans <<std::endl;
                    std::cout << "Emission:\t" << emm << std::endl;
#endif
                    double viterbiValue = emm + trans + trell[sequencePosition-1][previousState].viti;
                    double prelimValue = emm + trans + trell[sequencePosition-1][previousState].forw;
                    
                    if(viterbiValue>-INFINITY){
                        DefinedStates->insert(currentStatePtr);
                        _nextStates->insert(currentStatePtr->getToBegin(),currentStatePtr->getToEnd());
                        
                        setViterbi(viterbiValue);
                        setForward(prelimValue);
                    }
                }
            }
        }
#ifdef VERBOSE
        std::cout <<std::endl;
#endif
      
        return;
    }
    
    //!Calculate the forward score for a given cell
    void simpleTrellis::calcForward(){
        
        for(currentStatesIterator = currentStates->begin(); currentStatesIterator!=currentStates->end();currentStatesIterator++){
            currentStatePtr = (*currentStatesIterator);
            currentState = currentStatePtr->getIterator();
            double emm;
            
            //If emission is not yet calculated for the cell
            if (!trell[sequencePosition][currentState].emmCalculated){
                
                //Calculate emission
                emm=currentStatePtr->get_emission(*seq, sequencePosition);
                
                //Check for external definition
                if (externalDefinitionSet){
                    emm+=seq->getWeight(sequencePosition, currentState);
                }
                
                //Store emission value in cell
                trell[sequencePosition][currentState].emmCalculated=true;
                trell[sequencePosition][currentState].emm=emm;
            }
            else{  //emm previously calculated
                
                //Get previously calculated emm value;
                emm = trell[sequencePosition][currentState].emm;
            }
            
            std::vector<state*>::iterator it;
            
            for(it=currentStatePtr->getFromBegin();it!=currentStatePtr->getFromEnd();it++){
                previousStatePtr=(*it);
                
                
                if (previousStates->count(previousStatePtr)){
                    previousState = previousStatePtr->getIterator();
                    
                    
                    double trans;
                    if (trell[sequencePosition][currentState].trans.count(previousStatePtr)){
                        trans = trell[sequencePosition][currentState].trans[previousStatePtr];
                    }
                    else{
                        trans = getTransition();
                        trell[sequencePosition][currentState].trans[previousStatePtr]=trans;
                    }
                    
#ifdef VERBOSE
                    std::cout << "Position:\t" << sequencePosition <<std::endl;
                    std::cout << "State:\t" << currentState <<std::endl;
                    std::cout << "Previous State:\t"<<previousState << std::endl;
                    std::cout << "Transition:\t" << trans <<std::endl;
                    std::cout << "Emission:\t" << emm << std::endl;
#endif
                    double prelimValue = emm + trans + trell[sequencePosition-1][previousState].forw;
                    
                    if(prelimValue>-INFINITY){
                        DefinedStates->insert(currentStatePtr);
                        _nextStates->insert(currentStatePtr->getToBegin(),currentStatePtr->getToEnd());

                        setForward(prelimValue);
                    }
                }
            }
        }
#ifdef VERBOSE
        std::cout <<std::endl;
#endif
        return;
    }
    
    
    //!Calculate viterbi score for give cell from previous cell                                                             
    void simpleTrellis::calcViterbi(){
        
        //For all the possible states at this position
        for(currentStatesIterator = currentStates->begin(); currentStatesIterator!=currentStates->end();currentStatesIterator++){
            //Get pointer to current state
            currentStatePtr = (*currentStatesIterator);
            currentState = currentStatePtr->getIterator();
            
            //Calculate emission for the current cell
            double emm;
            //If emission is not yet calculated for the cell
            if (!trell[sequencePosition][currentState].emmCalculated){
                
                //Calculate emission
                emm=currentStatePtr->get_emission(*seq, sequencePosition);
                
                //Check for external definition
                if (externalDefinitionSet){
                    emm+=seq->getWeight(sequencePosition, currentState);
                }
                
                //Store emission value in cell
                trell[sequencePosition][currentState].emmCalculated=true;
                trell[sequencePosition][currentState].emm=emm;
            }
            else{  //emm was previously calculated
                
                //Get previously calculated emm value;
                emm = trell[sequencePosition][currentState].emm;
            }
            
            std::vector<state*>::iterator it;
            
            //Get the possible transition to this state
            for(it=currentStatePtr->getFromBegin();it!=currentStatePtr->getFromEnd();it++){
                //Get pointer to previous state
                previousStatePtr=(*it);
                previousState=previousStatePtr->getIterator();
                
                //Only need to calculate transition if the viterbi > -INFINITY
                if (trell[sequencePosition-1][previousState].viti!=-INFINITY){
                    
                    //Calculate transition
                    double trans;
                    //If transition is already calculated
                    if (trell[sequencePosition][currentState].trans.count(previousStatePtr)){
                        trans = trell[sequencePosition][currentState].trans[previousStatePtr];
                    }
                    else{
                        trans = getTransition();
                        trell[sequencePosition][currentState].trans[previousStatePtr]=trans;
                    }
#ifdef VERBOSE
                    std::cout << "Position:\t" << sequencePosition <<std::endl;
                    std::cout << "State:\t" << currentState <<std::endl;
                    std::cout << "Previous State:\t"<<previousState << std::endl;
                    std::cout << "Transition:\t" << trans <<std::endl;
                    std::cout << "Emission:\t" << emm << std::endl;
                    std::cout <<std::endl;
#endif
                    //Calculate viterbi value for cell
                    double viterbiValue = emm + trans + trell[sequencePosition-1][previousState].viti;
                    
                    //If current viterbi value > -Infinity then it's possible to transition from this state to next
                    //Set the next possible states
                    if(viterbiValue>-INFINITY){
                        _nextStates->insert(currentStatePtr->getToBegin(),currentStatePtr->getToEnd());
                        setViterbi(viterbiValue);
                    }
                }
            }
        }
        return;
    }
    
    
    
    //!Calculate viterbi score for give cell from previous cell
    void simpleTrellis::calcNthViterbi(size_t n){
        
        //For all the possible states at this position
        for(currentStatesIterator = currentStates->begin(); currentStatesIterator!=currentStates->end();currentStatesIterator++){
            //Get pointer to current state
            currentStatePtr = (*currentStatesIterator);
            currentState = currentStatePtr->getIterator();
            
            //Calculate emission for the current cell
            double emm;
            //If emission is not yet calculated for the cell
            if (!trell[sequencePosition][currentState].emmCalculated){
                
                //Calculate emission
                emm=currentStatePtr->get_emission(*seq, sequencePosition);
                
                //Check for external definition
                if (externalDefinitionSet){
                    emm+=seq->getWeight(sequencePosition, currentState);
                }
                
                //Store emission value in cell
                trell[sequencePosition][currentState].emmCalculated=true;
                trell[sequencePosition][currentState].emm=emm;
            }
            else{  //emm was previously calculated
                
                //Get previously calculated emm value;
                emm = trell[sequencePosition][currentState].emm;
            }
            
            std::vector<scores>* temp_scores = new std::vector<scores>;  //Store all nth_viterbi scores for current cell
            
            std::vector<state*>::iterator it;            
            //Get the possible transition to this state
            for(it=currentStatePtr->getFromBegin();it!=currentStatePtr->getFromEnd();it++){
                //Get pointer to previous state
                previousStatePtr=(*it);
                previousState=previousStatePtr->getIterator();
                
                            
                //Only need to calculate transition if the viterbi > -INFINITY
                if (trell[sequencePosition-1][previousState].viti!=-INFINITY){
                    
                    //Calculate transition
                    double trans;
                    //If transition is already calculated
                    if (trell[sequencePosition][currentState].trans.count(previousStatePtr)){
                        trans = trell[sequencePosition][currentState].trans[previousStatePtr];
                    }
                    else{
                        trans = getTransition();
                        trell[sequencePosition][currentState].trans[previousStatePtr]=trans;
                    }
#ifdef VERBOSE
                    std::cout << "Position:\t" << sequencePosition <<std::endl;
                    std::cout << "State:\t" << currentState <<std::endl;
                    std::cout << "Previous State:\t"<<previousState << std::endl;
                    std::cout << "Transition:\t" << trans <<std::endl;
                    std::cout << "Emission:\t" << emm << std::endl;
                    std::cout <<std::endl;
#endif
                    //Calculate viterbi value for cell
                    double viterbiValue = emm + trans + trell[sequencePosition-1][previousState].viti;
                    
                    //If current viterbi value > -Infinity then it's possible to transition from this state to next
                    //Set the next possible states
                    if(viterbiValue>-INFINITY){
                        _nextStates->insert(currentStatePtr->getToBegin(),currentStatePtr->getToEnd());
                        setViterbi(viterbiValue);
                    }
                    
                    for(size_t nth_iterator=0;nth_iterator<trell[sequencePosition-1][previousState].nth_viterbi_scores->size();++nth_iterator){
                        scores tmp;
                        tmp.viterbi_score = emm + trans + (*trell[sequencePosition-1][previousState].nth_viterbi_scores)[nth_iterator].viterbi_score;
                        tmp.traceback_state = previousState;
                        tmp.traceback_state_score = nth_iterator;
                        temp_scores->push_back(tmp);
                    }
                    
                }
            }
            //Sort scores and assign to current cell
            
            sort_scores(*temp_scores);
            if (temp_scores->size()>n){
                temp_scores->resize(n);
            }
            
            trell[sequencePosition][currentState].nth_viterbi_scores = temp_scores;
        }
        return;
    }

    
    
    
    //!Calculate backward score for the given cell
    void simpleTrellis::calcBackward(){
        
        //For all the possible states at this position
        for(currentStatesIterator = currentStates->begin(); currentStatesIterator!=currentStates->end();currentStatesIterator++){
            //Get pointer to current state
            currentStatePtr = (*currentStatesIterator);
            currentState = currentStatePtr->getIterator();
            
            
            //Calculate emission for the current cell
            double emm;
            //If emission is not yet calculated for the cell
            if (!trell[sequencePosition][currentState].emmCalculated){
                
                //Calculate emission
                emm=currentStatePtr->get_emission(*seq, sequencePosition);
                
                //Check for external definition
                if (externalDefinitionSet){
                    emm+=seq->getWeight(sequencePosition, currentState);
                }
                
                //Store emission value in cell
                trell[sequencePosition][currentState].emmCalculated=true;
                trell[sequencePosition][currentState].emm=emm;
            }
            else{  //emm was previously calculated
                
                //Get previously calculated emm value;
                emm = trell[sequencePosition][currentState].emm;
            }
            
            std::vector<state*>::iterator it;
            
            //Get the possible transition to this state
            for(it=currentStatePtr->getFromBegin();it!=currentStatePtr->getFromEnd();it++){
                //Get pointer to previous state
                previousStatePtr=(*it);
                previousState=previousStatePtr->getIterator();
                
//                std::cout <<"Current " << currentStatePtr->getName() << " FROM ";
//                std::cout << previousStatePtr->getName() << std::endl;
                //std::cout << previousStatePtr->getIterator() << std::endl;
                
                //Only need to calculate transition if the viterbi > -INFINITY
                if (trell[sequencePosition][currentState].back!=-INFINITY){
                    
                    //Calculate transition
                    double trans;
                    
                    //If transition is already calculated
                    if (trell[sequencePosition-1][currentState].trans.count(previousStatePtr)){
                        trans = trell[sequencePosition-1][currentState].trans[previousStatePtr];
                    }
                    else{
                        trans = getTransition();
                        trell[sequencePosition-1][currentState].trans[previousStatePtr]=trans;
                    }
//                    std::cout << "Position:\t" << sequencePosition <<std::endl;
//                    std::cout << "State:\t" << currentState <<std::endl;
//                    std::cout << "Previous State:\t"<<previousState << std::endl;
//                    std::cout << "Transition:\t" << exp(trans) <<std::endl;
//                    std::cout << "Emission:\t" << exp(emm) << std::endl;
//                    std::cout <<std::endl;

                    
                    //Calculate viterbi value for cell
                    double backward_value = emm + trans + trell[sequencePosition][currentState].back;
                    
                    //If current viterbi value > -Infinity then it's possible to transition from this state to next
                    //Set the next possible states
                    if(backward_value>-INFINITY){
                        _nextStates->insert(previousStatePtr->getToBegin(),previousStatePtr->getToEnd());
                        setBackward(sequencePosition-1,previousState,backward_value);
                        //std::cout << exp(backward_value) << std::endl;
                    }
                }
            }
        }
        
        ///Old backward functions
//        double prelim = getTransition() + previousStatePtr->get_emission(*seq,sequencePosition);
//        
//        if (externalDefinitionSet){
//            prelim+=seq->getWeight(sequencePosition, previousStatePtr->getIterator());
//        }
//        
//        //setBackward(  prelim + getBackward( sequencePosition , previousStatePtr->getIterator() ) );
        
        return;
    }
    

    //!Calculate the ending cell Viterbi scores
    void simpleTrellis::calcEndViterbi(){
        size_t stateSize=hmm->state_size();
        size_t seqSize=seq->getLength();
        
        for(int previousSt=0;previousSt<stateSize;previousSt++){
            
            double finalViterbiScore=getEndingTransition(previousSt)+trell[seqSize-1][previousSt].viti;
            
            if (finalViterbiScore>-INFINITY){
                if (finalViterbiScore>ending.viti){
                    ending.viti=finalViterbiScore;
                    ending.ptr=previousSt;
                }
            }
        }
        return;
    }
    
    //!Calculate the ending cell Viterbi scores
    void simpleTrellis::calcNthEndViterbi(size_t n){
        size_t stateSize=hmm->state_size();
        size_t seqSize=seq->getLength();
        
        std::vector<scores>* temp_nth = new std::vector<scores>;
        
        //Calculate the ending viterbi scores
        for(size_t previousSt=0;previousSt<stateSize;previousSt++){
            
            double trans = getEndingTransition(previousSt);
            double finalViterbiScore=trans+trell[seqSize-1][previousSt].viti;
            
            if (finalViterbiScore>-INFINITY){
                if (finalViterbiScore>ending.viti){
                    ending.viti=finalViterbiScore;
                    ending.ptr=previousSt;
                }
            }
            
            //Calculate the nth viterbi scores
            for(size_t i = 0; i<trell[seqSize-1][previousState].nth_viterbi_scores->size();++i){
                
                scores temp_score;
                temp_score.viterbi_score = trans + (*trell[seqSize-1][previousState].nth_viterbi_scores)[i].viterbi_score;
                temp_score.traceback_state = previousState;
                temp_score.traceback_state_score = i;
                
                temp_nth->push_back(temp_score);
            }
        }
        
        //Sort the top scores
        sort_scores(*temp_nth);
        
        //Resize the list to nth top
        if (temp_nth->size()>n){
            temp_nth->resize(n);
        }
        
        //Assign the Nth to ending cell
        ending.nth_viterbi_scores = temp_nth;
        return;
    }

    //!Calculate the ending cell Forward score
    void simpleTrellis::calcEndForward(){
        size_t stateSize=hmm->state_size();
        size_t seqSize=seq->getLength();
        
        for(int previousSt=0;previousSt<stateSize;previousSt++){
            
            double finalForwardScore=getEndingTransition(previousSt)+trell[seqSize-1][previousSt].forw;
            
            if (finalForwardScore>-INFINITY){
                ending.forw=addLog(ending.forw, finalForwardScore);
            }
        }
        
        probabilityOfSequence=ending.forw;
        
        return;
        
    }

    //!Calculate the ending cell Viterbi and Forward scores
    void simpleTrellis::calcEndForwardViterbi(){
        size_t stateSize=hmm->state_size();
        size_t seqSize=seq->getLength();
        
        for(int previousSt=0;previousSt<stateSize;previousSt++){
            double transitionValue=getEndingTransition(previousSt);
            
            double finalViterbiScore=transitionValue+trell[seqSize-1][previousSt].viti;
            double finalForwardScore=transitionValue+trell[seqSize-1][previousSt].forw;
            
            if (finalViterbiScore>-INFINITY){
                if (finalViterbiScore>ending.viti){
                    ending.viti=finalViterbiScore;
                    ending.ptr=previousSt;
                }
            }
            
            if (finalForwardScore>-INFINITY){
                ending.forw=addLog(ending.forw, finalForwardScore);
            }
        }
        
        probabilityOfSequence=ending.forw;
        
        return;
    }

    //!Calculate the posterior probability for the trellis cells
    void simpleTrellis::calcPosterior(){
        size_t seq_length = seq->getLength();
        size_t hmm_size = hmm->state_size();
        
        std::vector<double> temp (hmm_size,0);
        posterior = new std::vector< std::vector< double > > (seq_length,temp);
        posterior_pointer = new std::vector<int> (seq_length,-1);
        
        
        for(int i=0;i<seq->getLength();i++){
            double max(0.0);
            int max_ptr(-1);
            
            for (int j=0;j<hmm->state_size();j++){
                (*posterior)[i][j]= exp((trell[i][j].forw + trell[i][j].back) - probabilityOfSequence);
                
                if ((*posterior)[i][j]>max){
                    max = (*posterior)[i][j];
                    max_ptr = j;
                }
            }
            (*posterior_pointer)[i] = max_ptr;
        }
        return;
    }


    //!Get Traceback Pointer for the current cell
    int simpleTrellis::getPtr(){
        return trell[sequencePosition][currentState].ptr;
    }
    
    //!Get the Viterbi score for the current cell
    double& simpleTrellis::getViterbi(){
        return trell[sequencePosition][currentState].viti;
    }
    
    //!Get the Forward score for the current cell
    double& simpleTrellis::getForward(){
        return trell[sequencePosition][currentState].forw;
    }
    
    //!Get the Backward score for the current cell
    double& simpleTrellis::getBackward(){
        return trell[sequencePosition][currentState].back;
    }

    //!Get traceback pointer for given state and position
    //!\param seqPos Position within the sequence
    //!\param currentState State of index
    //!\return int Pointer to state in previous position (n-1)
    int simpleTrellis::getPtr(size_t seqPos, size_t currState){
        return trell[seqPos][currState].ptr;
    }
    
    //!Get Viterbi score for given state and position
    //!\param seqPos Position within the sequence
    //!\param currentState State of index
    //!\return double Viterbi score
    double& simpleTrellis::getViterbi(size_t seqPos, size_t currState){
        return trell[seqPos][currState].viti;
    }
    
    //!Get Forward score for given state and position
    //!\param seqPos Position within the sequence
    //!\param currentState State of index
    //!\return double Forward score
    double& simpleTrellis::getForward(size_t seqPos, size_t currState){
        return trell[seqPos][currState].forw;
    }
    
    //!Get Backward score for given state and position
    //!\param seqPos Position within the sequence
    //!\param currentState State of index
    //!\return double Backward score
    double& simpleTrellis::getBackward(size_t seqPos, size_t currState){
        return trell[seqPos][currState].back;
    }

    //!Set the viterbi score for the given cell
    //!\param val Viterbi score
    void simpleTrellis::setViterbi(double& value){
#ifdef VERBOSE
        std::cout << "Temp Viterbi Value:\t" << value <<std::endl;
        std::cout << "Previous Viterbi Value:\t" << trell[sequencePosition][currentState].viti << std::endl;
        
#endif 
        
        if (value > trell[sequencePosition][currentState].viti){
            trell[sequencePosition][currentState].viti=value;
            
            if (previousStatePtr==NULL){
                trell[sequencePosition][currentState].ptr = -1;
            }
            else{
                trell[sequencePosition][currentState].ptr = previousState;
            }
#ifdef VERBOSE 
            std::cout << "Pointer Set:\t" << previousState <<std::endl;
#endif
        }
        
        return;
    }
    
    
    //!Set the viterbi score for the given cell
    //!\param val Viterbi score
    void simpleTrellis::setNthViterbi(double& value){
#ifdef VERBOSE
        std::cout << "Temp Viterbi Value:\t" << value <<std::endl;
        std::cout << "Previous Viterbi Value:\t" << trell[sequencePosition][currentState].viti << std::endl;
        
#endif
        scores temp;
        temp.viterbi_score = value;
        temp.traceback_state = -1;
        
        trell[sequencePosition][currentState].nth_viterbi_scores = new std::vector<scores>(1,temp);
                
        if (value > trell[sequencePosition][currentState].viti){
            trell[sequencePosition][currentState].viti=value;
            
            if (previousStatePtr==NULL){
                trell[sequencePosition][currentState].ptr = -1;
            }
            else{
                trell[sequencePosition][currentState].ptr = previousState;
            }
#ifdef VERBOSE
            std::cout << "Pointer Set:\t" << previousState <<std::endl;
#endif
        }
        return;
    }
    
    
    
    //!Set the viterbi score for the given cell
    //!\param seqPos Sequence position
    //!\param currState State index
    //!\param value Viterbi score
    void simpleTrellis::setViterbi(size_t seqPos, size_t currState, double& value){
        trell[seqPos][currState].viti=value;
        
        return;
    }
    
    
    
    //!Set the Forward score for the given cell
    //!\param val Forward score
    void simpleTrellis::setForward(double& value){
        
        if (getForward()==-INFINITY){
            trell[sequencePosition][currentState].forw=value;
        }
        else{
            trell[sequencePosition][currentState].forw=addLog(trell[sequencePosition][currentState].forw, value);
        }
        
        return;
    }
    
    //!Set the Forward score for the given cell
    //!\param seqPos Sequence position
    //!\param currState State index
    //!\param value Forward score
    void simpleTrellis::setForward(size_t seqPos, size_t currState,double& value){
        if (getForward()==-INFINITY){
            trell[seqPos][currState].forw=value;
        }
        else{
            trell[seqPos][currState].forw=addLog(trell[seqPos][currState].forw, value);
        }
        return;
    }
    
    //!Set the Backward score for the given cell
    //!\param val Backward score
    void simpleTrellis::setBackward(double& value){
        
        //sequencePosition-1, currentState
        
        if (getBackward()==-INFINITY){
            trell[sequencePosition-1][currentState].back=value;
        }
        else{
            trell[sequencePosition-1][currentState].back=addLog(trell[sequencePosition-1][currentState].back, value);
        }
        
        if (value>-INFINITY){
            _nextStates->insert(currentStatePtr);
        }
        
        return;
    }
    
    
    //!Set the Backward score for the given cell
    //!\param seqPos Sequence position
    //!\param currState State index
    //!\param value Forward score
    void simpleTrellis::setBackward(size_t seqPos, size_t currState,double& value){
        if (getBackward(seqPos,currState)==-INFINITY){
            trell[seqPos][currState].back=value;
        }
        else{
            trell[seqPos][currState].back=addLog(trell[seqPos][currState].back, value);
        }
        return;
    }

    //!Set the traceback pointer for the given cell
    //!\param pointer Traceback pointer of given cell
    void simpleTrellis::setPtr(int pointer){
        trell[sequencePosition][currentState].ptr=pointer;
        return;
    }
    
    //!Set the Traceback pointer for the given cell
    //!\param seqPos Sequence position
    //!\param currState State index
    //!\param pointer Traceback Pointer
    void simpleTrellis::setPtr(size_t seqPos, size_t currState,int pointer){
        trell[seqPos][currState].ptr=pointer;
        return;
    }

    
    //!Perform traceback through trellis 
    //!\return [out] path trackback_path
    void simpleTrellis::traceback(traceback_path& path){
        
        if (ending.viti!=-INFINITY){
            path.setScore(ending.viti);
            int pointer=ending.ptr;
            path.push_back(ending.ptr);
            
            size_t seq_length = seq->getLength();
            if (seq_length==0){
                return;
            }
            
            for(size_t position=seq_length -1 ; position>0;position--){
                pointer=trell[position][pointer].ptr;
                path.push_back(pointer);
                
                if (pointer==-2){
                    path.clear();
                    return;
                }
            }
        }
        
        return;
    }
    
    //!Perform traceback through trellis
    //!\return [out] path trackback_path
    void simpleTrellis::traceback(traceback_path& path,size_t n){
        if (n<=ending.nth_viterbi_scores->size()){
            if (ending.viti!=-INFINITY){
                path.setScore((*ending.nth_viterbi_scores)[n].viterbi_score);
                
                int state_pointer=(*ending.nth_viterbi_scores)[n].traceback_state;
                size_t score_pointer = (*ending.nth_viterbi_scores)[n].traceback_state_score;
                
                path.push_back(state_pointer);
                
                size_t seq_length = seq->getLength();
                if (seq_length==0){
                    return;
                }
                
                for(size_t position=seq_length -1 ; position>0;position--){
                    state_pointer = (*trell[position][state_pointer].nth_viterbi_scores)[score_pointer].traceback_state;
                    score_pointer = (*trell[position][state_pointer].nth_viterbi_scores)[score_pointer].traceback_state_score;
                    
                    path.push_back(state_pointer);
                    
                    if (state_pointer==-2){
                        path.clear();
                        return;
                    }
                }
            }
        }
        return;
    }
    
    //!Perform stochastic Traceback on trellis (Not supported on simpleTrellis) 
    void simpleTrellis::traceStochViterbi(traceback_path& path){
        std::cerr << "Can't perform stochastic traceback on SimpleTrellis.  You must use a Stochastic Trellis\n";
        exit(1);
    }

    //!Perform stochastic Traceback on trellis (Not supported on simpleTrellis)
    void simpleTrellis::traceStochForward(traceback_path& path){
        std::cerr << "Can't perform stochastic traceback on SimpleTrellis.  You must use a Stochastic Trellis\n";
        exit(1);
    }

}
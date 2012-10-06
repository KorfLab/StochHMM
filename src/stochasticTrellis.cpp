//
//  stochasticTrellis.cpp

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

#include "stochasticTrellis.h"
namespace StochHMM{

    //!Create stochTrellis
    //!\param modl  Model file to use for decoding
    //!\param seqs  Sequences to use for decoding
    stochTrellis::stochTrellis(model* modl,sequences* seqs): basicTrellis(modl,seqs){
        
        size_t stateSize=hmm->state_size();
        size_t seqSize=seq->getLength();
        
        stochCell initial(stateSize);
        trell.assign(seqSize,std::vector<stochCell>(stateSize,initial));
        
        ending.posteriorProbability.assign(stateSize,-INFINITY);
        ending.viterbiProbability.assign(stateSize,-INFINITY);
        ending.log_trans.assign(stateSize,-INFINITY);
        
        return;
    }
    
    //TODO: Create a function to save and import trellis file
    
    
    //!Print the trellis to stdout
    void stochTrellis::print(){
        for(int i=0;i<seq->getLength();i++){
            for(int j=0;j<hmm->state_size();j++){
                std::cout << "Position: " << i << "\t" << "State: " << j <<std::endl;
                trell[i][j].print();
                std::cout << std::endl;
            }
        }
        return;
    }
    
    
    //!Calculate the viterbi and forward scores for a given cell                                                            
    void stochTrellis::calcForwardViterbi(){
        
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
    void stochTrellis::calcForward(){
        
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
    void stochTrellis::calcViterbi(){
        
        
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
                    double viterbiValue = emm + trans + trell[sequencePosition-1][previousState].viti;
                    
                    if(viterbiValue>-INFINITY){
                        DefinedStates->insert(currentStatePtr);
                        _nextStates->insert(currentStatePtr->getToBegin(),currentStatePtr->getToEnd());
                        
                        setViterbi(viterbiValue);
                    }
                }
            }
        }
#ifdef VERBOSE
        std::cout <<std::endl;
#endif
        return;
    }
    
    
    void stochTrellis::calcNthViterbi(size_t n){
        return;
    }
    
    void stochTrellis::setNthViterbi(double&){
        return;
    }
    
    void stochTrellis::calcNthEndViterbi(size_t n){
        return;
    }
    
    void stochTrellis::traceback(traceback_path& path, size_t n){
        return;
    }
    
    
    
    
    //!Calculate backward score for the given cell                                                         
    void stochTrellis::calcBackward(){
        
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
        return;
    }
    


    //!Calculate the ending cell Viterbi scores
    void stochTrellis::calcEndViterbi(){
        size_t stateSize=hmm->state_size();
        size_t seqSize=seq->getLength();
        
        for(int previousSt=0;previousSt<stateSize;previousSt++){
            
            double finalViterbiScore=getEndingTransition(previousSt)+trell[seqSize-1][previousSt].viti;
            
            if (finalViterbiScore>-INFINITY){
                ending.viterbiSum=addLog(ending.viterbiSum, finalViterbiScore);
                ending.viterbiProbability[previousSt]=finalViterbiScore;
                if (finalViterbiScore>ending.viti){
                    ending.viti=finalViterbiScore;
                    ending.ptr=previousSt;
                }
            }
        }
        
        ending.calcViterbiProb();
        
        for(int i=0;i<seqSize;i++){
            for(int j=0;j<stateSize;j++){
                trell[i][j].calcViterbiProb();
            }
        }
        return;
    }

    //!Calculate the ending cell Forward score
    void stochTrellis::calcEndForward(){
        size_t stateSize=hmm->state_size();
        size_t seqSize=seq->getLength();
        
        for(int previousSt=0;previousSt<stateSize;previousSt++){
            
            double finalForwardScore=getEndingTransition(previousSt)+trell[seqSize-1][previousSt].forw;
            
            if (finalForwardScore>-INFINITY){
                ending.forw=addLog(ending.forw, finalForwardScore);
                ending.posteriorProbability[previousSt]=finalForwardScore;
            }
        }
        
        probabilityOfSequence=ending.forw;
        ending.calcForwardProb();
        
        for(int i=0;i<seqSize;i++){
            for(int j=0;j<stateSize;j++){
                trell[i][j].calcForwardProb();
            }
        }
        
        return;

    }

    //!Calculate the ending cell Viterbi and Forward scores
    void stochTrellis::calcEndForwardViterbi(){
        size_t stateSize=hmm->state_size();
        size_t seqSize=seq->getLength();
        
        for(int previousSt=0;previousSt<stateSize;previousSt++){
            double transitionValue=getEndingTransition(previousSt);
            
            double finalViterbiScore=transitionValue+trell[seqSize-1][previousSt].viti;
            double finalForwardScore=transitionValue+trell[seqSize-1][previousSt].forw;
            
            if (finalViterbiScore>-INFINITY){
                ending.viterbiSum=addLog(ending.viterbiSum, finalViterbiScore);
                ending.viterbiProbability[previousSt]=finalViterbiScore;
                if (finalViterbiScore>ending.viti){
                    ending.viti=finalViterbiScore;
                    ending.ptr=previousSt;
                }
            }
            
            if (finalForwardScore>-INFINITY){
                ending.forw=addLog(ending.forw, finalForwardScore);
                ending.posteriorProbability[previousSt]=finalForwardScore;
            }
        }
        
        ending.calcForwardViterbiProb();
        probabilityOfSequence=ending.forw;
        
        for(int i=0;i<seqSize;i++){
            for(int j=0;j<stateSize;j++){
                trell[i][j].calcForwardViterbiProb();
            }
        }
        
        
        return;
    }

    //!Get Traceback Pointer for the current cell
    //!\return int index of state
    int stochTrellis::getPtr(){
        return trell[sequencePosition][currentState].ptr;
    }
    
    //!Get the Viterbi score for the current cell
    double& stochTrellis::getViterbi(){
        return trell[sequencePosition][currentState].viti;
    }
    
    //!Get the Forward score for the current cell
    double& stochTrellis::getForward(){
        return trell[sequencePosition][currentState].forw;
    }
    
    //!Get the Backward score for the current cell
    double& stochTrellis::getBackward(){
        return trell[sequencePosition][currentState].back;
    }
    
    //!Get traceback pointer for given state and position
    //!\param seqPos Position within the sequence
    //!\param currentState State of index
    //!\return int Pointer to state in previous position (n-1)
    int stochTrellis::getPtr(size_t seqPos, size_t currState){
        return trell[seqPos][currState].ptr;
    }
    
    //!Get Viterbi score for given state and position
    //!\param seqPos Position within the sequence
    //!\param currentState State of index
    //!\return double Viterbi score
    double& stochTrellis::getViterbi(size_t seqPos, size_t currState){
        return trell[seqPos][currState].viti;
    }
    
    //!Get Forward score for given state and position
    //!\param seqPos Position within the sequence
    //!\param currentState State of index
    //!\return double Forward score
    double& stochTrellis::getForward(size_t seqPos, size_t currState){
        return trell[seqPos][currState].forw;
    }

    //!Get Backward score for given state and position
    //!\param seqPos Position within the sequence
    //!\param currentState State of index
    //!\return double Backward score
    double& stochTrellis::getBackward(size_t seqPos, size_t currState){
        return trell[seqPos][currState].back;
    }

    
    //!Set the viterbi score for the given cell
    //!\param val Viterbi score
    void stochTrellis::setViterbi(double& value){
        
        if (value>trell[sequencePosition][currentState].viti){  //Take Maximum and assign pointer to it
            trell[sequencePosition][currentState].viti=value;
            setPtr(previousState);
        }
        
        if (previousState>=0){
            //Store the viterbi score from previous
            trell[sequencePosition][currentState].viterbiProbability[previousState]=value;
            
            //Store the new viterbi sum
            trell[sequencePosition][currentState].viterbiSum=addLog(trell[sequencePosition][currentState].viterbiSum, value);
        }
            
        if (value>-INFINITY){
//            nextStates->insert(currentStatePtr);
        }
        
        return;
    }


    //!Set the Forward score for the given cell
    //!\param val Forward score
    void stochTrellis::setForward(double& value){
        //Store forward contributions for each state
        if (previousState>=0){
            trell[sequencePosition][currentState].posteriorProbability[previousState]=value;
        }
        
        
        if (getForward()==-INFINITY){
            trell[sequencePosition][currentState].forw=value;
        }
        else{
            trell[sequencePosition][currentState].forw=addLog(trell[sequencePosition][currentState].forw, value);
        }
        
        if (value>-INFINITY){
//            nextStates->insert(currentStatePtr);
        }
        
        return;
    }

    //TODO:  Compare to the backward function and make standard 
    //!Set the Backward score for the given cell
    //!\param val Backward score
    void stochTrellis::setBackward(double& value){
        
        //sequencePosition-1, currentState
        
        if (getBackward()==-INFINITY){
            trell[sequencePosition-1][currentState].back=value;
        }
        else{
            trell[sequencePosition-1][currentState].back=addLog(trell[sequencePosition-1][currentState].back, value);
        }
        
        if (value>-INFINITY){
//            nextStates->insert(currentStatePtr);
        }
        
        return;
    }

    //!Set the traceback pointer for the given cell
    //!\param pointer Traceback pointer of given cell
    void stochTrellis::setPtr(int pointer){
        trell[sequencePosition][currentState].ptr=pointer;
        return;
    }

    //!Set the viterbi score for the given cell
    //!\param seqPos Sequence position
    //!\param currState State index
    //!\param value Viterbi score
    void stochTrellis::setViterbi(size_t seqPos, size_t currState, double& value){
        trell[seqPos][currState].viti=value;
        return;
    }

    //!Set the Forward score for the given cell
    //!\param seqPos Sequence position
    //!\param currState State index
    //!\param value Forward score
    void stochTrellis::setForward(size_t seqPos, size_t currState,double& value){
        if (getForward()==-INFINITY){
            trell[seqPos][currState].forw=value;
        }
        else{
            trell[seqPos][currState].forw=addLog(trell[seqPos][currState].forw, value);
        }
        return;
    }


    //!Set the Backward score for the given cell
    //!\param seqPos Sequence position
    //!\param currState State index
    //!\param value Forward score
    void stochTrellis::setBackward(size_t seqPos, size_t currState,double& value){
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
    void stochTrellis::setPtr(size_t seqPos, size_t currState,int pointer){
        trell[seqPos][currState].ptr=pointer;
        return;
    }

    //!Perform traceback through trellis 
    //!\return [out] path trackback_path
    void stochTrellis::traceback(traceback_path& path){
        if (ending.viti!=-INFINITY){
            int pointer=ending.ptr;
            path.push_back(ending.ptr);
            path.setScore(ending.viti);
            
            size_t seq_length = seq->getLength();
            if (seq_length==0){
                return;
            }
            
            for(size_t position=seq_length - 1 ; position>0;position--){
                pointer=trell[position][pointer].ptr;
                path.push_back(pointer);
                
                if (pointer==-2){  //If the pointer is the pre-initialized value (-2), then the traceback is invalid
                    path.clear();
                    return;
                }
            }
        }
        return;
        
    }

    //!Calculate the posterior probability for the trellis cells
    void stochTrellis::calcPosterior(){
        for(int i=0;i<seq->getLength();i++){
            for (int j=0;j<hmm->state_size();j++){
                std::cout << i << "\t" << j << "\t"  << exp((trell[i][j].forw + trell[i][j].back) - probabilityOfSequence) <<std::endl;
            }
        }
        return;
    }

    //!Perform stochastic traceback using viterbi scores
    //!\param path [out] traceback_path to store traceback in
    void stochTrellis::traceStochViterbi(traceback_path& path){
        
        if (ending.viti==-INFINITY){
            return;
        }
        
        size_t stateSize=hmm->state_size();
        size_t seqSize=seq->getLength();
        
        
        //Select ending cell
        double random=((double)rand()/((double)(RAND_MAX)+(double)(1)));
        double vprob_sum=0;
        int y=-1;
        
        //select last cell;
        for(int i=0;i<stateSize;i++){
            vprob_sum+=ending.viterbiProbability[i];		
            if (random<vprob_sum){ 
                y=i;
                break;
            }
        }
        
        if (y==-1){
            std::cerr << "Improper viterbi stochastic ending cell selection in trace_stoch_viterbi\n";
            exit(1);
        }
        
        int previous=-1;
        
        for(int seqPosition=seqSize-1;seqPosition>=1;seqPosition--){
            path.push_back(y);
            //previous=y;
            random=((double)rand()/((double)(RAND_MAX)+(double)(1)));
            vprob_sum=0;	
            
            int selection=-1;		
            
            for(int k=0;k<trell[seqPosition][y].viterbiProbability.size();k++){
                vprob_sum+=trell[seqPosition][y].viterbiProbability[k];
                
                if(random<=vprob_sum){
                    selection=k;
                    break;
                }
            }
            
            if (selection==-1){  //Invalid traceback
                path.clear();
                return;
            }
            
            y=selection;
        }
        path.push_back(y);
        
        return;
    }


    //!Perform stochastic traceback using the Forward scores
    //!\param path [out] traceback_path to store traceback in
    void stochTrellis::traceStochForward(traceback_path& path){
        
        if (ending.forw==-INFINITY){
            return;
        }
            
        size_t stateSize=hmm->state_size();
        size_t seqSize=seq->getLength();
        
        
        //Select ending cell
        double random=((double)rand()/((double)(RAND_MAX)+(double)(1)));
        double fprob_sum=0;
        int y=-1;
        
        //select last cell;
        for(int i=0;i<stateSize;i++){
            fprob_sum+=ending.posteriorProbability[i];		
            if (random<fprob_sum){ 
                y=i;
                break;
            }
        }
        
        if (y==-1){
            std::cerr << "Improper viterbi stochastic ending cell selection in trace_stoch_viterbi\n";
            exit(1);
        }
        
        
        for(int seqPosition=seqSize-1;seqPosition>=1;seqPosition--){
            path.push_back(y);
            //previous=y;
            random=((double)rand()/((double)(RAND_MAX)+(double)(1)));
            fprob_sum=0;	
            
            int selection=-1;		
            
            for(int k=0;k<trell[seqPosition][y].posteriorProbability.size();k++){
                fprob_sum+=trell[seqPosition][y].posteriorProbability[k];
                
                if(random<=fprob_sum){
                    selection=k;
                    break;
                }
            }
            
            if (selection==-1){  //Invalid traceback
                path.clear();
                return;
            }
            
            y=selection;
        }
        path.push_back(y);
        
        return;
    }

//    //!Perform traceback using the Posterior scores
//    void stochTrellis::tracePosterior(traceback_path& path){
//        if (ending.forw==-INFINITY){
//            return;
//        }
//        
//        size_t stateSize=hmm->state_size();
//        size_t seqSize=seq->getLength();
//        
//        
//        //Select ending cell
//        double fprob_max=0;
//        int y=-1;
//        
//        //select last cell;
//        for(int i=0;i<stateSize;i++){
//            if (fprob_max<ending.posteriorProbability[i]){ 
//                y=i;
//            }
//        }
//        
//        if (y==-1){
//            std::cerr << "Improper viterbi stochastic ending cell selection in trace_stoch_viterbi\n";
//            exit(1);
//        }
//        
//        
//        for(int seqPosition=seqSize-1;seqPosition>=1;seqPosition--){
//            path.push_back(y);
//            //previous=y;
//            fprob_max=0;	
//            
//            int selection=-1;		
//            
//            for(int k=0;k<trell[seqPosition][y].posteriorProbability.size();k++){
//                
//                if(fprob_max<trell[seqPosition][y].posteriorProbability[k]){
//                    selection=k;
//                }
//            }
//            
//            if (selection==-1){  //Invalid traceback
//                path.clear();
//                return;
//            }
//            
//            y=selection;
//        }
//        path.push_back(y);
//        
//        return;
//    }

}
//
//  nthTrellis.cpp

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

#include "nthTrellis.h"
namespace StochHMM{


    //!Create nthTrellis
    //!\param modl  Model file to use for decoding
    //!\param seqs  Sequences to use for decoding
    nthTrellis::nthTrellis(HMM* modl, sequences* seqs): basicTrellis(modl,seqs){
        size_t stateSize=model->state_size();
        size_t seqSize=seq->getLength();
        
        nthCell initial(stateSize);
        trell.assign(seqSize,vector<nthCell>(stateSize,initial));
        
        return;
    }
    
    
    nthTrellis::void calcViterbi(){
        
    }
    
    
    

    //!Print the trellis to stdout
    void nthTrellis::print(){
        for(int i=0;i<seq->getLength();i++){
            for(int j=0;j<model->state_size();j++){
                cout << "Position: " << i << "\t" << "State: " << j <<endl;
                trell[i][j].print();
                cout << endl;
            }
        }
        return;
    }

    //!Calculate the ending cell Viterbi scores
    void nthTrellis::calcEndViterbi(){
        int stateSize=model->state_size();
        int seqSize=seq->getLength();
        
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

    //!Calculate the ending cell Forward score
    void nthTrellis::calcEndForward(){
        int stateSize=model->state_size();
        int seqSize=seq->getLength();
        
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
    void nthTrellis::calcEndForwardViterbi(){
        int stateSize=model->state_size();
        int seqSize=seq->getLength();
        
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
    void nthTrellis::calcPosterior(){
        for(int i=0;i<seq->getLength();i++){
            for (int j=0;j<model->state_size();j++){
                cout << i << "\t" << j << "\t"  << exp((trell[i][j].forw + trell[i][j].back) - probabilityOfSequence) <<endl;
            }
        }
    }



    //!Get Traceback Pointer for the current cell
    //!\return int index of state
    int nthTrellis::getPtr(){
        return trell[sequencePosition][currentState].ptr;
    }

    //!Get the Viterbi score for the current cell
    double& nthTrellis::getViterbi(){
        return trell[sequencePosition][currentState].viti;
    }
    
    //!Get the Forward score for the current cell
    double& nthTrellis::getForward(){
        return trell[sequencePosition][currentState].forw;
    }

    //!Get the Backward score for the current cell
    double& nthTrellis::getBackward(){
        return trell[sequencePosition][currentState].back;
    }


    //!Get traceback pointer for given state and position
    //!\param seqPos Position within the sequence
    //!\param currentState State of index
    //!\return int Pointer to state in previous position (n-1)
    int nthTrellis::getPtr(int seqPos, int currState){
        return trell[seqPos][currState].ptr;
    }

    //!Get Viterbi score for given state and position
    //!\param seqPos Position within the sequence
    //!\param currentState State of index
    //!\return double Viterbi score
    double& nthTrellis::getViterbi(int seqPos, int currState){
        return trell[seqPos][currState].viti;
    }

    //!Get Forward score for given state and position
    //!\param seqPos Position within the sequence
    //!\param currentState State of index
    //!\return double Forward score
    double& nthTrellis::getForward(int seqPos, int currState){
        return trell[seqPos][currState].forw;
    }

    //!Get Backward score for given state and position
    //!\param seqPos Position within the sequence
    //!\param currentState State of index
    //!\return double Backward score
    double& nthTrellis::getBackward(size_t seqPos, size_t currState){
        return trell[seqPos][currState].back;
    }


    //!Set the viterbi score for the given cell
    //!\param val Viterbi score
    void nthTrellis::setViterbi(double value){
        
        if (value>trell[sequencePosition][currentState].viti){
            trell[sequencePosition][currentState].viti=value;
            setPtr(previousState);
        }
        
        if (value>-INFINITY){
            trell[sequencePosition][currentState].viterbiScores[previousState].first=value;
            trell[sequencePosition][currentState].viterbiScores[previousState].second=previousState;
            next->insert(currentState);
        }
        
        
        return;
    }

    //!Set the Forward score for the given cell
    //!\param val Forward score
    void nthTrellis::setForward(double value){
        
        if (getForward()==-INFINITY){
            trell[sequencePosition][currentState].forw=value;
        }
        else{
            trell[sequencePosition][currentState].forw=addLog(trell[sequencePosition][currentState].forw, value);
        }
        
        if (value>-INFINITY){
            next->insert(currentState);
        }
        
        return;
    }

    //!Set the Backward score for the given cell
    //!\param val Backward score
    void nthTrellis::setBackward(double value){
        
        //sequencePosition-1, currentState
        
        if (getBackward()==-INFINITY){
            trell[sequencePosition-1][currentState].back=value;
        }
        else{
            trell[sequencePosition-1][currentState].back=addLog(trell[sequencePosition-1][currentState].back, value);
        }
        
        if (value>-INFINITY){
            next->insert(currentState);
        }
        
        return;
    }

    //!Set the traceback pointer for the given cell
    //!\param pointer Traceback pointer of given cell
    void nthTrellis::setPtr(int pointer){
        trell[sequencePosition][currentState].ptr=pointer;
        return;
    }


    //!Set the viterbi score for the given cell
    //!\param seqPos Sequence position
    //!\param currState State index
    //!\param value Viterbi score
    void nthTrellis::setViterbi(int seqPos, int currState, double value){
        trell[seqPos][currState].viti=value;
        
        return;
    }

    //!Set the Forward score for the given cell
    //!\param seqPos Sequence position
    //!\param currState State index
    //!\param value Forward score
    void nthTrellis::setForward(int seqPos, int currState,double value){
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
    void nthTrellis::setBackward(int seqPos, int currState,double value){
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
    void nthTrellis::setPtr(int seqPos, int currState,int pointer){
        trell[seqPos][currState].ptr=pointer;
        return;
    }


    //!Perform traceback through trellis 
    //!\return [out] path trackback_path
    traceback_path* nthTrellis::traceback(){
        if (ending.viti!=-INFINITY){
            traceback_path* path= new(std::nothrow) traceback_path;
            if (path==NULL){
                std::cerr << "OUT OF MEMORY\nFile" << __FILE__ << "Line:\t"<< __LINE__ << std::endl;
                exit(1);
            }
            
            int pointer=ending.ptr;
            path->push_back(ending.ptr);
    #ifdef DEBUG_VERBOSE
            cout << ending.ptr << " ";
    #endif 
            for(int position=seq->getLength() -1 ; position>0;position--){
                pointer=trell[position][pointer].ptr;
                path->push_back(pointer);
    #ifdef DEBUG_VERBOSE
                cout << pointer << " ";
    #endif
            }
            cout << endl;
            return path;
        }
        else{
            return NULL;
        }
        
    }


}
//
//  trellisCells.cpp


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

#include "trellisCells.h"
namespace StochHMM{
    scores::scores():viterbi_score(-INFINITY),traceback_state(-2),traceback_state_score(0){
        return;
    }
    
    //!Create a simple Trellis Cell
    simpleCell::simpleCell():viti(-INFINITY), forw(-INFINITY),back(-INFINITY),ptr(-1),emm(-INFINITY),emmCalculated(false),statePtr(NULL), nth_viterbi_scores(NULL){
        return;
    }

    //! Print the trellis cell to stdout
    void simpleCell::print(){
        std::cout << "Viterbi: " << exp(viti) << "\n" <<"Forward: " << exp(forw) << "\n" <<"Backward: " << exp(back) << "\n" << "Vitebi Pointer: " << ptr <<std::endl;
        return;
    }
    
    simpleCell::~simpleCell(){
        if (nth_viterbi_scores!=NULL){
            delete nth_viterbi_scores;
            nth_viterbi_scores=NULL;
        }
        return;
    }
    
//    void simpleCell::operator=(simpleCell& rhs){
//        viti=rhs.viti;
//        forw=rhs.forw;
//        emmCalculated=rhs.emmCalculated;
//        emm=rhs.emm;
//        //trans = rhs.trans;
//        statePtr = rhs.statePtr;
//        //position = rhs.position;
//        //traceback_ptr = rhs.traceback_ptr;
//    }
    
    
    void simpleCell::clear(){
        viti=-INFINITY;
        forw=-INFINITY;
        
        emmCalculated=false;
        emm=-INFINITY;
        ptr=-1;
        // trans.clear();
        statePtr = NULL;
        //position = -1;
        //traceback_ptr=NULL;
        return;
    }
    
    
    

    //!Create a stochastic trellis cell
    //!\param state Number of states in the model
    stochCell::stochCell(size_t states): simpleCell(){
        posteriorProbability.assign(states,-INFINITY);
        viterbiProbability.assign(states,-INFINITY);
        log_trans.assign(states,-INFINITY);
        viterbiSum=-INFINITY;
        return;
    }
    
    //!Create a stochastic trellis cell 
    stochCell::stochCell(): simpleCell(){
        viterbiSum=-INFINITY;
    }
    
    //!Print the stochastic trellis cell to stdout
    void stochCell::print(){
        simpleCell::print();
        std::cout << "Viterbi Sum:\t" << exp(viterbiSum) << std::endl;
        std::cout << "Posterior Probability:";
        for (int i=0;i<posteriorProbability.size();i++){
            std::cout << "\t" << posteriorProbability[i];
        }
        std::cout << std::endl << "Viterbi Probability:";
        for (int i=0;i<viterbiProbability.size();i++){
            std::cout << "\t" << viterbiProbability[i];
        }
        std::cout <<std::endl;
        return;
    }
    
    //!Calculate the Viterbi and posterior probabilities for the cell
    //! Based on the Forward and Viterbi scores
    void stochCell::calcForwardViterbiProb(){
        if (forw==-INFINITY && viterbiSum==-INFINITY){
            return;
        }
        else if (forw>-INFINITY && viterbiSum>-INFINITY){
            for(int i=0;i<posteriorProbability.size();i++){
                posteriorProbability[i]=exp(posteriorProbability[i]-forw);
                viterbiProbability[i]=exp(viterbiProbability[i]-viterbiSum);
            }
        }
        else if (forw>-INFINITY){
            for(int i=0;i<posteriorProbability.size();i++){
                posteriorProbability[i]=exp(posteriorProbability[i]-forw);
            }
        }
        else if (viterbiSum>-INFINITY){
            for(int i=0;i<viterbiProbability.size();i++){
                viterbiProbability[i]=exp(viterbiProbability[i]-viterbiSum);
            }
        }
    }
    
    //!Calculates the forward probabilities of the cell
    void stochCell::calcForwardProb(){
        if (forw>-INFINITY){
            for(int i=0;i<posteriorProbability.size();i++){
                posteriorProbability[i]=exp(posteriorProbability[i]-forw);
            }
        }
        return;
    }
    
    
    //!Calculate the Viterbi Probabilities
    void stochCell::calcViterbiProb(){
        if (viterbiSum>-INFINITY){
            for(int i=0;i<posteriorProbability.size();i++){
                viterbiProbability[i]=exp(viterbiProbability[i]-viterbiSum);
            }
        }
        return;
    }

//    //!Create an nth trellis cell
//    nthCell::nthCell() : simpleCell(){
//        return;
//    }
//    
//    
//    //!Create an nth trellis cell with number of states
//    //!\param stateSize  Number of states in the model
//    nthCell::nthCell(int stateSize) : simpleCell(){
//        scores temp;
//        temp.viterbi_score = -INFINITY;
//        temp.traceback_state = -2;
//        temp.traceback_state_score = -2;
//        viterbiScores.assign(stateSize,temp);
//        return;
//    }
//
//    //!Print the nth trellis cell to stdout
//    void nthCell::print(){
//        for(int i=0;i<viterbiScores.size();i++){
//            std::cout << "Viterbi Score: " << viterbiScores[i].viterbi_score << std::endl;
//            std::cout << "Traceback State: " << viterbiScores[i].traceback_state << std::endl;
//            std::cout << "Traceback State Score: " << viterbiScores[i].traceback_state_score << std::endl;
//        }
//        return;
//    }

    //!Sort the viterbi scores in the nth trellis cells
    void sort_scores(std::vector<scores>& nth_scores){
        sort(nth_scores.begin(), nth_scores.end(), _vec_sort );
        return;
    }

    //Sort vector of pairs using the first value in the pair
    bool _vec_sort(const scores& i, const scores& j){
        return (i.viterbi_score > j.viterbi_score);
    }

}
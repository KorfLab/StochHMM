//
//  trellisCells.h


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


#ifndef TRELLISCELLS_H
#define TRELLISCELLS_H

#include <vector>
#include <limits>
#include <math.h>
#include <iostream>
#include <algorithm>
#include <map>
#include "stochTypes.h"
#include "state.h"
namespace StochHMM{

/*! \file  trellisCells.h
    \brief Definitions of trellis cell structure
    Trellis cells save different information depending on what kind of analysis
    is being done. For example, if a simple viterbi is being performed.
 You only need the viterbi score and no the viterbi scores coming from other states.
 */
    class simpleCell;
    
    
    class scores{
    public:
        scores();
        double viterbi_score;
        int traceback_state;
        size_t traceback_state_score;
    };

    bool _vec_sort(const scores& , const scores&);
    void sort_scores(std::vector<scores>&);

    //! \class simpleCell
    //! A class for that describe each cell in the trellis of for the simple 
    //! decoding algorithms
    
    class simpleCell{
    public:
        simpleCell();
        ~simpleCell();
        //void operator= (simpleCell&);
        void clear();
        void print(); //print cell to stdout;
        
        double viti;  //stores viterbi score
        double forw;  //stores forward 
        double back;  //stores backward algorithm scores
        
        int ptr;
        
        //Experimental
        bool emmCalculated;
        double emm; //Store emission probability for cell
        std::map<state*,double > trans;  //Store transition probabilities for cell
        
        
        state* statePtr; //State that current cell represents
        
        std::vector<scores>* nth_viterbi_scores;
    };
    
    
    //! \class stochCell
    //! Stochastic cell for trellis
    //! Stores viterbi and forward scores for all states
    class stochCell : public simpleCell{
    public:
        stochCell();
        stochCell(size_t);
        std::vector <double> posteriorProbability;  //stores traceback probability using forward
        std::vector <double> viterbiProbability; //stores traceback probability using viterbi(stochastic)
        std::vector <double> log_trans; //store transition probability for the backward algorithm
        double viterbiSum;  //sum of all viterbi scores. used for figuring vprob
        void print();
        void calcForwardViterbiProb();
        void calcViterbiProb();
        void calcForwardProb();
    };
    
   

    
//    //! \class nthCell
//    //! Stores n viterbi scores from each state
//    class nthCell : public simpleCell{
//    public:
//        nthCell();
//        nthCell(int);
//        
//                
//        std::vector<scores> viterbiScores;
//        
//        
//        void print();
//        void sortScores();
//    private:
//    };

}
#endif /*TRELLISCELLS_H*/
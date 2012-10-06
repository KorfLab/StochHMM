//
//  Counter.h
//  StochHMM
//
//  Created by Paul Lott on 1/19/12.
//  Copyright 2012 University of California, Davis.
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

#ifndef Counter_H
#define Counter_H

#include <iostream>
#include <vector>
#include <string>
#include "stochMath.h"
#include "track.h"
#include "sequence.h"
#include "sequences.h"
#include "sequenceStream.h"
#include "lexicalTable.h"
#include "index.h"

namespace StochHMM {

    extern bool quiet;

    typedef std::vector<std::vector<int> > table;

    enum countType {GENERAL, PERIODIC, PWM, MASK, NONE};

    class Counter{
    public:
        Counter();
        Counter(int, int, int, int, int);
        Counter(int, int, int, int);
        Counter(int, int, int);
        Counter(int);
        Counter(sequence&, countType, track*, int, int, int, int, int);
        Counter(track*, int, int, int, int, int);
        Counter(track*, int, int, int, int);
        Counter(track*, int, int, int);
        Counter(track*, int);
        Counter(std::vector<sequence*>&, countType, track*, int, int, int, int, int);
        Counter(sequenceStream&, countType, track*, int, int, int, int, int);

        //~Counter();
        
        bool countGeneral(sequence&, int, int);
        bool countGeneral(std::vector<sequence*>&, int, int);
        bool countGeneral(sequenceStream&, int, int);
        
        bool countPeriodic(sequence&, int, int);
        bool countPeriodic(std::vector<sequence*>&, int, int);
        bool countPeriodic(sequenceStream&, int, int);
        
        bool countPWM(sequence&, int, int);
        bool countPWM(std::vector<sequence*>&, int, int);
        bool countPWM(sequenceStream&, int, int);
        
        bool countMask(sequence&, int, int);
        bool countMask(std::vector<sequence*>&, int, int);
        //bool countMask(sequenceStream&, int, int);

        //!Working copy of countGeneral
        bool _count(sequence&, std::vector<lexicalTable*>&);
        
        void clear();
        
        //ACCESSOR
        void printTable();
        
        //MUTATOR
        void setOrder(size_t);
        void setPeriod(size_t);
        void setPseudoCount(size_t);
        void setTrack(track*);
        
    private:
        
        //!Total number of sequences counted
        int sequencesCounted;
        
        //!Bp of all sequences counted
        int lengthOfAllSequences;

        countType type;
        
        //!The track defined for our sequence(s)
        track* seqtrk;
        
        //!Creates x number of lexical tables and appends them to the vector of lexical tables
        void initializeTable(size_t);
        
        //!Used only in countPWM; stores the length of the first input sequence
        //!for making sure each sequence being counted is the same length
        int seqLength; 
        
        size_t period;
        size_t order;
        size_t pseudoCount;
        
        //!Indicates whether or not the variables are set
        //!When true, counting has begun, changing the variables would result in an error; must use clear() first to change
        //!When false, counting has not yet begun, variables can be changed
        bool variableSet; 
        
        //!The number of bp upstream to ignore
        int upstream;
        
        //!The number of bp downstream to ignore
        int downstream;
        
        std::vector<lexicalTable*> tables;
    };


}

#endif /*StochHMM_Counter_H*/

//
//  Counter.cpp
//  StochHMM
//
//  Created by Paul Lott on 1/19/12.
//  Copyright 2012 University of California, Davis.
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

#include "Counter.h"

namespace StochHMM {

    
    Counter::Counter(){
        period=0;
        upstream=0;
        downstream=0;
        order=0;
        pseudoCount=0;
        type = NONE;
        variableSet=false;
    }
    
    Counter::Counter(int per, int up, int down, int ord, int pseud) {
        period=per;
        upstream=up;
        downstream=down;
        order=ord;
        pseudoCount=pseud;
        type = NONE;
        variableSet=false;
        
    }
    
    Counter::Counter(int per, int up, int down, int ord) {
        period=per;
        upstream=up;
        downstream=down;
        order=ord;
        pseudoCount=0;
        type = NONE;
        variableSet=false;
        
    }
    
    Counter::Counter(int per, int up, int down) {
        period=per;
        upstream=up;
        downstream=down;
        order=0;
        pseudoCount=0;
        type = NONE;
        variableSet=false;
        
    }
    
    Counter::Counter(int per) {
        period=per;
        upstream=0;
        downstream=0;
        order=0;
        pseudoCount=0;
        type = NONE;
        variableSet=false;
    }
    
    Counter::Counter(sequence &seq, countType typ, track* trk, int per, int up, int down, int ord, int pseud){
        
        period = per;
        upstream = up;
        downstream = down;
        order = ord;
        pseudoCount = pseud;
        seqtrk = trk;
        variableSet = false;
        
        if (typ == GENERAL) {
            countGeneral(seq, up, down);
        }
        else if (typ == PERIODIC) {
            countPeriodic(seq, up, down);
        }
        else if (typ == PWM) {
            countPWM(seq, up, down);
        }
        else if (typ == MASK) {
            countMask(seq, up, down);
        }
        else if (typ == NONE) {
            std::cerr << "Need to define count type" << std::endl; 
            exit(1);
        }
        
    }
    
    Counter::Counter(track* trk, int per, int up, int down, int ord, int pseud) {
        seqtrk = trk;
        period=per;
        upstream=up;
        downstream=down;
        order=ord;
        pseudoCount=pseud;
        type = NONE;
        variableSet=false;
        
    }
    
    Counter::Counter(track* trk, int per, int up, int down, int ord) {
        seqtrk = trk;
        period=per;
        upstream=up;
        downstream=down;
        order=ord;
        pseudoCount=0;
        type = NONE;
        variableSet=false;
        
    }
    
    Counter::Counter(track* trk, int per, int up, int down) {
        seqtrk = trk;
        period=per;
        upstream=up;
        downstream=down;
        order=0;
        pseudoCount=0;
        type = NONE;
        variableSet=false;
        
    }
    
    Counter::Counter(track* trk, int per) {
        seqtrk = trk;
        period=per;
        upstream=0;
        downstream=0;
        order=0;
        pseudoCount=0;
        type = NONE;
        variableSet=false;
    }
    
    
    Counter::Counter(std::vector<sequence*> &seqs, countType typ, track* trk, int per, int up, int down, int ord, int pseud){
        
        period = per;
        upstream = up;
        downstream = down;
        order = ord;
        pseudoCount = pseud;
        seqtrk = trk;
        variableSet = false;
        
        if (typ == GENERAL) {
            countGeneral(seqs, up, down);
        }
        else if (typ == PERIODIC) {
            countPeriodic(seqs, up, down);
        }
        else if (typ == PWM) {
            countPWM(seqs, up, down);
        }
        else if (typ == MASK) {
            countMask(seqs, up, down);
        }
        else if (typ == NONE) {
            std::cerr << "Need to define count type" << std::endl; 
            exit(1);
        }
        
    }
    
    Counter::Counter(sequenceStream &seq, countType typ, track* trk, int per, int up, int down, int ord, int pseud){
    
        period = per;
        upstream = up;
        downstream = down;
        order = ord;
        pseudoCount = pseud;
        seqtrk = trk;
        variableSet = false;
        
        if (typ == GENERAL) {
            countGeneral(seq, up, down);
        }
        else if (typ == PERIODIC) {
            countPeriodic(seq, up, down);
        }
        else if (typ == PWM) {
            countPWM(seq, up, down);
        }
        else if (typ == MASK) {
            countMask(seq, up, down);
        }
        else if (typ == NONE) {
            std::cerr << "Need to define count type" << std::endl; 
            exit(1);
        }
        
    }
    
    //working copy of countGeneral
    bool Counter::_count(sequence &seq, std::vector<lexicalTable*>& tbl){

        bool counted = false;
        
        int seqln = seq.getLength();
        lengthOfAllSequences += seqln;
        
        for(int i=upstream; i<=seqln-downstream-1;i++){
            int per;
            
            if (type == MASK) {
                per = seq.getMask(i) - 1;
            }
            else {
                per = (i-upstream)%period;
            }
            
            std::pair<Index,Index> word_index;
            seq.get_index(i, order, word_index);
            int letters = word_index.first.size();
            int words = word_index.second.size();
            if (letters!=1 && words!=1){
                double val = 1.0f / static_cast<double> (letters * words);
                //increment all the applicable cells of the table by the value 'val'
                for (int i=0; i<words; i++) {
                    for (int j=0; j<letters; j++) {
                        (tbl[per])->incrementCountsDouble(word_index.second[i], word_index.first[j], val);
                    }
                }
                counted = true;
            }
            else{
                //increment the lexical table (letters,words)
                for (int i=0; i<words; i++) {
                    for (int j=0; j<letters; j++) {
                        (tbl[per])->incrementCounts(word_index.second[i], word_index.first[j]);
                    }
                }
                counted = true;
            }
        }
        
        sequencesCounted++;
        
        return counted;
    }
    
    bool Counter::countGeneral(sequence &seq, int up, int down) {
        
        
        if (variableSet == true) {
            if (up != upstream) {
                std::cerr << "Upstream: " << up << " is not initialized value: " << upstream << std::endl;
                exit(1);
            }
            else if (down != downstream) {
                std::cerr << "Downstream: " << down << " is not initialized value: " << downstream << std::endl;
                exit(1);
            }
        }
        else if (variableSet == false) {
            
            int seqln = seq.getLength();
            
            type = GENERAL;
            if (period != 1) {
                std::cout << "Period set to 1 for count general (was " << period << ")" << std::endl;
            }
            period = 1;
            initializeTable(1);
            
            //check if upstream is positive, not too small (smaller than order), 
            // and not too large (larger than seqLength)
            if (up < 0 || up < order || up > (seqln-1) ) {
                std::cerr << "Upstream must be non-negative, equal to or larger than the order, and smaller than the sequence length -1" << std::endl;
                exit(1);
            }
            
            //check if downstream is positive, and not too large (seqLength-downstream > upstream)
            if (down < 0 || (seqln-down-1) < up) {
                std::cerr << "Downstream must be non-negative, and at least as far in the sequence as upstream" << std::endl;
                exit(1);
            }
            upstream = up;
            downstream = down;
            variableSet = true;
        }
        
        bool counted;
        
        std::vector<lexicalTable*>& tbl=tables;
        counted = _count(seq,tbl);

        return counted;
    }
    
    bool Counter::countGeneral(std::vector<sequence*> &seqs, int up, int down){
        bool counted = true;
        for (size_t i=0; i<seqs.size(); i++) {
            if (countGeneral(*(seqs[i]), up, down) == false) {
                counted = false;
            }
        }
        return counted;
    }
    
    bool Counter::countGeneral(sequenceStream &seq, int up, int down) {
        
        if (variableSet == true) {
            if (up != upstream) {
                std::cerr << "Upstream: " << up << " is not initialized value: " << upstream << std::endl;
                exit(1);
            }
            else if (down != downstream) {
                std::cerr << "Downstream: " << down << " is not initialized value: " << downstream << std::endl;
                exit(1);
            }
        }
        else if (variableSet == false) {
            
            int seqln = seq.getLength();
            
            type = GENERAL;
            if (period != 1) {
                std::cout << "Period set to 1 for count general (was " << period << ")" << std::endl;
            }
            period = 1;
            initializeTable(1);
            
            //check if upstream is positive, not too small (smaller than order), 
            // and not too large (larger than seqLength)
            if (up < 0 || up < order || up > (seqln-1) ) {
                std::cerr << "Upstream must be non-negative, equal to or larger than the order, and smaller than the sequence length -1" << std::endl;
                exit(1);
            }
            
            //check if downstream is positive, and not too large (seqLength-downstream > upstream)
            if (down < 0 || (seqln-down-1) < up) {
                std::cerr << "Downstream must be non-negative, and at least as far in the sequence as upstream" << std::endl;
                exit(1);
            }
            upstream = up;
            downstream = down;
            variableSet = true;
        }
        
        bool counted;
        
        std::vector<lexicalTable*>& tbl=tables;
        counted = _count(seq,tbl);
        
        return counted;
    }
    
    bool Counter::countPWM(sequence &seq, int up, int down){
        
        int seqln = seq.getLength();
        
        if (variableSet == true) {
            if (up != upstream) {
                std::cerr << "Upstream: " << up << " is not initialized value: " << upstream << std::endl;
                exit(1);
            }
            else if (down != downstream) {
                std::cerr << "Downstream: " << down << " is not initialized value: " << downstream << std::endl;
                exit(1);
            }
            else if (seqln != seqLength) {
                std::cerr << "Seqlength: " << seqln << " is not initialized value: " << seqLength << std::endl;
                exit(1);
            }
        }
        //Setup variables on first call
        else if (variableSet == false) {
            type = PWM;
            initializeTable(period);
            seqLength = seqln;

            if (up < 0 || up < order || up > (seqln-1) ) {
                std::cerr << "Upstream must be non-negative, equal to or larger than the order, and smaller than the sequence length -1" << std::endl;
                exit(1);
            }
            
            if (down < 0 || (seqln-down-1) < up) {
                std::cerr << "Downstream must be non-negative, and at least as far in the sequence as upstream" << std::endl;
                exit(1);
            }
            
            upstream = up;
            downstream = down;
            variableSet = true;
        }
        
        bool counted;
        
        std::vector<lexicalTable*>& tbl=tables;
        counted = _count(seq,tbl);
        
        return counted;
    }
    
    bool Counter::countPWM(std::vector<sequence*>& seqs, int up, int down){
        bool counted = true;
        for (size_t i=0; i<seqs.size(); i++) {
            if (countPWM(*(seqs[i]), up, down) == false) {
                counted = false;
            }
        }
        return counted;
    }
    
    bool Counter::countPWM(sequenceStream &seq, int up, int down){
        
        int seqln = seq.getLength();
        
        if (variableSet == true) {
            if (up != upstream) {
                std::cerr << "Upstream: " << up << " is not initialized value: " << upstream << std::endl;
                exit(1);
            }
            else if (down != downstream) {
                std::cerr << "Downstream: " << down << " is not initialized value: " << downstream << std::endl;
                exit(1);
            }
            else if (seqln != seqLength) {
                std::cerr << "Seqlength: " << seqln << " is not initialized value: " << seqLength << std::endl;
                exit(1);
            }
        }
        //Setup variables on first call
        else if (variableSet == false) {
            type = PWM;
            initializeTable(period);
            seqLength = seqln;
            
            if (up < 0 || up < order || up > (seqln-1) ) {
                std::cerr << "Upstream must be non-negative, equal to or larger than the order, and smaller than the sequence length -1" << std::endl;
                exit(1);
            }
            
            if (down < 0 || (seqln-down-1) < up) {
                std::cerr << "Downstream must be non-negative, and at least as far in the sequence as upstream" << std::endl;
                exit(1);
            }
            
            upstream = up;
            downstream = down;
            variableSet = true;
        }
        
        bool counted;
        
        std::vector<lexicalTable*>& tbl=tables;
        counted = _count(seq,tbl);
        
        return counted;
    }
    
    bool Counter::countPeriodic(sequence& seq, int up, int down) {
        
        if (variableSet == true) {
            if (up != upstream) {
                std::cerr << "Upstream: " << up << " is not initialized value: " << upstream << std::endl;
                exit(1);
            }
            else if (down != downstream) {
                std::cerr << "Downstream: " << down << " is not initialized value: " << downstream << std::endl;
                exit(1);
            }
        }
        else if (variableSet == false) {
            
            int seqln = seq.getLength();
            
            type = PERIODIC;
            initializeTable(period);
            
            if (up < 0 || up < order || up > (seqln-1) ) {
                std::cerr << "Upstream must be non-negative, equal to or larger than the order, and smaller than the sequence length -1" << std::endl;
                exit(1);
            }
            
            if (down < 0 || (seqln-down-1) < up) {
                std::cerr << "Downstream must be non-negative, and at least as far in the sequence as upstream" << std::endl;
                exit(1);
            }
            
            upstream = up;
            downstream = down;
            variableSet = true;
        }
        
        bool counted;
        
        std::vector<lexicalTable*>& tbl=tables;
        counted = _count(seq,tbl);
        
        return counted;
        
    }
    
    bool Counter::countPeriodic(std::vector<sequence*> &seqs, int up, int down){
        bool counted = true;
        for (size_t i=0; i<seqs.size(); i++) {
            if (countPeriodic(*(seqs[i]), up, down) == false) {
                counted = false;
            }
        }
        return counted;
    }
    
    bool Counter::countPeriodic(sequenceStream& seq, int up, int down) {
        
        if (variableSet == true) {
            if (up != upstream) {
                std::cerr << "Upstream: " << up << " is not initialized value: " << upstream << std::endl;
                exit(1);
            }
            else if (down != downstream) {
                std::cerr << "Downstream: " << down << " is not initialized value: " << downstream << std::endl;
                exit(1);
            }
        }
        else if (variableSet == false) {
            
            int seqln = seq.getLength();
            
            type = PERIODIC;
            initializeTable(period);
            
            if (up < 0 || up < order || up > (seqln-1) ) {
                std::cerr << "Upstream must be non-negative, equal to or larger than the order, and smaller than the sequence length -1" << std::endl;
                exit(1);
            }
            
            if (down < 0 || (seqln-down-1) < up) {
                std::cerr << "Downstream must be non-negative, and at least as far in the sequence as upstream" << std::endl;
                exit(1);
            }
            
            upstream = up;
            downstream = down;
            variableSet = true;
        }
        
        bool counted;
        
        std::vector<lexicalTable*>& tbl=tables;
        counted = _count(seq,tbl);
        
        return counted;
        
    }
    
    
    bool Counter::countMask(sequence& seq, int up, int down){
        
        if (variableSet == true) {
            if (up != upstream) {
                std::cerr << "Upstream: " << up << " is not initialized value: " << upstream << std::endl;
                exit(1);
            }
            else if (down != downstream) {
                std::cerr << "Downstream: " << down << " is not initialized value: " << downstream << std::endl;
                exit(1);
            }
        }
        else if (variableSet == false) {
            
            int seqln = seq.getLength();
            
            type = MASK;
            initializeTable(period);
            
            if (up < 0 || up < order || up > (seqln-1) ) {
                std::cerr << "Upstream must be non-negative, equal to or larger than the order, and smaller than the sequence length -1" << std::endl;
                exit(1);
            }
            
            if (down < 0 || (seqln-down-1) < up) {
                std::cerr << "Downstream must be non-negative, and at least as far in the sequence as upstream" << std::endl;
                exit(1);
            }
            
            upstream = up;
            downstream = down;
            variableSet = true;
        }
        
        int max_mask_size = seq.getMaxMask();
        size_t num_tables = tables.size();    
        
        if (max_mask_size > num_tables){
            initializeTable(max_mask_size-num_tables);
        }

        bool counted;
        
        std::vector<lexicalTable*>& tbl=tables;
        counted = _count(seq,tbl);
        
        return counted;
        
    }
    
    bool Counter::countMask(std::vector<sequence*> &seqs, int up, int down){
        bool counted = true;
        for (size_t i=0; i<seqs.size(); i++) {
            //std::cout << seqs[i]->stringify() << std::endl;
            if (countMask(*(seqs[i]), up, down) == false) {
                counted = false;
            }
        }
        return counted;
    }
    
    void Counter::initializeTable(size_t table_count){
        size_t cols = seqtrk->getAlphaSize();
        size_t rows = integerPower<size_t>(cols,order);
        //size_t cols = pow(rows,order);
        
        for (int i = 0; i < table_count; i++) {
            lexicalTable *table = new lexicalTable();
            table->createTable(rows, cols, pseudoCount, COUNTS);
            //std::cout << order << std::endl;
            table->addTrack(seqtrk, order);
            //table->print();
            //exit(0);
            tables.push_back(table);
        }
        
        return;
    }
    
    void Counter::clear() {
        for(int i=0; i<tables.size(); i++) {
            delete tables[i];
        }
        order=0;
        period=0;
        upstream=0;
        downstream=0;
        pseudoCount=0;
        variableSet=false;
        type = NONE;
    }
    
    void Counter::printTable() {
        
        std::cout << "Number of sequences counted: " << sequencesCounted << std::endl;
        std::cout << "Number of bp counted: " << lengthOfAllSequences << std::endl;
        std::cout << "Upstream/downstream settings: " << upstream << "/" << downstream << std::endl;
        std::cout << "Order: " << order << std::endl;
        
        for(int i=0; i < tables.size(); i++) {
            tables[i]->print();
        }
    }
    
    void Counter::setOrder(size_t ord) {
        if (variableSet == false) {
            order = ord;
        }
        else {
            std::cerr << "Already started counting, order was set to " << order << ". Must use clear() to set order again." << std::endl;
            exit(1);
        }
    }
    
    void Counter::setPeriod(size_t per) {
        if (variableSet == false) {
            period = per;
        }
        else {
            std::cerr << "Already started counting, period was set to " << period << ". Must use clear() to set period again." << std::endl;
            exit(1);
        }
    }
    
    void Counter::setPseudoCount(size_t pseud) {
        if (variableSet == false) {
            pseudoCount = pseud;
        }
        else {
            std::cerr << "Already started counting, pseudocount was set to " << pseudoCount << 
            ". Must use clear() to set pseudocount again." << std::endl;
            exit(1);
        }
    }
    
    void Counter::setTrack(track* trk) {
        if (variableSet== false) {
            seqtrk = trk;
        }
        else {
            std::cerr << "Track is already set for counting. Must use clear() to set again." << std::endl;
        }
    }
}    
    
    //counts::counts(countDirection dir, countType typ, int per, int up, int down,int ord, int pseudo){
    //    direct=dir;
    //    type=typ;
    //    period=per;
    //    upstream=up;
    //    downstream=down;
    //    PWMSeqSize=0;
    //    order=ord;
    //    pseudoCount=pseudo;
    //    
    //    if (period<=0 && type==PERIODIC){
    //        std::cerr << "Period must be greater than or equal to one for Periodic Counting" << std::endl;
    //        exit(1);
    //    }
    //    else if (type==PERIODIC){
    //        for(int i=0;i<period;i++){
    //            if (direct==FORWARD || direct==BOTH){
    //                std::vector<int> temp (4,pseudoCount);
    //                forward.push_back(new table(POWER4[order][1],temp));
    //            }
    //            
    //            if (direct==REVERSE || direct==BOTH){
    //                std::vector<int> temp (4,pseudoCount);
    //                reverse.push_back(new table(POWER4[order][1],temp));
    //            }
    //        }
    //    }
    //}
    //
    //counts::~counts(){
    //    for(int i=0;i<forward.size();i++){
    //        delete forward[i];
    //    }
    //    for(int i=0;i<reverse.size();i++){
    //        delete forward[i];
    //    }
    //}
    //
    //
    //void counts::count(sequence& seq){
    //    
    //    size_t seqSize=seq.size();
    //    
    //    if (seqSize==0){
    //        return;
    //    }
    //    
    //    if (type==PERIODIC){
    //        if (direct == FORWARD){
    //            countPeriodic(seq,upstream,downstream,FORWARD);
    //        }
    //        else if (direct == REVERSE){
    //            seq.reverseComplement();
    //            countPeriodic(seq,downstream,upstream,REVERSE);
    //        }
    //        else {
    //            countPeriodic(seq,upstream,downstream,FORWARD);
    //            
    //            seq.reverseComplement();
    //            countPeriodic(seq,downstream,upstream,REVERSE);
    //        }
    //    }
    //    else{
    //        if (PWMSeqSize==0){
    //            PWMSeqSize=seqSize;
    //            
    //            initPWM(seqSize-(upstream+downstream));
    //        }
    //        else if (seqSize != PWMSeqSize){
    //            cerr << seq.getHeader() << " is not the same size as the first sequence.   When doing PWM all the sequence size must be the same." <<endl;
    //            return;
    //        }
    //        
    //        if (direct==FORWARD){
    //            countPWM(seq,upstream,downstream,FORWARD);
    //        }
    //        else if (direct==REVERSE){
    //            seq.reverseComplement();
    //            countPWM(seq,downstream,upstream,REVERSE);
    //        }
    //        else{
    //            countPWM(seq,upstream,downstream,FORWARD);
    //            seq.reverseComplement();
    //            countPWM(seq,downstream,upstream,REVERSE);
    //        }
    //    }
    //}
    //
    //
    //void counts::countPeriodic(sequence& seq,int up, int down,countDirection dir){
    //    
    //    std::vector<sequence*>Counted++;
    //    lengthOfAllSequences+=seq.size()-(upstream+downstream);
    //    
    //    vector<table*>& tbl=(dir==FORWARD) ? forward : reverse;
    //    
    //    for(int i=up;i<seq.size()-down;i++){
    //        
    //        int per;
    //        
    //        if (seq.isMasked()){
    //            per=seq.getMask(i);
    //            if (per<0){
    //                continue;
    //            }
    //        }
    //        else{
    //            per=(i-up)%period;
    //        }
    //        
    //        
    //        int xCoord=seq[i];
    //        if (xCoord==-1){
    //            return;
    //        }
    //        
    //        if (order==0){
    //            (*tbl[per])[0][xCoord]++;
    //        }
    //        else{
    //            int yCoord=0;
    //            if ((i-order)<0){  //If the order would have us look before the start of the sequence, we need to skip to the next position;
    //                continue;
    //            }
    //            for(int k=order;k>=1;k--){
    //                int prev_letter=seq[i-k];
    //                if (prev_letter==-1){
    //                    return;
    //                }
    //                yCoord+=prev_letter*POWER4[k-1][1];
    //            }
    //            
    //            
    //            (*tbl[per])[yCoord][xCoord]++;
    //        }
    //    }
    //}
    //
        //
    //
    //void counts::countPWM(sequence& seq,int up ,int down, countDirection dir){
    //    
    //    sequencesCounted++;
    //    lengthOfAllSequences+=seq.size()-(upstream+downstream);
    //    
    //    vector<table*>& tbl=(dir==FORWARD) ? forward : reverse;
    //    
    //    for(int i=up;i<seq.size()-down;i++){
    //        int xCoord=seq[i];
    //        (*tbl[period])[i-up][xCoord]++;
    //    }
    //    return;
    //}
    //
    //
    //void counts::print(bool Words){
    //    
    //    if (type==PERIODIC){
    //        if (!quiet){
    //            cout << "#Periodic:\t" << period << endl;
    //            cout << "#Order:\t" << order << endl;
    //            cout << "#PseudoCount:\t" << pseudoCount <<endl;
    //            cout << "#Upstream(bp):\t" << upstream << endl;
    //            cout << "#Downstream(bp):\t" << downstream << endl;
    //        }
    //        
    //        
    //        if (direct==FORWARD){
    //            if (!quiet){
    //                cout << "#Number of Sequences:\t" << sequencesCounted <<endl;
    //                cout << "#Average Size:\t" << (double) lengthOfAllSequences/ (double) sequencesCounted << endl;
    //                cout << "#Direction:\tForward" <<endl;
    //            }
    //            for(int i=0;i<forward.size();i++){
    //                cout << "PERIOD:\t" << i+1 <<endl;
    //                if (Words){
    //                    cout << "PREV\tA\tC\tG\tT\n";
    //                }
    //                for(int j=0;j<(*forward[i]).size();j++){
    //                    if (Words){
    //                        cout << indexToWord(j, order) << "\t";
    //                    }
    //                    cout << (*forward[i])[j][0] << "\t" ;
    //                    cout << (*forward[i])[j][1] << "\t" ;
    //                    cout << (*forward[i])[j][2] << "\t" ;
    //                    cout << (*forward[i])[j][3] << endl;
    //                }
    //                cout << endl;
    //            }
    //        }
    //        else if (direct==REVERSE){
    //            if (!quiet){
    //                cout << "#Number of Sequences:\t" << sequencesCounted <<endl;
    //                cout << "#Average Size:\t" << (double) lengthOfAllSequences/ (double) sequencesCounted << endl;
    //                cout << "#Direction:\tReverse" <<endl;
    //            }
    //            for(int i=0;i<reverse.size();i++){
    //                cout << "PERIOD:\t" << i+1 <<endl;
    //                if (Words){
    //                    cout << "PREV\tA\tC\tG\tT\n";
    //                }
    //                for(int j=0;j<(*reverse[i]).size();j++){
    //                    if (Words){
    //                        cout << indexToWord(j, order) << "\t";
    //                    }
    //                    cout << (*reverse[i])[j][0] << "\t" ;
    //                    cout << (*reverse[i])[j][1] << "\t" ;
    //                    cout << (*reverse[i])[j][2] << "\t" ;
    //                    cout << (*reverse[i])[j][3] << endl;
    //                }
    //                cout << endl;
    //            }
    //        }
    //        else if (direct==BOTH){
    //            if (!quiet){
    //                cout << "#Number of Sequences:\t" << sequencesCounted/2 <<endl;
    //                cout << "#Average Size:\t" << (double) lengthOfAllSequences/ (double) sequencesCounted << endl;
    //                cout << "#Direction:\BOTH" <<endl;
    //            }
    //            for(int i=0;i<reverse.size();i++){
    //                cout << "PERIOD:\t" << i+1 <<endl;
    //                if (Words){
    //                    cout << "PREV\tA\tC\tG\tT\n";
    //                }
    //                for(int j=0;j<(*reverse[i]).size();j++){
    //                    if (Words){
    //                        cout << indexToWord(j, order) << "\t";
    //                    }
    //                    cout << (*reverse[i])[j][0] + (*forward[i])[j][0]<< "\t";
    //                    cout << (*reverse[i])[j][1] + (*forward[i])[j][1]<< "\t";
    //                    cout << (*reverse[i])[j][2] + (*forward[i])[j][2]<< "\t";
    //                    cout << (*reverse[i])[j][3] + (*forward[i])[j][3]<< endl;
    //                }
    //                cout << endl;
    //            }
    //        }
    //    }
    //    else{
    //        if (!quiet){
    //            cout << "#PWM:\t" << PWMSeqSize-(upstream+downstream) << endl;
    //            cout << "#PseudoCount:\t" << pseudoCount <<endl;
    //            cout << "#Upstream(bp):\t" << upstream << endl;
    //            cout << "#Downstream(bp):\t" << downstream << endl;
    //        }
    //        if (direct==FORWARD){
    //            if (!quiet){
    //                cout << "#Number of Sequences:\t" << sequencesCounted <<endl;
    //                cout << "#Direction:\tForward" <<endl;
    //                cout << "POSITION\tA\tC\tG\tT\n";
    //            }
    //            for(int i=0;i<forward.size();i++){
    //                for(int j=0;j<(*forward[i]).size();j++){
    //                    cout << j+1 << "\t";
    //                    cout << (*forward[i])[j][0] << "\t" ;
    //                    cout << (*forward[i])[j][1] << "\t" ;
    //                    cout << (*forward[i])[j][2] << "\t" ;
    //                    cout << (*forward[i])[j][3] << endl;
    //                }
    //                cout << endl;
    //            }
    //            
    //        }
    //        else if (direct==REVERSE){
    //            if (!quiet){
    //                cout << "#Number of Sequences:\t" << sequencesCounted <<endl;
    //                cout << "#Direction:\tReverse" <<endl;
    //                cout << "#POSITION\tA\tC\tG\tT\n";
    //            }
    //            for(int i=0;i<forward.size();i++){
    //                for(int j=0;j<(*forward[i]).size();j++){
    //                    cout << j+1 << "\t";
    //                    cout << (*reverse[i])[j][0] << "\t" ;
    //                    cout << (*reverse[i])[j][1] << "\t" ;
    //                    cout << (*reverse[i])[j][2] << "\t" ;
    //                    cout << (*reverse[i])[j][3] << endl;
    //                }
    //                cout << endl;
    //            }
    //            
    //        }
    //        else if (direct==BOTH){
    //            if (!quiet){
    //                cout << "#Number of Sequences:\t" << sequencesCounted <<endl;
    //                cout << "#Direction:\tBOTH" <<endl;
    //                cout << "#POSITION\tA\tC\tG\tT\n";
    //            }
    //            for(int i=0;i<forward.size();i++){
    //                for(int j=0;j<(*forward[i]).size();j++){
    //                    cout << j+1 << "\t";
    //                    cout << (*reverse[i])[j][0] + (*forward[i])[j][0]<< "\t" ;
    //                    cout << (*reverse[i])[j][1] + (*forward[i])[j][1]<< "\t" ;
    //                    cout << (*reverse[i])[j][2] + (*forward[i])[j][2]<< "\t" ;
    //                    cout << (*reverse[i])[j][3] + (*forward[i])[j][3]<< endl;
    //                }
    //                cout << endl;
    //            }
    //            
    //        }
    //        
    //    }
    //    return;
    //    
    //}
    //
    ////Compute word from INDEX VALUE
    //string counts::indexToWord(int index, int length){
    //	index++;
    //	string output="";
    //	while (length>0){
    //		double dreg=POWER4[length-1][1];
    //		float value=index/dreg;
    //		if (value<=1){
    //			output+="A";}
    //		else if (value<=2){
    //			output+="C";
    //			index-=dreg;}
    //		else if (value<=3){
    //			output+="G";
    //			index-=2*dreg;}
    //		else {
    //			output+="T";
    //			index-=3*dreg;
    //		}
    //		length--;
    //	}
    //	return output;
    //}


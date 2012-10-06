//
//  alignment.h
//  StochHMM
//
//  Created by Paul Lott on 5/18/12.
//  Copyright (c) 2012 University of California, Davis. All rights reserved.
//

#ifndef StochHMM_alignment_h
#define StochHMM_alignment_h

#include <iostream>
#include <vector>
#include <string>
#include <map>
#include <algorithm>
#include <stdint.h>
#include "track.h"
#include "sequence.h"
#include "seqTools.h"

namespace StochHMM {
    
    enum tracebackDirection { NONE, DIAGONAL, LEFT, UP };
    enum alignment_type { cGlobal, cLocal, cGlocal};
    
    
    class cell {
    public:
        cell();
        inline void setTraceback(tracebackDirection dir){traceback = dir;}
        inline void setScore(double value){score=value;}
        inline void setDiag(double value){scores[0]=value;}
        inline void setLeft(double value){scores[1]=value;}
        inline void setUp(double value){scores[2]=value;}
        inline double getScore(){return score;}
        inline tracebackDirection getTraceback(){return traceback;}
    private:
        tracebackDirection traceback;
        double score;
        std::vector<double> scores;    //(DIAGONAL, LEFT, UP)
    };
    
    typedef std::vector<std::vector<double> > mmMatrix;
    
    class alignment{
    public:
        alignment();
        
        
        double align(alignment_type typ);
        double align(sequence& target, sequence& query, alignment_type typ);
        double align(sequence& target, sequence& query, alignment_type typ,  double match, double mismatch, double gap, double gapext);
        double align(sequence& target, sequence& query, alignment_type typ, const mmMatrix&, double gap, double gapext);
        
        double calcLambda(std::vector<double>&);
        double estimateLambda();
        double calc_pvalue();
        
        void setMatrix(mmMatrix&);
        void setMatch(double);
        void setMismatch(double);
        
        inline void setGap(double value){gap = value;}
        inline void setGapExt(double value){gap_ext = value;}
        
        void setTarget(sequence&);
        void setQuery (sequence&);
        
        inline int getHighScore(){return high_score;}
        
        inline size_t getAlignmentSize(){
            if (end_position!=0) return end_position-start_position+1; 
            return 0;
        }
        
        
        std::string getAlignment();
        void printAlignment();
        
        void printTrellis();
        
        std::string traceback();
        std::string stochasticTraceback();
        
    private:
        
        sequence* target;
        sequence* query;
        
        track* tr;
        
        std::vector<std::vector<cell> >* trellis;
        std::vector<std::vector<double> >* mMatrix;
                
        double gap;
        double gap_ext;
        
        double high_score;
        size_t start_position;
        size_t end_position;
        
        void _reset();
        void _resetSeqs();
        double _getGapScore(size_t,size_t,tracebackDirection);
        void _initTrellis(size_t, size_t);
        
    };
    
}



#endif

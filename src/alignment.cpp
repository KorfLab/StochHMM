//
//  alignment.cpp
//  StochHMM
//
//  Created by Paul Lott on 5/18/12.
//  Copyright (c) 2012 University of California, Davis. All rights reserved.
//

#include "alignment.h"

namespace StochHMM{
    
    
    cell::cell():traceback(NONE),score(0.0){
        scores.insert(scores.begin(),3,0.0);
    }
    
    
    alignment::alignment():gap(0),gap_ext(0),high_score(0),start_position(0),target(NULL),query(NULL),tr(NULL){
        
    
    }
    
    void alignment::_reset(){
        target=NULL;
        query=NULL;
        tr=NULL;
        delete mMatrix;
        mMatrix=NULL;
        gap=0;
        gap_ext=0;
        high_score=0;
        start_position=0;
        end_position=0;
        return;
    }
    
    void alignment::_resetSeqs(){
        target=NULL;
        query=NULL;
        delete mMatrix;
        mMatrix=NULL;
        high_score=0;
        start_position=0;
        end_position=0;
    }
    
    void alignment::setMatch(double value){
        if (target==NULL && query==NULL){
            std::cerr << "Either target or query sequence must be set before defining match matrix" << std::endl;
            return;
        }
        else{
            size_t alphaSize=(target!=NULL) ? target->size() : query->size();
            if (mMatrix==NULL){
                std::vector<double> row (alphaSize, 0.0);
                mMatrix = new mmMatrix;
                mMatrix->insert(mMatrix->begin(),alphaSize,row);
                //std::cout << mMatrix->size() <<std::endl;
                //std::cout << (*mMatrix)[0].size();
            }
            else if (mMatrix->size()!=alphaSize){
                delete mMatrix;
                mMatrix=NULL;
                std::vector<double> row (alphaSize,0.0);
                mMatrix = new mmMatrix;
                mMatrix->insert(mMatrix->begin(),alphaSize,row);
            }
            
            for(size_t i=0;i<alphaSize;++i){
                (*mMatrix)[i][i]=value;
            }
            
        }
        return;
    }
    
    void alignment::setMatrix(mmMatrix& matrix){
        if (target==NULL && query==NULL){
            std::cerr << "Either target or query sequence must be set before defining match matrix" << std::endl;
            return;
        }
        else{
            size_t alphaSize=(target!=NULL) ? target->size() : query->size();
            if (mMatrix==NULL){
                mMatrix=&matrix;
            }
            else{
                delete mMatrix;
                mMatrix=&matrix;
            }
        }
        return;
    }
    
    void alignment::setMismatch(double value){
        if (target==NULL && query==NULL){
            std::cerr << "Either target or query sequence must be set before defining match matrix" << std::endl;
            return;
        }
        else{
            size_t alphaSize=(target!=NULL) ? target->size() : query->size();
            if (mMatrix==NULL){
                std::vector<double> row (0.0,alphaSize);
                mMatrix = new mmMatrix;
                mMatrix->insert(mMatrix->begin(),alphaSize,row);
            }
            else if (mMatrix->size()!=alphaSize){
                delete mMatrix;
                mMatrix=NULL;
                std::vector<double> row (0.0,alphaSize);
                mMatrix = new mmMatrix;
                mMatrix->insert(mMatrix->begin(),alphaSize,row);
            }
            
            for(size_t row=0;row<alphaSize;++row){
                for(size_t col=0;col<alphaSize;++col){
                    if (row!=col){
                        (*mMatrix)[row][col]=value;
                    }
                }
            }
        }
        return;
    }
    
    void alignment::setTarget(sequence& seq){
        track* temp = seq.getTrack();
        if (tr==temp){
            target= &seq;
            high_score=0;
            start_position=0;
            end_position=0;
        }
        else if (tr==NULL){
            target = &seq;
            tr=temp;
        }
        else{
            _resetSeqs();
            target = &seq;
            tr=temp;
        }
        return;
    }
    
    void alignment::setQuery (sequence& seq){
        track* temp = seq.getTrack();
        if (tr==temp){
            query= &seq;
            high_score=0;
            start_position=0;
            end_position=0;
        }
        else if (tr==NULL){
            query = &seq;
            tr=temp;
        }
        else{
            _resetSeqs();
            query = &seq;
            tr=temp;
        }
        return;
    }
    
    double alignment::_getGapScore(size_t queryIndex, size_t targetIndex, tracebackDirection dir){
        if (gap == gap_ext){
            return gap;
        }
        else{
            if (dir == LEFT){
                if ( (*trellis)[queryIndex][targetIndex].getTraceback() == LEFT){
                    return gap_ext;
                }
                else{
                    return gap;
                }
            }
            else{
                if ( (*trellis)[queryIndex][targetIndex].getTraceback() == UP){
                    return gap_ext;
                }
                else{
                    return gap;
                }
            }
        }
    }
    
        
    void alignment::_initTrellis(size_t rows, size_t columns){
        if (trellis!=NULL){
            delete trellis;
            trellis=NULL;
        }
        
        cell x;
        x.setTraceback(NONE);
        x.setScore(0.0);
        x.setDiag(0.0);
        x.setLeft(0.0);
        x.setUp(0.0);
        std::vector<cell> row (columns+1,x);
        trellis= new std::vector<std::vector<cell> > (rows+1, row);
        
        (*trellis)[0][0].setTraceback(NONE);
        (*trellis)[0][1].setScore(gap);
        (*trellis)[0][1].setTraceback(LEFT);
        (*trellis)[1][0].setScore(gap);
        (*trellis)[1][0].setTraceback(UP);
        
        for(size_t i=2;i<=columns;i++){
            (*trellis)[0][i].setScore((*trellis)[0][i-1].getScore()+gap_ext);
            (*trellis)[0][i].setTraceback(LEFT);
        }
        
        for(size_t j=2;j<=rows;j++){
            (*trellis)[j][0].setScore((*trellis)[j-1][0].getScore()+gap_ext);
            (*trellis)[j][0].setTraceback(UP);
        }
        return;
    }
    
    
    void alignment::printTrellis(){
        for (size_t row=0;row<trellis->size();++row){
            for (size_t col = 0; col<(*trellis)[row].size();++col){
                std::cout << (*trellis)[row][col].getScore() << " " << (*trellis)[row][col].getTraceback() << "\t";
            }
            std::cout << std::endl;
        }
    }
    
    double alignment::align(sequence& tgt, sequence& qu, alignment_type typ, double mtch, double mismtch, double gp, double gp_ext){
        setTarget(tgt);
        setQuery(qu);
        setMatch(mtch);
        setMismatch(mismtch);
        setGap(gp);
        setGapExt(gp_ext);
        return align(typ);
    }


    
    double alignment::align (alignment_type typ){  //(sequence &top, sequence &side, double gap, double mismatch, double match){
        
        int top_size=target->size();
        int side_size=query->size();
        
        _initTrellis(side_size, top_size);
        double max_score(-INFINITY);
        size_t max_row_index(0);
        size_t max_col_index(0);
                
        //Fill
        for(size_t top_index=1;top_index<=top_size;top_index++){
            for(size_t side_index=1;side_index<=side_size;side_index++){
                //table[side_index-1][top_index-1].setDiag();
                
                short targetLetter = query->seqValue(side_index-1);
                short queryLetter =  target->seqValue(top_index-1);
                
                double diagonal_score = (*trellis)[side_index-1][top_index-1].getScore() + (*mMatrix)[queryLetter][targetLetter];
                
                
                
                double left_score=(*trellis)[side_index][top_index-1].getScore()+_getGapScore(side_index, top_index-1, LEFT);
                
                
                double up_score=(*trellis)[side_index-1][top_index].getScore()+_getGapScore(side_index-1, top_index, UP);
                
                
                
                if (typ == cLocal){
                    if (diagonal_score<0.0){
                        diagonal_score=0.0;
                    }
                    else if (diagonal_score>max_score){
                        max_score = diagonal_score;
                        max_row_index=side_index;
                        max_col_index=top_index;
                    }
                    
                    if (left_score<0.0){
                        left_score=0.0;
                    }
                    else if (left_score>max_score){
                        max_score = left_score;
                        max_row_index=side_index;
                        max_col_index=top_index;
                    }
                    
                    if (up_score<0.0){
                        up_score=0.0;
                    }
                    else if (up_score>max_score){
                        max_score = up_score;
                        max_row_index=side_index;
                        max_col_index=top_index;
                    }
                }
                
                
                if (diagonal_score>=up_score && diagonal_score >= left_score){
                    (*trellis)[side_index][top_index].setScore(diagonal_score);
                    (*trellis)[side_index][top_index].setTraceback(DIAGONAL);
                                    }
                else if(up_score>=left_score){
                    (*trellis)[side_index][top_index].setScore(up_score);
                    (*trellis)[side_index][top_index].setTraceback(UP);
                }
                else{
                    (*trellis)[side_index][top_index].setScore(left_score);
                    (*trellis)[side_index][top_index].setTraceback(LEFT);
                }
            }
        }
        
        std::string align_top;
        std::string align_side;
        std::string align_status;
        
        size_t top_iter;
        size_t side_iter;
        
        if (typ == cGlobal){
            top_iter=top_size;
            side_iter=side_size;
        }
        else if (typ == cLocal){
            top_iter=max_col_index;
            side_iter=max_row_index;
        }
        
        
        while (top_iter!=0 || side_iter!=0){
            
            std::cout << side_iter << "," << top_iter << "\t" << (*trellis)[side_iter][top_iter].getScore() << std::endl;
            
            switch ((*trellis)[side_iter][top_iter].getTraceback()){
                case DIAGONAL:
                    align_side+=query->getSymbol(--side_iter);
                    align_top+=target->getSymbol(--top_iter);
                    align_status+= (query->seqValue(side_iter) != target->seqValue(top_iter)) ? ' ' : '|';
                    break;
                case UP:                    
                    align_side+=query->getSymbol(--side_iter);
                    align_top+="-";
                    align_status+=" ";
                    break;
                case LEFT:
                    //cout << "Case LEFT:" << side_size << "\t" << top_size << endl;
                    
                    align_top+=target->getSymbol(--top_iter);
                    align_side+="-";
                    align_status+=" ";
                    break;
                default:
                    //cout << "We have a problem Houston\n" << endl;
                    break;
            }
            
            
            if (typ == cLocal && fabs((*trellis)[side_iter][top_iter].getScore())<0.001){
                break;
            }
        }
        
        std::reverse(align_top.begin(), align_top.end());
        std::reverse(align_side.begin(), align_side.end());
        std::reverse(align_status.begin(), align_status.end());
        
        std::cout << align_top << std::endl << align_status << std::endl << align_side << std::endl;

        printTrellis();
        return (*trellis)[side_size][top_size].getScore();
    }
    
    
    double alignment::calcLambda(std::vector<double>& frequencies){
        
        return 0;
    }
    
    double alignment::estimateLambda(){
        size_t alphaSize = mMatrix->size();
        double estimatedFreq= 1/static_cast<double>(alphaSize);
        
        double expectedScore(0);
        for ( size_t row=0; row<alphaSize; ++row){
            for ( size_t col=0; col<alphaSize; ++col){
                
            }
        }
    }
    
    
    double alignment::calc_pvalue(){
        std::vector<double> scores;
        
        for(int i=0;i<100;++i){
            sequence* oldQuery=query;
            sequence shuffled=shuffle(query);
            setQuery(shuffled);
            double score = align(cGlobal);
            scores.push_back(score);
            std::cout<< score << std::endl;
            query=oldQuery;
        }
        
        return 0.0;
    }

    
    
//    alignment local(std::string &top, std::string &side,int match,int mismatch, int gap, int min_score){
//        int top_size=top.size();
//        int side_size=side.size();
//        int max_top_index;
//        int max_side_index;
//        int max_score=0;
//        alignment results;
//        
//        //Initialize
//        cell x;
//        x.traceback=NONE;
//        x.score=0;
//        std::vector<cell> row (top_size+1,x);
//        std::vector<std::vector<cell> > table (side_size+1,row);
//        
//        //Fill
//        for(int top_index=1;top_index<=top_size;top_index++){
//            for(int side_index=1;side_index<=side_size;side_index++){
//                int diagonal_score;
//                int up_score;
//                int left_score;
//                
//                diagonal_score = (top[top_index-1] == side[side_index-1]) ? table[side_index-1][top_index-1].score+match : table[side_index-1][top_index-1].score-mismatch;
//                left_score=table[side_index][top_index-1].score-gap;
//                up_score=table[side_index-1][top_index].score-gap;
//                
//                
//                if (diagonal_score>=up_score && diagonal_score >= left_score){
//                    if (diagonal_score>0){
//                        table[side_index][top_index].score=diagonal_score;
//                        table[side_index][top_index].traceback=DIAGONAL;
//                    }
//                }
//                else if(up_score>=left_score){
//                    if (up_score>0){
//                        table[side_index][top_index].score=up_score;
//                        table[side_index][top_index].traceback=UP;
//                    }
//                }
//                else{
//                    if (left_score>0){
//                        table[side_index][top_index].score=left_score;
//                        table[side_index][top_index].traceback=LEFT;
//                    }
//                }
//                
//                if (table[side_index][top_index].score>max_score){
//                    max_side_index=side_index;
//                    max_top_index=top_index;
//                    max_score=table[side_index][top_index].score;
//                }
//            }
//        }
//        
//        /*
//         //Print table
//         for(int side_index=0;side_index<=side_size;side_index++){
//         for(int top_index=0;top_index<=top_size;top_index++){
//         cout << table[side_index][top_index].score << ",";
//         }
//         cout << endl;
//         }
//         */
//        
//        std::string align_top;
//        std::string align_side;
//        
//        if (max_score<min_score){
//            results.score=max_score;
//            return results;
//        }
//        
//        //cout << "Side_INDEX:\t" << max_side_index << endl;
//        //cout << "Top_INDEX:\t" << max_top_index << endl;
//        
//        if (max_side_index<side_size){
//            int difference=side_size-max_side_index;
//            results.end_position=max_top_index+difference;
//            results.begin_position=results.end_position-(side_size-1);
//        }
//        else{
//            results.end_position=max_top_index;
//            results.begin_position=max_top_index-(side_size-1);
//        }
//        
//        
//        top_size=max_top_index;
//        side_size=max_side_index;
//        int pointer=table[side_size][top_size].traceback;
//        
//        //cout << side_size << "\t" << top_size << endl;
//        while (pointer!=NONE){
//            //cout << "SIDE_SIZE\t" << side_size << endl;
//            //cout << "TOP_SIZE\t" << top_size << endl;
//            //cout << "TRACEBACK\t" << table[side_size][top_size].traceback << endl;
//            //cout << "SCORE\t" << table[side_size][top_size].score << endl;
//            
//            
//            switch (pointer){
//                case DIAGONAL:
//                    //cout << "Case Diagonal:" << side_size << "\t" << top_size << endl;
//                    align_side+=side[side_size-1];
//                    align_top+=top[top_size-1];
//                    top_size--;
//                    side_size--;
//                    pointer=table[side_size][top_size].traceback;
//                    break;
//                case UP:
//                    //cout << "Case UP:" << side_size << "\t" << top_size << endl;
//                    
//                    align_side+=side[side_size-1];
//                    align_top+="-";
//                    side_size--;
//                    pointer=table[side_size][top_size].traceback;
//                    break;
//                case LEFT:
//                    //cout << "Case LEFT:" << side_size << "\t" << top_size << endl;
//                    
//                    align_top+=top[top_size-1];
//                    align_side+="-";
//                    top_size--;
//                    pointer=table[side_size][top_size].traceback;
//                    break;
//                default:
//                    //cout << "We have a problem Houston\n" << endl;
//                    break;
//            }
//        }
//        
//        
//        results.target=align_top;
//        results.query=align_side;
//        results.score=max_score;
//        
//        //cout << align_top << endl << align_side <<endl<< max_score <<endl;
//        
//        return results;
//    }
}

